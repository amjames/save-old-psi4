/*
 * @BEGIN LICENSE
 *
 * Psi4: an open-source quantum chemistry software package
 *
 * Copyright (c) 2007-2016 The Psi4 Developers.
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *
 * @END LICENSE
 */

#include <boost/filesystem.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>
#include <cmath>
#include <utility>

#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <physconst.h>
#include <libmints/mints.h>
#include <psi4-dec.h>
#include <libpsi4util/libpsi4util.h>
#include <libmints/mints.h>


#include "denscubeprop.h"
#include "csg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

DensityCubeProperties::DensityCubeProperties(SharedWavefunction wfn) :
  Prop(wfn)
{
  common_init();
}
DensityCubeProperties::~DensityCubeProperties()
{
}
void DensityCubeProperties::common_init()
{
  Options &options = Process::environment.options;
  print_ = options.get_int("PRINT");
  cubepath_ =  options.get_str("CUBEPROP_FILEPATH");
  do_pop_analysis_ = options.get_bool("CUBEPROP_EUD_POP");

  for(size_t ind = 0; ind < options["CUBEPROP_TASKS"].size(); ind++){
    std::string task_name = options["CUBEPROP_TASKS"][ind].to_string();
    add(task_name);
  }

  grid_ = boost::shared_ptr<CubicScalarGrid>(new CubicScalarGrid(basisset_, options));
  grid_->set_filepath(options.get_str("CUBEPROP_FILEPATH"));

}
void DensityCubeProperties::print_header()
{
    outfile->Printf( "  ********************************************\n");
    outfile->Printf( "  ** ==> Density Cube Properties (v1.0) <== **\n");
    outfile->Printf( "  **      Effectively unpaired density      **\n" )
    outfile->Printf( "  **    analysis and CubeFile generation    **\n");
    outfile->Printf( "  **           by Andrew M. James           **\n");
    outfile->Printf( "  ********************************************\n\n");
    grid_->print_header();
    outfile->Flush();
}
void DensityCubeProperties::print_EUD_summary(
    std::vector<std::tuple<double,double,double,int,int>> info)
{
  std::sort(info.begin(),info.end());
  outfile->Printf("   Largest Contributions\n");
  outfile->Printf("   orb    NOON(a)    NOON(b)      u.p \n");
  outfile->Printf("-------------------------------------------\n");
  char** labels = basisset_->molecule()->irrep_labels();
  double total_up = 0.00;
  double total_near_occ=0.00;
  double total_near_uocc=0.00;
  for(auto tup: info)
  {
    double nu = std::get<0>(tup);
    total_up += nu;
    if(nu > 0.5){
      outfile->Printf(" %2d %4s   %9.6f  %9.6f   %9.6f\n",
        std::get<4>(tup),labels[std::get<3>(tup)],
        std::get<1>(tup),std::get<2>(tup));
    } else {
      double na = std::get<1>(tup);
      double nb = std::get<2>(tup);
      if((na+nb) > 1.5)
        total_near_occ+=nu;
      else
        if((na+nb) <0.5)
          total_near_uocc+=nu;
    }
  }
  outfile->Printf("-------------------------------------------\n");
  outfile->Printf(    "Total up e-             = %9.6f\n",total_up);
  outfile->Printf(    "Nearly DOCC             = %9.6f\n",total_near_occ);
  outfile->Printf(    "Nearly UOCC             = %9.6f\n",total_near_uocc);
}
SharedMatrix DensityCubeProperties::Unpaired_D_so2ao(SharedMatrix Dso)
{
  std::string new_name(Dso->name()+"->ao");
  double* temp = new double[AO2USO_->max_ncol() * AO2USO_->max_nrow()];
  //  Ca_so _
  SharedMatrix Dao(new Matrix(new_name,
      basisset_->nbf(),basisset_->nbf()));

  int symm = Dso->symmetry();
  for (int h = 0; h < AO2USO_->nirrep(); h++){

    int nao = AO2USO_->rowspi()[0];
    int nsol = AO2USO_->colspi()[h];
    int nsor = AO2USO_->colspi()[h^symm];

    if (!nsol || !nsor )continue;
    double** Ulp = AO2USO_->pointer(h);
    double** Urp = AO2USO_->pointer(h^symm);
    double** DSOp = Dso->pointer(h^symm);
    double** DAOp = Dao->pointer();

    C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
    C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
  }
  delete[] temp;
  return Dao;
}
SharedMatrix DensityCubeProperties::Du_s_mo()
{
  outfile->Printf("       Effectively unpaired-electron analysis\n");
  outfile->Printf("                 f_u(na,nb) = \n");
  outfile->Printf("            (na+nb)^2*[(1-na)+(1-nb)]^2\n");
  outfile->Printf("==================================================\n");
  std::pair<SharedMatrix,SharedVector> NaOa = Na_mo();
  std::pair<SharedMatrix,SharedVector> NbOb = Nb_mo();
  SharedMatrix Na = NaOa.first;
  SharedMatrix Nb = NbOb.first;
  SharedVector Oa = NaOa.second;
  SharedVector Ob = NbOb.second;
  SharedMatrix Du(new Matrix("Effectively Unpaired Density",Na->rowspi(),Na->colspi()));
  Du->zero();
  std::vector<std::tuple<double,double,double,int,int>> f_s_info;
  //totals for my own analysis
  for(int h =0; h < Oa->nirrep(); h++){
    for(int i =0; i < Oa->dimpi()[h]; i++){
      double na=Oa->get(h,i);
      double nb=Ob->get(h,i);
      double nt=na+nb;
      double del=(1-na)+(1-nb);
      double tot =(nt*nt)*(del*del);
      Du->set(h,i,i,tot);
      f_s_info.push_back(std::make_tuple(tot,na,nb,h,i));
    }
  }
  print_EUD_summary(f_s_info);
  outfile->Printf("Tr(Du(NO)                        = %9.3f\n",Du->trace());
  Du->back_transform(Na);
  outfile->Printf("Tr(Du(MO))                       = %9.3f\n",Du->trace());
  return Du;
}
SharedMatrix DensityCubeProperties::Du_s_so()
{
  SharedMatrix Du_mo = Du_s_mo();
  SharedMatrix Du_so(new Matrix("Du(SO)",Ca_so_->rowspi(),Ca_so_->colspi()));

  int symm = Du_mo->symmetry();
  double* temp = new double[Ca_so_->max_ncol() * Cb_so_->max_nrow()];
  for(int h=0; h < Du_mo->nirrep();h++){
    int nmol = Ca_so_->colspi()[h];
    int nmor = Ca_so_->colspi()[h^symm];
    int nsol = Ca_so_->rowspi()[h];
    int nsor = Ca_so_->rowspi()[h^symm];
    if(!nmol|| !nmor || !nsol || !nsor)continue;
    double** Clp = Ca_so_->pointer(h);
    double** Crp = Ca_so_->pointer(h^symm);
    double** Dmop = Du_mo->pointer(h^symm);
    double** Dsop = Du_so->pointer(h^symm);
    C_DGEMM('N','T',nmol,nsor,nmor,1.0,Dmop[0],nmor,Crp[0],nmor,0.0,temp,nsor);
    C_DGEMM('N','N',nsol,nsor,nmol,1.0,Clp[0],nmol,temp,nsor,0.0,Dsop[0],nsor);
  }
  delete[] temp;
  outfile->Printf("Tr(Du(SO))                       = %9.3f\n",Du_so->trace());
  return Du_so;
}
SharedMatrix DensityCubeProperties::Du_s_ao()
{
  return Unpaired_D_so2ao(Du_so);
}
void DensityCubeProperties::mulliken_EUD(SharedMatrix Du_ao)
{
    outfile->Printf("   Mulliken Pop analysis: (a.u.)\n");
    outfile->Printf("----------------------------------\n");
    boost::shared_ptr<Molecule> mol = basisset_->molecule();
    double* Qt = new double[mol->natom()];

    ::memset(Qt,'\0',mol->natom()*sizeof(double));

    SharedMatrix DuS(Du_ao->clone());
    boost::shared_ptr<OneBodyAOInt>overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);
    DuS->gemm(false,false,1.0,Du_ao,S,0.0);

    //compute D_up*S
    for(int mu =0; mu < basisset_->nbf(); mu++){
      int shellmu = basisset_->function_to_shell(mu);
      int A = basisset_->shell_to_center(shellmu);
      Qt[A] += DuS->get(0,mu,mu);
    }
    outfile->Printf( "Center  Symbol Unpaired Pop\n");
    double sumt = 0.00;
    for(int A =0; A < mol->natom(); A++){
      outfile->Printf("%5d    %2s    %8.5f\n", \
          A+1,mol->label(A).c_str(),Qt[A]);
      sumt+= Qt[A];
    }
    outfile->Printf("----------------------------------\n");
    outfile->Printf("Total unpaired = %8.5f\n",sumt);
    outfile->Printf("----------------------------------\n");
    delete[] Qt;
}
SharedMatrix DensityCubeProperties::compute_EUD_S()
{
  SharedMatrix Ds_ao = Du_s_ao();
  if(do_pop_analysis_){
    mulliken_EUD(Ds_ao);
  }
  return Ds_ao;
}
void DensityCubeProperties::compute_EUD(std::string type)
{
  if (same_dens_)
    throw PSIEXCEPTION("DENSCUBEPROP: a/b density is the same requesting EUD makes no sense");

  if(type == "S" || type=="ALL"){
    SharedMatrix Du= compute_EUD_S();
    grid_->compute_density(Du_ao,"EUD_S")
    SharedMatrix DuS(Du_ao->clone());
    boost::shared_ptr<OneBodyAOInt>overlap(integral_->ao_overlap());
    SharedMatrix S(new Matrix("S",basisset_->nbf(),basisset_->nbf()));
    overlap->compute(S);
    DuS->gemm(false,false,1.0,Du,S,0.0);
    grid_->compute_density(DuS,"EUD_SmOL");
  }
}
void DensityCubeProperties::compute()
{
  print_header();
  std::stringstream ss;
  ss << cubepath_ << "/" << "geom.xyz";
  boost::filesystem::path data_dir(cubepath_);
  if(not boost::filesystem::is_directory(data_dir)){
    outfile->Printf("CubeFile output path (\"%s\") is not found. Creating it...",cubepath_.c_str());
    auto ret = boost::filesystem::create_directory(cubepath_);
    if(not ret){
      throw PSIEXCEPTION(" A problem occurred, could not create path\n");
    }
    outfile->Printf(" Done!\n");
  }
  basisset_->molecule()->save_xyz_file(ss.str());
  if (tasks_.count("EUD")){
    compute_EUD("ALL");
  }
  if (tasks_.count("EUD_S")){
    compute_EUD("S");
  }
}
}//namespace

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

  for(size_t ind = 0; ind < options["CUBEPROP_TASKS"].size(); ind++){
    std::string task_name = options["CUBEPROP_TASKS"][ind].to_string();
    add(task_name);
  }

  boost::shared_ptr<Molecule> mol = basisset_->molecule();
  grid_ = boost::shared_ptr<CubicScalarGrid>(new CubicScalarGrid(basisset_, options));
  grid_->set_filepath(options.get_str("CUBEPROP_FILEPATH"));

}

void DensityCubeProperties::print_header()
{
    outfile->Printf( "  ==> Density based plotting tool (v2.0) <==\n\n");
    grid_->print_header();
    outfile->Flush();
}

SharedMatrix DensityCubeProperties::Gen_D_mo2so(SharedMatrix Dmo)
{
  std::string new_name(Dmo->name()+"->SO");
    SharedMatrix Dso = SharedMatrix(new Matrix(new_name,
          Cb_so_->nirrep(),Cb_so_->nrow(),Cb_so_->ncol()));
  for(int h =0; h < Dmo->nirrep(); h++){
    int nmo = Cb_so_->colspi()[h];
    int nso = Cb_so_->rowspi()[h];

    if (!nmo || !nso) continue;
    double** Dmop = Dmo->pointer(h);
    double** Cp  = Cb_so_->pointer(h);
    double** Dsop = Dso->pointer(h);
    C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Dmop[0],nmo,0.0,Dsop[0],nmo);
  }
  return Dso;

}

SharedMatrix DensityCubeProperties::Gen_D_so2ao(SharedMatrix Dso)
{
  std::string new_name(Dso->name()+"->ao");
  SharedMatrix Dao(new Matrix(new_name,
      Cb_so_->nrow(),Cb_so_->ncol()));

  int offset = 0;
  for (int h = 0; h < Cb_so_->nirrep(); h++){

    int ncol = Cb_so_->ncol();
    int nmo = Cb_so_->colspi()[h];
    int nso = AO2USO_->colspi()[h];
    int nao = AO2USO_->rowspi()[h];

    if (!nmo || !nso || !nao)continue;
    double** Dsop = Dso->pointer(h);
    double** Up = AO2USO_->pointer(h);
    double** Daop = Dao->pointer(h);

    C_DGEMM(
        'N','N',
        nao,nmo,nso,
        1.0,Up[0],nso,
        Dsop[0],nmo,0.0,
        &Daop[0][offset],ncol
        );
    offset += nmo;
  }
  return Dao;


}

SharedMatrix DensityCubeProperties::Gen_D_mo2ao(SharedMatrix Dmo)
{
  //first half mo->so
  SharedMatrix Dso = Gen_D_mo2so(Dmo);
  //then  so->ao
  SharedMatrix Dao = Gen_D_so2ao(Dso);
  std::string new_name (Dmo->name()+"->ao");
  Dao->set_name(new_name);
  return Dao;
}

SharedMatrix DensityCubeProperties::Dt_ao()
{
  double* temp = new double[AO2USO_->max_ncol()*AO2USO_->max_nrow()];
  SharedMatrix D = SharedMatrix(new Matrix("Dt (ao basis)", basisset_->nbf(),basisset_->nbf()));
  SharedMatrix Dt = Dt_so();
  int symm = Da_so_->symmetry();
  for( int h = 0; h < AO2USO_->nirrep(); ++h ){
    int nao = AO2USO_->rowspi()[0];
    int nsol = AO2USO_->colspi()[h];
    int nsor = AO2USO_->colspi()[h^symm];
    if (!nsol || !nsor) continue;
    double** Ulp = AO2USO_->pointer(h);
    double** Urp = AO2USO_->pointer(h^symm);
    double** DSOp = Dt->pointer(h^symm);
    double** DAOp = D->pointer();
    C_DGEMM('N','T',nsol,nao,nsor,1.0,DSOp[0],nsor,Urp[0],nsor,0.0,temp,nao);
    C_DGEMM('N','N',nao,nao,nsol,1.0,Ulp[0],nsol,temp,nao,1.0,DAOp[0],nao);
  }
  delete[] temp;
  return D;
}

void DensityCubeProperties::compute_densities(std::string key)
{
  SharedMatrix Dt = Dt_ao();
  std::string total = key+"_total";
  grid_->compute_density(Dt,total);
  if (!same_dens_){
    std::string alpha = key+"_alpha";
    SharedMatrix Da = Da_ao();
    grid_->compute_density(Da,alpha);
    std::string beta = key+"_beta";
    SharedMatrix Db = Db_ao();
    grid_->compute_density(Db,beta);
    SharedMatrix Ds = Da->clone();
    Ds->subtract(Db);
    std::string spin = key+"_spin";
    grid_->compute_density(Ds,spin);
  }
}

void DensityCubeProperties::compute_EUD(std::string type, bool atomic_contrib)
{
  if (same_dens_)
    throw PSIEXCEPTION("DENSCUBEPROP: a/b density is the same requesting EUD makes no sense");
  if (type == "D-mat"){
    std::string typekey("EUD_D-mat");
    /*
     * Evaluate
     * D(i,j) = 2(Da(i,j)+Db(i,j)) - sum_k[(Da(i,k)+Db(i,k))*(Db(k,j)+Db(k,j))]
     */
    SharedMatrix Dt = Da_mo();
    SharedMatrix Dmo = Dt->clone();
    Dmo->scale(2.0);
    SharedMatrix Dt2 =Dt->clone();
    Dt2->gemm(false,false,-1.0,Dt,Dt,0.0);
    Dmo->add(Dt2);

    // transform Dmo ->ao basis for plotting
    SharedMatrix Dao = Gen_D_mo2ao(Dmo);
    grid_->compute_density(Dao,typekey);
  }
  if(type == "D-NatOrb"){
    std::string typekey("EUD_D-NatOrb");
    std::pair<SharedMatrix, SharedVector> pair = Nt_mo();
    SharedMatrix N_mo = pair.first;
    SharedVector O = pair.second;
    SharedVector f_d(O->clone());
    outfile->Printf("        Number of effectively unpaired Electrons\n");
    outfile->Printf("========================================================\n");
    outfile->Printf(" index       No occ num         Nu\n");
    outfile->Printf("---------------------------------------\n");
    int index = 0;
    for(int h = 0; h < O->nirrep(); h++){
      for(int i =0; i < O->dim(h); i++){
        double nd = O->get(h,i);
        f_d->set(h,i,2*nd-(nd*nd));
        outfile->Printf(" %5d     %9.3f    %9.3f\n",index,nd,f_d->get(h,i));
      }
    }
    SharedMatrix Dmo(new Matrix(N_mo->nirrep(),N_mo->rowspi(),N_mo->colspi()));
    Dmo->set_diagonal(f_d);
    Dmo->back_transform(N_mo);
    SharedMatrix Dao = Gen_D_mo2ao(Dmo);

    grid_->compute_density(Dao,typekey);
  }
}

void DensityCubeProperties::compute()
{
  print_header();
  if (tasks_.count("DENSITY"))
    compute_densities("D");
  if (tasks_.count("EUD")){
    compute_EUD("D-mat");
    compute_EUD("D-NatOrb");
  }
}
}//namespace

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

SharedMatrix DensityCubeProperties::Unpaired_Da_mo2so(SharedMatrix Dmo)
{
  std::string new_name(Dmo->name()+"->SO");
  SharedMatrix Dso = SharedMatrix(new Matrix("temp",
      Ca_so_->nirrep(),Ca_so_->rowspi(),Cb_so_->colspi()));
  for(int h =0; h < Dmo->nirrep(); h++){
    int nmo = Ca_so_->colspi()[h];
    int nso = Ca_so_->rowspi()[h];

    if (!nmo || !nso) continue;
    double** Dmop = Dmo->pointer(h);
    double** Cp  = Ca_so_->pointer(h);
    double** Dsop = Dso->pointer(h);
    C_DGEMM('N','N',nso,nmo,nmo,1.0,Cp[0],nmo,Dmop[0],nmo,0.0,Dsop[0],nmo);
  }

  return Dso;

}

SharedMatrix DensityCubeProperties::Unpaired_Db_mo2so(SharedMatrix Dmo)
{
  std::string new_name(Dmo->name()+"->SO");
  SharedMatrix Dso = SharedMatrix(new Matrix("temp",
      Cb_so_->nirrep(),Cb_so_->rowspi(),Cb_so_->colspi()));
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

SharedMatrix DensityCubeProperties::Unpaired_Da_so2ao(SharedMatrix Dso)
{
  std::string new_name(Dso->name()+"->ao");
  SharedMatrix Dao(new Matrix(new_name,
      Ca_so_->nrow(),Ca_so_->ncol()));

  int offset = 0;
  for (int h = 0; h < Ca_so_->nirrep(); h++){

    int ncol = Ca_so_->ncol();
    int nmo = Ca_so_->colspi()[h];
    int nso = AO2USO_->colspi()[h];
    int nao = AO2USO_->rowspi()[h];

    if (!nmo || !nso || !nao)continue;
    double** Dsop = Dso->pointer(h);
    double** Up = AO2USO_->pointer(h);
    double** Daop = Dao->pointer(0);

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

SharedMatrix DensityCubeProperties::Unpaired_Db_so2ao(SharedMatrix Dso)
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
    double** Daop = Dao->pointer(0);

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

SharedMatrix DensityCubeProperties::Unpaired_Da_mo2ao(SharedMatrix Dmo)
{
  //first half mo->so
  SharedMatrix Dso = Unpaired_Da_mo2so(Dmo);
  //then  so->ao
  SharedMatrix Dao = Unpaired_Da_so2ao(Dso);
  std::string new_name (Dmo->name()+"->ao");
  Dao->set_name(new_name);
  return Dao;
}

SharedMatrix DensityCubeProperties::Unpaired_Db_mo2ao(SharedMatrix Dmo)
{
  //first half mo->so
  SharedMatrix Dso = Unpaired_Db_mo2so(Dmo);
  //then  so->ao
  SharedMatrix Dao = Unpaired_Db_so2ao(Dso);
  std::string new_name (Dmo->name()+"->ao");
  Dao->set_name(new_name);
  return Dao;
}

// prints table of NO occupation numbers, contributions to NUp with
// alpha/beta HONO-LUNO info
void DensityCubeProperties::print_Num_UP_info(
    std::vector<boost::tuple<int,int,double,int,int,double,double>>Upmetric,
    std::string fdef
  )
{
  //tuple holds
  //aNO(h),aNO(i)aNOON, ,bNOON, Nup-contribution
  char** labels = basisset_->molecule()->irrep_labels();
  //construct header: line indent 4sp.  +40 chars
  outfile->Printf("    Nonzero contributions: f_up(n) = \n");
  //line indent for function definition
  std::string func_def("    ");
  int indent = ((int)(40 - fdef.length())/2);
  for(int spno =0; spno < indent -1; spno++) func_def += " ";
  func_def += fdef;
  func_def += "\n";
  outfile->Printf(func_def.c_str());
  outfile->Printf("=================================================\n");
  outfile->Printf("                   Occupation          Unpaired \n");
  outfile->Printf("                     Number          contribution\n");
  outfile->Printf("-------------------------------------------------\n");
  int nalpha = wfn_->nalpha();
  int nbeta  = wfn_->nbeta();
  double Nup_total = 0.00;
  for(int index = 0; index < wfn_->nmo(); ++index){
    double Nup_cont = boost::get<6>(Upmetric[index]);
    if(Nup_cont != 100.00){
      int irr_a = boost::get<0>(Upmetric[index]);
      int idx_a = boost::get<1>(Upmetric[index]);
      double NOONa = boost::get<2>(Upmetric[index]);
      int irr_b = boost::get<3>(Upmetric[index]);
      int idx_b = boost::get<4>(Upmetric[index]);
      double NOONb = boost::get<5>(Upmetric[index]);
      if(index < nalpha){
        outfile->Printf("alpha-#%2d   %4d%3s %8.3f \n",
            index, idx_a,labels[irr_a],NOONa);
      }else{
        outfile->Printf("alpha-#%2d   %4d%3s %8.3f \n",
            index, idx_a,labels[irr_a],NOONa);
      }
      if(index < nalpha){
        outfile->Printf(" beta-#%2d   %4d%3s %8.3f   %8.3f \n",
            index, idx_b,labels[irr_b],NOONb,Nup_cont);
      }else{
        outfile->Printf(" beta-#%2d   %4d%3s %8.3f   %8.3f\n",
            index, idx_b,labels[irr_b],NOONb,Nup_cont);
      }
      Nup_total += Nup_cont;
    }
  }
  outfile->Printf("-------------------------------------------------\n");
  outfile->Printf("TOTAL Nu = %8.3f e-\n",Nup_total);
  outfile->Printf("-------------------------------------------------\n");
}



                  // --> EUD S Function  <--//
std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sa_mo()
{
  std::pair<SharedMatrix,SharedVector> pair = Na_mo();
  SharedMatrix N = pair.first;
  SharedVector O = pair.second;
  SharedVector Os(new Vector("Alpha Nup(S)",O->nirrep(),O->dimpi()));
  SharedMatrix Smo(new Matrix("Alpha unpaired Denisty(S)",
        N->nirrep(),N->rowspi(),N->colspi()));
  for(int h=0; h<O->nirrep(); ++h){
    for(int i =0; i < O->dimpi()[h]; ++i){
      double ni = O->get(h,i);
      double ns = ni*ni*(1-ni)*(1-ni);
      Os->set(h,i,ns);
    }
  }
  Smo->set_diagonal(Os);
  Smo->back_transform(N);
  return make_pair(Smo,Os);
}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sb_mo()
{
  std::pair<SharedMatrix,SharedVector> pair = Nb_mo();
  SharedMatrix N = pair.first;
  SharedVector O = pair.second;
  SharedVector Os(new Vector("Beta Nup(S)",O->nirrep(),O->dimpi()));
  SharedMatrix Smo(new Matrix("Beta unpaired Denisty(S)",
        N->nirrep(),N->rowspi(),N->colspi()));
  for(int h=0; h<O->nirrep(); ++h){
    for(int i =0; i < O->dimpi()[h]; ++i){
      double ni = O->get(h,i);
      double ns = ni*ni*(1-ni)*(1-ni);
      Os->set(h,i,ns);
    }
  }
  Smo->set_diagonal(Os);
  Smo->back_transform(N);
  return make_pair(Smo,Os);

}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sa_so()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Sa_mo();
  SharedMatrix Smo = pair.first;
  SharedVector Os = pair.second;
  SharedMatrix Sso = Unpaired_Da_mo2so(Smo);
  return make_pair(Sso,Os);

}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sb_so()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Sb_mo();
  SharedMatrix Smo = pair.first;
  SharedVector Os = pair.second;
  SharedMatrix Sso = Unpaired_Db_mo2so(Smo);
  return make_pair(Sso,Os);


}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sa_ao()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Sa_so();
  SharedMatrix Sso = pair.first;
  SharedVector Os = pair.second;
  SharedMatrix Sao = Unpaired_Da_so2ao(Sso);
  return make_pair(Sso,Os);
}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Sb_ao()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Sb_so();
  SharedMatrix Sso = pair.first;
  SharedVector Os = pair.second;
  //SharedMatrix Sao = Unpaired_Db_so2ao(Sso);
  SharedMatrix Sao = Unpaired_Db_so2ao(pair.first);
  return make_pair(Sso,Os);
}

std::pair<SharedMatrix, SharedVector> DensityCubeProperties::compute_EUD_S()
{
  std::pair<SharedMatrix,SharedVector> pair_a = EUD_Sa_ao();
  std::pair<SharedMatrix,SharedVector> pair_b = EUD_Sb_ao();
  outfile->Printf("Sa\n");
  SharedMatrix Sa = pair_a.first;
  Sa->print_out();
  outfile->Printf("Sb\n");
  SharedMatrix Sb = pair_b.first;
  Sb->print_out();
  SharedVector Os_a = pair_a.second;
  SharedVector Os_b = pair_b.second;
  SharedMatrix St = Sa->clone();
  St->add(Sb);
  outfile->Printf("St\n");
  SharedVector Os_t(Os_a->clone());
  Os_t->add(Os_b);
  std::pair<SharedMatrix,SharedVector> no_pair_a = Na_mo();
  std::pair<SharedMatrix,SharedVector> no_pair_b = Nb_mo();
  SharedVector Oa = no_pair_a.second;
  SharedVector Ob = no_pair_b.second;
  outfile->Printf("Oa\n");
  Oa->print_out();
  outfile->Printf("Ob\n");
  //aNO(h),aNO(i)aNOON, ,bNOON, Nup-contribution
  std::vector<boost::tuple<int,int,double,int,int,double,double> > S_metric;
  for(int h = 0; h < Os_t->nirrep(); h++){
    for(int i = 0; i < Os_t->dimpi()[h]; i++){
      S_metric.push_back(boost::tuple<int,int,double,int,int,double,double>(
            i,h,Oa->get(h,i),i,h,Ob->get(h,i),Os_t->get(h,i)));
    }
  }
  DensityCubeProperties::print_Num_UP_info(S_metric,std::string("n^2*(2-n)^2"));
  return make_pair(St,Os_t);
}


                  // --> EUD U Function  <--//
std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ua_mo()
{
  std::pair<SharedMatrix,SharedVector> pair = Na_mo();
  SharedMatrix N = pair.first;
  SharedVector O = pair.second;
  SharedVector Ou(new Vector("Alpha Nup(U)",O->nirrep(),O->dimpi()));
  SharedMatrix Umo(new Matrix("Alpha unpaired Denisty(U)",
        N->nirrep(),N->rowspi(),N->colspi()));
  for(int h=0; h<O->nirrep(); ++h){
    for(int i =0; i < O->dimpi()[h]; ++i){
      double ni = O->get(h,i);
      double nu = std::min(ni,1-ni);
      Ou->set(h,i,nu);
    }
  }
  Umo->set_diagonal(Ou);
  Umo->back_transform(N);
  return make_pair(Umo,Ou);
}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ub_mo()
{
  std::pair<SharedMatrix,SharedVector> pair = Nb_mo();
  SharedMatrix N = pair.first;
  SharedVector O = pair.second;
  SharedVector Ou(new Vector("Alpha Nup(U)",O->nirrep(),O->dimpi()));
  SharedMatrix Umo(new Matrix("Alpha unpaired Denisty(U)",
        N->nirrep(),N->rowspi(),N->colspi()));
  for(int h=0; h<O->nirrep(); ++h){
    for(int i =0; i < O->dimpi()[h]; ++i){
      double ni = O->get(h,i);
      double nu = std::min(ni,1-ni);
      Ou->set(h,i,nu);
    }
  }
  Umo->set_diagonal(Ou);
  Umo->back_transform(N);
  return make_pair(Umo,Ou);
}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ua_so()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Ua_mo();
  SharedMatrix Umo = pair.first;
  SharedVector Ou = pair.second;
  SharedMatrix Uso = Unpaired_Da_mo2so(Umo);
  return make_pair(Uso,Ou);

}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ub_so()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Ub_mo();
  SharedMatrix Umo = pair.first;
  SharedVector Ou = pair.second;
  SharedMatrix Uso = Unpaired_Db_mo2so(Umo);
  return make_pair(Uso,Ou);

}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ua_ao()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Ua_so();
  SharedMatrix Umo = pair.first;
  SharedVector Ou = pair.second;
  SharedMatrix Uso = Unpaired_Da_so2ao(Umo);
  return make_pair(Uso,Ou);

}

std::pair<SharedMatrix,SharedVector> DensityCubeProperties::EUD_Ub_ao()
{
  std::pair<SharedMatrix,SharedVector> pair = EUD_Ub_so();
  SharedMatrix Umo = pair.first;
  SharedVector Ou = pair.second;
  SharedMatrix Uso = Unpaired_Db_so2ao(Umo);
  return make_pair(Uso,Ou);

}
std::pair<SharedMatrix,SharedVector> DensityCubeProperties::compute_EUD_U()
{
  std::pair<SharedMatrix,SharedVector> pair_a = EUD_Ua_ao();
  std::pair<SharedMatrix,SharedVector> pair_b = EUD_Ub_ao();
  SharedMatrix Ua = pair_a.first;
  SharedMatrix Ub = pair_b.first;
  SharedVector Ou_a = pair_a.second;
  SharedVector Ou_b = pair_b.second;
  SharedMatrix Ut = Ua->clone();
  Ut->add(Ub);
  SharedVector Ou_t(Ou_a->clone());
  Ou_t->add(Ou_b);
  std::pair<SharedMatrix,SharedVector> no_pair_a = Na_mo();
  std::pair<SharedMatrix,SharedVector> no_pair_b = Nb_mo();
  SharedVector Oa = no_pair_a.second;
  SharedVector Ob = no_pair_b.second;
  //aNO(h),aNO(i)aNOON, ,bNOON, Nup-contribution
  std::vector<boost::tuple<int,int,double,int,int,double,double> > U_metric;
  for(int h = 0; h < Ou_t->nirrep(); h++){
    for(int i = 0; i < Ou_t->dimpi()[h]; i++){
      U_metric.push_back(boost::tuple<int,int,double,int,int,double,double>(
            h,i,Oa->get(h,i),h,i,Ob->get(h,i),Ou_t->get(h,i)));
    }
  }
  DensityCubeProperties::print_Num_UP_info(U_metric,std::string("min(n,2-n)"));
  return make_pair(Ut,Ou_t);
}


void DensityCubeProperties::compute_EUD(std::string type, bool atomic_contrib)
{
  if (same_dens_)
    throw PSIEXCEPTION("DENSCUBEPROP: a/b density is the same requesting EUD makes no sense");
  /* std::string Skey("EUD_S"); */
  /* std::pair<SharedMatrix,SharedVector> S_pair = compute_EUD_S(); */
  /* std::string Ukey("EUD_U"); */
  /* std::pair<SharedMatrix,SharedVector> U_pair = compute_EUD_U(); */
  /* grid_->compute_density(S_pair.first,Skey); */
  /* grid_->compute_density(U_pair.first,Ukey); */
  outfile->Printf("Ca_so_\n");
  Ca_so_->print_out();

  outfile->Printf("S_so = so basis overlap\n");
  SharedMatrix S = overlap_so();
  S->print_out();

  if(wfn_->name() == "CCEnergyWavefunction"){
    set_Da_mo(wfn_->get_OPDM("A"));
    set_Db_mo(wfn_->get_OPDM("B"));
  }

  std::pair<SharedMatrix,SharedVector> NaOa = Na_mo();
  outfile->Printf("Na_mo: \n");
  NaOa.first->print_out();

  outfile->Printf("Oa:\n");
  NaOa.second->print_out();

  std::pair<SharedMatrix,SharedVector> NbOb = Nb_mo();
  outfile->Printf("Nb_mo: \n");
  NaOa.first->print_out();

  outfile->Printf("Oa: \n");
  NaOa.second->print_out();

  std::pair<SharedMatrix,SharedVector> Spair = compute_EUD_S();
  SharedVector NupS = Spair.second;
  double num_unpaired = 0.00;
  for(int h =0; h < NupS->nirrep(); ++h){
    for(int i=0; i <NupS->dimpi()[h];++i){
      double ns = NupS->get(h,i);
      oufile->Printf("NupS[%2d][%2d] = %9.3f\n",h,i,ns);
      num_unpaired+=ns;
    }
  }
  outfile->Printf("# effectively unpaired e- = %9.3f\n",num_upaired);

}

void DensityCubeProperties::compute()
{
  print_header();
  if (tasks_.count("EUD")){
    compute_EUD("ALL");
  }
}
}//namespace

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

#include <psi4-dec.h>

#include <libpsi4util/libpsi4util.h>
#include <libmints/mints.h>

#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>

#include "cubeprop.h"
#include "csg.h"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace psi;
using namespace boost;
using namespace std;

namespace psi {

CubeProperties::CubeProperties(SharedWavefunction wfn) :
    options_(Process::environment.options)
{
    basisset_ = wfn->basisset();
    boost::shared_ptr<PetiteList> PL(new PetiteList(basisset_,wfn->integral()));
    SharedMatrix SO2AO = PL->sotoao();
    S_ = wfn->S();

    Ca_ = wfn->Ca_subset("AO", "ALL");
    Da_ = wfn->Da_subset("AO");

    if (wfn->same_a_b_orbs()) {
        Cb_ = Ca_;
    } else {
        Cb_ = wfn->Cb_subset("AO", "ALL");
    }

    if (wfn->same_a_b_dens()) {
        Db_ = Da_;
    } else {
        Db_ = wfn->Db_subset("AO");
    }

    std::pair<SharedMatrix,SharedVector> alpha_NOpair =
      compute_NOpair(SO2AO,wfn);
    Na_ = alpha_NOpair.first;
    NOONa_ = alpha_NOpair.second;

    int nirrep = wfn->nirrep();
    Dimension nmopi = wfn->nmopi();
    // Gather orbital information
    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_a_.push_back(boost::tuple<double,int,int>(wfn->epsilon_a()->get(h,i),i,h));
        }
    }
    std::sort(info_a_.begin(), info_a_.end(), std::less<boost::tuple<double,int,int> >()); // Sort as in wfn
    for (int h = 0; h < nirrep; h++) {
        for (int i = 0; i < (int)nmopi[h]; i++) {
            info_b_.push_back(boost::tuple<double,int,int>(wfn->epsilon_b()->get(h,i),i,h));
        }
    }
    std::sort(info_b_.begin(), info_b_.end(), std::less<boost::tuple<double,int,int> >()); // Sort as in wfn

    common_init();
}
CubeProperties::~CubeProperties()
{
}
void CubeProperties::common_init()
{
    grid_ = boost::shared_ptr<CubicScalarGrid>(new CubicScalarGrid(basisset_, options_));
    grid_->set_filepath(options_.get_str("CUBEPROP_FILEPATH"));
}
void CubeProperties::print_header()
{
    outfile->Printf( "  ==> One Electron Grid Properties (v2.0) <==\n\n");
    grid_->print_header();
    outfile->Flush();
}
std::pair<SharedMatrix,SharedVector> CubeProperties::compute_NOpair(
    boost::shared_ptr<Matrix>so2ao,boost::shared_ptr<Wavefunction> wfn){
  // total density AO basis
  //SharedMatrix Da = wfn->Da_subset("MO");
  SharedMatrix Da = wfn->Da();
  //SharedMatrix Db = wfn->Db_subset("MO");
  SharedMatrix Db = wfn->Db();
  //SharedMatrix X(new Matrix(Dt->rowspi(),Dt->colspi()));
  SharedMatrix X = S_->clone();
  outfile->Printf("S matrix");
  S_->print_out();

  //transform overlap to AO basis
  //X->remove_symmetry(S_,so2ao);
  //create sym orthogonalization matrix
  outfile->Printf("S^{1/2}");
  X->power(0.5);
  X->print_out();
  //Da->transform(X);
  //Da->back_transform(X);
  SharedMatrix Dt = Da->clone();
  //Db->transform(X);
  //Db->back_transform(X);
  Dt->add(Db);
  //Transform Dt(ao) to orthogonal ao basis
  Dt->transform(X);
  outfile->Printf("Dt 'orthogonalized' so basis\n");
  Dt->print_out();

  SharedVector vals(new Vector(Dt->rowspi()));
  SharedMatrix N(new Matrix(Dt->rowspi(),Dt->colspi()));
  SharedMatrix NO(new Matrix(Dt->rowspi(),Dt->colspi()));
  Dt->diagonalize(N,vals,descending);
  outfile->Printf("Dt eigenvectors\n");
  X = S_->clone();
  X->power(-0.5);
  NO->gemm(false,false,1.0,X,N,0.0);
  NO->print_out();
  outfile->Printf("NOON \n");
  vals->print_out();
  std::pair<SharedMatrix,SharedVector> ret(NO,vals);
  return ret;

}
void CubeProperties::compute_properties()
{
    print_header();

    std::string filepath = options_.get_str("CUBEPROP_FILEPATH");
    std::stringstream ss;
    ss << filepath << "/" << "geom.xyz";

    // Is filepath a valid directory?
    boost::filesystem::path data_dir(filepath);
    if(not boost::filesystem::is_directory(data_dir)){
        printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath.c_str());
        outfile->Printf("Filepath \"%s\" is not valid.  Please create this directory.\n",filepath.c_str());
        outfile->Flush();
        exit(Failure);
    }

    basisset_->molecule()->save_xyz_file(ss.str());

    for (size_t ind = 0; ind < options_["CUBEPROP_TASKS"].size(); ind++) {
        std::string task = options_["CUBEPROP_TASKS"][ind].to_string();

        if (task == "DENSITY") {
            boost::shared_ptr<Matrix> Dt(Da_->clone());
            boost::shared_ptr<Matrix> Ds(Da_->clone());
            Dt->copy(Da_);
            Ds->copy(Da_);
            Dt->add(Db_);
            Ds->subtract(Db_);
            compute_density(Dt, "Dt");
            compute_density(Ds, "Ds");
            compute_density(Da_, "Da");
            compute_density(Db_, "Db");
       } else if (task == "EUD") {
          boost::shared_ptr<Vector> Nu(NOONa_->clone());
          boost::shared_ptr<Vector> Nd(NOONa_->clone());
          boost::shared_ptr<Vector> Ns(NOONa_->clone());
          Nu->set_name("Nu = min(1,n)");
          Ns->set_name("Ns = n^2(2-n)^2");
          Nd->set_name("Nd = 2*n - n^2");
          Nu->zero();
          Nd->zero();
          Ns->zero();
          double sumNs= 0.00;
          double sumNu= 0.00;
          double sumNd= 0.00;
          for(int h=0; h < NOONa_->nirrep(); ++h){
            for(int i =0; i <NOONa_->dim(h); ++i){
              double ni = NOONa_->get(h,i);
              double ns = (ni*ni)*((2-ni)*(2-ni));
              sumNs+= ns;
              double nu = 1-std::abs(1-ni);
              sumNu +=nu;
              double nd = 2*ni-(ni*ni);
              sumNd += nd;
              Nu->set(h,i,nu);
              Ns->set(h,i,ns);
              Nd->set(h,i,nd);
            }
          }
          outfile->Printf("         number(s) of effectively unpaired electrons\n");
          outfile->Printf("=========================================================\n");
          Nu->print_out();
          Ns->print_out();
          Nd->print_out();
          SharedMatrix U = Da_->clone();
          SharedMatrix D = Da_->clone();
          SharedMatrix S = Da_->clone();
          U->zero();
          U->set_diagonal(Nu);
          U->back_transform(Na_);
          D->zero();
          D->set_diagonal(Nd);
          D->back_transform(Na_);
          S->zero();
          S->set_diagonal(Ns);
          S->back_transform(Na_);
          compute_density(U,"updens_U");
          compute_density(D,"updens_D");
          compute_density(S,"updens_S");
        } else if (task == "ESP") {
            boost::shared_ptr<Matrix> Dt(Da_->clone());
            Dt->copy(Da_);
            Dt->add(Db_);
            compute_esp(Dt);
        } else if (task == "ORBITALS") {
            std::vector<int> indsa0;
            std::vector<int> indsb0;

            if (options_["CUBEPROP_ORBITALS"].size() == 0) {
                for (int ind = 0; ind < Ca_->colspi()[0]; ind++) {
                    indsa0.push_back(ind);
                }
                for (int ind = 0; ind < Cb_->colspi()[0]; ind++) {
                    indsb0.push_back(ind);
                }
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_ORBITALS"].size(); ind++) {
                    int val = options_["CUBEPROP_ORBITALS"][ind].to_integer();
                    if (val > 0) {
                        indsa0.push_back(abs(val) - 1);
                    } else {
                        indsb0.push_back(abs(val) - 1);
                    }
                }
            }
            std::vector<string> labelsa;
            std::vector<string> labelsb;
            CharacterTable ct = basisset_->molecule()->point_group()->char_table();
            for (size_t ind = 0; ind < indsa0.size(); ++ind){
                int i = get<1>(info_a_[indsa0[ind]]);
                int h = get<2>(info_a_[indsa0[ind]]);
                labelsa.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            for (size_t ind = 0; ind < indsb0.size(); ++ind){
                int i = get<1>(info_b_[indsb0[ind]]);
                int h = get<2>(info_b_[indsb0[ind]]);
                labelsb.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            }
            if (indsa0.size()) compute_orbitals(Ca_, indsa0,labelsa, "Psi_a");
            if (indsb0.size()) compute_orbitals(Cb_, indsb0,labelsb, "Psi_b");
        //} else if (task == "NATURAL_ORBITALS"){
            //std::vector<int> indsa0;
            //std::vector<int> indsb0;

            //if (options_["CUBEPROP_NATURAL_ORBITALS"].size() == 0) {
            //    for (int ind = 0; ind < Na_->colspi()[0]; ind++) {
            //        indsa0.push_back(ind);
            //    }
            //    for (int ind = 0; ind < Nb_->colspi()[0]; ind++) {
            //        indsb0.push_back(ind);
            //    }
            //} else {
            //    for (size_t ind = 0; ind < options_["CUBEPROP_NATURAL_ORBITALS"].size(); ind++) {
            //        int val = options_["CUBEPROP_NATURAL_ORBITALS"][ind].to_integer();
            //        if (val > 0) {
            //            indsa0.push_back(abs(val) - 1);
            //        } else {
            //            indsb0.push_back(abs(val) - 1);
            //        }
            //    }
            //}
            //std::vector<string> labelsa;
            //std::vector<string> labelsb;
            //CharacterTable ct = basisset_->molecule()->point_group()->char_table();
            //for (size_t ind = 0; ind < indsa0.size(); ++ind){
            //    int i = get<1>(info_a_[indsa0[ind]]);
            //    int h = get<2>(info_a_[indsa0[ind]]);
            //    labelsa.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            //}
            //for (size_t ind = 0; ind < indsb0.size(); ++ind){
            //    int i = get<1>(info_b_[indsb0[ind]]);
            //    int h = get<2>(info_b_[indsb0[ind]]);
            //    labelsb.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
            //}
            //if (indsa0.size()) compute_orbitals(Na_, indsa0,labelsa, "NO_a");
            //if (indsb0.size()) compute_orbitals(Nb_, indsb0,labelsb, "NO_b");

        } else if (task == "BASIS_FUNCTIONS") {
            std::vector<int> inds0;
            if (options_["CUBEPROP_BASIS_FUNCTIONS"].size() == 0) {
                for (int ind = 0; ind < basisset_->nbf(); ind++) {
                    inds0.push_back(ind);
                }
            } else {
                for (size_t ind = 0; ind < options_["CUBEPROP_BASIS_FUNCTIONS"].size(); ind++) {
                    inds0.push_back(options_["CUBEPROP_BASIS_FUNCTIONS"][ind].to_integer() - 1);
                }
            }
            compute_basis_functions(inds0, "Phi");
        } else if (task == "LOL") {
            compute_LOL(Da_, "LOLa");
            compute_LOL(Db_, "LOLb");
        } else if (task == "ELF") {
            compute_ELF(Da_, "ELFa");
            compute_ELF(Db_, "ELFb");
        } else {
            throw PSIEXCEPTION(task + "is an unrecognized PROPERTY_TASKS value");
        }
    }
}
void CubeProperties::compute_density(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_density(D, key);
}
void CubeProperties::compute_esp(boost::shared_ptr<Matrix> Dt, const std::vector<double>& w)
{
    grid_->compute_density(Dt, "Dt");
    grid_->compute_esp(Dt, w, "ESP");
}
void CubeProperties::compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::vector<std::string>& labels, const std::string& key)
{
    grid_->compute_orbitals(C, indices, labels, key);
}
void CubeProperties::compute_basis_functions(const std::vector<int>& indices, const std::string& key)
{
    grid_->compute_basis_functions(indices, key);
}
void CubeProperties::compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_LOL(D, key);
}
void CubeProperties::compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key)
{
    grid_->compute_ELF(D, key);
}

void CubeProperties::compute_EUD(boost::shared_ptr<Matrix> Na,boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> n)
{
  //boost::shared_ptr<IntegralFactory> integral(basisset_);
  SharedMatrix Na_ao(new Matrix(Ca_->colspi(),Ca_->rowspi()));
  SharedMatrix Nb_ao(new Matrix(Cb_->colspi(),Cb_->rowspi()));
  //boost::shared_ptr<PetiteList> pt(basisset_,integral);
  //SharedMatrix so2ao = pt->sotoao();

  Na_ao->gemm(false,false,1.0,Ca_,Na,0.0);
  Nb_ao->gemm(false,false,1.0,Cb_,Nb,0.0);
  Nb->transform(Cb_);
  Na->transform(Ca_);
  SharedMatrix Nt(Na_ao->clone());
  Nt->copy(Na_ao);
  Nt->add(Nb_ao);
  SharedMatrix EUD(Da_->clone());

  SharedVector nu(n->clone());
  nu->zero();
  if( options_["CUBEPROP_EUD_NUMBER_FUNCTION"].size() == 0){
    for(int h =0; h < n->nirrep(); ++h){
      for(int i =0; i< n->dim(h); ++i){
        double nu_i;
        double ni = n->get(h,i);
        nu->set(h,i,(ni*ni)*((2-ni)*(2-ni)));
      }
    }
    outfile->Printf("         number of effectively unpaired electrons\n");
    outfile->Printf("=========================================================\n");
    nu->print_out();
    EUD->zero();
    EUD->set_diagonal(nu);
    EUD->transform(Nt);
    outfile->Printf("       Density of effectively unpaired electrons\n");
    outfile->Printf("=========================================================\n");
    EUD->print_out();
    grid_->compute_density(EUD,"EUD");
  }else{
    return;
    /*
    for (size_t f_ind = 0;
      f_ind < options_["CUBEPROP_EUD_NUMBER_FUNCTION"].size();
      ++f_ind
      ){
        for(int h =0; h < n->nirrep(); ++h){
          for(int i =0; i< n->dim(h); ++i){
          }
        }
    }
    */
  }
}
void CubeProperties::compute_natural_orbitals(boost::shared_ptr<Matrix>NO){
  std::vector<int> indsa0;

  if (options_["CUBEPROP_ORBITALS"].size() == 0) {
      for (int ind = 0; ind < NO->colspi()[0]; ind++) {
          indsa0.push_back(ind);
      }
  } else {
      for (size_t ind = 0; ind < options_["CUBEPROP_ORBITALS"].size(); ind++) {
          int val = options_["CUBEPROP_ORBITALS"][ind].to_integer();
          if (val > 0) {
              indsa0.push_back(abs(val) - 1);
          }
      }
  }
  std::vector<string> labelsa;
  CharacterTable ct = basisset_->molecule()->point_group()->char_table();
  for (size_t ind = 0; ind < indsa0.size(); ++ind){
      int i = get<1>(info_a_[indsa0[ind]]);
      int h = get<2>(info_a_[indsa0[ind]]);
      labelsa.push_back(to_string(i + 1) + "-" + ct.gamma(h).symbol());
  }

  if (indsa0.size()) compute_orbitals(NO, indsa0,labelsa, "NO_");

}
void CubeProperties::compute_EUD_direct(boost::shared_ptr<Matrix> Na,boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> n)
{
  outfile->Printf(" DIRECT Computation of EUD is not availing \n");
  return;
}
}//namespace psi

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

#ifndef _psi_src_lib_libcubeprop_dcubeprop_h_
#define _psi_src_lib_libcubeprop_dcubeprop_h_

#include <map>

#include <libmints/typedefs.h>
#include <libmints/oeprop.h>
#include <libmints/wavefunction.h>
#include <bin/ccenergy/ccwave.h>
namespace boost {
  template<class T> class shared_ptr;
}
namespace psi {

class CubicScalarGrid;
class Matrix;
class Wavefunction;
class IntegralFactory;
class BasisSet;

class  DensityCubeProperties : public Prop{

protected:
        // => Required O/L from Base Prop //
  // show title and input options
  void print_header();
  // initialize important variables
  void common_init();

  // Grid-based property computer
  boost::shared_ptr<CubicScalarGrid> grid_;


  // ==>helpers that have no purpose except to prevent repeating code <==//

  // Compute Dt_so and transform to ao basis for plotting
  SharedMatrix Dt_ao();

  // Transform some kind of density from mo to so basis
  SharedMatrix Unpaired_Da_mo2so(SharedMatrix Dmo);
  SharedMatrix Unpaired_Db_mo2so(SharedMatrix Dmo);
  // Transform some kind of density from so to ao basis
  SharedMatrix Unpaired_Da_so2ao(SharedMatrix Dso);
  SharedMatrix Unpaired_Db_so2ao(SharedMatrix Dso);
  // Transform some kind of density from mo to ao basis
  // This just chains the two above together
  SharedMatrix Unpaired_Da_mo2ao(SharedMatrix Dmo);
  SharedMatrix Unpaired_Db_mo2ao(SharedMatrix Dmo);

  std::pair<SharedMatrix,SharedVector> EUD_Sa_mo();
  std::pair<SharedMatrix,SharedVector> EUD_Sb_mo();
  std::pair<SharedMatrix,SharedVector> EUD_Sa_so();
  std::pair<SharedMatrix,SharedVector> EUD_Sb_so();
  std::pair<SharedMatrix,SharedVector> EUD_Sa_ao();
  std::pair<SharedMatrix,SharedVector> EUD_Sb_ao();
  std::pair<SharedMatrix,SharedVector> compute_EUD_S();

  std::pair<SharedMatrix,SharedVector> EUD_Ua_mo();
  std::pair<SharedMatrix,SharedVector> EUD_Ub_mo();
  std::pair<SharedMatrix,SharedVector> EUD_Ua_so();
  std::pair<SharedMatrix,SharedVector> EUD_Ub_so();
  std::pair<SharedMatrix,SharedVector> EUD_Ua_ao();
  std::pair<SharedMatrix,SharedVector> EUD_Ub_ao();
  std::pair<SharedMatrix,SharedVector> compute_EUD_U();

  void print_Num_UP_info(
      std::vector<boost::tuple<int,int,double,int,int,double,double>> Upmetric,
      std::string fdef);

public:
    // => Constructors <= //

    /// Construct a DensityCubeProperties object from a Wavefunction (possibly with symmetry in wfn)
    DensityCubeProperties(SharedWavefunction wfn);

    /// Common Destructor
    virtual ~DensityCubeProperties();

    // => High-Level Property Computers <= //

    /// Compute all relevant properties from options object specifications
    void compute();

    // => Low-Level Property Computers (Do not use unless you are an advanced client code) <= //

    // Compute a density grid task (key.cube)
    //void compute_densities(const std::string key="D");
    // Compute an ESP grid task (Dt.cube and ESP.cube)
    //void compute_natural_orbitals(const std::vector<int>& indices, const std::vector<std::string>& labels, const std::string& key);
    // Compute effectively unpaired electron density (EUD) grid task
    //  type.cube, also prints number of effectively unpaired electrons to
    //  output, do compute atomic contributions?
    void compute_EUD(std::string type="U", bool atomic_contrib=false);

};

}

#endif

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

#ifndef _psi_src_lib_libcubeprop_cubeprop_h_
#define _psi_src_lib_libcubeprop_cubeprop_h_

#include <map>

#include <libmints/typedefs.h>
#include <libmints/wavefunction.h>

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
  // Global options object
  Options& options_;

  // Grid-based property computer
  boost::shared_ptr<CubicScalarGrid> grid_;

        // => Natural orbital info <= //
  // (Note: SO basis stored to keep the same conventions as base class.) //
  // alpha Natural Orbitals
  SharedMatrix NOa_so_;
  // alpha occupation numbers
  SharedVector Oca_;
  // beta Natural Orbitals
  SharedMatrix NOb_so_;
  // beta occupation numbers
  SharedVector Ocb_;

  // print requested tasks, and other input info //
  void print_header();


public:
    // => Constructors <= //

    /// Construct a CubeProperties object from a Wavefunction (possibly with symmetry in wfn)
    CubeProperties(SharedWavefunction wfn);
    CubeProperties(SharedWavefunction wfn, SharedMatrix OPDM);
    CubeProperties(SharedWavefunction wfn, SharedMatrix OPDM_a, SharedMatrix OPDM_b);

    /// Common Destructor
    virtual ~CubeProperties();

    // => High-Level Property Computers <= //

    /// Compute all relevant properties from options object specifications
    void compute_properties();

    // => Low-Level Property Computers (Do not use unless you are an advanced client code) <= //

    /// Obligatory title info
    void print_header();
    /// Compute a density grid task (key.cube)
    void compute_density(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ESP grid task (Dt.cube and ESP.cube)
    void compute_esp(boost::shared_ptr<Matrix> Dt, const std::vector<double>& nuc_weights = std::vector<double>());
    /// Compute an orbital task (key_N.cube, for 0-based indices of C)
    void compute_orbitals(boost::shared_ptr<Matrix> C, const std::vector<int>& indices, const std::vector<std::string>& labels, const std::string& key);
    /// Compute a basis function task (key_N.cube, for 0-based indices of basisset_)
    void compute_basis_functions(const std::vector<int>& indices, const std::string& key);
    /// Compute a LOL grid task (key.cube)
    void compute_LOL(boost::shared_ptr<Matrix> D, const std::string& key);
    /// Compute an ELF grid task (key.cube)
    void compute_ELF(boost::shared_ptr<Matrix> D, const std::string& key);

    /// Compute effectively unpaired electron density (EUD) grid task (key.cube)
    void compute_EUD(boost::shared_ptr<Matrix> Na, boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> n);
    void compute_EUD_direct(boost::shared_ptr<Matrix> Na, boost::shared_ptr<Matrix> Nb, boost::shared_ptr<Vector> n);
    void compute_natural_orbitals(boost::shared_ptr<Matrix> NO);

};

}

#endif

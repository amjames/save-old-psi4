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
  bool do_pop_analysis_;

  // Grid-based property computer
  boost::shared_ptr<CubicScalarGrid> grid_;


  // Transform some kind of density from so to ao basis
  SharedMatrix Unpaired_D_so2ao(SharedMatrix Dso);

  SharedMatrix compute_EUD_S(bool mulliken);
  SharedMatrix Du_s_mo();
  SharedMatrix Du_s_so();
  SharedMatrix Du_s_ao();
  void mulliken_EUD(SharedMatrix Du_ao);


  std::string cubepath_;

  void print_EUD_summary(
      std::vector<std::tuple<double,double,double,int,int>>info
      );

public:
    // => Constructor <= //
    DensityCubeProperties(SharedWavefunction wfn);
    /// Common Destructor
    virtual ~DensityCubeProperties();

    /// Compute all relevant properties from options object specifications
    void compute();

    // => Low-Level Property Computers <= //
    // =>         (advanced)           <= //
    /// compute effectively unpaired electron density

    /*   available function "types":
     *    Du_s -> type="S"-> Head-Gordon,CPL,2003 eq(18)
     *
     *    option:
     *      atomic_contrib:
     *        do compute atomic pop analysis of UP density?
     *        This will give the number of UP electrons associated with each
     *        atom, useful for deducing regions of possible chemical activity.
     *        Gives an idea of how well localized a radical is for example
     */
    void compute_EUD(std::string type="ALL", bool atomic_contrib=true);

};

}

#endif

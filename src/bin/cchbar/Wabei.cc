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

/*! \file
    \ingroup CCHBAR
    \brief Enter brief description of file here
*/
#include <cstdio>
#include <string>
#include <libdpd/dpd.h>
#include <libqt/qt.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {

/* Wabei(): Computes all contributions to the Wabei HBAR matrix
** elements, whose spin-orbital definition is:
**
** Wabei = <ab||ei> - Fme t(mi,ab) + t(i,f) <ab||ef>
**            (I)         (II)             (IIIa)
**   - P(ab) t(i,f) t(m,b) <am||ef> + 1/2 t(i,f) t(mn,ab) <mn||ef>
**               (IIIb)                            (IIIc)
**   + 1/2 P(ab) t(i,f) t(m,a) t(n,b) <mn||ef>
**                   (IIId)
**   + 1/2 t(mn,ab) <mn||ei> + 1/2 P(ab) t(m,a) t(n,b) <mn||ei>
**             (IVa)                         (IVb)
**   + P(ab) t(mi,fb) <am||ef> - P(ab) t(m,a) <mb||ei>
**               (V)                    (VI)
**   - P(ab) t(m,a) t(ni,fb) <mn||ef>
**                 (VII)
**
**  [cf. Gauss and Stanton, JCP 103, 3561-3577 (1995)]
*/

void Wabei_RHF(void);
void Wabei_ROHF(void);
void WABEI_UHF(void);
void NEW_WABEI_UHF(void);
void Wabei_UHF(void);
void NEW_Wabei_UHF(void);
void WAbEi_UHF(void);
void NEW_WAbEi_UHF(void);
void WaBeI_UHF(void);
void NEW_WaBeI_UHF(void);

void Wabei_build(void)
{
  outfile->Printf("in Wabei_build() ");
  if(params.ref == 0) Wabei_RHF();
  else if(params.ref == 1) Wabei_ROHF();
  else if(params.ref == 2) {
    if (params.new_Wabei_AAAA){
      outfile->Printf("\n\tUsing new Wabei_AAAA_UHF");
      NEW_WABEI_UHF();
    }else{
      outfile->Printf("\n\tusing old Wabei_AAAA_UHF");
      WABEI_UHF();
    }
    if (params.new_Wabei_ABAB){
      outfile->Printf("\n\t Using new Wabei_ABAB_UHF");
      NEW_WAbEi_UHF();
    }else{
      outfile->Printf("\n\t Using old Wabei_ABAB_UHF");
      WAbEi_UHF();
    }
    if (params.new_Wabei_BABA){
      outfile->Printf("\n\t Using new Wabei_BABA_UHF");
      NEW_WaBeI_UHF();
    }else{
      outfile->Printf("\n\t Using old Wabei_BABA_UHF");
      WaBeI_UHF();
    }
    if (params.new_Wabei_BBBB){
      outfile->Printf("\n\t Using new Wabei_BBBB_UHF");
      NEW_Wabei_UHF();
    }else{
      outfile->Printf("\n\t Using old Wabei_BBBB_UHF");
      Wabei_UHF();
    }
  }
}

}} // namespace psi::cchbar

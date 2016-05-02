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
#include <libdpd/dpd.h>
#include <libciomr/libciomr.h>
#include <libqt/qt.h>
#include <libpsio/psio.h>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {
void debug_check(void);
void debug_check(void)
{
  dpdbuf4 W;
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_copy(&W, PSIF_CC_TMP0, "WEiAb-check");
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_TMP0, rspq, 25, 29, "WeIaB-check", 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 25, 29, 25, 29, 0, "WeIaB-check");
  global_dpd_->buf4_print(&W,"outfile",1 );
  global_dpd_->buf4_close(&W);
}

/* WaBeI_UHF(): Computes all contributions to the aBeI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (eI,aB) ordering and is referred to on disk as "WaBeI".
**
** The spin-orbital expression for the Wabei elements is:
**
** Wabei = <ab||ei> - Fme t_mi^ab + t_i^f <ab||ef>
**         - P(ab) t_m^b <am||ef>t_i^f + 1/2 tau_mn^ab <mn||ef> t_i^f
**         + 1/2 <mn||ei> tau_mn^ab - P(ab) <mb||ef> t_mi^af
**         - P(ab) t_m^a { <mb||ei> - t_ni^bf <mn||ef> }
**
** (cf. Gauss and Stanton, JCP 103, 3561-3577 (1995).)
**
** For the aBeI spin case, we evaluate these contractions with two
** target orderings, (aB,eI) and (eI,aB), depending on the term.
** After all terms have been evaluated, the (aB,eI) terms are sorted
** into (eI,aB) ordering and both groups arer added together.
**
** TDC, June 2002
*/
void build_Z1A_BABA();
void build_Z1B_BABA();
void Bints_read_Fe_block(dpdbuf4* Buf, int eFirrep, int Firrep, int eIndex);

void WaBeI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;
  psio_tocprint(PSIF_CC_HBAR);

  /**** Term I ****/

  /** W(eI,aB) <--- <eI|aB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WeIaB");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(eI,aB) <--- - F_me t_mI^aB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  //global_dpd_->buf4_print(&W,"outfile",1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term III ****/  /** This will require special out-of-core code **/

  /** Z(Ie,Ba) <--- t_I^F <Fe|Ba> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 28, 24, 28, 0, "Z(Ie,Ba)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &B, &Z, 1, 0, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  /** Z(Ie,Ba) --> W'(aB,eI) **/
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, srqp, 29, 25, "W'(aB,eI)");
  global_dpd_->buf4_close(&Z);
 // debug_check();


  /**** Term IV ****/

  /** WaBeI <-- - t_M^B <Ma|Fe> t_I^F - t_m^a <mB|eF> t_I^F
      Evaluate in three steps:
          (1) Z_MaeI = - <aM|eF> t_I^F [stored (aM,eI)]
	  (2) Z_mBeI = <mB|eF> t_I^F   [stored (mB,eI)]
          (3) WaBeI <-- t_M^B Z_MaeI - t_m^a Z_mBeI
       Store targets in:  W(eI,aB)  and   W'(aB,eI)
  **/

  /** Z(aM,eI) <-- - <aM|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_M^B Z(aM,eI) --> W(eI,aB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(mB,eI) <-- <mB|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0, "F <iA|bC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** - t_m^a Z(mB,eI) --> W'(aB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Term V ****/

  /** WaBeI <-- tau_mN^aB <mN|eF> t_I^F
      Evaluate in two steps:
         (1) Z_mNeI = <mN|eF> t_I^F
         (2) WaBeI <-- tau_mN^aB Z_mNeI
      Store target in W'(aB,eI)
  **/

  /** Z(mN,eI) <-- <mN|eF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 23, 29, 23, 29, 0, "D <iJ|aB>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);

  /** tau_mN^aB Z(mN,eI) --> W'(aB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 23, 25, 23, 25, 0, "Z(mN,eI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term VI ****/

  /** tau_mN^aB <mN|eI> --> W'(aB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tauiJaB");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 23, 25, 23, 25, 0, "E <iJ|aK>");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);

  /**** Term VII ****/

  /** WaBeI <-- <Bm|Fe> t_Im^Fa - <am||ef> t_Im^Bf - <aM|eF> t_IM^BF
      Evaluate in five steps:
        (1) Sort <Bm|Fe> to F(Be,mF) ordering.  (** Note that we assume a sort has already
            been done for <am||ef> and <aM|eF> in the WABEI and Wabei terms. **)
        (2) Z'(Be,Ia) = F(Be,mF) T(Ia,mF)
        (3) Sort Z'(Be,Ia) --> W(eI,aB)
        (4) Z''(ae,IB) = - F(ae,mf) T(IB,mf) - F(ae,MF) T(IB,MF)
	(5) Sort Z''(ae,IB) --> W(eI,aB)

      NB: The storage for the sorts is expensive and will eventually require out-of-core
          codes.
  **/

  /** <Bm|Fe> --> F(Be,mF) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, psqr, 28, 27, "F <Ai|Bc> (Ac,iB)");
  global_dpd_->buf4_close(&F);

  /** <Bm|Fe> t_Im^Fa --> Z(Be,Ia) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 28, 27, 28, 27, 0, "F <Ai|Bc> (Ac,iB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(Be,Ia) --> W(eI,aB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "Z(Be,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrsp, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&Z);

  /** Z''(ae,IB) <-- - <am||ef> t_Im^Bf **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(ae,IB) <-- -<aM|eF> t_IM^BF **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 20, 15, 20, 0, "F <aI|bC> (ab,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z''(ai,IB) --> W(eI,aB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0, "Z(ae,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_HBAR, qrps, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&Z);


  /**** Terms VIII and IX ****/

  /** WaBeI <-- - t_m^a { <mB|eI> + t_In^Bf <mn||ef> + t_IN^BF <mN|eF> }
                + t_M^B {-<Ma|Ie> + t_In^Fa <Mn|Fe> }
      Evaluate in three steps:
         (1) Z_mBeI =  <mB|eI> + t_In^Bf <mn||ef> + tIN^BF <mN|eF>  [stored (mB,eI)]
         (2) Z_MaeI = -<Ma|Ie> + t_In^Fa <Mn|Fe>                    [stored (aM,eI)]
         (3) WaBeI <-- - t_m^a Z_mBeI + t_M^B Z_MaeI
      Store targets in     W'(aB,eI) and  W(eI,aB)
  **/

  /** Z(mB,eI) <-- <mB|eI> **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "Z(mB,eI)");
  global_dpd_->buf4_close(&D);

  /** <mn||ef> t_In^Bf --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_IN^BF --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(me,IB) --> Z(mB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 27, 25, "Z(mB,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W'(aB,eI) <-- - t_m^a Z(mB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(aM,eI) <-- - <Ma|Ie> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_TMP0, qpsr, 25, 25, "Z(aM,eI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) <-- t_In^Fa <Mn|Fe> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) --> Z(aM,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 25, 25, "Z(aM,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(eI,aB) <-- t_M^B Z_MaeI **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(aB,eI) and W(eI,aB) terms into WeIaB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&W);

  //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  //global_dpd_->buf4_print(&W, "outfile",1);
  //global_dpd_->buf4_close(&W);
}
void NEW_WaBeI_UHF(void)
{
  timer_on("UHF_WaBeI(NEW)");
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E,C, F1, F2, W1, W2, Tau;
  double value, alpha, beta;
  double core_total;
  int maxrows;
  int coltot, rowtot, incore;
  int Gef, Gei, Gab, Ge, Gi, Gf,Gfe, Gmi, Gm, nrows, ncols, nlinks, EE, e, row, Gnm, eidx;
  int Gma, ma, m, a, Ga, Gb, I, i, mi, BA,BM, ei, ab, ba, b, BB, fb, bf, fe, ef, mb, am;
  int Gam, Gmb, h;
  //psio_tocprint(PSIF_CC_HBAR);

  /***** Term I *****/
  if(params.print == 2) outfile->Printf("F<eI|aB> -> WaBeI ... ");
  /** W(eI,aB) <--- <eI|aB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 25, 29, 25, 29, 0, "F <aI|bC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WeIaB");
  global_dpd_->buf4_close(&F);
  if(params.print == 2) outfile->Printf("done\n");

  /**** Term II ****/
  //This is wrong
  /** W(eI,aB) <--- - F_ME t_Mi^Ab **/
  if(params.print == 2) outfile->Printf("\t-F_me t_mI^aB -> WaBeI ... ");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 23, 29, 23, 29, 0, "tiJaB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 2, 3, "Fme");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1, 1);
  global_dpd_->file2_close(&Fme);
  //global_dpd_->buf4_print(&W,"outfile",1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&T2);
  if(params.print == 2) outfile->Printf("done\n");
  /**** Term IIIa ****/

  /**t_I^F <Fe|aB>  **/
  if(params.print == 2) outfile->Printf("\tB*T1 -> WaBeI ... ");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 28, 28, 28, 28, 0, "B <Ab|Cd>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gfe = Gab = Gei; /* W and B are totally symmetric */
    for(Ge=0; Ge<moinfo.nirreps; Ge++){
      Gf= Ge ^ Gfe;
      Gi= Gf;
      B.matrix[Gfe] = global_dpd_->dpd_block_matrix(moinfo.avirtpi[Gf],B.params->coltot[Gfe]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      nrows = moinfo.aoccpi[Gi];
      ncols = W.params->coltot[Gei];
      nlinks = moinfo.avirtpi[Gf];
      if(nrows && ncols && nlinks){
        for(int nume=0; nume< moinfo.bvirtpi[Ge]; nume++){
          eidx = moinfo.bvir_off[Ge] + nume;
          Bints_read_Fe_block(&B, Gfe, Gf, eidx);
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][eidx],moinfo.aoccpi[Gi]);
          C_DGEMM('n','n',nrows,ncols,nlinks,1.0,T1.matrix[Gi][0],nlinks,B.matrix[Gfe][0],ncols,
              1.0,W.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&W, Gei,W.row_offset[Gei][eidx],moinfo.aoccpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gfe],moinfo.avirtpi[Gf],B.params->coltot[Gfe]);
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.aoccpi[Gi],W.params->coltot[Gei]);
    }
  }
  //global_dpd_->buf4_print(&W,"outfile",1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  if(params.print & 2) outfile->Printf( "done.\n");

  if(params.print & 2) outfile->Printf("\t\tD*T1*Tau+ E*Tau -> WaBeI ...");

  /**** Terms IIIc + IIId + IVa+IVb ****/
  /*
   * 4 terms can be expressed as - (Tau_Mn^Ab W_mNeI)
   * Notes:
   *      1. W_MnIe intermediate is read from disk (Mn,eI)order to temp buffer Z
   *      2. W_MnIe is sorted to (eI,mN) order, Saved to disk, Re-Read into buffer Z
   *            in (eI, mN) order
   *      3. Tau_ijab (mN,aB) is read from disk.
   *      5. Read W_aBeI (eI, aB) into buffer W.
   *      4. Loop over eI(row index) of W_eIaB target:
   * --AMJ 1/16
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 22, 25, 22,25, 0, "WMnIe (Mn,eI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rsqp, 25, 23, "WMnIe (eI,nM)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR,  0, 25, 23, 25, 23, 0, "WMnIe (eI,nM)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  23, 29,  23, 29, 0, "tauiJaB");
  for(Gei=0; Gei< moinfo.nirreps; Gei++) {
    Gab = Gnm = Gei; /* Everything is totally symmetric */
    nrows = T.params->rowtot[Gnm];
    ncols = T.params->coltot[Gab];
    if (nrows && ncols) {
      global_dpd_->buf4_mat_irrep_init(&Z, Gei);
      global_dpd_->buf4_mat_irrep_rd(&Z, Gei);
      global_dpd_->buf4_mat_irrep_init(&T, Gnm);
      global_dpd_->buf4_mat_irrep_rd(&T, Gnm);
      global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
      for(ei=0; ei< W.params->rowtot[Gei]; ei++){
         global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
         C_DGEMV('t',nrows,ncols,1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
             1.0,W.matrix[Gei][0],1);
         global_dpd_->buf4_mat_irrep_row_wrt(&W, Gei, ei);
      }
      global_dpd_->buf4_mat_irrep_row_close(&W, Gei);
      global_dpd_->buf4_mat_irrep_close(&T, Gnm);
      global_dpd_->buf4_mat_irrep_close(&Z, Gei);
    }
  }
  global_dpd_->buf4_close(&T);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);
  if(params.print & 2) outfile->Printf("done\n");

  /**** Term IIIb + V ****/
  /**
   * WAbEi <-- (T2+T1*T1) *F
   *
   * Z1a(Ia,mF) = t(I,F)t(m,a) - T(ma,IF)
   * <mB|eF>(Be,mF)Z1b(Ia,mF) = W1(Be,Ia)           [Contract 444]
   * Z1b(IB,MF) = t(MB,IF) - t(I,F)t(M,B)
   * <aM|eF>(ae,MF)Z1b(IB,MF) = W2(ae,IB)           [Contract 444]
   * <am||ef>(ae,mf)T2(IB,mf) =>W2(ae,IB)           [Contract 444]
   * W2(ae,IB) sort axpy(qrps) WAbEi (eI,aB)
   * W1(Be,Ia) sort axpy(qrsp) WAbEi (eI,aB)
   *
   */
  if(params.print & 2) outfile->Printf("\t\tWaBeI <-- (T2+T1*T1) * F ...");
  if(!params.wabei_lowdisk){
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 27, 29, 27, 29, 0 , "F <iA|bC>");
    global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, qrps, 28,27, "F <iA|bC> (Ab,iC)");
    incore =1;
    core_total=0;
    for(h=0; h<moinfo.nirreps; h++) {
       coltot= F.params->coltot[h];
       if(coltot)
        maxrows=DPD_BIGNUM/coltot;
       if(maxrows<1){
         outfile->Printf("\n Wabei_UHF(AAAA) Error: A single row of OVVV > 2 GW.\n");
         exit(PSI_RETURN_FAILURE);
       }
       else maxrows = DPD_BIGNUM;
       rowtot = F.params->rowtot[h];
       for(; rowtot > maxrows; rowtot -= maxrows) {
          if (core_total > (core_total +2*maxrows*coltot)) incore =0;
          else core_total += 2*maxrows*coltot;
       }
       if(core_total > (core_total + 2*rowtot*coltot)) incore =0;
       core_total += 2*rowtot*coltot;
    }
    if(core_total > dpd_memfree()) incore = 0;
    if(!incore && (params.print == 1)){
       outfile->Printf("\n Wabei_UHF(AAAA) Error: no out-of-core algorithim for(T2+T1*T1)*F -> Wabei.\n");
       outfile->Printf("core required: %d, DPD_MEMFREE: %d",core_total, dpd_memfree());
       exit(PSI_RETURN_FAILURE);
    }
    incore =1;
    global_dpd_->buf4_close(&F);

    /** Z1a = t_I^F t_m^a + t_mI^aF **/
    outfile->Printf("\n /** Z1a = t_I^F t_m^a + t_mI^aF **/\n " );
    build_Z1A_BABA();

    /** W1(Be,Ia)<--  <mB|eF>Z1a(Ia,mF) **/
    outfile->Printf(" /** W1(Be,Ia)<--  <mB|eF>Z1a(Ia,mF) **/\n ");
    global_dpd_->buf4_init(&F,PSIF_CC_FINTS, 0, 28, 27, 28, 27, 0, "F <iA|bC> (Ab,iC)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 27, 24 ,27, 0, "Z1a(Ia,mF)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "W1(Be,Ia)");
    if(incore)global_dpd_->contract444(&F, &Z, &W, 0, 0, -1, 0);
    //global_dpd_->buf4_print(&W,"outfile",1);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&F);

    /** Z1b = +t_MI^BF - t_I^F t_M^B **/
    outfile->Printf( " /** Z1b = +t_MI^BF - t_I^F t_M^B **/ \n");
    build_Z1B_BABA();

    /** W2(IB,ae)<-- <aM|eF>Z1b(IB,MF) **/
    outfile->Printf(" /** W2(IB,ae)<-- <aM|eF>Z1b(IB,MF) **/\n ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 20, 15, 20, 0, "F <aI|bC> (ab,IC)");
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0 , 0, 20, 20, 20, 20, 0, "Z1b(IB,MF)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0,  0, 15, 20, 15, 20, 0, "W2(ae,IB)");
    if(incore)global_dpd_->contract444(&F, &Z, &W, 0, 0, 1, 0);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    /** W2(IB,ae)<--t_Im^Bf<am||ef> **/
    outfile->Printf(" /** W2(IB,ae)<--t_Im^Bf<am||ef> **/\n ");
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 15, 30, 15, 30, 0, "F <ai||bc> (ab,ic)");
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0,  "tIAjb");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0,  0, 15, 20, 15, 20, 0, "W2(ae,IB)");
    if(incore)global_dpd_->contract444(&F, &T, &W, 0,0, 1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    /** Add W2 and W1 to target **/
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 28, 24, 28, 24, 0, "W1(Be,Ia)");
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrsp, 25, 29, "WeIaB",1 );
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 15, 20, 15, 20, 0,"W2(ae,IB)");
    /* replace with ooc or smater code above */
    global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, qrps, 25, 29, "WeIaB",1 );
    global_dpd_->buf4_close(&W);
  }else{
   // Once I get this working correctly I will worry about the low-disk case
    outfile->Printf("\nWABEI_UHF(BABA) Error: No low-disk algorithim for (T2+T1*T1)*F ->WaBeI\n");
    exit(PSI_RETURN_FAILURE);
  }
  if(params.print & 2)outfile->Printf("done!\n");

  /**** Terms VIII and IX ****/

  /** WaBeI <-- - t_m^a { <mB|eI> + t_In^Bf <mn||ef> + t_IN^BF <mN|eF> }
                + t_M^B {-<Ma|Ie> + t_In^Fa <Mn|Fe> }
      Evaluate in three steps:
         (1) Z_mBeI =  <mB|eI> + t_In^Bf <mn||ef> + tIN^BF <mN|eF>  [stored (mB,eI)]
         (2) Z_MaeI = -<Ma|Ie> + t_In^Fa <Mn|Fe>                    [stored (aM,eI)]
         (3) WaBeI <-- - t_m^a Z_mBeI + t_M^B Z_MaeI
      Store targets in     W'(aB,eI) and  W(eI,aB)
  **/

  /** Z(mB,eI) <-- <mB|eI> **/
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 27, 25, 27, 25, 0, "D <iJ|aB> (iB,aJ)");
  global_dpd_->buf4_copy(&D, PSIF_CC_TMP0, "Z(mB,eI)");
  global_dpd_->buf4_close(&D);

  /** <mn||ef> t_In^Bf --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 30, 30, 30, 0, "D <ij||ab> (ia,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <mN|eF> t_IN^BF --> Z(me,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 30, 20, 30, 20, 0, "D <Ij|Ab> (ia,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(me,IB) --> Z(mB,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 30, 20, 30, 20, 0, "Z(me,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 27, 25, "Z(mB,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W'(aB,eI) <-- - t_m^a Z(mB,eI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 27, 25, 27, 25, 0, "Z(mB,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->contract244(&T1, &Z, &W, 0, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /** Z(aM,eI) <-- - <Ma|Ie> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 24, 24, 24, 24, 0, "C <Ia|Jb>");
  global_dpd_->buf4_sort(&C, PSIF_CC_TMP0, qpsr, 25, 25, "Z(aM,eI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) <-- t_In^Fa <Mn|Fe> **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 24, 27, 24, 27, 0, "D <Ij|Ab> (Ib,jA)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(Me,Ia) --> Z(aM,eI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 24, 24, 24, 0, "Z(Me,Ia)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, spqr, 25, 25, "Z(aM,eI)", 1);
  global_dpd_->buf4_close(&Z);

  /** W(eI,aB) <-- t_M^B Z_MaeI **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 25, 25, 25, 25, 0, "Z(aM,eI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&Z, &T1, &W, 1, 0, 0, 1, 1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);

  /**** Combine accumulated W'(aB,eI) and W(eI,aB) terms into WeIaB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 29, 25, 29, 25, 0, "W'(aB,eI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 25, 29, "WeIaB", 1);
  global_dpd_->buf4_close(&W);

  //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 25, 29, 25, 29, 0, "WeIaB");
  //global_dpd_->buf4_print(&W, "outfile",1);
  //global_dpd_->buf4_close(&W);
  timer_off("UHF_WaBeI(NEW)");

}

/*
 ** Bints_read_Fe_block
 ** To avoid sorting the B integrals <Ab|Cd> to <aB|cD>
 ** I will read the rows [allF][eval] into a block matrix
 ** for every e. this may require some out of core fanciness
 ** or may be more expensive than the sort, profiling will tell
 **
 ** The matrix for the buffer should have been appropriately sized
 ** previously with block_matrix()
*/
void Bints_read_Fe_block(dpdbuf4* Buf,int GeF, int GF, int eIdx)
{
  double value;               //Matrix element
  int Fe;                     //row index (buffer)
  int Ba;                     //col index (file)
  int aB;                     //col index (file)
  int col;                    //col index (buffer matrix)
  int FIdx, aIdx, BIdx;       //orbital indices
  int nF, na, nB;             //orbital counters
  int Ga, GB;                 //Irrep counters

  /* I am using these for readability */
  /*   and to avoid looking up things in moinfo */
  /*   over and over */
  int FMax,aMax,BMax;         //orbital counter limits
  int FStart, aStart,BStart;  //orbital offset holders
  dpdparams4 spoofedParams = global_dpd_->params4[29][29];

  //for debuging checks
  int coltot = Buf->params->coltot[GeF];
  FMax = moinfo.avirtpi[GF];
   /* remove this timer after profiling */
  timer_on("Noncontiguous Bints Read (WaBeI build)");
  /* prepare the input buffer from the file */
  global_dpd_->file4_mat_irrep_row_init(&(Buf->file), GeF);

  /* loop over F values */
  FStart = moinfo.avir_off[GF];

  for( nF =0; nF < FMax; nF++){
    FIdx = FStart+nF;
    Fe = Buf->params->rowidx[FIdx][eIdx];
    global_dpd_->file4_mat_irrep_row_rd(&(Buf->file), GeF,Fe);

    /* count the cols of the target matrix */
    col =0;

    /* Loop over irreps of the beta virtual space */
    for(Ga=0; Ga<moinfo.nirreps; Ga++){
      /* Compute irrep of B */
      GB=GeF^Ga;
      /* orbital indices will start at the offset */
      /* each irrep in each space */
      aStart = moinfo.bvir_off[Ga];
      BStart = moinfo.avir_off[GB];
      aMax = moinfo.bvirtpi[Ga];
      BMax = moinfo.avirtpi[GB];

      /* loop over the a's */
      for(na=0; na< aMax; na++){
        /* compute a orbital index */
        aIdx = aStart+na;

        /* loop over b's */
        for(nB=0; nB<BMax; nB++,col++){
          /* compute B orbital index */
          BIdx = BStart +nB;
          /* Get the file index for this value */
          Ba = Buf->file.params->colidx[BIdx][aIdx];
          aB = spoofedParams.colidx[aIdx][BIdx];
          /* get the value */
          //value = Buf->file.matrix[GeF][0][Ba];
          value = Buf->file.matrix[GeF][0][aB];
          /* shove it into the Buff's matrix */
          if (col > coltot) outfile->Printf("\nError: Colls have run amok\n");
          Buf->matrix[GeF][nF][col] = value;
        }
      }
    }
  }
  global_dpd_->file4_mat_irrep_row_close(&(Buf->file), GeF);
  timer_off("Noncontiguous Bints Read (WaBeI build)");
}


void build_Z1A_BABA()
{
  dpdfile2 TIF, Tma;
  dpdbuf4 Z, T2;
  int I, F, m, a, Iorb, Forb, aorb, morb, GI, Gm, Ga, GF;
  int h, row, col;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 24, 27, 24, 27, 0, "tIbjA" );
  //global_dpd_->buf4_scmcopy(&T2, PSIF_CC_TMP0, "Z1a(iA,Mf)",-1);
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1a(Ia,mF)");
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 24, 27, 24, 27, 0, "Z1a(Ia,mF)");
  global_dpd_->file2_init(&TIF, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_init(&Tma, PSIF_CC_OEI, 0, 2, 3, "tia");
  global_dpd_->file2_mat_init(&TIF);
  global_dpd_->file2_mat_init(&Tma);
  global_dpd_->file2_mat_rd(&Tma);
  global_dpd_->file2_mat_rd(&TIF);

  for(h = 0; h<moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row = 0; row< Z.params->rowtot[h]; row++){
      Iorb = Z.params->roworb[h][row][0];
      aorb = Z.params->roworb[h][row][1];
      I = TIF.params->rowidx[Iorb];
      a = Tma.params->colidx[aorb];
      GI = TIF.params->psym[Iorb];
      Ga = Tma.params->qsym[aorb];
      for(col = 0; col < Z.params->coltot[h]; col++ ){
        morb = Z.params->colorb[h][col][0];
        Forb = Z.params->colorb[h][col][1];
        m = Tma.params->rowidx[morb];
        F = TIF.params->colidx[Forb];
        Gm = Tma.params->psym[morb];
        GF = TIF.params->qsym[Forb];

        if(GI == GF && Ga == Gm){
          Z.matrix[h][row][col] += TIF.matrix[GI][I][F] * Tma.matrix[Gm][m][a];
        }

      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&TIF);
  global_dpd_->file2_mat_close(&Tma);
  global_dpd_->file2_close(&TIF);
  global_dpd_->file2_close(&Tma);
  global_dpd_->buf4_close(&Z);
}

void build_Z1B_BABA()
{
  dpdfile2 T1;
  dpdbuf4 Z, T2;
  int I, F, B, M, Iorb, Forb, Borb, Morb, GI, GM, GB, GF;
  int h, row, col;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB" );
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1b(IB,MF)");
  global_dpd_->buf4_close(&T2);

  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1b(IB,MF)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  for(h = 0; h < moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z, h);
    global_dpd_->buf4_mat_irrep_rd(&Z, h);
    for(row = 0; row < Z.params->rowtot[h]; row++){
      Iorb = Z.params->roworb[h][row][0];
      Borb = Z.params->roworb[h][row][1];
      I = T1.params->rowidx[Iorb];
      B = T1.params->colidx[Borb];
      GI = T1.params->psym[Iorb];
      GB = T1.params->qsym[Borb];
      for(col = 0; col < Z.params->coltot[h]; col++){
        Morb = Z.params->colorb[h][col][0];
        Forb = Z.params->colorb[h][col][1];
        M = T1.params->rowidx[Morb];
        F = T1.params->colidx[Forb];
        GM = T1.params->psym[Morb];
        GF = T1.params->qsym[Forb];

        if(GI == GF && GB == GM) {
           Z.matrix[h][row][col] -= T1.matrix[GI][I][F] *T1.matrix[GM][M][B];
        }
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z, h);
    global_dpd_->buf4_mat_irrep_close(&Z, h);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
}
}} // namespace psi::cchbar

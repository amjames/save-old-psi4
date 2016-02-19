/*
 *@BEGIN LICENSE
 *
 * PSI4: an ab initio quantum chemistry software package
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
 *@END LICENSE
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
#include <cmath>
#include "MOInfo.h"
#include "Params.h"
#define EXTERN
#include "globals.h"

namespace psi { namespace cchbar {
void build_UHF_Z1(void);
void debug_check(void);
void debug_break(void);
void debug_break(void){
  return;
}
/* WABEI_UHF(): Computes all contributions to the ABEI spin case of
** the Wabei HBAR matrix elements.  The final product is stored in
** (EI,AB) ordering and is referred to on disk as "WEIAB".
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
** For the ABEI spin case, we evaluate these contractions with two
** target orderings, (AB,EI) and (EI,AB), depending on the term.
** After all terms have been evaluated, the (AB,EI) terms are sorted
** into (EI,AB) ordering and both groups arer added together.
**
** TDC, June 2002
*/

void debug_check(void)
{
  dpdbuf4 W;
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_copy(&W, PSIF_CC_TMP0, "WEIAB-check");
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_TMP0, rspq, 21, 7, "WEIAB-check", 1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 21, 7, 21, 7, 0, "WEIAB-check");
  global_dpd_->buf4_print(&W,"outfile",1 );
  global_dpd_->buf4_close(&W);
}
void WABEI_UHF(void)
{
  dpdfile2 Fme, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C;

  /**** Term I ****/

  /** W(EI,AB) <--- <EI||AB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEIAB");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  /** W(EI,AB) <--- - F_ME t_MI^AB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&Fme, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->contract244(&Fme, &T2, &W, 0, 0, 0, -1.0, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_close(&Fme);
  global_dpd_->buf4_close(&T2);

  /**** Term IIIa ****/

  /** W'(AB,EI) <--- <AB||EF> t_I^F **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 7, 5, 5, 5, 1, "B <AB|CD>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&B, &T1, &W, 3, 1, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  global_dpd_->buf4_close(&W);


  /**** Term IIIc+IIId ****/

  /** WABEI <-- 1/2 tau_MN^AB <MN||EF> t_I^F
      Evaluate in two steps:
         (1) Z_MNEI = <MN||EF> t_I^F
         (2) WABEI <-- 1/2 tau_MN^AB Z_MNEI
      Store target in W'(AB,EI)
  **/

  /** Z(MN,EI) <-- <MN||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 2, 5, 2, 5, 0, "D <IJ||AB> (I>J,AB)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&D, &T1, &Z, 3, 1, 0, 1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** tau_MN^AB Z(MN,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 2, 21, 2, 21, 0, "Z(MN,EI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->contract444(&T2, &Z, &W, 1, 1, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);

  /**** Term IVa + IVb****/

  /** tau_MN^AB <MN||EI> --> W'(AB,EI) **/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 2, 7, 2, 7, 0, "tauIJAB");
  global_dpd_->buf4_init(&E, PSIF_CC_EINTS, 0, 2, 21, 2, 21, 0, "E <IJ||KA> (I>J,AK)");
  global_dpd_->contract444(&T2, &E, &W, 1, 1, -1, 1);
  global_dpd_->buf4_close(&E);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  //debug_check();

  /**** Term IIIb ****/

  /** WABEI <-- t_M^B <MA||EF> t_I^F - t_M^A <MB||EF> t_I^F
      Evaluate in two steps:
          (1) Z_MBEI = <MB||EF> t_I^F
          (2) WABEI <-- t_M^B Z_MAEI - t_M^A Z_MBEI
  **/

  /** Z(MB,EI) <-- - <MB||EF> t_I^F **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 20, 5, 20, 5, 1, "F <IA|BC>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract424(&F, &T1, &Z, 3, 1, 0, -1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** t_M^A Z(MB,EI) --> Z1(AB,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z, &Z1, 0, 0, 0, 1.0, 0.0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qprs, 5, 21, "Z2(BA,EI)");
  global_dpd_->buf4_close(&Z1);

  /** Z1(AB,EI) - Z2(BA,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z1(AB,EI)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z2(BA,EI)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1.0);
  global_dpd_->buf4_close(&Z2);

  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_axpy(&Z1, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z1);
  //debug_check();
  /**** Term V ****/

  /** WABEI <-- <BM||EF> t_IM^AF + <Bm|Ef> t_Im^Af - <AM||EF> t_IM^BF - <Am|Ef> t_Im^Bf
      Evaluate in six steps:
        (1) Sort <BM||EF> and <Bm|Ef> to F(BE,MF) and F(BE,mf) ordering.
        (2) Z(BE,IA) = F(BE,MF) T(IA,MF) + F(BE,mf) T(IA,mf)
        (3) Sort Z(BE,IA) --> Z'(EI,AB)
	(4) Sort Z'(EI,AB) --> Z''(EI,BA)
        (5) AXPY: Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA)
        (6) AXPY: W(EI,AB) <-- Z'(EI,AB)
      NB: The storage for the sorts is expensive and will eventually require out-of-core
          codes.
  **/

  /** <BM||EF> --> F(BE,MF) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 5, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5, 20, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_close(&F);

  /** <Bm|Ef> --> (BE,mf) **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
  global_dpd_->buf4_sort(&F, PSIF_CC_FINTS, prqs, 5, 30, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_close(&F);

  /** <BM||EF> t_IM^AF --> Z(BE,IA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** <Bm|Ef> t_Im^Af --> Z(BE,IA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 5, 30, 5, 30, 0, "F <Ai|Bc> (AB,ic)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&F, &T2, &Z, 0, 0, -1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&F);
  global_dpd_->buf4_close(&Z);

  /** Z(BE,IA) --> Z'(EI,AB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "Z(BE,IA)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrsp, 21, 5, "Z'(EI,AB)");
  global_dpd_->buf4_close(&Z);

  /** Z'(EI,AB) --> Z''(EI,BA) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, pqsr, 21, 5, "Z''(EI,BA)");
  global_dpd_->buf4_close(&Z);

  /** Z'(EI,AB) = Z'(EI,AB) - Z''(EI,BA) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z''(EI,BA)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** W(EI,AB) <-- Z'(EI,AB) **/
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "Z'(EI,AB)");
  global_dpd_->buf4_axpy(&Z, &W, 1);
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_close(&W);
  //debug_check();

  /**** Terms VI+ VII ****/

  /** WABEI <-- -P(AB) t_M^A { <MB||EI> + t_IN^BF <MN||EF> + t_In^Bf <Mn|Ef> }
      Evaluate in two steps:
         (1) Z_MBEI = <MB||EI> + t_IN^BF <MN||EF> + tIn^Bf <Mn|Ef>
         (2) WABEI <-- - t_M^A Z_MBEI + t_M^B Z_MAEI
      Store target in W'(AB,EI)
  **/

  /** Z(MB,EI) <-- <MB||EI> **/
  global_dpd_->buf4_init(&C, PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
  global_dpd_->buf4_copy(&C, PSIF_CC_TMP0, "Z(MB,EI)");
  global_dpd_->buf4_close(&C);
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->buf4_scm(&Z, -1);
  global_dpd_->buf4_close(&Z);

  /** <MN||EF> t_IN^BF --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 0);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** <Mn|Ef> t_In^Bf --> Z(ME,IB) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_init(&D, PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1, 1);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&D);
  global_dpd_->buf4_close(&Z);

  /** Z(ME,IB) --> Z(MB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z(ME,IB)");
  global_dpd_->buf4_sort_axpy(&Z, PSIF_CC_TMP0, psqr, 20, 21, "Z(MB,EI)", 1);
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) <-- -t_M^A Z(MB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 21, 20, 21, 0, "Z(MB,EI)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->contract244(&T1, &Z1, &Z, 0, 0, 0, -1, 0);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) --> Z'(BA,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qprs, 5, 21, "Z'(BA,EI)");
  global_dpd_->buf4_close(&Z);

  /** Z(AB,EI) = Z(AB,EI) - Z'(BA,EI) **/
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z'(BA,EI)");
  global_dpd_->buf4_axpy(&Z2, &Z1, -1);
  global_dpd_->buf4_close(&Z2);
  global_dpd_->buf4_close(&Z1);

  /** Z(AB,EI) --> W'(AB,EI) **/
  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 5, 21, 5, 21, 0, "Z(AB,EI)");
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_axpy(&Z, &W, 1.0);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);
  //debug_check();

  /**** Combine accumulated W'(AB,EI) and W(EI,AB) terms into WEIAB ****/
  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 7, 21, 7, 21, 0, "W'(AB,EI)");
  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, rspq, 21, 7, "WEIAB", 1);
  global_dpd_->buf4_close(&W);
}

void NEW_WABEI_UHF(void)
{
  timer_on("UHF_WABEI(NEW)");
  dpdfile2 FME, T1;
  dpdbuf4 F, W, T2, B, Z, Z1, Z2, D, T, E, C, F1, F2, W1, W2, Tau;
  double value, alpha, beta;
  int Gef, Gei, Gab, Ge, Gi, Gf, Gmi, Gm, nrows, ncols, nlinks, EE, e, row, Gnm;
  int Gma, ma, m, a, Ga, Gb, I, i, mi,  ei, ab, ba, b, BB, fb, bf, fe, ef, mb, am;
  int Gam, Gmb;
  double ***WW1, ***WW2;
  int h, incore, core_total, rowtot, coltot, maxrows;

  /**** Term I ****/
  if(params.print == 2) {
    outfile->Printf( "\n\t\tF<AI|BC> -> Wabei...");

  }
  if(params.print & 2) outfile->Printf( "done.\n");
  /** W(EI,AB) <--- <EI||AB> **/
  global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 7, 21, 5, 1, "F <AI|BC>");
  global_dpd_->buf4_copy(&F, PSIF_CC_HBAR, "WEIAB");
  global_dpd_->buf4_close(&F);

  /**** Term II ****/

  if(params.print == 2) {
    outfile->Printf( "\t\t FME*T2 -> Wabei...");
  }
  /** W(EI,AB) <--- - F_ME t_MI^AB **/
  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 0, 7, 2, 7, 0, "tIJAB");
  global_dpd_->file2_init(&FME, PSIF_CC_OEI, 0, 0, 1, "FME");
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  /* global_dpd_->contract244(&FME,&T2,&W, 0, 0, 0,-1.0,1.0); */
  global_dpd_->file2_mat_init(&FME);
  global_dpd_->file2_mat_rd(&FME);
  for(Gei=0; Gei < moinfo.nirreps; Gei++) {
    Gmi = Gab = Gei;
    global_dpd_->buf4_mat_irrep_init(&T2,Gmi);
    global_dpd_->buf4_mat_irrep_rd(&T2,Gmi);
    row=0;
    for(Ge=0; Ge<moinfo.nirreps; Ge++){
      Gm= Ge;
      Gi= Gm ^ Gmi;
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      nrows = moinfo.aoccpi[Gm];
      ncols = moinfo.aoccpi[Gi] * W.params->coltot[Gei];
      if(nrows && ncols){
        for(EE=0; EE< moinfo.avirtpi[Ge]; EE++){
          e = moinfo.avir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
          C_DGEMV('t',nrows,ncols, -1.0,&T2.matrix[Gmi][row][0],ncols,
              &FME.matrix[Gm][0][EE],moinfo.avirtpi[Ge], 1.0,W.matrix[Gei][0],1);
          global_dpd_->buf4_mat_irrep_wrt_block(&W,Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
        }
      }
      row+= moinfo.aoccpi[Gm]* moinfo.aoccpi[Gi];
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.aoccpi[Gi],W.params->coltot[Gei]);
    }
    global_dpd_->buf4_mat_irrep_close(&T2,Gmi);

  }
  global_dpd_->file2_mat_close(&FME);
  global_dpd_->file2_close(&FME);
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_close(&W);
  if(params.print & 2) outfile->Printf( "done.\n");

  /**** Term IIIa ****/

  /** W(EI,AB) <--- <AB||EF> t_I^F **/
  if(params.print & 2) {
    outfile->Printf( "\t\tB*T1 -> Wabei...");

  }
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&B, PSIF_CC_BINTS, 0, 5, 7, 5, 5, 1, "B <AB|CD>");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);
  for(Gef=0; Gef < moinfo.nirreps; Gef++) {
    Gei = Gab = Gef; /* W and B are totally symmetric */
    for(Ge=0; Ge<moinfo.nirreps; Ge++){
      Gf= Ge ^ Gef;
      Gi= Gf;
      B.matrix[Gef] = global_dpd_->dpd_block_matrix(moinfo.avirtpi[Gf],B.params->coltot[Gef]);
      W.matrix[Gei] = global_dpd_->dpd_block_matrix(moinfo.aoccpi[Gi],W.params->coltot[Gei]);
      nrows = moinfo.aoccpi[Gi];
      ncols = W.params->coltot[Gei];
      nlinks = moinfo.avirtpi[Gf];
      if(nrows && ncols){
        for(EE=0; EE< moinfo.avirtpi[Ge]; EE++){
          e = moinfo.avir_off[Ge] + EE;
          global_dpd_->buf4_mat_irrep_rd_block(&B, Gef, B.row_offset[Gef][e],moinfo.avirtpi[Gf]);
          global_dpd_->buf4_mat_irrep_rd_block(&W, Gei, W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
          C_DGEMM('n','n',nrows,ncols,nlinks,1.0,T1.matrix[Gi][0],nlinks,B.matrix[Gef][0],ncols,
              1.0,W.matrix[Gei][0],ncols);
          global_dpd_->buf4_mat_irrep_wrt_block(&W, Gei,W.row_offset[Gei][e],moinfo.aoccpi[Gi]);
        }
      }
      global_dpd_->free_dpd_block(B.matrix[Gef],moinfo.avirtpi[Gf],W.params->coltot[Gef]);
      global_dpd_->free_dpd_block(W.matrix[Gei],moinfo.aoccpi[Gi],W.params->coltot[Gei]);
    }
  }
  //global_dpd_->buf4_print(&W,"outfile",1);
  global_dpd_->buf4_close(&W);
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&B);
  if(params.print & 2) outfile->Printf( "done.\n");


  if(params.print & 2) outfile->Printf("\t\tD*T1*Tau+ E*Tau ...");

  /**** Terms IIIc + IIId + IVa+IVb ****/
  /*
   * 4 terms can be expressed as - (Tau_mn^ab W_MNEI)
   * Notes:
   *      1. W_MNIE intermediate is read from disk (M>N-,EI)order to temp buffer Z
   *      2. W_MNIE is sorted to (EI,M>N-) order, Saved to disk, Re-Read into buffer Z
   *            in (EI, M>N-) order
   *      3. Tau_IJAB (MN,AB) is read from disk.
   *      5. Read W_ABEI (EI, A>B-) into buffer W.
   *      4. Loop over EI(row index) of W_EIAB target:
   * --AMJ 1/16
   */
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 2, 21, 2,21, 0, "WMNIE (M>N,EI)");
  global_dpd_->buf4_sort(&Z, PSIF_CC_HBAR, rspq, 21, 2, "WMNIE (EI,M>N)");
  global_dpd_->buf4_close(&Z);
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR,  0, 21, 7, 21, 7, 0, "WEIAB");
  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR,  0, 21, 2, 21, 2, 0, "WMNIE (EI,M>N)");
  global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0,  2, 7,  2, 7, 0, "tauIJAB");
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
         C_DGEMV('t',nrows,ncols,-1,T.matrix[Gei][0],ncols,Z.matrix[Gei][ei],1,
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
  //global_dpd_->buf4_print(&W,"outfile",1);
  global_dpd_->buf4_close(&W);
  if(params.print & 2) outfile->Printf("done\n");

  if(params.print & 2) outfile->Printf("\t\t (T2+T1*T1)*F... ");
  /**** Term IIIb + V  ****/
  /*
   * WEIAB <-- W2(EI,AB)-W2(EI,BA)
   * W2(EI,AB ) = W1(IB,EA) (sorted )
   * W1(IB,EA) = Z1(IB,MF) <MA||FE> (MF,EA) + t_Ij^Ab(IB,mf) <mA|fE> (mf,EA)
   * Z1(IB,MF) = t_IJ^AB(IB,MF) - t_I^F*t_M^B
   *
   * NB:
   * W2(EI,AB)-W2(EI,BA) achieved by reading in with A/S flag
   *
   * --AMJ 02/12/2016
   **/
  build_UHF_Z1();//Z(IA,MF)
  if(!params.wabei_lowdisk){
    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 21, 5, 21, 5, 1 , "F <AI|BC>");
    //                                                        BM||EF   BE MF
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 5, 20, "F <AI||BC> (AB,IC)");
     //can we run contractions fully in core?
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
    if(!incore && (params.print & 2)){
       outfile->Printf("\n Wabei_UHF(AAAA) Error: no out-of-core algorithim for(T2+T1*T1)*F -> Wabei.\n");
       outfile->Printf("core required: %d, DPD_MEMFREE: %d",core_total, dpd_memfree);
       exit(PSI_RETURN_FAILURE);
    }
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1(IA,MF)");
    //global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
                                                                //F <BM||EF> (BE,MF)
    global_dpd_->buf4_init(&F, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "F <AI||BC> (AB,IC)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 20, 5, 20,  0, "W1(BE,IA)");
    if(incore)global_dpd_->contract444(&F,&Z, &W, 0, 0, -1, 0);
    global_dpd_->buf4_close(&W);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&Z);

    global_dpd_->buf4_init(&F, PSIF_CC_FINTS, 0, 26, 28, 26, 28, 0, "F <Ai|Bc>");
    global_dpd_->buf4_sort(&F, PSIF_CC_TMP0, prqs, 5 ,30,"F <Ai|Bc> (AB,ic)");
    global_dpd_->buf4_close(&F);                                  //tIAmf
    global_dpd_->buf4_init(&T, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0,"tIAjb");
    global_dpd_->buf4_init(&F, PSIF_CC_TMP0,  0, 5, 30, 5, 30,  0,"F <Ai|Bc> (AB,ic)");
    global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "W1(BE,IA)");
    if(incore) global_dpd_->contract444(&F, &T, &W, 0, 0, -1, 1);
    global_dpd_->buf4_close(&F);
    global_dpd_->buf4_close(&T);
    global_dpd_->buf4_close(&W);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 5, 20, 5, 20, 0, "W1(BE,IA)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, qrsp, 21, 5, "W2(EI,AB)");
    global_dpd_->buf4_close(&Z1);

    global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "W2(EI,AB)");
    global_dpd_->buf4_sort(&Z1, PSIF_CC_TMP0, pqsr, 21, 5, "W2'(EI,BA)");
    global_dpd_->buf4_close(&Z1);

    /* global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "W2(EI,AB)"); */
    /* global_dpd_->buf4_init(&Z2, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 0, "W2'(EI,BA)"); */
    /* global_dpd_->buf4_axpy(&Z2,&Z1, -1); */
    /* global_dpd_->buf4_close(&Z1); */
    /* global_dpd_->buf4_close(&Z2); */

    global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 5, 21, 7, 0, "WEIAB" );
    global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 5, 21, 5, 1, "W2(EI,AB)");
    global_dpd_->buf4_axpy(&Z,&W,1.0);
    //global_dpd_->buf4_print(&W,"outfile",1);
    global_dpd_->buf4_close(&Z);
    global_dpd_->buf4_close(&W);

    //global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7, 0, "WEIAB" );
    //global_dpd_->buf4_print(&W,"outfile",1);
    //global_dpd_->buf4_close(&W);
    psio_close(PSIF_CC_TMP0,0);  //Z1, sorted Fints removed from disk
    psio_open(PSIF_CC_TMP0, PSIO_OPEN_NEW);
  }
  else{
   // Once I get this working correctly I will worry about the low-disk case
    outfile->Printf("\nWABEI_UHF(AAAA) Error: No low-disk algorithim for (T2+T1*T1)*F ->Wabei\n");
    exit(PSI_RETURN_FAILURE);
  }

  if(params.print & 2 ) outfile->Printf("done\n");
  /** Term VI + VII **/

  /*- Pab t_M^A {<MB||EI> +t_IN^BF<MN||EF> -t_nI^fB<Mn|Ef>}
   * 1.   t_IN^BF * <MN||EF> + t_nI^fB * <Mn|Ef> --> Z(ME,IB)
   * 2.   Z(ME,IB) --sort--> Z(EI,MB)
   * 3. - t_M^A( <MB||EI> + Z(EI,MB) ) --> W'(EI,AB)
   * 4. WABEI <-- W'(EI,AB)- W'(EI,AB)
   */
//  if(params.print & 2 ) outfile->Printf("\t\tT1*(C+D*T2)-->WEIAB");
//        /* <MN||EF>(ME,NF)*t_IN^BF(NF,IB) --> Z(ME,IB) */
//  global_dpd_->buf4_init(&D,  PSIF_CC_DINTS, 0, 20, 20, 20, 20, 0, "D <IJ||AB> (IA,JB)");
//  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
//  global_dpd_->buf4_init(&Z,  PSIF_CC_TMP0,  0, 20, 20, 20, 20, 0, "Z(ME,IB)");
//  global_dpd_->contract444(&D, &T2, &Z, 0, 0, 1.0, 0.0);
//  global_dpd_->buf4_close(&D);
//  global_dpd_->buf4_close(&T2);
//        /* - <Mn|Ef> (ME,nf)*T_nI^fB(nf,IB) --> Z(ME,IB) */
//  global_dpd_->buf4_init(&D,  PSIF_CC_DINTS, 0, 20, 30, 20, 30, 0, "D <Ij|Ab> (IA,jb)");
//  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 30, 20, 30, 0, "tIAjb");
//  global_dpd_->contract444(&D,&T2, &Z, 0, 0, -1.0,1.0);
//  global_dpd_->buf4_close(&D);
//  global_dpd_->buf4_close(&T2);
//        /* Z(ME,IB) -- sort --> Z(EI,MB) */
//        /*   pq,rs                qr,ps */
//  global_dpd_->buf4_sort(&Z, PSIF_CC_TMP0, qrps, 21, 20, "Z(EI,MB)");
//  global_dpd_->buf4_close(&Z);
//
//      /* - t_M^A ( <MB||EI> + Z(EI,MB) ) --> W(EI,AB) */                      /*MB,EI*/
//  global_dpd_->buf4_init(&C,  PSIF_CC_CINTS, 0, 20, 21, 20, 21, 0, "C <IA||JB> (IA,BJ)");
//  global_dpd_->buf4_sort_axpy(&C,PSIF_CC_TMP0, rspq, 21, 20, "Z(EI,MB)",  -1.0);
//  global_dpd_->buf4_close(&C);
//  global_dpd_->buf4_init(&Z, PSIF_CC_TMP0, 0, 21, 20, 21, 20, 0, "Z(EI,MB)");
//  global_dpd_->buf4_init(&W, PSIF_CC_TMP0, 0, 21, 5, 21, 7, 0, "W'EIAB");
//  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
//  global_dpd_->contract424(&Z, &T1 ,&W,2,0,0,-1,0);
//  global_dpd_->file2_close(&T1);
//  global_dpd_->buf4_sort_axpy(&W, PSIF_CC_HBAR, pqsr, 21, 7, "WEIAB",-1);
//  global_dpd_->buf4_close(&Z);
//  global_dpd_->buf4_init(&Z, PSIF_CC_HBAR, 0, 21, 5, 21, 7,0, "WEIAB");
//  global_dpd_->buf4_axpy(&W, &Z,1);
//  global_dpd_->buf4_close(&W);
//  global_dpd_->buf4_close(&Z);
//  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21, 7, 21, 7,0, "WEIAB");
//  //global_dpd_->buf4_print(&W,"outfile",1);
//  global_dpd_->buf4_close(&W);
//  //global_dpd_->file2_mat_init(&T1);
  //global_dpd_->file2_mat_rd(&T1);

  /*
  for(Gei=0; Gei < moinfo.nirreps; Gei++){
    Gmb = Gei; // Z is totally symmetric
    global_dpd_->buf4_mat_irrep_row_init(&Z, Gei);
    global_dpd_->buf4_mat_irrep_row_init(&W, Gei);
    for(ei=0; ei< Z.params->rowtot[Gei];ei++){
      global_dpd_->buf4_mat_irrep_row_rd(&W, Gei, ei);
      global_dpd_->buf4_mat_irrep_row_rd(&Z, Gei, ei);
      for(Gm =0; Gm < moinfo.nirreps; Gm++){
        Ga = Gm; //T1 is totally symmetric
        Gb = Gm ^ Gmb;
        nrows = moinfo.avirtpi[Ga];
        ncols = moinfo.avirtpi[Gb];
        nlinks = moinfo.aoccpi[Gm];
        mb=Z.col_offset[Gei][Gm];
        ab=W.col_offset[Gei][Ga];
        if(nrows && ncols && nlinks){
          C_DGEMM('t','n',nrows,ncols,nlinks,1.0,T1.matrix[Gm][0],nrows,
              &(Z.matrix[Gei][0][mb]),ncols,1,&(W.matrix[Gei][0][ab]),ncols);
          }
        nrows = moinfo.avirtpi[Gb];
        ncols = moinfo.avirtpi[Ga];
        nlinks = moinfo.aoccpi[Gm];
        am = Z.col_offset[Gei][Ga];
        ba = W.col_offset[Gei][Gb];
        if(nrows&& ncols && nlinks){
          C_DGEMM('t','n',nrows,ncols,nlinks,-1.0,&(Z.matrix[Gei][0][am]),nrows,
              T1.matrix[Gm][0],ncols,1,&(W.matrix[Gei][0][ba]),ncols);
        }
      }
      global_dpd_->buf4_mat_irrep_row_wrt(&W,Gei,ei);
    }
    global_dpd_->buf4_mat_irrep_row_close(&W,Gei);
    global_dpd_->buf4_mat_irrep_row_close(&Z, Gei);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&W);
  global_dpd_->buf4_close(&Z);
  */
  global_dpd_->buf4_init(&W, PSIF_CC_HBAR, 0, 21,7, 21, 7, 0, "WEIAB" );
  global_dpd_->buf4_print(&W, "outfile", 1);
  global_dpd_->buf4_close(&W);

  if(params.print & 2) outfile->Printf("done.\n");
  timer_off("UHF_WABEI(NEW)");
}


/*
**
** AMJ 1/2016
*/
void build_UHF_Z1(void)
{
  dpdbuf4 T2, Z1, Fint;
  dpdfile2 T1;
  int row,col,h;
  int GI,GB,GM,GF,GA;
  int i,b,m,f,a;
  int I,B,M,F,A;

  global_dpd_->buf4_init(&T2, PSIF_CC_TAMPS, 0, 20, 20, 20, 20, 0, "tIAJB");
  global_dpd_->buf4_copy(&T2, PSIF_CC_TMP0, "Z1(IA,MF)");
  global_dpd_->buf4_close(&T2);
  global_dpd_->buf4_init(&Z1, PSIF_CC_TMP0, 0, 20, 20, 20, 20, 0, "Z1(IA,MF)");
  global_dpd_->file2_init(&T1, PSIF_CC_OEI, 0, 0, 1, "tIA");
  global_dpd_->file2_mat_init(&T1);
  global_dpd_->file2_mat_rd(&T1);

  for(h = 0; h < moinfo.nirreps; h++){
    global_dpd_->buf4_mat_irrep_init(&Z1, h);
    global_dpd_->buf4_mat_irrep_rd(&Z1, h);
    for(row = 0; row< Z1.params->rowtot[h]; row ++){
      i  = Z1.params->roworb[h][row][0];
      a  = Z1.params->roworb[h][row][1];
      I  = T1.params->rowidx[i];
      A  = T1.params->colidx[a];
      GI = T1.params->psym[i];
      GA = T1.params->qsym[a];
      for(col =0; col < Z1.params->coltot[h]; col ++){
        m = Z1.params->colorb[h][col][0];
        f = Z1.params->colorb[h][col][1];
        M = T1.params->rowidx[m];
        F = T1.params->colidx[f];
        GM = T1.params->psym[m];
        GF = T1.params->qsym[f];

        if( GI == GF && GA == GM ){
          Z1.matrix[h][row][col] -= (T1.matrix[GI][I][F] * T1.matrix[GM][M][A]);
        }
      }
    }
    global_dpd_->buf4_mat_irrep_wrt(&Z1, h);
    global_dpd_->buf4_mat_irrep_close(&Z1, h);
  }
  global_dpd_->file2_mat_close(&T1);
  global_dpd_->file2_close(&T1);
  global_dpd_->buf4_close(&Z1);
}//build_Z1

}} // namespace psi::cchbar


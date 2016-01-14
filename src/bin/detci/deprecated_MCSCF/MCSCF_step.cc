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
    \ingroup DETCAS
    \brief Enter brief description of file here
*/

#include <libqt/qt.h>
#include <cmath>
#include <psifiles.h>
#include <libpsio/psio.h>
#include <libciomr/libciomr.h>
#include <libparallel/ParallelPrinter.h>
#include "MCSCF.h"
#define EXTERN
#include "globaldefs.h"
#include "structs.h"
#include "globals.h"

namespace psi { namespace detci {

#define MO_HESS_MIN 1.0E-2


/*
** calc_orb_step()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and an approximate orbital Hessian
**
** C. David Sherrill
** April 1998
*/
void MCSCF::calc_orb_step(int npairs, double *grad, double *hess_diag, double *theta)
{

  int pair;
  double numer, denom;

  for (pair=0; pair<npairs; pair++) {
    numer = grad[pair];
    denom = hess_diag[pair];
    if (denom < 0.0) {
      outfile->Printf("Warning: MO Hessian denominator negative\n");
      denom = -denom;
    }
    if (denom < MO_HESS_MIN) {
      outfile->Printf("Warning: MO Hessian denominator too small\n");
      denom = MO_HESS_MIN;
    }
    theta[pair] =  - numer / denom;
  }

}


/*
** calc_orb_step_full()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and a square matrix orbital Hessian
**
** C. David Sherrill
** September 2003
*/
void MCSCF::calc_orb_step_full(int npairs, double *grad, double **hess, double *theta)
{
  double **hess_inv;
  double **hess_copy; /* for testing! */
  int i,j;
  double tval;
  int solved;
  double *BVector;
  int *pivots;
  double hess_det = 1.0;
  int *indx;
  double biggest_step;

  hess_copy = block_matrix(npairs, npairs);
  indx = init_int_array(npairs);

  for (i=0; i<npairs; i++) {
    for (j=0; j<npairs; j++) {
      hess_copy[i][j] = hess[i][j];
    }
  }

  ludcmp(hess_copy,npairs,indx,&hess_det);
  for (j=0;j<npairs;j++){
    hess_det *= hess_copy[j][j];
  }
  outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);

  /*
     if the orbital Hessian is not positive definite, we may have some
     trouble converging the orbitals.  Guarantee it's positive definite
     by levelshifting
  */
  if (MCSCF_Parameters.level_shift) {
    while (hess_det < MCSCF_Parameters.determ_min) {
      outfile->Printf("Level shifting the hessian by %8.3E\n",MCSCF_Parameters.shift);
      for (i=0;i<npairs;i++) {
        hess[i][i] += MCSCF_Parameters.shift;
      }
      for (i=0;i<npairs;i++) {
        for (j=0;j<npairs;j++) {
          hess_copy[i][j] = hess[i][j];
        }
      }
      ludcmp(hess_copy,npairs,indx,&hess_det);
      for (j=0;j<npairs;j++){
        hess_det *= hess_copy[j][j];
      }
      outfile->Printf("The determinant of the hessian is %8.3E\n",hess_det);
    }
    outfile->Printf("Determinant of the hessian is greater than %8.3E\n",
      MCSCF_Parameters.determ_min);
  }


  /* let's re-copy hess into hess_copy because we ludcmp'd hess_copy */
  for (i=0;i<npairs;i++) {
    for (j=0;j<npairs;j++) {
      hess_copy[i][j] = hess[i][j];
    }
  }

  if (!MCSCF_Parameters.invert_hessian) { /* solve H delta = - g */
    outfile->Printf("Solving system of linear equations for orbital step...");
    BVector = init_array(npairs);
    pivots = init_int_array((npairs * (npairs - 1))/2);
    for(i=0;i<npairs;i++){
      BVector[i] = -grad[i];
      theta[i] = 0.0;
    }
    solved = C_DGESV(npairs,1,&(hess_copy[0][0]),npairs,pivots,BVector,npairs);
    if (solved == 0) {
      outfile->Printf("equations solved!\n");
      for(i=0;i<npairs;i++) {
        theta[i] = BVector[i];
      }
    }
    else {
      throw PsiException("FAILED TO SOLVE FOR THETA VALUES\n", __FILE__, __LINE__) ;
    }
    free(BVector);
    free(pivots);
  } /* end solution of linear equations H delta = -g */

  else { /* direct inversion of orbital Hessian */
    outfile->Printf("Attempting to directly invert the Hessian matrix\n");
    hess_inv = block_matrix(npairs,npairs);

    /* note: this will destroy hessian matrix; don't use it again later! */
    invert_matrix(hess_copy,hess_inv,npairs,"outfile");

    /* debug check */
    mmult(hess_inv,0,hess,0,hess_copy,0,npairs,npairs,npairs,0);
    outfile->Printf("Hessian * Hessian inverse = \n");
    print_mat(hess_copy,npairs,npairs,"outfile");
    outfile->Printf("\n");

    /* step = - B^{-1} * g */
    zero_arr(theta,npairs);
    /* the line below seems to have trouble unless I take out the -1
       which should be there, and even then it's not really working */
    /*
    C_DGEMV('n',npairs,npairs,-1.0,hess_inv[0],npairs,grad,1,0.0,theta,1);
    */

    for (i=0; i<npairs; i++) {
      tval = 0.0;
      for (j=0; j<npairs; j++) {
        tval += hess_inv[i][j] * grad[j];
      }
      theta[i] = - tval;
    }
    free_block(hess_inv);
  } /* end direct inversion of Hessian */

  /* make sure the step is not too big */
  biggest_step = 0.0;
  for (i=0; i<npairs; i++) {
    tval = theta[i];
    if (fabs(tval) > biggest_step) biggest_step = fabs(tval);
  }
  outfile->Printf("\nLargest step in theta space is %12.6lf \n", biggest_step);
  if (biggest_step > MCSCF_Parameters.step_max) {
    outfile->Printf("Scaling the step\n");
    for (i=0;i<npairs;i++) {
      theta[i] = theta[i] * MCSCF_Parameters.step_max / biggest_step;
    }
  }

  free_block(hess_copy);
  free(indx);
}


/*
** calc_orb_step_bfgs()
**
** This function calculates the step in theta space for the orbitals
** given the orbital gradient and a square matrix orbital Hessian INVERSE.
** With the inverse already available, this is very straightforward.
**
** C. David Sherrill
** March 2004
*/
void MCSCF::calc_orb_step_bfgs(int npairs, double *grad, double **hess, double *theta)
{

  int i, j;
  double tval, biggest_step;

  for (i=0; i<npairs; i++) {
    tval = 0.0;
    for (j=0; j<npairs; j++) {
      tval += hess[i][j] * grad[j];
    }
    theta[i] = - tval;
  }

  /* make sure the step is not too big */
  biggest_step = 0.0;
  for (i=0; i<npairs; i++) {
    tval = theta[i];
    if (fabs(tval) > biggest_step) biggest_step = fabs(tval);
  }
  outfile->Printf("\nLargest step in theta space is %12.6lf \n", biggest_step);
  if (biggest_step > MCSCF_Parameters.step_max) {
    outfile->Printf("Largest allowed step %12.6lf --- scaling the step\n",
      MCSCF_Parameters.step_max);
    for (i=0;i<npairs;i++) {
      theta[i] = theta[i] * MCSCF_Parameters.step_max / biggest_step;
    }
  }

}


/*
** print_step
**
** This function prints out the information for a given orbital iteration
*/

void MCSCF::print_step(int iter, int npairs, int steptype, OutFile& IterSummaryOut)
{

   IterSummaryOut.Printf("%5d %5d %14.9lf %14.9lf %20.12lf", iter+1,
     npairs, MCSCF_CalcInfo.scaled_mo_grad_rms, MCSCF_CalcInfo.mo_grad_rms,
     MCSCF_CalcInfo.energy);

  if (steptype == 0)
    IterSummaryOut.Printf(" %9s\n", "CONV");
  else if (steptype == 1)
    IterSummaryOut.Printf(" %9s\n", "NR");
  else if (steptype == 2)
    IterSummaryOut.Printf(" %9s\n", "DIIS");
  else {
    IterSummaryOut.Printf(" ?\n");
    outfile->Printf("(print_step): Unrecognized steptype %d\n", steptype);
  }

}


}} // end namespace psi::detci


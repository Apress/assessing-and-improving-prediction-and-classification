/******************************************************************************/
/*                                                                            */
/*  MutInf_C - Mutual information for continuous data                         */
/*                                                                            */
/******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "info.h"

#define DEBUG 0

/*
--------------------------------------------------------------------------------

   MutualInformationParzen

--------------------------------------------------------------------------------
*/

double outer_crit ( double t ) ;
double inner_crit ( double t ) ;

// These are for passing to the integrands inner_crit() and outer_crit()
static ParzDens_1 *this_dens_dep = NULL ;
static ParzDens_1 *this_dens_trial = NULL ;
static ParzDens_2 *this_dens_bivar = NULL ;

static double accuracy ;   // Set integration accuracy per n

MutualInformationParzen::MutualInformationParzen (
   int nn ,              // Number of cases
   double *dep_vals ,    // They are here
   int ndiv )            // Number of divisions of range, typically 5-10
{
   n = nn ;
   n_div = ndiv ;
   depvals = NULL ;
   dens_dep = NULL ;

   MEMTEXT ( "MutualInformationParzen constructor" ) ;
   depvals = (double *) MALLOC ( n * sizeof(double) ) ;
   assert (depvals != NULL) ;

   memcpy ( depvals , dep_vals , n * sizeof(double) ) ;

   dens_dep = new ParzDens_1 ( n , depvals , n_div ) ;
   assert (dens_dep != NULL) ;
}

MutualInformationParzen::~MutualInformationParzen ()
{
   MEMTEXT ( "MutualInformationParzen destructor" ) ;
   FREE ( depvals ) ;
   delete dens_dep ;
}

double MutualInformationParzen::mut_inf ( double *x )
{
   double criterion ;

   MEMTEXT ( "MutualInformationParzen::compute()" ) ;

   this_dens_dep = dens_dep ;

   this_dens_trial = new ParzDens_1 ( n , x , n_div ) ;
   assert (this_dens_trial != NULL) ;

   this_dens_bivar = new ParzDens_2 ( n , depvals , x , n_div ) ;
   assert (this_dens_bivar != NULL) ;

   accuracy = (n > 200)  ?  1.e-5 : 1.e-6 ;

   criterion = integrate ( this_dens_trial->low , this_dens_trial->high ,
                  (this_dens_trial->high - this_dens_trial->low) / 10.0 ,
                  1.e-6 , accuracy , outer_crit ) ;

   delete this_dens_trial ;
   delete this_dens_bivar ;

   return criterion ;
}

/*
   This pair of routines are called by integrate() to return the integrand.
   inner_crit() does the actual work of defining the function being integrated.
   outer_crit is the wrapper for 2-D integration.

   NOTE... The ParzDens_1 and ParzDens_2 classes in PARZDENS.CPP convert the
   data to a normal distribution, which is usually good.  In this special case,
   it is silly to use a Parzen window to estimate the one-dimensional density,
   because we know what it is: normal!  So in this code I have used #if 0 to
   remove the Parzen computation, and replaced it with a straight evaluation
   of the normal density.  I left the code here for users who may wish to
   remove the normalization from the Parzen window to experiment.  If you do so,
   be sure to change the #if 0 to #if 1 here.
*/

static double this_x, this_px ;  // Needed for two-dimensional integration

double inner_crit ( double t )
{
   double py, pxy, term ;
#if 0
   py = this_dens_dep->density ( t ) ; // General case
#else
   py = exp ( -0.5 * t * t ) / sqrt ( 2.0 * PI ) ; // Only if Parzen normalized
#endif
   pxy = this_dens_bivar->density ( t , this_x ) ;
   term = this_px * py ;
   if (term < 1.e-30)
      term = 1.e-30 ;
   term = pxy / term ;
   if (term < 1.e-30)
      term = 1.e-30 ;
   return pxy * log ( term ) ;
}

double outer_crit ( double t )
{
   double val, high, low ;

   high = this_dens_dep->high ;
   low = this_dens_dep->low ;
   this_x = t ;
#if 0
   this_px = this_dens_trial->density ( this_x ) ;
#else
   this_px = exp ( -0.5 * this_x * this_x ) / sqrt ( 2.0 * PI ) ;
#endif
   val = integrate ( low , high , (high - low) / 10.0 , 1.e-7 , 0.1 * accuracy ,
                     inner_crit ) ;
   return val ;
}

/*
--------------------------------------------------------------------------------

   Adaptive partitioning method of Darbellay and Vajda
   IEEE Transactions on Information Theory Vol. 45 No. 4  May 1999

--------------------------------------------------------------------------------
*/

MutualInformationAdaptive::MutualInformationAdaptive (
   int nn ,              // Number of cases
   double *dep_vals ,    // They are here
   int respect_ties ,    // Treat ties as if discrete classes?
   double crit )         // Chi-square test criterion, typically 6.0
{
   int i, *indices ;
   double *work ;

   MEMTEXT ( "MutualInformationAdaptive constructor" ) ;

   n = nn ;
   chi_crit = crit ;

/*
   Convert the 'dependent' variable to ranks
*/

   indices = (int *) MALLOC ( n * sizeof(int) ) ;
   assert ( indices != NULL ) ;

   work = (double *) MALLOC ( n * sizeof(double) ) ;
   assert ( work != NULL ) ;

   y = (int *) MALLOC ( n * sizeof(int) ) ;
   assert ( y != NULL ) ;

   if (respect_ties) {
      y_tied = (int *) MALLOC ( n * sizeof(int) ) ;
      assert ( y_tied != NULL ) ;
      }
   else
      y_tied = NULL ;

   for (i=0 ; i<n ; i++) {
      work[i] = dep_vals[i] ;
      indices[i] = i ;
      }

   qsortdsi ( 0 , n-1 , work , indices ) ;

   for (i=0 ; i<n ; i++) {
      y[indices[i]] = i ;  // We now have ranks
      if (! respect_ties)
         continue ;
      if (i < n-1  &&
            work[i+1] - work[i] < 1.e-12 * (1.0+fabs(work[i])+fabs(work[i+1])))
         y_tied[i] = 1 ;
      else
         y_tied[i] = 0 ;
      }

#if DEBUG
   if (respect_ties) {
      int k = 0 ;
      printf ( "\nConstructor ties..." ) ;
      for (i=0 ; i<n-1 ; i++) {
#if 0
         printf ( "\n%16.12lf  %16.12le  %16.12le  %16.12lf  %d",
            work[i], fabs(work[i]-work[i+1]),
            1.e-12 * (1.0+fabs(work[i])+fabs(work[i+1])),
            fabs(work[i]-work[i+1]), y_tied[i] ) ;
#endif
         if (y_tied[i])
            ++k ;
         }
      printf ( "\nY tied count = %d of %d", k, n ) ;
      }
#endif

   FREE ( indices ) ;
   FREE ( work ) ;
}

MutualInformationAdaptive::~MutualInformationAdaptive ()
{
   MEMTEXT ( "MutualInformationAdaptive destructor" ) ;
   FREE ( y ) ;
   if (y_tied != NULL)
      FREE ( y_tied ) ;
}

double MutualInformationAdaptive::mut_inf ( double *xraw , int respect_ties )
{
   int i, k, ix, iy, *indices, *current_indices, *x, nstack, splittable ;
   int fullXstart, fullXstop, fullYstart, fullYstop, ipos ;
   int trialXstart[4], trialXstop[4], trialYstart[4], trialYstop[4] ;
   int ipx, ipy, xcut[4], ycut[4], iSubRec, *x_tied, ioff ;
   int X_AllTied, Y_AllTied ;
   int centerX, centerY, currentDataStart, currentDataStop ;
   int actual[4], actual44[16] ;
   double *work, expected[16], diff, testval, xfrac[4], yfrac[4] ;
   double px, py, pxy, MI ;

struct {
   int Xstart ;     // X value (rank) at which this rectangle starts
   int Xstop ;      // And stops
   int Ystart ;     // Ditto for Y
   int Ystop ;
   int DataStart ;  // Starting index into indices for the cases in this
   int DataStop ;   // rectangle, and the (inclusive) ending index
} stack[256] ;  // Wildly conservative

   MEMTEXT ( "MutualInformationAdaptive::compute()" ) ;

   indices = (int *) MALLOC ( n * sizeof(int) ) ;
   assert ( indices != NULL ) ;

   current_indices = (int *) MALLOC ( n * sizeof(int) ) ;
   assert ( current_indices != NULL ) ;

   work = (double *) MALLOC ( n * sizeof(double) ) ;
   assert ( work != NULL ) ;

   x = (int *) MALLOC ( n * sizeof(int) ) ;
   assert ( x != NULL ) ;

   if (respect_ties) {
      x_tied = (int *) MALLOC ( n * sizeof(int) ) ;
      assert ( x_tied != NULL ) ;
      }
   else
      x_tied = NULL ;

/*
   Convert this 'independent' variable to ranks
*/

   for (i=0 ; i<n ; i++) {
      work[i] = xraw[i] ;
      indices[i] = i ;
      }

   qsortdsi ( 0 , n-1 , work , indices ) ;

   for (i=0 ; i<n ; i++) {
      x[indices[i]] = i ;  // We now have ranks
      if (! respect_ties)
         continue ;
      if (i < n-1  &&
            work[i+1] - work[i] < 1.e-12 * (1.0+fabs(work[i])+fabs(work[i+1])))
         x_tied[i] = 1 ;
      else
         x_tied[i] = 0 ;
      }

#if DEBUG
   if (respect_ties) {
      k = 0 ;
      printf ( "\nCompute ties..." ) ;
      for (i=0 ; i<n-1 ; i++) {
#if 0
         printf ( "\n%16.12lf  %16.12le  %16.12le  %16.12lf  %d",
            work[i], fabs(work[i]-work[i+1]),
            1.e-12 * (1.0+fabs(work[i])+fabs(work[i+1])),
            fabs(work[i]-work[i+1]), x_tied[i] ) ;
#endif
         if (x_tied[i])
            ++k ;
         }
      printf ( "\nX tied count = %d of %d", k, n ) ;
   }
#endif

/*
   The array 'indices' indexes the cases.
   The contents of a rectangle will always be defined by starting and stopping
   (inclusive) positions in this array.  These case indices will be shuffled
   as needed to place the correct cases in a rectangle.  Note that this
   method implies automatic nesting of rectangles, as defined by nested subsets
   of this array.
   Initialize indices to be the entire input rectangle.
*/

   for (i=0 ; i<n ; i++)
      indices[i] = i ;

/*
   Initialize the rectangle stack to have one entry, the entire rectangle
*/

   stack[0].Xstart = 0 ;
   stack[0].Xstop = n-1 ;
   stack[0].Ystart = 0 ;
   stack[0].Ystop = n-1 ;
   stack[0].DataStart = 0 ;
   stack[0].DataStop = n-1 ;
   nstack = 1 ;   // This is the top-of-stack pointer: One item in the stack

/*
   Main loop
   As long as there is at least one rectangle in the stack awaiting possible
   splitting, we process it.
*/

   MI = 0.0 ;   // Will cumulate mutual information here
   while (nstack > 0) {

      // Get the rectangle pushed onto the stack most recently
      --nstack ;                                    // Pop the rectangle
      fullXstart = stack[nstack].Xstart ;           // Starting X rank
      fullXstop = stack[nstack].Xstop ;             // And ending
      fullYstart = stack[nstack].Ystart ;           // Ditto for Y
      fullYstop = stack[nstack].Ystop ;
      currentDataStart = stack[nstack].DataStart ;  // The cases start here
      currentDataStop = stack[nstack].DataStop ;    // And end here

#if DEBUG
      printf ( "\n\nstack=%d   data: %d %d", nstack, currentDataStart, currentDataStop ) ;
      printf ( "\nX: %d %d    Y: %d %d", fullXstart, fullXstop, fullYstart, fullYstop ) ;
#endif

      // Do a trial 2x2 split.  Adjust the split so it does not split ties.

      centerX = (fullXstart + fullXstop) / 2 ;
      X_AllTied = (x_tied != NULL)  &&  (x_tied[centerX] != 0) ;
      if (X_AllTied) {
         for (ioff=1 ; centerX-ioff >= fullXstart ; ioff++) {
            if (! x_tied[centerX-ioff]) {
               X_AllTied = 0 ;
               centerX -= ioff ;
#if DEBUG
               printf ( "\n  X offset = -%d", ioff ) ;
#endif
               break ;
               }
            if (centerX + ioff == fullXstop)
               break ;
            if (! x_tied[centerX+ioff]) {
               X_AllTied = 0 ;
               centerX += ioff ;
#if DEBUG
               printf ( "\n  X offset = %d", ioff ) ;
#endif
               break ;
               }
            }
         }

      centerY = (fullYstart + fullYstop) / 2 ;
      Y_AllTied = (y_tied != NULL)  &&  (y_tied[centerY] != 0) ;
      if (Y_AllTied) {
         for (ioff=1 ; centerY-ioff >= fullYstart ; ioff++) {
            if (! y_tied[centerY-ioff]) {
               Y_AllTied = 0 ;
               centerY -= ioff ;
#if DEBUG
               printf ( "\n  Y offset = -%d", ioff ) ;
#endif
               break ;
               }
            if (centerY + ioff == fullYstop)
               break ;
            if (! y_tied[centerY+ioff]) {
               Y_AllTied = 0 ;
               centerY += ioff ;
#if DEBUG
               printf ( "\n  Y offset = %d", ioff ) ;
#endif
               break ;
               }
            }
         }

#if DEBUG
      printf ( "\n  X_AllTied=%d  Y_AllTied=%d", X_AllTied, Y_AllTied ) ;
#endif

#if 0 /*!!!!!!*/
      X_AllTied = Y_AllTied = 0 ;
      centerX = (fullXstart + fullXstop) / 2 ;
      centerY = (fullYstart + fullYstop) / 2 ;
#endif

      if (X_AllTied  ||  Y_AllTied)
         splittable = 0 ;
      else {
         trialXstart[0] = trialXstart[1] = fullXstart ;
         trialXstop[0] = trialXstop[1] = centerX ;
         trialXstart[2] = trialXstart[3] = centerX+1 ;
         trialXstop[2] = trialXstop[3] = fullXstop ;
         trialYstart[0] = trialYstart[2] = fullYstart ;
         trialYstop[0] = trialYstop[2] = centerY ;
         trialYstart[1] = trialYstart[3] = centerY+1 ;
         trialYstop[1] = trialYstop[3] = fullYstop ;

         // Compute the expected count in each of the four sub-rectangles
         for (i=0 ; i<4 ; i++)
            expected[i] = (currentDataStop - currentDataStart + 1) *
                 (trialXstop[i]-trialXstart[i]+1.0) / (fullXstop-fullXstart+1.0) *
                 (trialYstop[i]-trialYstart[i]+1.0) / (fullYstop-fullYstart+1.0) ;

         // Compute the actual count in each of the four sub-rectangles
         actual[0] = actual[1] = actual[2] = actual[3] = 0 ;
         for (i=currentDataStart ; i<=currentDataStop ; i++) {
            k = indices[i] ;
            if (x[k] <= centerX) {
               if (y[k] <= centerY)
                  ++actual[0] ;
               else 
                  ++actual[1] ;
               }
            else {
               if (y[k] <= centerY)
                  ++actual[2] ;
               else 
                  ++actual[3] ;
               }
            }

         // Do the chi-square test for this 2x2 split
         testval = 0.0 ;
         for (i=0 ; i<4 ; i++) {
            diff = fabs ( actual[i] - expected[i] ) - 0.5 ;
            testval += diff * diff / expected[i] ;
            }

#if DEBUG
         printf ( "\nCenter: %d %d", centerX, centerY ) ;
         printf ( "\nExpected: %.2lf %.2lf %.2lf %.2lf",
                  expected[0], expected[1], expected[2], expected[3] ) ;
         printf ( "\nActual: %d %d %d %d    Chisq=%.4lf",
                  actual[0], actual[1], actual[2], actual[3], testval ) ;
#endif

         splittable = (testval > chi_crit)  ?  1 : 0 ;

         // If the 2x2 splitting test said that the full rectangle is uniform,
         // and hence should not be split, but the rectangle is big, we should
         // do a finer test.  It may be that the full rectangle really is not
         // uniform, but the four sub-rectangles just happen to have the expected
         // number of cases.  Go to a 4x4 split.

         if (! splittable && fullXstop-fullXstart > 30 && fullYstop-fullYstart > 30) {
            ipx = fullXstart - 1 ;
            ipy = fullYstart - 1 ;
            for (i=0 ; i<4 ; i++) {
               xcut[i] = (fullXstop - fullXstart + 1) * (i+1) / 4 + fullXstart - 1 ;
               xfrac[i] = (xcut[i] - ipx) / (fullXstop - fullXstart + 1.0) ;
               ipx = xcut[i] ;
               ycut[i] = (fullYstop - fullYstart + 1) * (i+1) / 4 + fullYstart - 1 ;
               yfrac[i] = (ycut[i] - ipy) / (fullYstop - fullYstart + 1.0) ;
               ipy = ycut[i] ;
               }
            // Compute expected counts
            for (ix=0 ; ix<4 ; ix++) {
               for (iy=0 ; iy<4 ; iy++) {
                  expected[ix*4+iy] = xfrac[ix] * yfrac[iy] *
                                      (currentDataStop - currentDataStart + 1) ;
                  actual44[ix*4+iy] = 0 ;
                  }
               }
            // Compute actual counts
            for (i=currentDataStart ; i<=currentDataStop ; i++) {
               k = indices[i] ;
               for (ix=0 ; ix<3 ; ix++) {
                  if (x[k] <= xcut[ix])
                     break ;
                  }
               for (iy=0 ; iy<3 ; iy++) {
                  if (y[k] <= ycut[iy])
                     break ;
                  }
               ++actual44[ix*4+iy] ;
               }
            // Compute the chi-square test
            testval = 0.0 ;
            for (ix=0 ; ix<4 ; ix++) {
               for (iy=0 ; iy<4 ; iy++) {
                  diff = fabs ( actual44[ix*4+iy] - expected[ix*4+iy] ) - 0.5 ;
                  testval += diff * diff / expected[ix*4+iy] ;
                  }
               }
#if DEBUG
            printf ( "\n2x2 test failed.  4x4 = %.4lf", testval ) ;
#endif
            splittable = (testval > 3 * chi_crit)  ?  1 : 0 ; // Hard-coded as 22 in original algorithm
            } // If trying 4x4 split
         } // Else not all tied

      // If splittable, splint this full rectangle into 2x2 sub-rectangles

      if (splittable) {
         
         // Make a copy of the case indices of this rectangle
         // We will rearrange the cases in the block comprising this rectangle
         // so that the new order reflects the partitioning.

         for (i=currentDataStart ; i<=currentDataStop ; i++)
            current_indices[i] = indices[i] ;

         ipos = currentDataStart ;  // Will rearrange indices starting here
         for (iSubRec=0 ; iSubRec<4 ; iSubRec++) { // Check all 4 sub-rectangles
#if DEBUG
            printf ( "\nSubrec %d  n=%d", iSubRec, actual[iSubRec] ) ;
#endif

            // Is this sub-rectangle big enough to be worth pushing onto stack?
            if (actual[iSubRec] >= 3) {
#if DEBUG
               printf ( "   Pushing to %d (%d %d)", nstack+1, ipos, ipos + actual[iSubRec] - 1 ) ;
#endif
               // Push this sub-rectangle onto the stack
               stack[nstack].Xstart = trialXstart[iSubRec] ;
               stack[nstack].Xstop = trialXstop[iSubRec] ;
               stack[nstack].Ystart = trialYstart[iSubRec] ;
               stack[nstack].Ystop = trialYstop[iSubRec] ;
               stack[nstack].DataStart = ipos ;
               stack[nstack].DataStop = ipos + actual[iSubRec] - 1 ;
               ++nstack ;
               // Act according to which sub-rectangle this is
               // We rearrange to bring this sub-rectangle to the correct place
               if (iSubRec == 0) {
                  for (i=currentDataStart ; i<=currentDataStop ; i++) {
                     k = current_indices[i] ;
                     if (x[k] <= centerX  &&  y[k] <= centerY)
                        indices[ipos++] = current_indices[i] ;
                     }
                  }
               else if (iSubRec == 1) {
                  for (i=currentDataStart ; i<=currentDataStop ; i++) {
                     k = current_indices[i] ;
                     if (x[k] <= centerX  &&  y[k] > centerY)
                        indices[ipos++] = current_indices[i] ;
                     }
                  }
               else if (iSubRec == 2) {
                  for (i=currentDataStart ; i<=currentDataStop ; i++) {
                     k = current_indices[i] ;
                     if (x[k] > centerX  &&  y[k] <= centerY)
                        indices[ipos++] = current_indices[i] ;
                     }
                  }
               else {  // iSubRec == 3
                  for (i=currentDataStart ; i<=currentDataStop ; i++) {
                     k = current_indices[i] ;
                     if (x[k] > centerX  &&  y[k] > centerY)
                        indices[ipos++] = current_indices[i] ;
                     }
                  }
               } // If this sub-rectangle is large enough to be worth pushing

            else {  // This sub-rectangle is small, so get its contribution now
               if (actual[iSubRec] > 0) { // It only contributes if it has cases
                  px = (trialXstop[iSubRec] - trialXstart[iSubRec] + 1.0) / n ;
                  py = (trialYstop[iSubRec] - trialYstart[iSubRec] + 1.0) / n ;
                  pxy = (double) actual[iSubRec] / n ;
                  MI += pxy * log ( pxy / (px * py) ) ;
#if DEBUG
                  printf ( "   Inc=%.4lf  MI=%.4lf  (px=%.4lf py=%.4lf pxy=%.4lf)",
                             pxy * log ( pxy / (px * py) ), MI, px, py, pxy ) ;
#endif
                  }
               } // Else this sub-rectangle is too small to push, so process it
            } // For all 4 sub-rectangles
         } // If splitting
      else {  // Else the chi-square tests failed, so we do not split
         px = (fullXstop - fullXstart + 1.0) / n ;
         py = (fullYstop - fullYstart + 1.0) / n ;
         pxy = (currentDataStop - currentDataStart + 1.0) / n ;
         MI += pxy * log ( pxy / (px * py) ) ;
#if DEBUG
         printf ( "\nNo split: Inc=%.4lf  MI=%.4lf  (px=%.4lf py=%.4lf pxy=%.4lf)",
                             pxy * log ( pxy / (px * py) ), MI, px, py, pxy ) ;
#endif
         }
      } // While rectangles in the stack

   FREE ( indices ) ;
   FREE ( current_indices ) ;
   FREE ( work ) ;
   FREE ( x ) ;
   if (x_tied != NULL)
      FREE ( x_tied ) ;

   return MI ;
}

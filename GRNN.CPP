/******************************************************************************/
/*                                                                            */
/*  GRNN - General Regression Neural Network                                  */
/*                                                                            */
/*         This implementation uses a primitive annealing training method     */
/*         It is not good enough for serious work.                            */
/*         It could be brought up to decent standards by using the included   */
/*         annealing just as a starting point, following it with good         */
/*         refinement.                                                        */
/*         Also, a user friendly version would have provision for progress    */
/*         reports and user interruption.  And last but not least, error      */
/*         checks like failure to allocate sufficient memory should be        */
/*         included in any serious inplementation.                            */
/*                                                                            */
/*                                                                            */
/*  To use this class:                                                        */
/*    1) Construct a new instance of the class                                */
/*    2) Call add_case() exactly ncases times, each time providing the        */
/*       nin+nout vector of inputs and outputs.                               */
/*    3) Call train()                                                         */
/*    4) Call predict() as many times as desired                              */
/*    5) Optionally, call reset() and go to step 2                            */
/*                                                                            */
/*  This does not include any checks for insufficient memory.                 */
/*  It also assumes that the user calls add_case exactly ncases times         */
/*  and does not check for failure to do so.                                  */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "grnn.h"

double normal () ;
#define EPS1 1.e-180

/*
--------------------------------------------------------------------------------

   Constructor, destructor, reset(), add_case()

--------------------------------------------------------------------------------
*/

GRNN::GRNN ( int ncase , int nin , int nout )
{
   ncases = ncase ;
   ninputs = nin ;
   noutputs = nout ;
   tset = (double *) malloc ( ncases * (ninputs + noutputs) * sizeof(double) ) ;
   sigma = (double *) malloc ( ninputs * sizeof(double) ) ;
   outwork = (double *) malloc ( noutputs * sizeof(double) ) ;
   reset () ;
}


GRNN::~GRNN ()
{
   if (tset != NULL)
      free ( tset ) ;
   if (sigma != NULL)
      free ( sigma ) ;
   if (outwork != NULL)
      free ( outwork ) ;
}

/*
   This discards any existing training data.
   It does not need to be called after construction, but it must
   be called if the user wants to reuse the GRNN object for a new dataset.
*/

void GRNN::reset ()
{
   nrows = 0 ;      // No rows (via add_case()) yet present
   trained = 0 ;    // Training not done yet
}

/*
   Build the training set one case at a time.
   The user must call this member EXACTLY ncases times after construction
   or a call to reset(), and before a call to train().
*/

void GRNN::add_case ( double *newcase )
{
   if (nrows >= ncases)  // Careful user never lets this happen
      return ;           // But cheap insurance

   memcpy ( tset + nrows * (ninputs + noutputs) , newcase ,
            (ninputs + noutputs) * sizeof(double) ) ;
   ++nrows ;
}


/*
--------------------------------------------------------------------------------

   predict() - Given an input vector, compute output using trained model

--------------------------------------------------------------------------------
*/

void GRNN::predict (
   double *input ,     // Input vector
   double *output      // Returned output
   )
{
   int icase, iout, ivar ;
   double *dptr, diff, dist, psum ;

   for (iout=0 ; iout<noutputs ; iout++) // For each output
      output[iout] = 0.0 ;               // Will sum kernels here
   psum = 0.0 ;                          // Denominator sum

   for (icase=0 ; icase<ncases ; icase++) {  // Do all training cases

      dptr = tset + (ninputs + noutputs) * icase ; // Point to this case
      dist = 0.0 ;                           // Will sum distance here
      for (ivar=0 ; ivar<ninputs ; ivar++) { // All variables in this case
         diff = input[ivar] - dptr[ivar] ;   // Input minus case
         diff /= sigma[ivar] ;               // Scale per sigma
         dist += diff * diff ;               // Cumulate Euclidean distance
         }

      dist = exp ( -dist ) ;                // Apply the Gaussian kernel

      if (dist < EPS1)                      // If this case is far from all
         dist = EPS1 ;                      // prevent zero density

      dptr += ninputs ;                     // Outputs stored after inputs
      for (ivar=0 ; ivar<noutputs ; ivar++) // For every output variable
         output[ivar] += dist * dptr[ivar] ;// Cumulate numerator
      psum += dist ;                        // Cumulate denominator
      } // For all training cases

   for (ivar=0 ; ivar<noutputs ; ivar++)
      output[ivar] /= psum ;
}

/*
--------------------------------------------------------------------------------

   execute() - Given sigma weights, pass through the training set, return MSE.

--------------------------------------------------------------------------------
*/

double GRNN::execute ()
{
   int itest, icase, iout, ivar ;
   double *dptr, *tptr, diff, dist, psum, err ;

   err = 0.0 ;

   for (itest=0 ; itest<ncases ; itest++) {
      tptr = tset + (ninputs + noutputs) * itest ; // Test case

      for (iout=0 ; iout<noutputs ; iout++) // For each output
         outwork[iout] = 0.0 ;              // Will sum kernels here
      psum = 0.0 ;                          // Denominator sum

      for (icase=0 ; icase<ncases ; icase++) {  // Do all training cases

         if (icase == itest)                    // Do not include test case
            continue ;                          // In trial kernel

         dptr = tset + (ninputs + noutputs) * icase ; // Point to this case
         dist = 0.0 ;                           // Will sum distance here
         for (ivar=0 ; ivar<ninputs ; ivar++) { // All variables in this case
            diff = tptr[ivar] - dptr[ivar] ;    // Test case minus kernel case
            diff /= sigma[ivar] ;               // Scale per sigma
            dist += diff * diff ;               // Cumulate Euclidean distance
            }

         dist = exp ( -dist ) ;                // Apply the Gaussian kernel

         if (dist < EPS1)                      // If this case is far from all
            dist = EPS1 ;                      // prevent zero density

         dptr += ninputs ;                     // Outputs stored after inputs
         for (ivar=0 ; ivar<noutputs ; ivar++) // For every output variable
            outwork[ivar] += dist * dptr[ivar] ;// Cumulate numerator
         psum += dist ;                        // Cumulate denominator
         } // For all training cases

      tptr += ninputs ;                        // Outputs stored after inputs
      for (ivar=0 ; ivar<noutputs ; ivar++) {
         outwork[ivar] /= psum ;               // Predicted output
         diff = outwork[ivar] - tptr[ivar] ;   // Predicted minus actual
         err += diff * diff ;                  // Cumulate squared error
         }

      } // For all test cases

   err /= ncases * noutputs ;                  // MSE

   return err ;
}

/*
--------------------------------------------------------------------------------

   train() - Train the model

   After add_case has been called exactly ncases times, this must be called
   to train the model.

   This routine is the weak point in this GRNN class.  The training algorithm
   is relatively slow and inaccurate.
   It is an excellent starting point for refinement, having a high probability
   of finding a solution near a global minimum.  Addition of a good refinement
   algorithm would make this class usable in real applications.

--------------------------------------------------------------------------------
*/

void GRNN::anneal_train (
   int n_outer ,      // Number of outer loop iterations, perhaps 10-20
   int n_inner ,      // Number of inner loop iterations, perhaps 100-10000
   double start_std   // Starting standard deviation of log weights, about 3.0
   )
{
   int i, inner, outer ;
   double error, best_error, std, *best_wts, *test_wts, *center ;

/*
   Best_wts keeps track of the best (log) sigma weights.
   Center is the center around which perturbation is done.
   It starts at zero.  After completion of each pass through the inner loop
   it is changed to best_wts.
*/

   best_wts = (double *) malloc ( ninputs * sizeof(double) ) ;
   test_wts = (double *) malloc ( ninputs * sizeof(double) ) ;
   center = (double *) malloc ( ninputs * sizeof(double) ) ;

   for (i=0 ; i<ninputs ; i++)
      center[i] = 0.0 ;

   best_error = -1.0 ;
   std = start_std ;

   for (outer=0 ; outer<n_outer ; outer++) {
      for (inner=0 ; inner<n_inner ; inner++) {

         for (i=0 ; i<ninputs ; i++) {
            test_wts[i] = center[i] + std * normal() ;
            sigma[i] = exp ( test_wts[i] ) ;
            }

         error = execute () ;
         if ((best_error < 0.0)  ||  (error < best_error)) {
            best_error = error ;
            memcpy ( best_wts , test_wts , ninputs * sizeof(double) ) ;
            }
         } // For inner loop iterations
      memcpy ( center , best_wts , ninputs * sizeof(double) ) ;
      std *= 0.7 ;
      } // For outer loop iterations

   for (i=0 ; i<ninputs ; i++)
      sigma[i] = exp ( best_wts[i] ) ;

   trained = 1 ;    // Training complete
   free ( best_wts ) ;
   free ( test_wts ) ;
   free ( center ) ;
}

/*
   This is customized for this demonstration
*/

void GRNN::train ()
{
   anneal_train ( 10 , 100 , 3.0 ) ;
}

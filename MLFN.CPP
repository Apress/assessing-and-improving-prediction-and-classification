/******************************************************************************/
/*                                                                            */
/*  MLFN - Multiple Layer Feedforward Network                                 */
/*                                                                            */
/*         This implementation uses a primitive annealing training method     */
/*         It is not good enough for serious work.                            */
/*         It could be brought up to decent standards by using the included   */
/*         annealing just as a starting point, following it with good         */
/*         refinement like LM.                                                */
/*         Also, a user friendly version would have provision for progress    */
/*         reports and user interruption.  And last but not least, error      */
/*         checks like failure to allocate sufficient memory should be        */
/*         included in any serious implementation.                            */
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
#include "mlfn.h"

double normal () ;

/*
--------------------------------------------------------------------------------

   Constructor, destructor, reset(), add_case()

   The temporary work areas tset (the training set)
   and svd (the SingularValueDecomp) are used only for training.
   (Except add_case() cumulates the training set in tset.)
   At the end of train() they could be deleted to free up system memory.
   But they are not deleted here because the MLFN object will be reused
   many times. So we avoid repeated allocations and deletions,
   which slow the system.

--------------------------------------------------------------------------------
*/

MLFN::MLFN ( int ncase , int nin , int nout , int nhid )
{
   ncases = ncase ;
   ninputs = nin ;
   noutputs = nout ;
   nhidden = nhid ;
   probs = NULL ;
   svd = new SingularValueDecomp ( ncase , nhid+1 ) ;
   tset = (double *) malloc ( ncases * (ninputs + noutputs) * sizeof(double) ) ;
   inwts = (double *) malloc ( nhidden * (ninputs + 1) * sizeof(double) ) ;
   outwts = (double *) malloc ( noutputs * (nhidden + 1) * sizeof(double) ) ;
   hid = (double *) malloc ( nhidden * sizeof(double) ) ;
   reset () ;
}


MLFN::~MLFN ()
{
   if (svd != NULL)    // This and tset could be freed and set to NULL by train()
      delete svd ;     // But they are not here.
   if (tset != NULL)   // Nevertheless, demonstrate checking for this.
      free ( tset ) ;
   if (inwts != NULL)
      free ( inwts ) ;
   if (outwts != NULL)
      free ( outwts ) ;
   if (hid != NULL)
      free ( hid ) ;
   if (probs != NULL)
      free ( probs ) ;
}

/*
   This discards any existing training data.
   It does not need to be called after construction, but it must
   be called if the user wants to reuse the MLFN object for a new dataset.
*/

void MLFN::reset ()
{
   nrows = 0 ;      // No rows (via add_case()) yet present
   trained = 0 ;    // Training not done yet
}

/*
   Build the training set one case at a time.
   The user must call this member EXACTLY ncases times after construction
   or a call to reset(), and before a call to train().
   The first version assumes the usual situation that all cases are to be
   weighted identically.  The second version allows the user to specify
   a probability for relative weighting of the case.
   These versions must NOT be mixed after construction!
   Once an MLFN object has been constructed, all calls to add_case() must
   be one version or the other, even if reset() is called.
*/

void MLFN::add_case ( double *newcase )
{
   if (nrows >= ncases)  // Careful user never lets this happen
      return ;           // But cheap insurance

   memcpy ( tset + nrows * (ninputs + noutputs) , newcase ,
            (ninputs + noutputs) * sizeof(double) ) ;
   ++nrows ;
}

void MLFN::add_case ( double *newcase , double prob )
{
   if (nrows >= ncases)  // Careful user never lets this happen
      return ;           // But cheap insurance

   memcpy ( tset + nrows * (ninputs + noutputs) , newcase ,
            (ninputs + noutputs) * sizeof(double) ) ;

   if (probs == NULL)
      probs = (double *) malloc ( ncases * sizeof(double) ) ;

   probs[nrows] = prob ;
   ++nrows ;
}

/*
--------------------------------------------------------------------------------

   predict() - Given an input vector, compute output using trained model

--------------------------------------------------------------------------------
*/

void MLFN::predict (
   double *input ,     // Input vector
   double *output      // Returned output
   )
{
   int i, j ;
   double sum, *dptr, temp1, temp2 ;

/*
   Compute hidden layer activations
*/

   for (i=0 ; i<nhidden ; i++) {
      dptr = inwts + i * (ninputs + 1) ;  // Weights for this neuron
      sum = dptr[ninputs] ;       // Constant
      for (j=0 ; j<ninputs ; j++)
         sum += dptr[j] * input[j] ;
      if (sum > 150.0)
         sum = 150.0 ;
      if (sum < -150.0)
         sum = -150.0 ;
      temp1 = exp ( sum ) ;
      temp2 = 1.0 / temp1 ;
      hid[i] = (temp1 - temp2) / (temp1 + temp2) ; // Hyperbolic tangent
      }

/*
   Compute output layer activations
*/

   for (i=0 ; i<noutputs ; i++) {
      dptr = outwts + i * (nhidden + 1) ;  // Weights for this neuron
      sum = dptr[nhidden] ;       // Constant
      for (j=0 ; j<nhidden ; j++)
         sum += dptr[j] * hid[j] ;
      output[i] = sum ;
      }
}

/*
--------------------------------------------------------------------------------

   execute() - Given input weights, pass through the training set.
               Compute optimal output weights and return MSE.
               Note that when mean squared error is being optimized,
               there is no need to treat the output weights as optimizable
               parameters.  This is because for any set of input weights,
               the MSE optimal output weights can be computed via singular
               value decomposition.  This remarkable fact is ignored by
               the majority of MLFN training packages!

--------------------------------------------------------------------------------
*/

double MLFN::execute ()
{
   int i, j, iout, icase ;
   double err, sum, diff, *aptr, *bptr, *dptr, *tptr, temp1, temp2 ;
   double prob ;

   aptr = svd->a ;   // Singular value decomposition design matrix

   for (icase=0 ; icase<nrows ; icase++) {
      tptr = tset + icase * (ninputs + noutputs) ;  // This training case
      if (probs == NULL)         // If all cases equally weighted
         prob = 1.0 ;            // Do nothing
      else                       // But if user specifies a distribution
         prob = sqrt ( probs[icase] ) ;   // Weight case per its probability

/*
   Compute hidden layer activations.
   Put them in design matrix, with constant last.
*/

      for (i=0 ; i<nhidden ; i++) {
         dptr = inwts + i * (ninputs + 1) ;  // Weights for this neuron
         sum = dptr[ninputs] ;       // Constant
         for (j=0 ; j<ninputs ; j++)
            sum += dptr[j] * tptr[j] ;
         if (sum > 150.0)
            sum = 150.0 ;
         if (sum < -150.0)
            sum = -150.0 ;
         temp1 = exp ( sum ) ;
         temp2 = 1.0 / temp1 ;
         *aptr++ = prob * (temp1 - temp2) / (temp1 + temp2) ; // Hyperbolic tangent
         } // For computing each hidden activation

      *aptr++ = prob ;   // Constant
      } // For all training cases

   svd->svdcmp () ;

/*
   For each output, solve for optimal output weights.
*/

   for (iout=0 ; iout<noutputs ; iout++) {
      bptr = svd->b ;
      for (icase=0 ; icase<nrows ; icase++) {
         tptr = tset + icase * (ninputs + noutputs) + ninputs ; // Case's outputs
         *bptr++ = tptr[iout] ;         // This output forms RHS
         if (probs != NULL)             // If user specified distribution
            *(bptr-1) *= sqrt ( probs[icase] ) ; // Weight this case
         }

      dptr = outwts + iout * (nhidden + 1) ;  // Weights for this neuron
      svd->backsub ( 1.e-4 , dptr ) ;         // Find those weights.
      } // For each output

/*
   The weights for the outputs are now in place in outwts.
   Pass through the training set, computing the activations
   of the hidden layer, then the activations of the output neurons.
   Compare this attained activation to the desired in the training sample,
   and cumulate the mean squared error.
   With a little work, we could avoid computing the hidden layer activations,
   because we could have instructed the svd constructor to preserve the
   design matrix which is them in the uniform probability case.  But this
   is overkill here.  Just bite the bullet and recompute.
*/

   err = 0.0 ;

   for (icase=0 ; icase<nrows ; icase++) {
      tptr = tset + icase * (ninputs + noutputs) ;  // This training case
      for (i=0 ; i<nhidden ; i++) {
         dptr = inwts + i * (ninputs + 1) ;  // Weights for this neuron
         sum = dptr[ninputs] ;       // Constant
         for (j=0 ; j<ninputs ; j++)
            sum += dptr[j] * tptr[j] ;
         if (sum > 150.0)
            sum = 150.0 ;
         if (sum < -150.0)
            sum = -150.0 ;
         temp1 = exp ( sum ) ;
         temp2 = 1.0 / temp1 ;
         hid[i] = (temp1 - temp2) / (temp1 + temp2) ; // Hyperbolic tangent
         } // For computing each hidden activation

      for (i=0 ; i<noutputs ; i++) {
         dptr = outwts + i * (nhidden + 1) ;  // Weights for this neuron
         sum = dptr[nhidden] ;       // Constant
         for (j=0 ; j<nhidden ; j++)
            sum += dptr[j] * hid[j] ;
         diff = sum - tptr[ninputs+i] ;
         if (probs != NULL)                     // If user specified distribution
            err += diff * diff * probs[icase] ; // Weight this case
         else 
            err += diff * diff ;
         }
      } // For all training cases

   err /= nrows * noutputs ;

   return err ;
}

/*
--------------------------------------------------------------------------------

   train() - Train the model

   After add_case has been called exactly ncases times, this must be called
   to train the model.

   This routine is the weak point in this MLFN class.  The training algorithm
   is relatively slow and inaccurate.
   It is an excellent starting point for refinement, having a high probability
   of finding a solution near a global minimum.  Addition of a good refinement
   algorithm, such as LM, would make this class usable in real applications.

--------------------------------------------------------------------------------
*/

void MLFN::anneal_train (
   int n_outer ,      // Number of outer loop iterations, perhaps 10-20
   int n_inner ,      // Number of inner loop iterations, perhaps 100-10000
   double start_std   // Starting standard deviation of weights, about 2.0
   )
{
   int i, inner, outer ;
   double error, best_error, std, *best_wts, *center ;

/*
   Best_wts keeps track of the best input weights.
   (The corresponding best output weights are computed by SVD.)
   Center is the center around which perturbation is done.
   It starts at zero.  After completion of each pass through the inner loop
   it is changed to best_wts.
*/

   best_wts = (double *) malloc ( nhidden * (ninputs + 1) * sizeof(double) ) ;
   center = (double *) malloc ( nhidden * (ninputs + 1) * sizeof(double) ) ;

   for (i=0 ; i<nhidden*(ninputs+1) ; i++)
      center[i] = 0.0 ;

   best_error = -1.0 ;
   std = start_std ;

   for (outer=0 ; outer<n_outer ; outer++) {
      for (inner=0 ; inner<n_inner ; inner++) {

         for (i=0 ; i<nhidden*(ninputs+1) ; i++)
            inwts[i] = center[i] + std * normal() ;
         error = execute () ;
         if ((best_error < 0.0)  ||  (error < best_error)) {
            best_error = error ;
            memcpy ( best_wts , inwts , nhidden*(ninputs+1) * sizeof(double) ) ;
            }
         } // For inner loop iterations
      memcpy ( center , best_wts , nhidden*(ninputs+1) * sizeof(double) ) ;
      std *= 0.8 ;
      } // For outer loop iterations

   memcpy ( inwts , best_wts , nhidden*(ninputs+1) * sizeof(double) ) ;
   execute () ;     // Computes output weights for predict() use later
   trained = 1 ;    // Training complete
   free ( best_wts ) ;
   free ( center ) ;
#if 0               // Allow reuse
   free ( tset ) ;  
   tset = NULL ;
   delete svd ;
   svd = NULL ;
#endif
}

/*
   This is customized for this demonstration
*/

void MLFN::train ()
{
   anneal_train ( 10 , 100 , 1.0 ) ;
}

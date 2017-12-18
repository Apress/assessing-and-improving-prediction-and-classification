/******************************************************************************/
/*                                                                            */
/*  ARCING_M - Compare bagging and AdaBoost for multiple classification       */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "mlfn.h"
#include "minimize.h"

double unifrand () ;
double normal () ;

/*
--------------------------------------------------------------------------------

   Notes on the component models...

   For clarity and to avoid the need for fancy parameter passing,
   these routines assume that anything needed by the model is a global or
   static in this module.  These things are assumed to have been already
   constructed and ready for use.  In particular, the following model
   characteristics and functions are assumed:
      ---> An array 'models' of nmodels pointers to instances of the prediction
           model is constructed and ready to use.  These instances must be
           totally independent of one another.  Each instance has a single
           output.  A separate model will be used for each class.
      ---> Model::reset() resets the model in preparation for a new training set
      ---> Model::add_case() adds a new case to the model's training set
           This includes both straight and implied importance versions.
      ---> Model::train() trains the model.
      ---> Model::predict() predicts using the model.
           Numeric predictions need not be bounded by Model::predict because
           they will be limited to [-1,1] here.  The model for a given class
           should be trained to predict +1, and all other models (classes)
           trained to predict -1.

--------------------------------------------------------------------------------
*/

static MLFN **model ;      // Created and deleted in main, for actual error
static MLFN **models ;     // Ditto, for arcing classifiers

/*
--------------------------------------------------------------------------------

   Bagging

--------------------------------------------------------------------------------
*/

class Bagging {

public:

   Bagging ( int n , int nin , int nout , double *tset , int nmods ) ;
   ~Bagging () ;
   void numeric_predict ( double *input , double *output ) ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   int nout ;         // Number of outputs
   double *work ;     // Work area for holding output and training case
   int *count ;       // Classification voting counter
} ;

Bagging::Bagging (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nouts ,        // Number of outputs
   double *tset ,     // Training cases, ntrain by (nin+nout)
   int nmods          // Number of bootstrap replications
   )
{
   int i, k, iboot, iout ;
   double *tptr ;

   nout = nouts ;
   nmodels = nmods ;
   work = (double *) malloc ( (nin + nout) * sizeof(double) ) ;
   count = (int *) malloc ( nout * sizeof(int) ) ;

/*
   Build the bootstrap training sets and train each model
*/

   for (iboot=0 ; iboot<nmodels ; iboot++) {

      for (iout=0 ; iout<nout ; iout++)
         models[iboot*nout+iout]->reset() ;  // Prepares the reusable model

      for (i=0 ; i<n ; i++) {           // Build this bootstrap training set
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         tptr = tset + k * (nin + nout) ;    // Point to this case
         memcpy ( work , tptr , nin * sizeof(double) ) ;  // This case's input
         for (iout=0 ; iout<nout ; iout++) { // Each output handled by separate model
            work[nin] = tptr[nin+iout] ;     // This output
            models[iboot*nout+iout]->add_case ( work ) ; // Add it to the model's training set
            }
         }

      for (iout=0 ; iout<nout ; iout++)
         models[iboot*nout+iout]->train() ;  // Train the model
      }
}

Bagging::~Bagging ()
{
   if (work != NULL)
      free ( work ) ;
   if (count != NULL)
      free ( count ) ;
}

void Bagging::numeric_predict ( double *input , double *output )
{
   int i, imodel ;
   double out ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      for (i=0 ; i<nout ; i++) {
         models[imodel*nout+i]->predict ( input , &out ) ;
         if (out > 1.0)      // This hard limiting aids stability
            out = 1.0 ;      // Design your model so that [-1,1] is its 
         if (out < -1.0)     // natural range, so this has minimal impact
            out = -1.0 ;
         output[i] += out ;
         }
      }

   for (i=0 ; i<nout ; i++)
      output[i] /= nmodels ;
}

int Bagging::class_predict ( double *input )
{
   int i, imodel, ibest, bestcount ;
   double best, out ;

   for (i=0 ; i<nout ; i++)
      count[i] = 0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      for (i=0 ; i<nout ; i++) {
         models[imodel*nout+i]->predict ( input , &out ) ;
         if ((i == 0)  ||  (out > best)) {
            best = out ;
            ibest = i ;
            }
         }
      ++count[ibest] ;   // Count winners
      }

/*
   At this time, count[i] is the number of times output i was the greatest
*/

   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (count[i] > bestcount)) {
         bestcount = count[i] ;
         ibest = i ;
         }
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   AdaBoostMH - AdaBoost algorithm for many-class classification in which
      the underlying model provides a class plus numerical confidence.
      Each class has its own underlying model with output dedicated to this class.
      These outputs must have been trained to attempt to produce +1 for the
      output corresponding to the correct class, and -1 for all other outputs.
      For compatibility with this model as well as the published literature,
      each output value is hard limited here to +/-1.  It is good if this is
      the natural range of your model, so that the hard limiting has little
      or no impact.  This limiting aids stability when the underlying model
      has a chance of occasionally producing wild output values.
      Also, having a known natural range helps keep the search for an
      optimal alpha efficient.

--------------------------------------------------------------------------------
*/

class AdaBoostMH {

public:

   AdaBoostMH ( int n , int nin , int nout , double *tset , int nmods ) ;
   ~AdaBoostMH () ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   int nout ;         // Number of outputs (classes)
   double *work ;     // Work area for holding output and training case
   double *alpha ;    // Nmods long alpha constant for each model
   double *dist ;     // N by nout long probability distribution
   double *u ;        // N by nout long work area for saving model's error products
} ;

/*
   This local routine is the criterion function that is passed to the
   minimization routines that compute the optimal alpha.
*/

static int local_n, local_nout ;
static double *local_dist, *local_u ;

static double alpha_crit ( double trial_alpha )
{
   int i, k ;
   double sum ;

   k = local_n * local_nout ;
   sum = 0.0 ;
   for (i=0 ; i<k ; i++)
      sum += local_dist[i] * exp ( -trial_alpha * local_u[i] ) ;
   return sum ;
}

/*
   Constructor
*/

AdaBoostMH::AdaBoostMH (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nouts ,        // Number of outputs
   double *tset ,     // Training cases, ntrain by (nin+nout)
   int nmods          // Number of models
   )
{
   int i, iout, imodel, ngood, nbad, n_nout ;
   double *tptr, *dptr, temp, sum, h, uu ;
   double x1, y1, x2, y2, x3, y3 ;

   nout = nouts ;
   nmodels = nmods ;
   n_nout = n * nout ;

   work = (double *) malloc ( (nin + nout) * sizeof(double) ) ;
   alpha = (double *) malloc ( nmodels * sizeof(double) ) ;
   dist = (double *) malloc ( n_nout * sizeof(double) ) ;
   u = (double *) malloc ( n_nout * sizeof(double) ) ;

/*
   Initialize distribution to be uniform
*/

   temp = 1.0 / n_nout ;
   for (i=0 ; i<n_nout ; i++)
      dist[i] = temp ;

/*
   Main training loop trains sequence of models
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {

      for (i=0 ; i<nout ; i++)
         models[imodel*nout+i]->reset() ;     // Prepares the reusable model

      for (i=0 ; i<n ; i++) {                 // Build this training set
         tptr = tset + i * (nin + nout) ;     // Point to this case
         dptr = dist + i * nout ;             // Its output probability weights
         memcpy ( work , tptr , nin * sizeof(double) ) ;  // This case's input
         for (iout=0 ; iout<nout ; iout++) {  // Each output handled by separate model
            work[nin] = tptr[nin+iout] ;      // This output
            models[imodel*nout+iout]->add_case ( work , dptr[iout] ) ;
            }
         }

      for (iout=0 ; iout<nout ; iout++)
         models[imodel*nout+iout]->train() ;  // Train the model

/*
   Compute the optimal alpha.
   In the standard algorithm, the sign of an output 'h' is the class prediction
   (+ for this class, - for all others) and the magnitude is the model's
   confidence in the prediction.  Theoretically, h need not be bounded because
   alpha can always compensate for the scaling of h.
   However, things seem to go best when the model is trained to produce
   predictions (h) in a fixed known range such as the [-1,1] used here.
   Among other things, this lets the search for alpha be made efficiently.

   The first step is to compute and save u, the error array.
*/

      ngood = nbad = 0 ;                // Degenerate if all cases good or bad
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + nout) ;  // Point to this case
         for (iout=0 ; iout<nout ; iout++) {
            models[imodel*nout+iout]->predict ( tptr , &h ) ;
            if (h > 1.0)       // This hard limiting aids stability
               h = 1.0 ;       // Design your model so that [-1,1] is its 
            if (h < -1.0)      // natural range, so this has minimal impact
               h = -1.0 ;
            uu = h * tptr[nin+iout] ; // Error is predicted times true
            u[i*nout+iout] = uu ;    // Save it for optimization of alpha
            if (uu > 0.0)            // Class prediction is correct
               ++ngood ;             // This lets us detect degenerate situation
            if (uu < 0.0)            // Class prediction is incorrect
               ++nbad ;              // If degenerate, optimal alpha does not exist
            }
         }

      if (nbad == 0) {           // Unusual situation of model never failing
         nmodels = imodel + 1 ;  // No sense going any further
         alpha[imodel] = 0.5 * log ( (double) n ) ; // Heuristic big value
         break ;                 // So exit training loop early
         }

      if (ngood == 0) {          // Unusual situation of model being worthless
         nmodels = imodel ;      // No sense going any further
         break ;                 // So exit training loop early
         }

      local_n = n ;
      local_nout = nout ;
      local_dist = dist ;
      local_u = u ;

      glob_min ( -1.0 , 1.0 , 3 , 0 , 0.0 , alpha_crit , &x1 , &y1 ,
                 &x2 , &y2 , &x3 , &y3 ) ;

      brentmin ( 20 , 0.0 , 1.e-6 , 1.e-4 , alpha_crit , &x1 , &x2 , &x3 , y2 ) ;

      alpha[imodel] = x2 ;

/*
   Adjust the relative weighting, then rescale to make it a distribution
*/

      sum = 0.0 ;
      for (i=0 ; i<n_nout ; i++) {
         dist[i] *= exp ( -alpha[imodel] * u[i] ) ;
         sum += dist[i] ;        // Cumulate total probability for rescale
         }

      for (i=0 ; i<n_nout ; i++) // A probability distribution must sum to 1.0
         dist[i] /= sum ;        // Make it so

      } // For all models

   free ( dist ) ;
   free ( u ) ;
}

AdaBoostMH::~AdaBoostMH ()
{
   if (alpha != NULL)
      free ( alpha ) ;
   if (work != NULL)
      free ( work ) ;
}

/*
   Make a class prediction.
   It returns 0 for the first class, 1 for the second class, et cetera.
*/

int AdaBoostMH::class_predict ( double *input )
{
   int iout, imodel, ibest ;
   double sum, best, h ;

   if (nmodels == 0)   // Abnormal condition of no decent models
      return -1 ;      // Return an error flag

   for (iout=0 ; iout<nout ; iout++) {
      sum = 0.0 ;
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel*nout+iout]->predict ( input , &h ) ;
         if (h > 1.0)       // This hard limiting aids stability
            h = 1.0 ;       // Design your model so that [-1,1] is its 
         if (h < -1.0)      // natural range, so this has minimal impact
            h = -1.0 ;
         sum += alpha[imodel] * h ; // Weighted sum of models' predictions for iout
         }
      if ((iout == 0)  ||  (sum > best)) {
         best = sum ;
         ibest = iout ;
         }
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   AdaBoostOC - AdaBoost algorithm for many-class classification in which
      the underlying model provides a class only; no numerical confidence.
      These models will each be strictly binary, being trained to produce
      +1 for one state and -1 for the other state.
      Only the sign of the underlying model's prediction will be used.
      The magnitude is ignored.

      This algorithm uses output coding for classification.

      Note that the original work by Robert Schapire uses {0,1} for
      coloring and training.  This implementation uses {-1,1} because
      many models train better with symmetric output possibilities.
      The mathematics is identical.  The only difference is in the
      symbolic coding.

--------------------------------------------------------------------------------
*/

class AdaBoostOC {

public:

   AdaBoostOC ( int ncases , int nins , int nouts , double *tset , int nmods ) ;
   ~AdaBoostOC () ;
   int class_predict ( double *input ) ;

private:
   double colorcrit ( int *coloring ) ;// Compute coloring criterion
   int n ;              // Number of training cases
   int nmodels ;        // Number of models (nmods in constructor call)
   int nin ;            // Number of inputs
   int nout ;           // Number of outputs (classes)
   double *work ;       // Work area for holding output and training case
   double *alpha ;      // Nmods long alpha constant for each model
   double *dist ;       // N long case probability distribution
   double *err_dist ;   // N by nout long error probability distribution
   int *colormap ;      // Nclasses long coloring map
   int *bestcolor ;     // Nclasses long best coloring map
   int *tclass ;        // N long class membership of training cases
   int *h ;             // N long predictions for training set
   double *w ;          // Nout by nout 'w' matrix
} ;

/*
   Constructor
*/

AdaBoostOC::AdaBoostOC (
   int ncases ,       // Number of training cases
   int nins ,         // Number of inputs
   int nouts ,        // Number of outputs
   double *tset ,     // Training cases, ntrain by (nin+nout)
   int nmods          // Number of models
   )
{
   int i, k, kk, iout, jout, imodel, *h, icolor, *coloring ;
   double *tptr, *dptr, temp, best, sum, out, denom ;

   n = ncases ;
   nin = nins ;
   nout = nouts ;
   nmodels = nmods ;

   work = (double *) malloc ( (nin + nout) * sizeof(double) ) ;
   alpha = (double *) malloc ( nmodels * sizeof(double) ) ;
   dist = (double *) malloc ( n * sizeof(double) ) ;
   err_dist = (double *) malloc ( n * nout * sizeof(double) ) ;
   colormap = (int *) malloc ( nmodels * nout * sizeof(int) ) ;
   bestcolor = (int *) malloc ( nout * sizeof(int) ) ;
   tclass = (int *) malloc ( n * sizeof(int) ) ;
   w = (double *) malloc ( nout * nout * sizeof(double) ) ;
   h = (int *) malloc ( n * sizeof(int) ) ;

/*
   Compute and save class of each training case.
   This speeds repetitive operations later.
*/

   for (i=0 ; i<n ; i++) {                 // For entire training set
      tptr = tset + i * (nin + nout) ;     // Point to this case
      for (iout=0 ; iout<nout ; iout++) {
         if (tptr[nin+iout] > 0.0) {       // +1 for this class, -1 for others
            tclass[i] = iout ;
            break ;
            }
         }
      }

/*
   Initialize error distribution to be uniform over errors
*/

   temp = 1.0 / (n * (nout - 1)) ;
   for (i=0 ; i<n ; i++) {                 // For entire training set
      tptr = tset + i * (nin + nout) ;     // Point to this case
      dptr = err_dist + i * nout ;         // Its error probability weights
      for (iout=0 ; iout<nout ; iout++)
         dptr[iout] = (tptr[nin+iout] > 0.0)  ?  0.0  :  temp ;
      }

/*
   Main training loop trains sequence of models
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      coloring = colormap + imodel * nout ;  // Each model has its own coloring

/*
   Compute w
*/

      for (iout=0 ; iout<nout ; iout++) {
         for (jout=0 ; jout<nout ; jout++) {
            sum = 0.0 ;
            for (i=0 ; i<n ; i++) {
               if (jout == tclass[i])
                  sum += err_dist[i * nout + iout] ;
               }
            w[iout*nout+jout] = sum ;
            }
         }

/*
   Find the optimal coloring by testing ALL possibilities.
   This is VERY time consuming if nout is large!
   If this routine may be used for more than a half dozen or so classes,
   use better optimization.  Actually, selecting the best from among
   a modest number of random choices will do quite well in practice.

   This algorithm breaks out of the main loop the moment the highest bit
   is flipped.  This is because all subsequent colorings would represent
   partitions already tested, simply having signs reversed.
*/

      for (i=0 ; i<nout ; i++)
         coloring[i] = -1 ;

      best = -1.e80 ;
      for (;;) {
         k = nout ;               // Start flipping with bit nout-1
         while (k--) {            // Count bits
            if (coloring[k] > 0)  // If this bit is 1
               coloring[k] = -1 ; // Flip it to -1 and continue to next bit
            else {                // If this bit is -1
               coloring[k] = 1 ;  // Flip it to 1
               break ;            // And do not go to next bit
               }
            } // While flipping bits for next trial coloring
         if (k <= 0)              // Happens when all combinations tested
            break ;               // (Keeping highest bit at -1)
         temp = colorcrit ( coloring ) ; // Test this coloring
         if (temp > best) {       // If it is the best so far
            best = temp ;         // Update best and save this coloring
            memcpy ( bestcolor , coloring , nout * sizeof(int) ) ;
            }
         } // For all possible colorings

      memcpy ( coloring , bestcolor , nout * sizeof(int) ) ; // Get best

/*
   The optimal coloring has been found.
   Compute the casewise probability distribution.
*/

   denom = 0.0 ;
   for (i=0 ; i<n ; i++) {
      dptr = err_dist + i * nout ;         // Its error probability weights
      icolor = coloring [ tclass[i] ] ;
      sum = 0.0 ;
      for (iout=0 ; iout<nout ; iout++) {
         if (coloring[iout] != icolor)
            sum += dptr[iout] ;
         }
      dist[i] = sum ;
      denom += sum ;
      }

   for (i=0 ; i<n ; i++)
      dist[i] /= denom ;

/*
   Build the training set and train the model.
*/

      models[imodel]->reset() ;               // Prepares the reusable model

      for (i=0 ; i<n ; i++) {                 // Build this training set
         tptr = tset + i * (nin + nout) ;     // Point to this case
         memcpy ( work , tptr , nin * sizeof(double) ) ;  // This case's input
         work[nin] = coloring[tclass[i]] ;    // Target output
         models[imodel]->add_case ( work , dist[i] ) ;
         }

      models[imodel]->train() ;  // Train the model

/*
   Compute the error and alpha
*/

      sum = 0.0 ;
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + nout) ;  // Point to this case
         dptr = err_dist + i * nout ;      // Its error probability weights
         models[imodel]->predict ( tptr , &out ) ; // Make its predictions
         if (out > 0.0)         // Only the sign of the prediction is used
            h[i] = 1 ;          // To determine the binary decision
         else 
            h[i] = -1 ;
         kk = (coloring[tclass[i]] != h[i])  ?  1  :  0 ;
         for (iout=0 ; iout<nout ; iout++) {
            k = kk ;
            if (coloring[iout] == h[i])
               ++k ;
            sum += k * dptr[iout] ;
            }
         }

      sum *= 0.5 ;

      if (sum < 1.e-12)         // Prevent undefined fpt operations below
         sum = 1.e-12 ;

      if (sum > 1.0 - 1.e-12)
         sum = 1.0 - 1.e-12 ;

      alpha[imodel] = 0.5 * log ( (1.0 - sum) / sum ) ;

/*
   Update the error probability distribution
*/

      denom = 0.0 ;
      for (i=0 ; i<n ; i++) {
         dptr = err_dist + i * nout ;      // Its error probability weights
         kk = (coloring[tclass[i]] != h[i])  ?  1  :  0 ;
         for (iout=0 ; iout<nout ; iout++) {
            k = kk ;
            if (coloring[iout] == h[i])
               ++k ;
            dptr[iout] *= exp ( alpha[imodel] * k ) ;
            denom += dptr[iout] ;
            }
         }

      for (i=0 ; i<n ; i++) {          // Normalize to sum to 1.0
         dptr = err_dist + i * nout ;
         for (iout=0 ; iout<nout ; iout++)
            dptr[iout] /= denom ;
         }

      } // For all models

   free ( dist ) ;
   free ( err_dist ) ;
   free ( w ) ;
   free ( h ) ;
   free ( tclass ) ;
   free ( bestcolor ) ;
}

AdaBoostOC::~AdaBoostOC ()
{
   if (alpha != NULL)
      free ( alpha ) ;
   if (work != NULL)
      free ( work ) ;
   if (colormap != NULL)
      free ( colormap ) ;
}

/*
   Compute coloring criterion that is maximized
*/

double AdaBoostOC::colorcrit ( int *coloring )
{
   int i, j ;
   double sum ;

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      for (j=0 ; j<nout ; j++) {
         if (coloring[i] != coloring[j])
            sum += w[i*nout+j] ;
         }
      }

   return sum ;
}

/*
   Make a class prediction.
   It returns 0 for the first class, 1 for the second class, et cetera.
*/

int AdaBoostOC::class_predict ( double *input )
{
   int iout, imodel, ibest, hh, *coloring ;
   double best, out ;

   if (nmodels == 0)   // Abnormal condition of no decent models
      return -1 ;      // Return an error flag

   for (iout=0 ; iout<nout ; iout++)
      work[iout] = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      coloring = colormap + imodel * nout ;  // Each model has its own coloring
      models[imodel]->predict ( input , &out ) ;
      if (out > 0.0)   // Only the sign of the prediction is used
         hh = 1 ;      // The magnitude is ignored
      else 
         hh = -1 ;
      for (iout=0 ; iout<nout ; iout++) {
         if (coloring[iout] == hh)
            work[iout] += alpha[imodel] ;
         }
      }

   for (iout=0 ; iout<nout ; iout++) {
      if ((iout == 0)  ||  (work[iout] > best)) {
         best = work[iout] ;
         ibest = iout ;
         }
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   Optional main to test it

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, j, k, ntries, itry, nsamps, nclasses, nmodels, ndone, ibest ;
   double *x, *test, *work, *tptr, separation, out, diff, temp, temp2, best ;
   double sum_numeric_error, sum_class_error, sum_train_error, train_error ;
   double bagging_numeric_error, bagging_class_error, bagging_train_error ;
   double adaboost_mh_class_error ;
   double adaboost_mh_train_error ;
   double adaboost_oc_class_error ;
   double adaboost_oc_train_error ;
   Bagging *bagging ;
   AdaBoostMH *adaboost_mh ;
   AdaBoostOC *adaboost_oc ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 6) {
      printf (
         "\nUsage: ARCING_M  nsamples  nclasses  nmodels  ntries  separation" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nclasses = atoi ( argv[2] ) ;
   nmodels = atoi ( argv[3] ) ;
   ntries = atoi ( argv[4] ) ;
   separation = atof ( argv[5] ) ;
#else
   nsamps = 50 ;
   nclasses = 5 ;
   nmodels = 10 ;
   ntries = 100 ;
   separation = 0.7 ;
#endif

   if ((nsamps <= 0)  ||  (nclasses <= 1)  ||  (nmodels <= 0)
    || (ntries <= 0)  ||  (separation < 0.0)) {
      printf (
         "\nUsage: ARCING_M  nsamples  nclasses  nmodels  ntries  separation" ) ;
      exit ( 1 ) ;
      }

   separation *= 4.0 ;  // For compatibility with ARCING.CPP in binary case

/*
   Allocate memory and initialize
*/

   model = (MLFN **) malloc ( nclasses * sizeof(MLFN *) ) ;
   for (i=0 ; i<nclasses ; i++)
      model[i] = new MLFN ( nsamps , 2 , 1 , 2 ) ;

   models = (MLFN **) malloc ( nmodels * nclasses * sizeof(MLFN *) ) ;
   for (i=0 ; i<nmodels*nclasses ; i++)
      models[i] = new MLFN ( nsamps , 2 , 1 , 2 ) ;

   work = (double *) malloc ( (nclasses+2) * sizeof(double) ) ;
   x = (double *) malloc ( nsamps * (nclasses+2) * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * (nclasses+2) * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   sum_numeric_error = 0.0 ;   // For comparison purposes, real error
   sum_class_error = 0.0 ;     // For comparison purposes, class error
   sum_train_error = 0.0 ;     // For comparison purposes, training error

   bagging_numeric_error = 0.0 ;
   bagging_class_error = 0.0 ;
   bagging_train_error = 0.0 ;
   adaboost_mh_class_error = 0.0 ;
   adaboost_mh_train_error = 0.0 ;
   adaboost_oc_class_error = 0.0 ;
   adaboost_oc_train_error = 0.0 ;

   for (itry=0 ; itry<ntries ; itry++) {
      ndone = itry + 1 ;

/*
   Generate the data.
   It is bivariate clusters with moderate positive correlation.
   One class is shifted above and to the left of the other class.
   We use x as the dataset for all resampling algorithms.
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[(nclasses+2)*i] = normal () ;
         x[(nclasses+2)*i+1] = .7071 * x[(nclasses+2)*i]  +  .7071 * normal () ;
         k = (int) (unifrand() * nclasses) ;  // Pick a class for this case
         if (k >= nclasses)                   // Should never happen
            k = nclasses - 1 ;                // But cheap insurance
         x[(nclasses+2)*i] -= k * separation - 0.5 * (k % 2) * separation ;
         x[(nclasses+2)*i+1] += k * separation - 0.5 * (k % 2) * separation ;
         for (j=0 ; j<nclasses ; j++)
            x[(nclasses+2)*i+2+j] = (j == k)  ?  1.0  :  -1.0 ;
         }

      for (i=0 ; i<10*nsamps ; i++) {
         test[(nclasses+2)*i] = normal () ;
         test[(nclasses+2)*i+1] = .7071 * test[(nclasses+2)*i]  +  .7071 * normal () ;
         k = (int) (unifrand() * nclasses) ;  // Pick a class for this case
         if (k >= nclasses)                   // Should never happen
            k = nclasses - 1 ;                // But cheap insurance
         test[(nclasses+2)*i] -= k * separation - 0.5 * (k % 2) * separation ;
         test[(nclasses+2)*i+1] += k * separation - 0.5 * (k % 2) * separation ;
         for (j=0 ; j<nclasses ; j++)
            test[(nclasses+2)*i+2+j] = (j == k)  ?  1.0  :  -1.0 ;
         }

/*
   Train a model with this data and test it on an independent test set.
   This gives us a basis of comparison for the resampling methods.
   The numeric error is the ordinary mean squared error.
*/

      for (j=0 ; j<nclasses ; j++)     // Each output handled by separate model
         model[j]->reset () ;

      for (i=0 ; i<nsamps ; i++) {
         tptr = x + i * (nclasses+2) ;    // Point to this case
         memcpy ( work , tptr , 2 * sizeof(double) ) ;  // This case's input
         for (j=0 ; j<nclasses ; j++) {   // Each output handled by separate model
            work[2] = tptr[2+j] ;         // This output
            model[j]->add_case ( work ) ; // Add it to the model's training set
            }
         }

      for (j=0 ; j<nclasses ; j++)     // Each output handled by separate model
         model[j]->train() ;           // Train the model

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         tptr = x + i * (nclasses+2) ;    // Point to this case
         for (j=0 ; j<nclasses ; j++) {
            model[j]->predict ( tptr , &out ) ;
            if ((j == 0)  ||  (out > best)) {
               best = out ;
               ibest = j ;
               }
            }
         if (tptr[2+ibest] < 0.0)
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp = temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         tptr = test + i * (nclasses+2) ;    // Point to this case
         for (j=0 ; j<nclasses ; j++) {
            model[j]->predict ( tptr , &out ) ;
            if ((j == 0)  ||  (out > best)) {
               best = out ;
               ibest = j ;
               }
            if (out > 1.0)
               out = 1.0 ;
            if (out < -1.0)
               out = -1.0 ;
            diff = out - tptr[j+2] ;
            temp += diff * diff ;
            }
         if (tptr[2+ibest] < 0.0)
            temp2 += 1.0 ;
         }

      sum_train_error += train_error ;
      sum_numeric_error += temp / (10 * nsamps * nclasses) ;
      sum_class_error += temp2 / (10 * nsamps) ;

      printf ( "\n\n\nDid%5d    Observed error: Numeric = %8.4lf  Class =%7.4lf (%7.4lf)",
               ndone, sum_numeric_error / ndone, sum_class_error / ndone ,
               sum_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

/*
   Handle bagging
*/

      bagging = new Bagging ( nsamps , 2 , nclasses , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         tptr = x + i * (nclasses+2) ;    // Point to this case
         k = bagging->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp = temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         tptr = test + i * (nclasses+2) ;    // Point to this case
         bagging->numeric_predict ( tptr , work ) ;
         for (j=0 ; j<nclasses ; j++) {
            diff = work[j] - tptr[j+2] ;
            temp += diff * diff ;
            }
         k = bagging->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            temp2 += 1.0 ;
         }

      bagging_numeric_error += temp / (10 * nsamps * nclasses) ;
      bagging_class_error += temp2 / (10 * nsamps) ;
      bagging_train_error += train_error ;
      delete bagging ;

      printf ( "\n             Bagging error: Numeric = %8.4lf  Class =%7.4lf (%7.4lf)",
               bagging_numeric_error / ndone, bagging_class_error / ndone ,
               bagging_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

/*
   Handle adaboost_mh
*/

      adaboost_mh = new AdaBoostMH ( nsamps , 2 , nclasses , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         tptr = x + i * (nclasses+2) ;    // Point to this case
         k = adaboost_mh->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         tptr = test + i * (nclasses+2) ;    // Point to this case
         k = adaboost_mh->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            temp2 += 1.0 ;
         }

      adaboost_mh_class_error += temp2 / (10 * nsamps) ;
      adaboost_mh_train_error += train_error ;
      delete adaboost_mh ;

      printf ( "\nAdaBoostMH error:                               Class =%7.4lf (%7.4lf)",
               adaboost_mh_class_error / ndone ,
               adaboost_mh_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

/*
   Handle adaboost_oc
*/

      adaboost_oc = new AdaBoostOC ( nsamps , 2 , nclasses , x , nmodels * nclasses ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         tptr = x + i * (nclasses+2) ;    // Point to this case
         k = adaboost_oc->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         tptr = test + i * (nclasses+2) ;    // Point to this case
         k = adaboost_oc->class_predict ( tptr ) ;
         if (tptr[k+2] < 0.0)
            temp2 += 1.0 ;
         }

      adaboost_oc_class_error += temp2 / (10 * nsamps) ;
      adaboost_oc_train_error += train_error ;
      delete adaboost_oc ;

      printf ( "\nAdaBoostOC error:                               Class =%7.4lf (%7.4lf)",
               adaboost_oc_class_error / ndone ,
               adaboost_oc_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

     } // For all tries


   for (i=0 ; i<nclasses ; i++)
      delete model[i] ;
   free ( model ) ;
   for (i=0 ; i<nmodels*nclasses ; i++)
      delete models[i] ;
   free ( models ) ;
   free ( x ) ;
   free ( work ) ;
   free ( test ) ;

   return EXIT_SUCCESS ;
}

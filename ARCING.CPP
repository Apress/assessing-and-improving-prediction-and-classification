/******************************************************************************/
/*                                                                            */
/*  ARCING - Compare bagging and AdaBoost methods for binary classification   */
/*           All routines here are strictly binary.                           */
/*           Multiple class algorithms are in ARCING_M.CPP                    */
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
           totally independent of one another.
      ---> The training set must use a value of +1 to signify the first class,
           and -1 for the second class.
      ---> Model::reset() resets the model in preparation for a new training set
      ---> Model::add_case() adds a new case to the model's training set
           This includes both straight and implied importance versions.
      ---> Model::train() trains the model.
      ---> Model::predict() predicts using the model.
           This includes both numeric and classification versions.
           Numeric predictions need not be bounded by Model::predict because
           they will be hard limited to [-1,1] here.  It is strongly suggested
           that predictions have a natural domain of [-1,1] so that the
           limiting here will have minimal impact.

--------------------------------------------------------------------------------
*/

static MLFN *model ;      // Created and deleted in main, for actual error
static MLFN **models ;    // Ditto, for arcing classifiers

/*
--------------------------------------------------------------------------------

   Bagging

--------------------------------------------------------------------------------
*/

class Bagging {

public:

   Bagging ( int n , int nin , double *tset , int nmods ) ;
   ~Bagging () ;
   void numeric_predict ( double *input , double *output ) ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
} ;

Bagging::Bagging (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout)
   int nmods          // Number of models (bootstrap replications here)
   )
{
   int i, k, iboot ;
   double *tptr ;

   nmodels = nmods ;

/*
   Build the bootstrap training sets and train each model
*/

   for (iboot=0 ; iboot<nmodels ; iboot++) {

      models[iboot]->reset() ;          // Prepares the reusable model

      for (i=0 ; i<n ; i++) {           // Build this bootstrap training set
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         tptr = tset + k * (nin + 1) ;  // Point to this case
         models[iboot]->add_case ( tptr ) ; // Add it to the model's training set
         }

      models[iboot]->train () ; // Train this model
      }
}

Bagging::~Bagging ()
{
}

void Bagging::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double outwork ;

   *output = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &outwork ) ;
      if (outwork > 1.0)   // Impose hard limiting for stability
         outwork = 1.0 ;   // For a well designed model this will have
      if (outwork < -1.0)  // Little or no impact
         outwork = -1.0 ;
      *output += outwork ; // Bagging output is mean across all models
      }

   *output /= nmodels ;
}

int Bagging::class_predict ( double *input )
{
   int imodel, count0, count1 ;
   double outwork ;

   count0 = count1 = 0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &outwork ) ;
      if (outwork > 0.0)        // If the model predicts the first class
         ++count0 ;             // Tally this vote
      else if (outwork < 0.0)   // Other class
         ++count1 ;
      }

   if (count0 > count1)         // If the first class got the most votes
      return 0 ;                // Return this result
   return 1 ;                   // No provision for ties in this version
}

/*
--------------------------------------------------------------------------------

   AdaBoostBinaryNoConf - Simplest AdaBoost algorithm for strictly binary
      classification in which the underlying model provides only a class
      prediction, with no accompanying numerical confidence.
      This implementation explicitly discards the numeric value of the
      model's prediction, using only the sign as the class prediction.

--------------------------------------------------------------------------------
*/

class AdaBoostBinaryNoConf {

public:

   AdaBoostBinaryNoConf ( int n , int nin , double *tset , int nmods ) ;
   ~AdaBoostBinaryNoConf () ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *alpha ;    // Nmods long alpha constant for each model
   double *dist ;     // N long probability distribution
   double *h ;        // N long work area for saving model's predictions
} ;

AdaBoostBinaryNoConf::AdaBoostBinaryNoConf (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout, where nout=1)
   int nmods          // Number of models
   )
{
   int i, imodel ;
   double *tptr, temp, eps, beta, out ;

   nmodels = nmods ;
   alpha = (double *) malloc ( nmodels * sizeof(double) ) ;
   dist = (double *) malloc ( n * sizeof(double) ) ;
   h = (double *) malloc ( n * sizeof(double) ) ;

/*
   Initialize distribution to be uniform
*/

   temp = 1.0 / n ;
   for (i=0 ; i<n ; i++)
      dist[i] = temp ;

/*
   Main training loop trains sequence of models
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {

      models[imodel]->reset() ;         // Prepares the reusable model

      for (i=0 ; i<n ; i++) {           // Build this training set
         tptr = tset + i * (nin + 1) ;  // Point to this case
         models[imodel]->add_case ( tptr , dist[i] ) ; // Add it to the model's training set
         }

      models[imodel]->train () ; // Train this model

/*
   Compute eps as the probability of error.  Use this to compute optimal alpha.
   The 'if' statement in this loop must be true when the actual class (as
   stored in the training set) is not the predicted class (as defined by the
   model's prediction).  This particular implementation uses a model having
   positive prediction for the first class, negative for the second.
*/

      eps = 0.0 ;
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + 1) ;  // Point to this case
         models[imodel]->predict ( tptr , &out ) ;
         if (out > 0.0)                 // If it predicts first class
            h[i] = 1.0 ;                // Flag it this way
         else                           // But if it predicts second class
            h[i] = -1.0 ;               // Other flag value
         if ((tptr[nin] * h[i]) < 0.0)  // Error test: true & predicted opp sign?
            eps += dist[i] ;            // Model erred if signs opposite
         }

      if (eps <= 0.0) {          // Unusual situation of model being perfect
         nmodels = imodel + 1 ;  // No sense going any further
         alpha[imodel] = 0.5 * log ( (double) n ) ; // Arbitrary large value
         break ;                 // So exit training loop early
         }

      if (eps > 0.5) {           // Unusual situation of model being worthless
         nmodels = imodel ;      // No sense going any further
         break ;                 // So exit training loop early
         }

      alpha[imodel] = 0.5 * log ( (1.0 - eps) / eps ) ;
      beta = eps / (1.0 - eps) ;

/*
   Adjust the relative weighting, then rescale to make it a distribution
*/

      temp = 0.0 ;
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + 1) ;  // Point to this case
         if ((tptr[nin] * h[i]) > 0.0)  // If this case correctly classified
            dist[i] *= beta ;   // Reduce its probability weighting
         temp += dist[i] ;      // Cumulate total probability for rescale
         }

      for (i=0 ; i<n ; i++)     // A probability distribution must sum to 1.0
         dist[i] /= temp ;      // Make it so
#if DEBUG
      printf ( "  k=%lf", temp ) ;
#endif

      } // For all models

   free ( dist ) ;
   free ( h ) ;
}

AdaBoostBinaryNoConf::~AdaBoostBinaryNoConf ()
{
   if (alpha != NULL)
      free ( alpha ) ;
}

/*
   Make a class prediction.
   For compatibility with the other models, it returns 0 for the first class
   and 1 for the second class.
   Recall that by definition, AdaBoostBinaryNoConf ignores the
   prediction's magnitude, being concerned only with the sign.
*/

int AdaBoostBinaryNoConf::class_predict ( double *input )
{
   int i ;
   double sum, out ;

   if (nmodels == 0)   // Abnormal condition of no decent models
      return -1 ;      // Return an error flag

   sum = 0.0 ;
   for (i=0 ; i<nmodels ; i++) {
      models[i]->predict ( input , &out ) ;
      if (out > 0.0)          // If it predicts first class
         sum += alpha[i] ;
      else if (out < 0.0)     // But if it predicts second class
         sum -= alpha[i] ;
      }

   if (sum > 0.0)
      return 0 ;    // Flags first class
   else 
      return 1 ;    // Flags second class
}

/*
--------------------------------------------------------------------------------

   AdaBoostBinaryNoConfSampled - AdaBoostBinaryNoConf modified to build
      training set by probability sampling; for models that cannot be trained
      using a probability distribution.
      This implementation uses a resolution factor of 5 for computing idist,
      the subscript table.  Higher factors result in fewer search comparisons
      later but require more memory.  This value of 5 is already very high,
      and any more would be overkill for any but the most critical apps.

--------------------------------------------------------------------------------
*/

class AdaBoostBinaryNoConfSampled {

public:

   AdaBoostBinaryNoConfSampled ( int n , int nin , double *tset , int nmods ) ;
   ~AdaBoostBinaryNoConfSampled () ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *alpha ;    // Nmods long alpha constant for each model
   double *dist ;     // N long probability distribution
   double *cdf ;      // N long work area for cumulative distribution function
   int *idist ;       // 5N long work area for distribution sampling subscripts
   double *h ;        // N long work area for saving model's predictions
} ;

AdaBoostBinaryNoConfSampled::AdaBoostBinaryNoConfSampled (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout, where nout=2)
   int nmods          // Number of models
   )
{
   int i, j, imodel, m ;
   double *tptr, temp, eps, beta, out ;

   nmodels = nmods ;
   m = 5 * n ;       // Resolution factor = 5 ;

   alpha = (double *) malloc ( nmodels * sizeof(double) ) ;
   dist = (double *) malloc ( n * sizeof(double) ) ;
   cdf = (double *) malloc ( n * sizeof(double) ) ;
   h = (double *) malloc ( n * sizeof(double) ) ;
   idist = (int *) malloc ( m * sizeof(int) ) ;

/*
   Initialize distribution to be uniform
*/

   temp = 1.0 / n ;
   for (i=0 ; i<n ; i++)
      dist[i] = temp ;

/*
   Main training loop trains sequence of models
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {

      models[imodel]->reset() ;         // Prepares the reusable model

/*
   Build the table for sampling from tset with distribution 'dist'
*/

      cdf[0] = dist[0] ;       // Compute cumulated distribution function
      for (i=1 ; i<n ; i++)
         cdf[i] = cdf[i-1] + dist[i] ;
      cdf[n-1] = 1.0 + 1.e-8 ; // Avoid fpt roundoff error causing overrun later

      j = -1 ;
      temp = 0.0 ;
      for (i=0 ; i<m ; i++) {
         while (temp <= (double) i) {
            ++j ;
            temp = m * cdf[j] ;
            }
         idist[i] = j ;
         }

/*
   Use this table to build this training set by randomly sampling from tset
*/

      for (i=0 ; i<n ; i++) {           // Build this training set
         temp = unifrand() ;            // Uniform [0,1)
         j = (int) (m * temp + 0.999999) - 1 ;
         if (j < 0)                     // Happens very rarely
            j = 0 ;                     // So must be prepared
         j = idist[j] ;                 // Table guaranteed <= true index
         while (temp > cdf[j])          // This refinement loop should almost
            ++j ;                       // always end the first time or quickly
         tptr = tset + j * (nin + 1) ;  // Point to this case
         models[imodel]->add_case ( tptr ) ; // Add it to the model's training set
         }

      models[imodel]->train () ; // Train this model

/*
   Compute eps as the probability of error.  Use this to compute optimal alpha.
   The 'if' statement in this loop must be true when the actual class (as
   stored in the training set) is not the predicted class (as defined by the
   model's prediction).  This particular implementation uses a model having
   positive output for the first class, and negative for the second class.
   It considers the fact that the model may be returning any real number.
   Because our class AdaBoostBinaryNoConfSampled by definition ignores the
   numeric value, the test here involves whether the sign of the
   true and predicted is the same.
*/

      eps = 0.0 ;
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + 1) ;  // Point to this case
         models[imodel]->predict ( tptr , &out ) ;
         if (out > 0.0)      // If it predicts first class
            h[i] = 1.0 ;     // Flag it this way
         else                // But if it predicts second class
            h[i] = -1.0 ;    // Other flag value
         if ((tptr[nin] * h[i]) < 0.0)  // Error test: Signs the same?
            eps += dist[i] ; // Model erred if signs opposite
         }

      if (eps <= 0.0) {          // Unusual situation of model being perfect
         nmodels = imodel + 1 ;  // No sense going any further
         alpha[imodel] = 0.5 * log ( (double) n ) ; // Arbitrary large value
         break ;                 // So exit training loop early
         }

      if (eps > 0.5) {           // Unusual situation of model being worthless
         nmodels = imodel ;      // No sense going any further
         break ;                 // So exit training loop early
         }

      alpha[imodel] = 0.5 * log ( (1.0 - eps) / eps ) ;
      beta = eps / (1.0 - eps) ;

/*
   Adjust the relative weighting, then rescale to make it a distribution
*/

      temp = 0.0 ;
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + 1) ;  // Point to this case
         if ((tptr[nin] * h[i]) > 0.0)  // If this case correctly classified
            dist[i] *= beta ;           // Reduce its probability weighting
         temp += dist[i] ;              // Cumulate total probability for rescale
         }

      for (i=0 ; i<n ; i++)     // A probability distribution must sum to 1.0
         dist[i] /= temp ;      // Make it so

      } // For all models

   free ( dist ) ;
   free ( h ) ;
   free ( cdf ) ;
   free ( idist ) ;
}

AdaBoostBinaryNoConfSampled::~AdaBoostBinaryNoConfSampled ()
{
   if (alpha != NULL)
      free ( alpha ) ;
}

/*
   Make a class prediction.
   For compatibility with the other models, it returns 0 for the first class
   and 1 for the second class.
   Recall that by definition, AdaBoostBinaryNoConfSampled ignores the
   prediction's magnitude, being concerned only with the sign.
*/

int AdaBoostBinaryNoConfSampled::class_predict ( double *input )
{
   int i ;
   double sum, out ;

   if (nmodels == 0)   // Abnormal condition of no decent models
      return -1 ;      // Return an error flag

   sum = 0.0 ;
   for (i=0 ; i<nmodels ; i++) {
      models[i]->predict ( input , &out ) ;
      if (out > 0.0)          // If it predicts first class
         sum += alpha[i] ;
      else if (out < 0.0)     // But if it predicts second class
         sum -= alpha[i] ;
      }

   if (sum > 0.0)
      return 0 ;    // Flags first class
   else 
      return 1 ;    // Flags second class
}

/*
--------------------------------------------------------------------------------

   AdaBoostBinary - AdaBoost algorithm for strictly binary classification
      in which the underlying model provides a class plus numerical confidence.
      The first class is signified by +1 in the training set, and the second
      by -1.  The underlying model should attempt to make similar predicitons.
      For prediction, the model is presumed to have predicted the first class
      if its prediction is positive, and the second class if negative.
      The magnitude of the prediction is interpreted as confidence.
      For compatibility with this model as well as the published literature,
      the 'h' value is hard limited here to +/-1.  It is good if this is
      the natural range of your model, so that the hard limiting has little
      or no impact.  This limiting aids stability when the underlying model
      has a chance of occasionally producing wild output values.
      Also, having a known natural range helps keep the search for an
      optimal alpha efficient.

--------------------------------------------------------------------------------
*/

class AdaBoostBinary {

public:

   AdaBoostBinary ( int n , int nin , double *tset , int nmods ) ;
   ~AdaBoostBinary () ;
   int class_predict ( double *input ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *alpha ;    // Nmods long alpha constant for each model
   double *dist ;     // N long probability distribution
   double *u ;        // N long work area for saving model's error products
} ;

/*
   This local routine is the criterion function that is passed to the
   minimization routines that compute the optimal alpha.
*/

static int local_n ;
static double *local_dist, *local_u ;

static double alpha_crit ( double trial_alpha )
{
   int i ;
   double sum ;

   sum = 0.0 ;
   for (i=0 ; i<local_n ; i++)
      sum += local_dist[i] * exp ( -trial_alpha * local_u[i] ) ;
   return sum ;
}

/*
   Constructor
*/

AdaBoostBinary::AdaBoostBinary (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout, where nout=1)
   int nmods          // Number of models
   )
{
   int i, imodel, ngood, nbad ;
   double *tptr, temp, sum, h ;
   double x1, y1, x2, y2, x3, y3 ;

   nmodels = nmods ;
   alpha = (double *) malloc ( nmodels * sizeof(double) ) ;
   dist = (double *) malloc ( n * sizeof(double) ) ;
   u = (double *) malloc ( n * sizeof(double) ) ;

/*
   Initialize distribution to be uniform
*/

   temp = 1.0 / n ;
   for (i=0 ; i<n ; i++)
      dist[i] = temp ;

/*
   Main training loop trains sequence of models
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {

      models[imodel]->reset() ;         // Prepares the reusable model

      for (i=0 ; i<n ; i++) {           // Build this training set
         tptr = tset + i * (nin + 1) ;  // Point to this case
         models[imodel]->add_case ( tptr , dist[i] ) ; // Add it to the model's training set
         }

      models[imodel]->train () ; // Train this model

/*
   Compute the optimal alpha.
   In the standard algorithm, the sign of h is the model's class prediction
   (+ for first class, - for second) and the magnitude of h is the model's
   confidence in the prediction.  Theoretically, h need not be bounded because
   alpha can always compensate for the scaling of h.
   However, things seem to go best when the model is trained to produce
   predictions (h) in a fixed known range.  Among other things, this lets
   the search for alpha be made efficiently.
*/

      sum = 0.0 ;
      ngood = nbad = 0 ;                // Degenerate if all cases good or bad
      for (i=0 ; i<n ; i++) {
         tptr = tset + i * (nin + 1) ;  // Point to this case
         models[imodel]->predict ( tptr , &h ) ;
         if (h > 1.0)  // Hard limiting for a potentially wild model
            h = 1.0 ;  // like a neural net helps stability
         if (h < -1.0)
            h = -1.0 ;
         u[i] = h * tptr[nin] ; // Error indicator is predicted times true
         if (u[i] > 0.0)        // Class prediction is correct
            ++ngood ;           // This lets us detect degenerate situation
         if (u[i] < 0.0)        // Class prediction is incorrect
            ++nbad ;            // If degenerate, optimal alpha does not exist
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
      for (i=0 ; i<n ; i++) {
         dist[i] *= exp ( -alpha[imodel] * u[i] ) ;
         sum += dist[i] ;      // Cumulate total probability for rescale
         }

      for (i=0 ; i<n ; i++)     // A probability distribution must sum to 1.0
         dist[i] /= sum ;       // Make it so

      } // For all models

   free ( dist ) ;
   free ( u ) ;
}

AdaBoostBinary::~AdaBoostBinary ()
{
   if (alpha != NULL)
      free ( alpha ) ;
}

/*
   Make a class prediction.
   For compatibility with the other models, it returns 0 for the first class
   and 1 for the second class.
*/

int AdaBoostBinary::class_predict ( double *input )
{
   int i ;
   double h, sum ;

   if (nmodels == 0)   // Abnormal condition of no decent models
      return -1 ;      // Return an error flag

   sum = 0.0 ;
   for (i=0 ; i<nmodels ; i++) {
      models[i]->predict ( input , &h ) ;
      if (h > 1.0)  // Hard limiting for a potentially wild model
         h = 1.0 ;  // like a neural net helps stability
      if (h < -1.0)
         h = -1.0 ;
      sum += alpha[i] * h ;
      }

   if (sum > 0.0)
      return 0 ;    // Flags first class
   else 
      return 1 ;    // Flags second class
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
   int i, k, ntries, itry, nsamps, nmodels, ndone, nhid ;
   double *x, *test, separation, out, diff, temp, temp2 ;
   double sum_numeric_error, sum_class_error, sum_train_error, train_error ;
   double bagging_numeric_error, bagging_class_error, bagging_train_error ;
   double adaboost_binary_noconf_class_error ;
   double adaboost_binary_noconf_train_error ;
   double adaboost_binary_noconf_sampled_class_error ;
   double adaboost_binary_noconf_sampled_train_error ;
   double adaboost_binary_class_error ;
   double adaboost_binary_train_error ;
   Bagging *bagging ;
   AdaBoostBinaryNoConf *adaboost_binary_noconf ;
   AdaBoostBinaryNoConfSampled *adaboost_binary_noconf_sampled ;
   AdaBoostBinary *adaboost_binary ;

   nhid = 1 ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf (
         "\nUsage: ARCING  nsamples  nmodels  ntries  separation" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nmodels = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   separation = atof ( argv[4] ) ;
#else
   nsamps = 100 ;
   nmodels = 2 ;
   ntries = 10 ;
   separation = 0.0 ;
#endif

   if ((nsamps <= 0)  ||  (nmodels <= 0)  ||  (ntries <= 0)  ||  (separation < 0.0)) {
      printf ( "\nUsage: ARCING  nsamples  nmodels  ntries  separation" ) ;
      exit ( 1 ) ;
      }

/*
   Allocate memory and initialize
*/

   model = new MLFN ( nsamps , 2 , 1 , nhid ) ;
   models = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;
   for (i=0 ; i<nmodels ; i++)
      models[i] = new MLFN ( nsamps , 2 , 1 , nhid ) ;

   x = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * 3 * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   sum_numeric_error = 0.0 ;   // For comparison purposes, real error
   sum_class_error = 0.0 ;     // For comparison purposes, class error
   sum_train_error = 0.0 ;     // For comparison purposes, training error

   bagging_numeric_error = 0.0 ;
   bagging_class_error = 0.0 ;
   bagging_train_error = 0.0 ;
   adaboost_binary_noconf_class_error = 0.0 ;
   adaboost_binary_noconf_train_error = 0.0 ;
   adaboost_binary_noconf_sampled_class_error = 0.0 ;
   adaboost_binary_noconf_sampled_train_error = 0.0 ;
   adaboost_binary_class_error = 0.0 ;
   adaboost_binary_train_error = 0.0 ;

   for (itry=0 ; itry<ntries ; itry++) {
      ndone = itry + 1 ;

/*
   Generate the data.
   It is bivariate clusters with moderate positive correlation.
   One class is shifted above and to one (random) side of of the other class.
   We use x as the dataset for all resampling algorithms.
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[3*i] = normal () ;
         x[3*i+1] = .7071 * x[3*i]  +  .7071 * normal () ;
         if (unifrand() > 0.5) {
            x[3*i] -= separation ;
            if (unifrand() > 0.8)
               x[3*i+1] += 5 * separation ;
            else 
               x[3*i+1] -= 5 * separation ;
            x[3*i+2] = 1.0 ;
            }
         else {
            x[3*i] += separation ;
            x[3*i+2] = -1.0 ;
            }
         }

      for (i=0 ; i<10*nsamps ; i++) {
         test[3*i] = normal () ;
         test[3*i+1] = .7071 * test[3*i]  +  .7071 * normal () ;
         if (unifrand() > 0.5) {
            test[3*i] -= separation ;
            if (unifrand() > 0.8)
               test[3*i+1] += 5 * separation ;
            else 
               test[3*i+1] -= 5 * separation ;
            test[3*i+2] = 1.0 ;
            }
         else {
            test[3*i] += separation ;
            test[3*i+2] = -1.0 ;
            }
         }

/*
   Train a model with this data and test it on an independent test set.
   This gives us a basis of comparison for the resampling methods.
   The numeric error is the ordinary mean squared error.
   The class error is zero if the order relationship of the predicted outputs
   is the same as the order relationship for the true values.  Else one.
*/

      model->reset () ;
      for (i=0 ; i<nsamps ; i++)
         model->add_case ( x + 3 * i ) ;
      model->train () ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         model->predict ( x + 3 * i , &out ) ;
         if (x[3*i+2] * out < 0.0)
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp = temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         model->predict ( test + 3 * i , &out ) ;
         if (test[3*i+2] * out < 0.0)
            temp2 += 1.0 ;
         if (out > 1.0)
            out = 1.0 ;
         if (out < -1.0)
            out = -1.0 ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }

      sum_train_error += train_error ;
      sum_numeric_error += temp / (10 * nsamps) ;
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

      bagging = new Bagging ( nsamps , 2 , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         k = bagging->class_predict ( x + 3 * i ) ;
         if ((x[3*i+2] > 0.0)  &&  (k != 0))
            train_error += 1.0 ;
         if ((x[3*i+2] < 0.0)  &&  (k != 1))
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp = temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         bagging->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         k = bagging->class_predict ( test + 3 * i ) ;
         if ((test[3*i+2] > 0.0)  &&  (k != 0))
            temp2 += 1.0 ;
         if ((test[3*i+2] < 0.0)  &&  (k != 1))
            temp2 += 1.0 ;
         }

      bagging_numeric_error += temp / (10 * nsamps) ;
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
   Handle adaboost_binary_noconf
*/

      adaboost_binary_noconf = new AdaBoostBinaryNoConf ( nsamps , 2 , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         k = adaboost_binary_noconf->class_predict ( x + 3 * i ) ;
         if ((x[3*i+2] > 0.0)  &&  (k != 0))
            train_error += 1.0 ;
         if ((x[3*i+2] < 0.0)  &&  (k != 1))
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         k = adaboost_binary_noconf->class_predict ( test + 3 * i ) ;
         if ((test[3*i+2] > 0.0)  &&  (k != 0))
            temp2 += 1.0 ;
         if ((test[3*i+2] < 0.0)  &&  (k != 1))
            temp2 += 1.0 ;
         }

      adaboost_binary_noconf_class_error += temp2 / (10 * nsamps) ;
      adaboost_binary_noconf_train_error += train_error ;
      delete adaboost_binary_noconf ;

      printf ( "\nAdaBoostBinaryNoconf error:                     Class =%7.4lf (%7.4lf)",
               adaboost_binary_noconf_class_error / ndone ,
               adaboost_binary_noconf_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

/*
   Handle adaboost_binary_noconf_sampled
*/

      adaboost_binary_noconf_sampled = new AdaBoostBinaryNoConfSampled ( nsamps , 2 , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         k = adaboost_binary_noconf_sampled->class_predict ( x + 3 * i ) ;
         if ((x[3*i+2] > 0.0)  &&  (k != 0))
            train_error += 1.0 ;
         if ((x[3*i+2] < 0.0)  &&  (k != 1))
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         k = adaboost_binary_noconf_sampled->class_predict ( test + 3 * i ) ;
         if ((test[3*i+2] > 0.0)  &&  (k != 0))
            temp2 += 1.0 ;
         if ((test[3*i+2] < 0.0)  &&  (k != 1))
            temp2 += 1.0 ;
         }

      adaboost_binary_noconf_sampled_class_error += temp2 / (10 * nsamps) ;
      adaboost_binary_noconf_sampled_train_error += train_error ;
      delete adaboost_binary_noconf_sampled ;

      printf ( "\nAdaBoostBinaryNoconfSampled error:              Class =%7.4lf (%7.4lf)",
               adaboost_binary_noconf_sampled_class_error / ndone ,
               adaboost_binary_noconf_sampled_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

/*
   Handle adaboost_binary
*/

      adaboost_binary = new AdaBoostBinary ( nsamps , 2 , x , nmodels ) ;

      train_error = 0.0 ;
      for (i=0 ; i<nsamps ; i++) {
         k = adaboost_binary->class_predict ( x + 3 * i ) ;
         if ((x[3*i+2] > 0.0)  &&  (k != 0))
            train_error += 1.0 ;
         if ((x[3*i+2] < 0.0)  &&  (k != 1))
            train_error += 1.0 ;
         }
      train_error /= nsamps ;

      temp2 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         k = adaboost_binary->class_predict ( test + 3 * i ) ;
         if ((test[3*i+2] > 0.0)  &&  (k != 0))
            temp2 += 1.0 ;
         if ((test[3*i+2] < 0.0)  &&  (k != 1))
            temp2 += 1.0 ;
         }

      adaboost_binary_class_error += temp2 / (10 * nsamps) ;
      adaboost_binary_train_error += train_error ;
      delete adaboost_binary ;

      printf ( "\nAdaBoostBinary error:                           Class =%7.4lf (%7.4lf)",
               adaboost_binary_class_error / ndone ,
               adaboost_binary_train_error / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

     } // For all tries


   delete model ;
   for (i=0 ; i<nmodels ; i++)
      delete models[i] ;
   free ( models ) ;
   free ( x ) ;
   free ( test ) ;

   return EXIT_SUCCESS ;
}

/******************************************************************************/
/*                                                                            */
/*  MULTCLAS - Compare methods for combining multiple class predictors        */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "..\mlfn.h"
#include "..\logistic.h"

#define DEBUG 0

double unifrand () ;
double normal () ;
void qsortd ( int first , int last , double *data ) ;
void qsortdsi ( int first , int last , double *data , int *slave ) ;

static MLFN **models ;
static MLFN **model_pairs ;

/*
--------------------------------------------------------------------------------

   All of these routines let the model be an external reference.
   This is generally regarded as sloppy coding, but there is a reason here.
   By making the models external, we avoid the need to pass a typed identifier
   in the parameter lists.  An alternative method would be to define a parent
   "Model" class, but this would mean redefining the specific model to reflect
   its parentage.  The method shown here fits in best with the other code.
   Feel free to modify it as desired.

   Also note that many of these model combiners include variables in the
   constructor call that are not needed!  For example, the Average class
   is so simple that it does not need the number of training cases, the
   number of inputs, or the training set.  Yet this information is passed
   to the constructor.  The only reason I did this was to make all of the
   combiner objects share a common interface.  When you code these algorithms
   for your application, you will obvious change this.
--------------------------------------------------------------------------------
*/

/*
--------------------------------------------------------------------------------

   Average - Compute the simple average of the predictions

   This is a 'full resolution' version of majority rule.

--------------------------------------------------------------------------------
*/

class Average {

public:

   Average ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Average () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
} ;

Average::Average (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
}

Average::~Average ()
{
   free ( outwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int Average::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         output[i] += outwork[i] ; // Average output across all models
      }

/*
   Output[i] now contains the sum of all models' predictions for class i.
   Find the largest, and normalize to unit sum.
   This normalization only makes sense if the models' predictions are 0-1
   probabilities or some such thing.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   Median - Compute the median of the predictions

--------------------------------------------------------------------------------
*/

class Median {

public:

   Median ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Median () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
   double *sortwork ; // Work vector nout * nmodels long
} ;

Median::Median (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   sortwork = (double *) malloc ( nout * nmodels * sizeof(double) ) ;
}

Median::~Median ()
{
   free ( outwork ) ;
   free ( sortwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int Median::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum, *rptr ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         sortwork[i*nmodels+imodel] = outwork[i] ;
      }

   for (i=0 ; i<nout ; i++) {           // Sort each class' outputs
      rptr = sortwork + i * nmodels ;   // Point to this class's output vector
      qsortd ( 0 , nmodels-1 , rptr ) ; // Median sort is faster
      if (nmodels % 2)                  // If odd, median is center
         output[i] = rptr[nmodels/2] ;
      else                              // Else mean of two around center
         output[i] = 0.5 * (rptr[nmodels/2-1] + rptr[nmodels/2]) ;
      }

/*
   Output[i] now contains the median of all models' predictions for class i.
   Find the largest, and normalize to unit sum.
   This normalization only makes sense if the model's predictions are 0-1
   probabilities or some such thing.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   MaxMax - Compute the maximum of the predictions

   This is best when every class has a specialist model good at finding it.

--------------------------------------------------------------------------------
*/

class MaxMax {

public:

   MaxMax ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~MaxMax () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
} ;

MaxMax::MaxMax (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
}

MaxMax::~MaxMax ()
{
   free ( outwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int MaxMax::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++) {
         if ((imodel == 0)  ||  (outwork[i] > output[i]))
            output[i] = outwork[i] ; // Max output across all models
         }
      }

/*
   Output[i] now contains the max of all models' predictions for class i.
   Find the largest, and normalize to unit sum.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   MaxMin - Compute the minimum of the predictions

   This is best when some models are good at culling classes.

--------------------------------------------------------------------------------
*/

class MaxMin {

public:

   MaxMin ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~MaxMin () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
} ;

MaxMin::MaxMin (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
}

MaxMin::~MaxMin ()
{
   free ( outwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int MaxMin::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++) {
         if ((imodel == 0)  ||  (outwork[i] < output[i]))
            output[i] = outwork[i] ; // Min output across all models
         }
      }

/*
   Output[i] now contains the min of all models' predictions for class i.
   Find the largest, and normalize to unit sum.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   Majority - Compute the winner via simple majority.

   Note that no tie breaking is done for the models's outputs, as in most
   applications the model will generate continuous outputs related to
   confidence.  If this is not the case, add tie breaking if desired.

--------------------------------------------------------------------------------
*/

class Majority {

public:

   Majority ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Majority () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
} ;

Majority::Majority (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
}

Majority::~Majority ()
{
   free ( outwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int Majority::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum, temp ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++) {     // Find which class wins for this model
         if ((i == 0)  ||  (outwork[i] > best)) {
            best = outwork[i] ;
            ibest = i ;
            }
         }
      output[ibest] += 1.0 ;  // Tally this winning class
      }

/*
   Output[i] now contains the number of models that chose class i.
   Find the winner, and normalize to unit sum.
   Note that we employ an interesting tie-breaking scheme.  By adding a
   uniform random number less than 1.0 to each count, we ensure that
   no true winner is ignored, but ties are randomly broken.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      temp = output[i] + 0.999 * unifrand() ;
      if ((i == 0)  ||  (temp > best)) {
         best = temp ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   Borda - Compute the winner via Borda count

   Note that no tie breaking is done for the models's outputs, as in most
   applications the model will generate continuous outputs related to
   confidence.  If this is not the case, add tie breaking if desired.

--------------------------------------------------------------------------------
*/

class Borda {

public:

   Borda ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Borda () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   int *iwork ;       // Work vector nout long
   double *outwork ;  // Work vector nout long
} ;

Borda::Borda (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   nout = nclasses ;
   nmodels = nmods ;
   iwork = (int *) malloc ( nout * sizeof(int) ) ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
}

Borda::~Borda ()
{
   free ( outwork ) ;
   free ( iwork ) ;
}

/*
   This returns the integer 0 through nout-1 for its class decision.
   It also returns class responses, normalized to sum to one.
*/

int Borda::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum, temp ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

/*
   Compute the class output vectors for each model, then sort them.
   The Borda count for each class is the sum across models of the number
   of classes in that model ranked worse.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         iwork[i] = i ;               // Initialize for sorted indices
      qsortdsi ( 0 , nout-1 , outwork , iwork ) ;
      for (i=0 ; i<nout ; i++)        // Could skip i=0
         output[iwork[i]] += i ;
      }

/*
   Output[i] now contains the Borda count of class i.
   Find the winner, and normalize to unit sum.
   Note that we employ an interesting tie-breaking scheme.  By adding a
   uniform random number less than 1.0 to each count, we ensure that
   no true winner is ignored, but ties are randomly broken.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      temp = output[i] + 0.999 * unifrand() ;
      if ((i == 0)  ||  (temp > best)) {
         best = temp ;
         ibest = i ;
         }
      sum += output[i] ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   Intersection - Use intersection rule to compute minimal class set

--------------------------------------------------------------------------------
*/

class Intersection {

public:

   Intersection ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Intersection () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
   int *iwork ;       // Work vector nout long
   int *rank_cuts ;   // Work vector nmodels long
} ;

Intersection::Intersection (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int i, j, k, imodel, icase, nbad ;
   double *case_ptr, best ;

   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   rank_cuts = (int *) malloc ( nmodels * sizeof(int) ) ;
   iwork = (int *) malloc ( nout * sizeof(int) ) ;

/*
   Pass through the training set, invoking all models for each case.
   For each model, keep track of the worst rank across the training set.
   In other words, for each model, find the minimum number of highest-ranked
   classes we would need to keep in order for that subset to contain the
   correct class for every case in the training set.  Naturally, the
   holy grail is one, meaning that for this outstanding model its highest
   output is always the correct class.  Dream on.  The worst situation
   is when this threshold is nout, meaning that we would need to keep the
   entire set of all classes in order to get the entire training set right.
   This disastrous situation means the model is worthless as far as this
   algorithm is concerned.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++)
      rank_cuts[imodel] = 0 ;  // Will keep track of minimal size needed

   for (icase=0 ; icase<n ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Invoke each model.
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->predict ( case_ptr , outwork ) ;
         // Count how many of this model's outputs equal or exceed the output
         // corresponding to the true class.
         best = outwork[k] ;  // Output of true class
         nbad = 1 ;  // Counts i=k, skipped below so we don't have to trust fpt
         for (i=0 ; i<nout ; i++) {
            if (i == k)
               continue ;
            if (outwork[i] >= best)  // This should be true for i=k
               ++nbad ;       // Count number needed in subset
            }
         if (nbad > rank_cuts[imodel])  // Keep track of the worst performance
            rank_cuts[imodel] = nbad ;  // of this model across training set
         } // For all models
      } // For all cases
}

Intersection::~Intersection ()
{
   free ( outwork ) ;
   free ( iwork ) ;
   free ( rank_cuts ) ;
}

/*
   This returns the number of classes in the minimal set.  This number
   ranges from zero (failure to find a common subset, indicating possible
   outlier) to nout (all classes are a possibility, which happens when the
   models are all weak).
   Ideally, it will return one, but this only happens regularly for easy tasks.
   It also returns class flags in 'output'.  Whatever number this routine
   returns as the subset count, this many of the flags will be set to one,
   indicating their membership in the subset.  The rest will be zero.
*/

int Intersection::classify ( double *input , double *output )
{
   int i, imodel, n ;

/*
   Compute the class output vectors for each model, then sort them.
*/

   for (i=0 ; i<nout ; i++)
      output[i] = 1.0 ;      // Will kill these off below

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         iwork[i] = i ;               // Initialize for sorted indices
      qsortdsi ( 0 , nout-1 , outwork , iwork ) ;
      n = nout - rank_cuts[imodel] ;  // This many classes in models bad part
      for (i=0 ; i<n ; i++)           // Pass through this bad part
         output[iwork[i]] = 0.0 ;     // And nix every class there
      }

   n = 0 ;                            // Will count survivors
   for (i=0 ; i<nout ; i++) {         // Check each class
      if (output[i] > 0.5)            // 1.0 if in intersection, else 0.0
         ++n ;
      }

   return n ;
}

/*
--------------------------------------------------------------------------------

   Union - Use union rule to compute minimal class set

--------------------------------------------------------------------------------
*/

class Union {

public:

   Union ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Union () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
   int *iwork ;       // Work vector nout long
   int *rank_cuts ;   // Work vector nmodels long
} ;

Union::Union (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int i, j, k, imodel, icase, nbad, bestrank, ibestrank ;
   double *case_ptr, best ;

   nout = nclasses ;
   nmodels = nmods ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   rank_cuts = (int *) malloc ( nmodels * sizeof(int) ) ;
   iwork = (int *) malloc ( nout * sizeof(int) ) ;

/*
   Pass through the training set, invoking all models for each case.
   For each case, find the model that provides the best rank, then keep track
   (separately, for each model) of the worst of these best cases.
   Naturally, the holy grail is ones across the board, meaning that for every
   case there was at least one model that ranked its true class number one.
   Any model that is never the best ranker is superfluous as far as this
   algorithm is concerned.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++)
      rank_cuts[imodel] = 0 ;  // Will keep track of minimal size needed

   for (icase=0 ; icase<n ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Invoke each model.
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->predict ( case_ptr , outwork ) ;
         // Count how many of this model's outputs equal or exceed the output
         // corresponding to the true class.
         best = outwork[k] ;  // Output of true class
         nbad = 1 ;  // Counts i=k, skipped below so we don't have to trust fpt
         for (i=0 ; i<nout ; i++) {
            if (i == k)
               continue ;
            if (outwork[i] >= best)  // This should be true for i=k
               ++nbad ;       // Count number needed in subset
            }
         // Nbad is now the rank of this case's true class according to imodel.
         if ((imodel == 0)  ||  (nbad < bestrank)) {
            bestrank = nbad ;     // Keep track of the best model's rank
            ibestrank = imodel ;  // And know which model was best for this case
            }
         } // For all models
      // The best rank was 'bestrank' and it was achieved by model 'ibestrank'
      if (bestrank > rank_cuts[ibestrank]) // Keep track of worst best performance
         rank_cuts[ibestrank] = bestrank ; // of this model across training set
      } // For all cases
}

Union::~Union ()
{
   free ( outwork ) ;
   free ( iwork ) ;
   free ( rank_cuts ) ;
}

/*
   This returns the number of classes in the minimal set.  This number
   ranges from one (no model had rank_cuts greater than one, and all models
   agreed on their choice --- ideal!) to nout (at least one model had
   rank_cuts=nout, or there was massive disagreement among the models).
   Ideally, it will return one, but this only happens regularly for easy tasks.
   It also returns class flags in 'output'.  Whatever number this routine
   returns as the subset count, this many of the flags will be set to one,
   indicating their membership in the subset.  The rest will be zero.
*/

int Union::classify ( double *input , double *output )
{
   int i, imodel, n ;

/*
   Compute the class output vectors for each model, then sort them.
*/

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;      // Will activate these below

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         iwork[i] = i ;               // Initialize for sorted indices
      qsortdsi ( 0 , nout-1 , outwork , iwork ) ;
      for (i=nout-rank_cuts[imodel] ; i<nout ; i++) // Pass through good part
         output[iwork[i]] = 1.0 ;     // And select every class there
      }

   n = 0 ;                            // Will count number in union
   for (i=0 ; i<nout ; i++) {         // Check each class
      if (output[i] > 0.5)            // 1.0 if in union, else 0.0
         ++n ;
      }

   return n ;
}

/*
--------------------------------------------------------------------------------

   Logit - Use logistic regression to find best class
           This uses one common weight vector for all classes.

--------------------------------------------------------------------------------
*/

class Logit {

public:

   Logit ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~Logit () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *inwork ;   // Work vector nmodels+1 long
   double *outwork ;  // Work vector nout long
   int *iwork ;       // Work vector nout long
   double *rankwork ; // Work vector nmodels * nout long
   Logistic *logit ;  // Logistic regression object
} ;

Logit::Logit (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int j, k, nbelow, iclass, imodel, icase ;
   double *case_ptr, best ;

   nout = nclasses ;
   nmodels = nmods ;
   inwork = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   iwork = (int *) malloc ( nout * sizeof(int) ) ;
   rankwork = (double *) malloc ( nmodels * nout * sizeof(double) ) ;
   logit = new Logistic ( n * nclasses , nmodels ) ;

/*
   Pass through the training set, building the logistic regression training
   set.  We add nout cases to that set for each original training case.
*/

   for (icase=0 ; icase<n ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Invoke each model
      // and save its outputs.
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( case_ptr , rankwork + imodel * nout ) ;
      // The nmodels by nout array contains the outputs for each model
      // Create the nout cases for logistic regression
      for (iclass=0 ; iclass<nout ; iclass++) {  // Add a case for each class
         for (imodel=0 ; imodel<nmodels ; imodel++) { // Each model is an input
            nbelow = 0 ;
            best = rankwork[imodel*nout+iclass] ;  // Output for this class
            for (j=0 ; j<nout ; j++) {   // Compute the number worse
               if (rankwork[imodel*nout+j] < best)
                  ++nbelow ;
               }
            inwork[imodel] = (double) nbelow / (double) nout ; // Mean # below
            }
         inwork[nmodels] = (iclass == k)  ?  1 : 0 ;
         logit->add_case ( inwork ) ;
         } // For all classes
      } // For all cases

   logit->train () ;
}

Logit::~Logit ()
{
   free ( inwork ) ;
   free ( outwork ) ;
   free ( iwork ) ;
   free ( rankwork ) ;
   delete logit ;
}

/*
   This is a weighted Borda count
*/

int Logit::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum, temp ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

/*
   Compute the class output vectors for each model, then sort them.
   The Borda count for each class is the sum across models of the number
   of classes in that model ranked worse.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         iwork[i] = i ;               // Initialize for sorted indices
      qsortdsi ( 0 , nout-1 , outwork , iwork ) ;
      temp = logit->coefs[imodel] ;   // Weight for this model
      for (i=0 ; i<nout ; i++)
         output[iwork[i]] += i * temp ;
      }

/*
   Output[i] now contains the (weighted) Borda count of class i.
   Weight them and find the winner, then normalize to unit sum.
   Note that tie-breaking is not needed because the weights almost guarantee
   no ties.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      temp = output[i] ;
      if ((i == 0)  ||  (temp > best)) {
         best = temp ;
         ibest = i ;
         }
      sum += temp ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   LogitSep - Use logistic regression to find best class.
              This uses separate weight vectors for each class.
              The training set should have at least ten times the number of
              models for EACH class represented!

--------------------------------------------------------------------------------
*/

class LogitSep {

public:

   LogitSep ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~LogitSep () ;
   int classify ( double *input , double *output ) ;

private:
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *inwork ;   // Work vector nmodels+1 long
   double *outwork ;  // Work vector nout long
   int *iwork ;       // Work vector nout long
   double *rankwork ; // Work vector nmodels * nout long
   Logistic **logit ; // Logistic regression objects (one for each class)
} ;

LogitSep::LogitSep (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int i, j, k, nbelow, iclass, imodel, icase ;
   double *case_ptr, best ;

   nout = nclasses ;
   nmodels = nmods ;
   inwork = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   iwork = (int *) malloc ( nout * sizeof(int) ) ;
   rankwork = (double *) malloc ( nmodels * nout * sizeof(double) ) ;
   logit = (Logistic **) malloc ( nout * sizeof(Logistic *) ) ;
   for (i=0 ; i<nout ; i++)
      logit[i] = new Logistic ( n , nmodels ) ;

/*
   Pass through the training set, building the logistic regression training
   sets, separately for each class.
*/

   for (icase=0 ; icase<n ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Invoke each model
      // and save its outputs.
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( case_ptr , rankwork + imodel * nout ) ;
      // The nmodels by nout array contains the outputs for each model
      // Create the nout cases for logistic regression
      for (iclass=0 ; iclass<nout ; iclass++) {  // Add a case for each class
         for (imodel=0 ; imodel<nmodels ; imodel++) { // Each model is an input
            nbelow = 0 ;
            best = rankwork[imodel*nout+iclass] ;  // Output for this class
            for (j=0 ; j<nout ; j++) {   // Compute the number worse
               if (rankwork[imodel*nout+j] < best)
                  ++nbelow ;
               }
            inwork[imodel] = (double) nbelow / (double) nout ; // Mean # below
            }
         inwork[nmodels] = (iclass == k)  ?  1 : 0 ;
         logit[iclass]->add_case ( inwork ) ;
         } // For all classes
      } // For all cases

   for (iclass=0 ; iclass<nout ; iclass++)
      logit[iclass]->train () ;
}

LogitSep::~LogitSep ()
{
   int i ;

   free ( inwork ) ;
   free ( outwork ) ;
   free ( iwork ) ;
   free ( rankwork ) ;
   for (i=0 ; i<nout ; i++)
      delete logit[i] ;
   free ( logit ) ;
}

/*
   This is a weighted Borda count
*/

int LogitSep::classify ( double *input , double *output )
{
   int i, imodel, ibest ;
   double best, sum, temp ;

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

/*
   Compute the class output vectors for each model, then sort them.
   The Borda count for each class is the sum across models of the number
   of classes in that model ranked worse.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      for (i=0 ; i<nout ; i++)
         iwork[i] = i ;               // Initialize for sorted indices
      qsortdsi ( 0 , nout-1 , outwork , iwork ) ;
      for (i=0 ; i<nout ; i++)
         output[iwork[i]] += i * logit[i]->coefs[imodel] ;
      }

/*
   Output[i] now contains the (weighted) Borda count of class i.
   Weight them and find the winner, then normalize to unit sum.
   Note that tie-breaking is not needed because the weights almost guarantee
   no ties.
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++) {
      temp = output[i] ;
      if ((i == 0)  ||  (temp > best)) {
         best = temp ;
         ibest = i ;
         }
      sum += temp ;
      }

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return ibest ;
}

/*
--------------------------------------------------------------------------------

   LocalAcc - Use local accuracy to choose the best model

--------------------------------------------------------------------------------
*/

class LocalAcc {

public:

   LocalAcc ( int n , int ninputs , int nclasses , double *tset , int nmods ) ;
   ~LocalAcc () ;
   int classify ( double *input , double *output ) ;

private:
   int knn ;          // K nearest neighbors will be used
   int ncases ;       // Number of cases
   int nin ;          // Number of original case inputs
   int nout ;         // Number of outputs (nclasses in constructor call)
   int nmodels ;      // Number of models (nmods in constructor call)
   double *outwork ;  // Work vector nout long
   int *iwork ;       // Work vector ncases long
   double *distwork ; // Work vector ncases long
   double *trnx ;     // Work vector ncases * nin long holds raw predictors
   int *trncls ;      // Work vector ncases * nmodels long holds model decisions
   int *trntrue ;     // Work vector ncases long holds true class of each case
   int classprep ;    // Very private flag so cross validation can skip prep
} ;

LocalAcc::LocalAcc (
   int n ,            // Number of training cases
   int ninputs ,      // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int i, j, k, iclass, imodel, icase, ibest, itemp ;
   int *clsptr1, *clsptr2, true_class, *knn_counts ;
   int knn_min, knn_max, knn_best ;
   double *case_ptr, *last_ptr, *testcase, *clswork, best ;

   ncases = n ;
   nin = ninputs ;
   nout = nclasses ;
   nmodels = nmods ;

   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   iwork = (int *) malloc ( ncases * sizeof(int) ) ;
   distwork = (double *) malloc ( ncases * sizeof(double) ) ;
   trnx = (double *) malloc ( ncases * nin * sizeof(double) ) ;
   trncls = (int *) malloc ( ncases * nmodels * sizeof(int) ) ;
   trntrue = (int *) malloc ( ncases * sizeof(int) ) ;

/*
   Pass through the training set, saving raw inputs in trnx.
   Invoke all models for each case.
   For each model, save its class decision (highest output) in trncls.
   Also save the true class of each case in trntrue.
*/

   for (icase=0 ; icase<ncases ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      memcpy ( trnx + icase * nin , case_ptr , nin * sizeof(double) ) ;
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Save it in trntrue.
      // Then invoke each model and save the class decision of each.
      trntrue[icase] = k ;     // Save true class of this case
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->predict ( case_ptr , outwork ) ;
         for (i=0 ; i<nout ; i++) {     // Find which class wins for this model
            if ((i == 0)  ||  (outwork[i] > best)) {
               best = outwork[i] ;
               ibest = i ;
               }
            } // For all outputs of this model
         trncls[icase*nmodels+imodel] = ibest ; // Save this model's decision
         } // For all models
      } // For all cases

   classprep = 1 ;   // Tell classify() that it must fully prepare

#if 0
// This option sets knn to a fixed arbitary (but probably decent) value
   knn = 3 ;
   return ;
#else

// This option finds optimal knn by (slow!) cross validation

   if (ncases < 20) {
      knn = 3 ;
      return ;
      }

   knn_min = 3 ;     // Require at least this size
   knn_max = 10 ;    // But at most this size

   testcase = (double *) malloc ( nin * sizeof(double) ) ;
   clswork = (double *) malloc ( nout * sizeof(double) ) ;
   knn_counts = (int *) malloc ( (knn_max - knn_min + 1) * sizeof(int) ) ;

   for (knn=knn_min ; knn<=knn_max ; knn++)
      knn_counts[knn-knn_min] = 0 ; // Will count correct decisions for each knn

   --ncases ;  // Cross validation uses reduced training set in classify().
   last_ptr = trnx + ncases * nin ; // Last case keeps getting swapped in

   for (icase=0 ; icase<=ncases ; icase++) { // Omitted-case loop

      // Copy case icase to testcase, then copy last case to its spot
      // Also swap data in trncls and trntrue
      case_ptr = trnx + icase * nin ; // Point to this case
      memcpy ( testcase , case_ptr , nin * sizeof(double) ) ;
      true_class = trntrue[icase] ;  // True class of test case
      if (icase < ncases) {  // Last case is already there, so no need to swap
         memcpy ( case_ptr , last_ptr , nin * sizeof(double) ) ;
         trntrue[icase] = trntrue[ncases] ;
         clsptr1 = trncls + icase * nmodels ;
         clsptr2 = trncls + ncases * nmodels ;
         for (i=0 ; i<nmodels ; i++) {        // Memcpy might be faster
            itemp = clsptr1[i] ;
            clsptr1[i] = clsptr2[i] ;
            clsptr2[i] = itemp ;
            }
         }

      classprep = 1 ;   // Tell classify() that it must fully prepare
      for (knn=knn_min ; knn<=knn_max ; knn++) {
         iclass = classify ( testcase , clswork ) ;
         if (iclass == true_class)        // Correct decision?
           ++knn_counts[knn-knn_min] ;    // Score this trial knn
         classprep = 0 ;   // Tell classify() that it does not need to prepare
         }

      // Done, so move original last case back to its slot and restore icase   
      if (icase < ncases) { // Last case does not need repair as nothing changed
         memcpy ( last_ptr , case_ptr , nin * sizeof(double) ) ;
         memcpy ( case_ptr , testcase , nin * sizeof(double) ) ;
         trntrue[ncases] = trntrue[icase] ;
         trntrue[icase] = true_class ;
         clsptr1 = trncls + icase * nmodels ;
         clsptr2 = trncls + ncases * nmodels ;
         for (i=0 ; i<nmodels ; i++) {        // Memcpy might be faster
            itemp = clsptr1[i] ;
            clsptr1[i] = clsptr2[i] ;
            clsptr2[i] = itemp ;
            }
         }
      }

   ++ncases ; // Restore size of training set to its correct value

   // See which trial value of knn had the best score.
   // Do not break ties.  Favor the lowest value for classification speed.
   for (knn=knn_min ; knn<=knn_max ; knn++) {
      if ((knn==knn_min)  ||   (knn_counts[knn-knn_min] > ibest)) {
         ibest = knn_counts[knn-knn_min] ;
         knn_best = knn ;
         }
      }

   knn = knn_best ;
   classprep = 1 ;   // Tell classify() that it must fully prepare

   free ( testcase ) ;
   free ( clswork ) ;
   free ( knn_counts ) ;
#endif
}

LocalAcc::~LocalAcc ()
{
   free ( outwork ) ;
   free ( iwork ) ;
   free ( distwork ) ;
   free ( trnx ) ;
   free ( trncls ) ;
   free ( trntrue ) ;
}

int LocalAcc::classify ( double *input , double *output )
{
   int i, k, icase, imodel, ibest, numer, denom, bestmodel, bestchoice ;
   double dist, *cptr, diff, best, crit, bestcrit, conf, bestconf, sum ;

/*
   Find the knn nearest neighbors, keeping track of their indices
   If knn is small, there are faster ways than sorting the whole thing.

   The classprep flag should always be true in normal use.  The only reason
   it is here is because cross validation in the constructor tries many
   values of knn, and this expensive preparation is only needed the
   first time.  The things computed in this block do not change if the
   input case remains the same.
*/

   if (classprep) {
      for (icase=0 ; icase<ncases ; icase++) {
         iwork[icase] = icase ;         // Save index for sorting later
         dist = 0.0 ;                   // Will cumulate Euclidean distance here
         cptr = trnx + icase * nin ;    // Point to this case
         for (i=0 ; i<nin ; i++) {      // For all original input variables
            diff = input[i] - cptr[i] ; // Input minus training case
            dist += diff * diff ;       // Cumulate Euclidean distance
            }
         distwork[icase] = dist ;       // Save distance
         }

      qsortdsi ( 0 , ncases-1 , distwork , iwork ) ;
      }

/*
   Invoke all models and find class decision of each.
   Search the knn nearest neighbors and identify those cases for
   which this model chose this class.
   Calculate the fraction of those cases in which the model was correct.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      sum = 0.0 ;
      for (i=0 ; i<nout ; i++) {     // Find which class wins for this model
         sum += outwork[i] ;         // For computing tie-breaking confidence
         if ((i == 0)  ||  (outwork[i] > best)) {
            best = outwork[i] ;
            ibest = i ;
            }
         } // For all outputs of this model
      conf = best / sum ;            // May be needed to break a tie later

      // This model chose class 'ibest' with confidence 'conf'

      denom = 0 ;    // Counts cases (in knn) in which this model chose ibest
      numer = 0 ;    // Counts cases (in denom) in which this model was correct

      for (icase=0 ; icase<knn ; icase++) {
         k = iwork[icase] ;  // Index of this case
         if (trncls[k*nmodels+imodel] == ibest) {
            ++denom ;        // This model chose this class
            if (ibest == trntrue[k])
               ++numer ;     // And it was a correct choice
            }
         } // For the knn nearest neighbors

      if (denom > 0)
         crit = (double) numer / (double) denom ;  // Model's worthiness
      else 
         crit = 0.0 ;

      if ((imodel == 0)  ||  (crit > bestcrit)) {
         bestcrit = crit ;
         bestmodel = imodel ;  // Not needed here, but may be useful some time
         bestchoice = ibest ;
         bestconf = conf ;
         memcpy ( output , outwork , nout * sizeof(double) ) ;
         }

      else if ( fabs (crit - bestcrit) < 1.e-10) {  // We must break a tie
         if (conf > bestconf) {      // Fairly arbitrary
            bestcrit = crit ;
            bestmodel = imodel ;
            bestchoice = ibest ;
            bestconf = conf ;
            memcpy ( output , outwork , nout * sizeof(double) ) ;
            }
         }

      } // For all models

/*
   Normalize the outputs of the best model for user's convenience
*/

   sum = 0.0 ;
   for (i=0 ; i<nout ; i++)
      sum += output[i] ;

   if (sum > 0.0) {
      for (i=0 ; i<nout ; i++)
         output[i] /= sum ;
      }

   return bestchoice ;
}

/*
--------------------------------------------------------------------------------

   FuzzyInt - Use fuzzy integral to combine decisions

--------------------------------------------------------------------------------
*/

class FuzzyInt {

public:

   FuzzyInt ( int n , int nin , int nclasses , double *tset , int nmods ) ;
   ~FuzzyInt () ;
   int classify ( double *input , double *output ) ;

private:
   double recurse ( double x ) ; // Recursively compute the final g(A)-1

   int nout ;          // Number of outputs (nclasses in constructor call)
   int nmodels ;       // Number of models (nmods in constructor call)
   int *iwork ;        // Work vector nmodels long
   double *outwork ;   // Work vector nout long
   double *sortwork ;  // Work vector nmodels * nout long
   double *g ;         // Model g-values, nmods long
   double lambda ;     // Overall lambda
} ;

FuzzyInt::FuzzyInt (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   int nclasses ,     // Number of outputs (classes)
   double *tset ,     // Training cases, n by (nin+nout)
   int nmods          // Number of models in 'models' array
   )
{
   int j, k, iclass, imodel, icase ;
   double *case_ptr, best, xlo, xhi, y, ylo, yhi, step ;

   nout = nclasses ;
   nmodels = nmods ;
   iwork = (int *) malloc ( nmodels * sizeof(int) ) ;
   outwork = (double *) malloc ( nout * sizeof(double) ) ;
   sortwork = (double *) malloc ( nmodels * nout * sizeof(double) ) ;
   g = (double *) malloc ( nmodels * sizeof(double) ) ;

/*
   Pass through the training set, invoking all models for each case.
   For each model, compute g, its success rate across the training set.
   This will range from zero for a model that always misclassifies, to one
   for a perfect model.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++)
      g[imodel] = 0.0 ;  // Will sum correct decisions for each model

   for (icase=0 ; icase<n ; icase++) {
      case_ptr = tset + icase * (nin + nout) ; // Point to this case
      // Find the true class of this case
      for (j=0 ; j<nout ; j++) {        // Scan case's output vector
         if ((j == 0)  ||  (case_ptr[nin+j] > best)) { // Outs after nin inputs
            best = case_ptr[nin+j] ;
            k = j ;   // Keep track of index of best so far
            }
         }
      // At this time, k is the true class of this case.  Invoke each model
      // and find its decision.
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->predict ( case_ptr , outwork ) ;
         for (j=0 ; j<nout ; j++) {  // Scan model's output vector
            if ((j == 0)  ||  (outwork[j] > best)) {
               best = outwork[j] ;
               iclass = j ;   // Keep track of index of best so far
               }
            }
         // If this model chose correctly, count it
         if (iclass == k)
            g[imodel] += 1.0 ;
         } // For all models
      } // For all cases

/*
   Divide by the number of cases to get a 0-1 success rate for each model.
   Then subtract the rate that would be obtained by simply guessing, and
   renormalize to get a 0-1 figure for the contribution of this model.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      g[imodel] /= n ;
      g[imodel] = (g[imodel] - 1.0 / nout) / (1.0 - 1.0 / nout) ;
      if (g[imodel] > 1.0)
         g[imodel] = 1.0 ;
      if (g[imodel] < 0.0)
         g[imodel] = 0.0 ;
      }

/*
   Compute lambda.  We know that it cannot be less than -1, so the first step
   is to test for that extreme condition, which will happen when all models
   are perfect.
*/

   xlo = lambda = -1.0 ;
   ylo = recurse ( xlo ) ;
   if (ylo >= 0.0)   // Theoretically should never exceed zero
      return ;       // But allow for pathological numerical problems

/*
   Now we must bound the root.  Step out until the function becomes positive.
   If all models are worthless, the root is infinite, so we must avoid that!
*/

   step = 1.0 ;

   for (;;) {
      xhi = xlo + step ;
      yhi = recurse ( xhi ) ;
      if (yhi >= 0.0)   // If we have just bracketed the root
         break ;        // We can quit the search
      if (xhi > 1.e5) { // In the unlikely case of extremely poor models
         lambda = xhi ; // Fudge a value
         return ;       // And quit
         }
      step *= 2.0 ;     // Keep increasing the step size to avoid many tries
      xlo = xhi ;       // Move onward
      ylo = yhi ;
      }

/*
   We have bracketed the root between (xlo, ylo) and (xhi, yhi).
   Bisection is a primitive way to refine, but it always succeeds.
*/

   for (;;) {
      lambda = 0.5 * (xlo + xhi) ;
      y = recurse ( lambda ) ;                  // Evaluate the function here
      if (fabs ( y ) < 1.e-8)                   // Primary convergence criterion
         break ;
      if (xhi - xlo < 1.e-10 * (lambda + 1.1))  // Backup criterion
         break ;
      if (y > 0.0) {
         xhi = lambda ;
         yhi = y ;
         }
      else {
         xlo = lambda ;
         ylo = y ;
         }
      }
}

FuzzyInt::~FuzzyInt ()
{
   free ( iwork ) ;
   free ( outwork ) ;
   free ( sortwork ) ;
   free ( g ) ;
}

/*
   This local routine recursively computes the final g(A) - 1, which must
   equal zero.  This is the function whose root, lambda, we seek.
*/

double FuzzyInt::recurse ( double x )
{
   int i ;
   double val ;

   val = g[0] ;
   for (i=1 ; i<nmodels ; i++)
      val += g[i] + x * g[i] * val ;

   return val - 1.0 ;
}

/*
   This returns the class decision, plus fuzzy class confidences.
*/

int FuzzyInt::classify ( double *input , double *output )
{
   int i, k, iclass, imodel ;
   double sum, gsum, *rptr, minval, maxmin, best ;

/*
   Invoke all models and store their (normalized) outputs.
   Note that if we can GUARANTEE that the models' outputs are already
   probabilities, this normalization is counterproductive.
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , outwork ) ;
      sum = 0.0 ;               // Will sum this model's outputs
      for (i=0 ; i<nout ; i++)  // So we can normalize them to
         sum += outwork[i] ;    // probability-like quantities
      for (i=0 ; i<nout ; i++)
         sortwork[i*nmodels+imodel] = outwork[i] / sum ;
      }

/*
   This main outer loop computes the fuzzy integral for each class
*/

   for (iclass=0 ; iclass<nout ; iclass++) { // Compute for each class

/*
   Sort the models according to their output for this class
*/

      for (imodel=0 ; imodel<nmodels ; imodel++)
         iwork[imodel] = imodel ;  // Initialize index identity vector
      rptr = sortwork + iclass * nmodels ;   // Point to this class's outputs
      qsortdsi ( 0 , nmodels-1 , rptr , iwork ) ; // Sort ascending

/*
   Run through the models from highest output (for this class) to lowest.
*/

      maxmin = 0.0 ;                    // Will keep track of max H here
      gsum = 0.0 ;                      // Will be cumulative g
      for (i=nmodels-1 ; i>=0 ; i--) {  // Sorted ascending, so max to min
         k = iwork[i] ;                 // Index of this model
         gsum += g[k] + lambda * g[k] * gsum ; // Cumulate g(A)
         if (gsum < rptr[k])            // Compare g so far to model's output
            minval = gsum ;             // We need the smaller of the two
         else 
            minval = rptr[k] ;
         if (minval > maxmin)           // Keep track of max of this min
            maxmin = minval ;           // (For the current class)
         }

      output[iclass] = maxmin ;
      } // For all classes

/*
   Return the class having maximum fuzzy integral.
   Don't worry about ties, which should be rare.
*/

   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         iclass = i ;
         }
      }

   return iclass ;
}

/*
--------------------------------------------------------------------------------

   Pairwise - Use pairwise coupling to combine decisions

--------------------------------------------------------------------------------
*/

class Pairwise {

public:

   Pairwise ( int nclasses , int *ntrain ) ;
   ~Pairwise () ;
   int classify ( double *input , double *output ) ;

private:

   int nout ;          // Number of outputs (nclasses in constructor call)
   int npairs ;        // Number of models (nclasses * (nclasses-1) / 2)
   int *nij ;          // Number of training cases used for each model (ntrain)
   double *rij ;       // Models' predicted Prob ( i given (i or j) ); i<j
   double *uij ;       // p-hat sub i / (p-hat sub i + p-hat sub j) ; i<j
                       // For rij and uij, just subtract from 1 if i>j
} ;

Pairwise::Pairwise (
   int nclasses ,     // Number of outputs (classes)
   int *ntrain        // Number of training cases used for each model
   )
{
   nout = nclasses ;
   npairs = nclasses * (nclasses-1) / 2 ;
   nij = ntrain ;

   rij = (double *) malloc ( npairs * sizeof(double) ) ;
   uij = (double *) malloc ( npairs * sizeof(double) ) ;
}

Pairwise::~Pairwise ()
{
   free ( rij ) ;
   free ( uij ) ;
}

/*
   This returns the class decision, plus class probabilities.
*/

int Pairwise::classify ( double *input , double *output )
{
   int i, j, k, iclass, iter ;
   double rr, best, numer, denom, sum, delta, oldval ;

/*
   Evaluate all models for this case
*/

   for (i=0 ; i<npairs ; i++) {
      model_pairs[i]->predict ( input , &rr ) ;
      if (rr > 0.999999)   // Prevent numerical difficulties later
         rr = 0.999999 ;
      if (rr < 0.000001)
         rr = 0.000001 ;
      rij[i] = rr ;
      }

#if DEBUG > 2
   rij[0] = 0.9 ;
   rij[1] = 0.4 ;
   rij[2] = 0.7 ;
#endif

/*
   Compute starting estimates as row means
*/

   for (i=0 ; i<nout ; i++)
      output[i] = 0.0 ;

   k = 0 ;
   for (i=0 ; i<nout-1 ; i++) {     // First class in pair
      for (j=i+1 ; j<nout ; j++) {  // Second class in pair
         rr = rij[k++] ;            // This model's Prob[i given i or j]; i<j
         output[i] += rr ;          // Term for i<j
         output[j] += 1.0 - rr ;    // And symmetric term
         }
      }

   for (i=0 ; i<nout ; i++)
      output[i] /= npairs ;

#if DEBUG
   printf ( "\n\nModel predictions: " ) ;
   for (i=0 ; i<npairs ; i++)
      printf ( " %lf", rij[i] ) ;
   printf ( "\nStarting estimates: " ) ;
   for (i=0 ; i<nout ; i++)
      printf ( " %lf", output[i] ) ;
   getch () ;
#endif

/*
   If all we need is the class (the index of the max output) then the next
   step of refinement is not needed.  The order will not change.
   This refinement provides more accurate class probabilities.
   They nearly always spread out compared to the initial estimates.
   We begin by computing uij, which approximates rij.
*/

   k = 0 ;
   for (i=0 ; i<nout-1 ; i++) {    // First class in pair
      for (j=i+1 ; j<nout ; j++)   // Second class in pair
         uij[k++] = output[i] / (output[i] + output[j]) ;
      }

   for (iter=0 ; iter<10000 ; iter++) {  // Cheap insurance against endless loop

      delta = 0.0 ; // Keeps track of max change in any output for convergence test
      for (i=0 ; i<nout ; i++) {      // Correct one output at a time

         numer = denom = 0.0 ;        // Will cumulate for correction factor
         for (j=0 ; j<nout ; j++) {   // Cumulate numer and denom
            if (i < j) {
               k = (i * (2 * nout - i - 3) - 2) / 2 + j ;
               numer += nij[k] * rij[k] ;
               denom += nij[k] * uij[k] ;
               }
            else if (i > j) {
               k = (j * (2 * nout - j - 3) - 2) / 2 + i ;
               numer += nij[k] * (1.0 - rij[k]) ;
               denom += nij[k] * (1.0 - uij[k]) ;
               }
            } // For j, cumulating numer and denom for correction factor

         // Correct this one output (trial p) and renormalize to sum to one
         oldval = output[i] ;  // For convergence test
         output[i] *= numer / denom ;
         sum = 0.0 ;
         for (j=0 ; j<nout ; j++)
            sum += output[j] ;
         for (j=0 ; j<nout ; j++)
            output[j] /= sum ;

         if (fabs(output[i]-oldval) > delta)  // How much did this output change?
            delta = fabs(output[i]-oldval) ;  // Keep track of max

         // Recompute uij from the modified output
         k = 0 ;
         for (i=0 ; i<nout-1 ; i++) {    // First class in pair
            for (j=i+1 ; j<nout ; j++)   // Second class in pair
               uij[k++] = output[i] / (output[i] + output[j]) ;
            }

#if DEBUG > 1
         printf ( "\nIteration %d: ", iter ) ;
         for (i=0 ; i<nout ; i++)
            printf ( " %lf", output[i] ) ;
#if DEBUG > 2
         getch () ;
#endif
#endif

         } // For i

      if (delta < 1.e-6)   // If max change of any output is small
         break ;           // No need to keep iterating; we have converged

      } // For all iterations

#if DEBUG
         printf ( "\nFinal: " ) ;
         for (i=0 ; i<nout ; i++)
            printf ( " %lf", output[i] ) ;
         getch () ;
#endif

/*
   We're done.  Output the class having maximum probability.
*/

   for (i=0 ; i<nout ; i++) {
      if ((i == 0)  ||  (output[i] > best)) {
         best = output[i] ;
         iclass = i ;
         }
      }

   return iclass ;
}


/*
--------------------------------------------------------------------------------

   Optional main to test it

   If there are four or more models, the fourth model is deliberately worthless.
   If there are five or more models, the fifth model has some wild outputs.

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, j, k, ntries, itry, nsamps, imodel, nmodels, divisor, ndone ;
   int iclass, ibest, nclasses, n, npairs, *ntrain_pair, nh ;
   double *x, *xbad, *xwild, *test, *input, best ;
   double spread, temp, *out, temp1, temp2, temp3 ;
   double *computed_err_raw ;
   double computed_err_average ;
   double computed_err_median ;
   double computed_err_maxmax ;
   double computed_err_maxmin ;
   double computed_err_intersection_1 ;
   double computed_err_intersection_2 ;
   double computed_err_intersection_3 ;
   double computed_err_union_1 ;
   double computed_err_union_2 ;
   double computed_err_union_3 ;
   double computed_err_majority ;
   double computed_err_borda ;
   double computed_err_logit ;
   double computed_err_logitsep ;
   double computed_err_localacc ;
   double computed_err_fuzzyint ;
   double computed_err_pairwise ;
   Average *average ;
   Median *median ;
   MaxMax *maxmax ;
   MaxMin *maxmin ;
   Intersection *intersection ;
   Union *union_rule ;
   Majority *majority ;
   Borda *borda ;
   Logit *logit ;
   LogitSep *logitsep ;
   LocalAcc *localacc ;
   FuzzyInt *fuzzyint ;
   Pairwise *pairwise ;

   int nhid = 4 ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 6) {
      printf (
         "\nUsage: MULTCLAS  nsamples  nclasses  nmodels  ntries  spread" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nclasses = atoi ( argv[2] ) ;
   nmodels = atoi ( argv[3] ) ;
   ntries = atoi ( argv[4] ) ;
   spread = atof ( argv[5] ) ;
#else
   nsamps = 10 ;
   nclasses = 3 ;
   nmodels = 3 ;
   ntries = 1000 ;
   spread = 1.0 ;
#endif

   if ((nsamps <= 0)  ||  (nclasses <= 1)  ||  (nmodels <= 0)  ||  (ntries <= 0)
    || (spread < 0.0)) {
      printf (
         "\nUsage: MULTCLAS  nsamples  nclasses  nmodels  ntries  spread" ) ;
      exit ( 1 ) ;
      }

   divisor = 1 ;

/*
   Allocate memory and initialize
*/

   npairs = nclasses * (nclasses-1) / 2 ;
   models = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;
   model_pairs = (MLFN **) malloc ( npairs * sizeof(MLFN *) ) ;

   for (i=0 ; i<nmodels ; i++)
      models[i] = new MLFN ( nsamps , 2 , nclasses , nhid ) ;

   x = (double *) malloc ( nsamps * (2+nclasses) * sizeof(double) ) ;
   xbad = (double *) malloc ( nsamps * (2+nclasses) * sizeof(double) ) ;
   xwild = (double *) malloc ( nsamps * (2+nclasses) * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * (2+nclasses) * sizeof(double) ) ;
   computed_err_raw = (double *) malloc ( nmodels * sizeof(double) ) ;
   input = (double *) malloc ( 3 * sizeof(double) ) ;
   out = (double *) malloc ( nclasses * sizeof(double) ) ;
   ntrain_pair = (int *) malloc ( npairs * sizeof(int) ) ;

   for (imodel=0 ; imodel<nmodels ; imodel++)
      computed_err_raw[imodel] = 0.0 ;

   computed_err_average = 0.0 ;
   computed_err_median = 0.0 ;
   computed_err_maxmax = 0.0 ;
   computed_err_maxmin = 0.0 ;
   computed_err_intersection_1 = 0.0 ;
   computed_err_intersection_2 = 0.0 ;
   computed_err_intersection_3 = 0.0 ;
   computed_err_union_1 = 0.0 ;
   computed_err_union_2 = 0.0 ;
   computed_err_union_3 = 0.0 ;
   computed_err_majority = 0.0 ;
   computed_err_borda = 0.0 ;
   computed_err_logit = 0.0 ;
   computed_err_logitsep = 0.0 ;
   computed_err_localacc = 0.0 ;
   computed_err_fuzzyint = 0.0 ;
   computed_err_pairwise = 0.0 ;

/*
   Main outer loop does all tries
*/

   for (itry=0 ; itry<ntries ; itry++) {
      ndone = itry + 1 ;

/*
   Generate the data.
   We use x as the dataset for all prediction algorithms.
   (Actually, for the fourth model (if any), x is modified to create xbad
   to provide useless training data for this one model.  And for the fifth model
   if any, the output is occasionally wild (far outside 0-1).
   Note that the class for each case is chosen randomly, EXCEPT for the first
   few cases.  This ensures that each class is represented by at least one case.
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[(2+nclasses)*i] = normal () ;
         x[(2+nclasses)*i+1] = normal () ;
         for (k=0 ; k<nclasses ; k++)      // Init all outputs to zero
            x[(2+nclasses)*i+2+k] = 0.0 ;  // (Only one, the true class will be 1)
         if (i < nclasses) // Ensure that each class has at least one case
            k = i ;
         else 
            k = (int) (unifrand() * nclasses) ; // Randomly choose the correct class
         if (k >= nclasses)                // Should never happen, but be safe
            k = nclasses - 1 ;
         x[(2+nclasses)*i+2+k] = 1.0 ;     // Flag the correct class
         x[(2+nclasses)*i] += k * spread ; // Set predictors appropriatly
         x[(2+nclasses)*i+1] -= k * spread ;
         }

      if (nmodels >= 4) {                  // Generate totally random data
         for (i=0 ; i<nsamps ; i++) {
            xbad[(2+nclasses)*i] = x[(2+nclasses)*i] ;
            xbad[(2+nclasses)*i+1] = x[(2+nclasses)*i+1] ;
            for (k=0 ; k<nclasses ; k++)
               xbad[(2+nclasses)*i+2+k] = (unifrand() < 1.0/nclasses) ? 1.0 : 0.0 ;
            }
         }

      if (nmodels >= 5) {                  // Generate wild data
         for (i=0 ; i<nsamps ; i++) {
            xwild[(2+nclasses)*i] = x[(2+nclasses)*i] ;
            xwild[(2+nclasses)*i+1] = x[(2+nclasses)*i+1] ;
            for (k=0 ; k<nclasses ; k++) {
               if (unifrand() < 0.1)
                  xwild[(2+nclasses)*i+2+k] =
                             x[(2+nclasses)*i+2+k] * 1000.0 - 500.0 ;
               }
            }
         }

      for (i=0 ; i<10*nsamps ; i++) {      // Build a test dataset
         test[(2+nclasses)*i] = normal () ;
         test[(2+nclasses)*i+1] = normal () ;
         for (k=0 ; k<nclasses ; k++)
            test[(2+nclasses)*i+2+k] = 0.0 ;
         k = (int) (unifrand() * nclasses) ;
         if (k >= nclasses)
            k = nclasses - 1 ;
         test[(2+nclasses)*i+2+k] = 1.0 ;
         test[(2+nclasses)*i] += k * spread ;
         test[(2+nclasses)*i+1] -= k * spread ;
         }

/*
   Train the individual models.
   Also test each individually after it is trained.
*/

      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->reset () ;
         if (imodel == 3) {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( xbad + (2+nclasses) * i ) ;
            }
         else if (imodel == 4) {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( xwild + (2+nclasses) * i ) ;
            }
         else {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( x + (2+nclasses) * i ) ;
            }

         models[imodel]->train () ;
   
         temp = 0.0 ;
         for (i=0 ; i<10*nsamps ; i++) {
            models[imodel]->predict ( test + (2+nclasses) * i , out ) ;
            for (j=0 ; j<nclasses ; j++) {
               if ((j == 0)  ||  (out[j] > best)) {
                  best = out[j] ;
                  iclass = j ;
                  }
               }
            for (j=0 ; j<nclasses ; j++) {
               if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
                  best = test[(2+nclasses)*i+2+j] ;
                  ibest = j ;
                  }
               }
            if (iclass != ibest)  // If predicted class not true class
               temp += 1.0 ;      // Count this misclassification
            }
         computed_err_raw[imodel] += temp / (10 * nsamps) ;
         } // For all models

/*
   Train the models for pairwise classification.
   First we need to pass through the training set and count how many
   training cases are in each model for a class pair.
   We create the models here and delete them later, inside the try loop.
*/

      imodel = 0 ;                         // Identifies model in pair
      for (i=0 ; i<nclasses-1 ; i++) {     // First class in a pair
         for (j=i+1 ; j<nclasses ; j++) {  // Second class in pair

            ntrain_pair[imodel] = 0 ;     // Will count training cases for this model
            for (k=0 ; k<nsamps ; k++) {  // For all training cases
               if ((x[(2+nclasses) * k + 2 + i] > 0.5)  // Case is first class?
                || (x[(2+nclasses) * k + 2 + j] > 0.5)) // Or second class?
                  ++ntrain_pair[imodel] ;
               }
            nh = (nhid < ntrain_pair[imodel] - 1) ? nhid : ntrain_pair[imodel] - 1 ; // Ensure legal
            model_pairs[imodel] = new MLFN (ntrain_pair[imodel] , 2 , 1 , nh ) ;
            model_pairs[imodel]->reset () ;

            for (k=0 ; k<nsamps ; k++) {  // For all training cases
               input[0] = x[(2+nclasses) * k] ;     // First predictor
               input[1] = x[(2+nclasses) * k + 1] ; // Second predictor
               if (x[(2+nclasses) * k + 2 + i] > 0.5) // Case is first class?
                  input[2] = 1.0 ;
               else if (x[(2+nclasses) * k + 2 + j] > 0.5) // Case is second class?
                  input[2] = 0.0 ;
               else            // This case is neither first or second class
                  continue ;   // So it doesn't belong in training set
               model_pairs[imodel]->add_case ( input ) ;
               } // For all training cases
            model_pairs[imodel]->train () ;
            ++imodel ;  // Next model in collection of pairs
            } // For j, which is second class in pair
         } // For i, which is first class in pair

/*
   Average
*/

      average = new Average ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = average->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_average += temp / (10 * nsamps) ;
      delete average ;

/*
   Median
*/

      median = new Median ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = median->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_median += temp / (10 * nsamps) ;
      delete median ;

/*
   MaxMax
*/

      maxmax = new MaxMax ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = maxmax->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_maxmax += temp / (10 * nsamps) ;
      delete maxmax ;

/*
   MaxMin
*/

      maxmin = new MaxMin ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = maxmin->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_maxmin += temp / (10 * nsamps) ;
      delete maxmin ;

/*
   Intersection
*/

      intersection = new Intersection ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp1 = temp2 = temp3 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         n = intersection->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (out[ibest] < 0.5) { // If true class not in subset
            temp1 += 1.0 ;       // Count this misclassification
            temp2 += 1.0 ;       // In all categories
            temp3 += 1.0 ;
            }
         else {                  // True class is in subset.  Consider its size
            if (n > 3)
               temp3 += 1.0 ;
            if (n > 2)
               temp2 += 1.0 ;
            if (n > 1)
               temp1 += 1.0 ;
            }
         } // For all test cases
      computed_err_intersection_1 += temp1 / (10 * nsamps) ;
      computed_err_intersection_2 += temp2 / (10 * nsamps) ;
      computed_err_intersection_3 += temp3 / (10 * nsamps) ;
      delete intersection ;

/*
   Union
*/

      union_rule = new Union ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp1 = temp2 = temp3 = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         n = union_rule->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (out[ibest] < 0.5) { // If true class not in subset
            temp1 += 1.0 ;       // Count this misclassification
            temp2 += 1.0 ;       // In all categories
            temp3 += 1.0 ;
            }
         else {                  // True class is in subset.  Consider its size
            if (n > 3)
               temp3 += 1.0 ;
            if (n > 2)
               temp2 += 1.0 ;
            if (n > 1)
               temp1 += 1.0 ;
            }
         } // For all test cases
      computed_err_union_1 += temp1 / (10 * nsamps) ;
      computed_err_union_2 += temp2 / (10 * nsamps) ;
      computed_err_union_3 += temp3 / (10 * nsamps) ;
      delete union_rule ;

/*
   Majority
*/

      majority = new Majority ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = majority->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_majority += temp / (10 * nsamps) ;
      delete majority ;

/*
   Borda
*/

      borda = new Borda ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = borda->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_borda += temp / (10 * nsamps) ;
      delete borda ;

/*
   Logit
*/

      logit = new Logit ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = logit->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_logit += temp / (10 * nsamps) ;
      delete logit ;

/*
   Logitsep
*/

      logitsep = new LogitSep ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = logitsep->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_logitsep += temp / (10 * nsamps) ;
      delete logitsep ;

/*
   LocalAcc
*/

      localacc = new LocalAcc ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = localacc->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_localacc += temp / (10 * nsamps) ;
      delete localacc ;

/*
   FuzzyInt
*/

      fuzzyint = new FuzzyInt ( nsamps , 2 , nclasses , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = fuzzyint->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_fuzzyint += temp / (10 * nsamps) ;
      delete fuzzyint ;

/*
   Pairwise
*/

      pairwise = new Pairwise ( nclasses , ntrain_pair ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         iclass = pairwise->classify ( test + (2+nclasses) * i , out ) ;
         for (j=0 ; j<nclasses ; j++) {
            if ((j == 0)  ||  (test[(2+nclasses)*i+2+j] > best)) {
               best = test[(2+nclasses)*i+2+j] ;
               ibest = j ;
               }
            }
         if (iclass != ibest)  // If predicted class not true class
            temp += 1.0 ;      // Count this misclassification
         }
      computed_err_pairwise += temp / (10 * nsamps) ;
      delete pairwise ;

      for (i=0 ; i<npairs ; i++)
         delete model_pairs[i] ;

/*
   Print results so far
*/

      temp = 0.0 ;
      printf ( "\n\n\nDid%5d    Raw errors:", ndone ) ;
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         printf ( "  %.4lf", computed_err_raw[imodel] / ndone ) ;
         temp += computed_err_raw[imodel] / ndone ;
         }
      printf ( "\n       Mean raw error = %8.5lf", temp / nmodels ) ;
      printf ( "\n        Average error = %8.5lf", computed_err_average / ndone ) ;
      printf ( "\n         Median error = %8.5lf", computed_err_median / ndone ) ;
      printf ( "\n         MaxMax error = %8.5lf", computed_err_maxmax / ndone ) ;
      printf ( "\n         MaxMin error = %8.5lf", computed_err_maxmin / ndone ) ;
      printf ( "\n       Majority error = %8.5lf", computed_err_majority / ndone ) ;
      printf ( "\n          Borda error = %8.5lf", computed_err_borda / ndone ) ;
      printf ( "\n          Logit error = %8.5lf", computed_err_logit / ndone ) ;
      printf ( "\n       Logitsep error = %8.5lf", computed_err_logitsep / ndone ) ;
      printf ( "\n       LocalAcc error = %8.5lf", computed_err_localacc / ndone ) ;
      printf ( "\n       FuzzyInt error = %8.5lf", computed_err_fuzzyint / ndone ) ;
      printf ( "\n       Pairwise error = %8.5lf", computed_err_pairwise / ndone ) ;
      printf ( "\n Intersection error 1 = %8.5lf", computed_err_intersection_1 / ndone ) ;
      printf ( "\n Intersection error 2 = %8.5lf", computed_err_intersection_2 / ndone ) ;
      printf ( "\n Intersection error 3 = %8.5lf", computed_err_intersection_3 / ndone ) ;
      printf ( "\n        Union error 1 = %8.5lf", computed_err_union_1 / ndone ) ;
      printf ( "\n        Union error 2 = %8.5lf", computed_err_union_2 / ndone ) ;
      printf ( "\n        Union error 3 = %8.5lf", computed_err_union_3 / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

     } // For all tries

   _getch () ;
   return EXIT_SUCCESS ;
}

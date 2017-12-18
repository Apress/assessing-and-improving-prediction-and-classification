/******************************************************************************/
/*                                                                            */
/*  GRNNGATE - Test GRNN gating                                               */
/*                                                                            */
/*  This uses an arbitrary number of univariate real predictors               */
/*  to make a combined univariate prediction.                                 */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "..\mlfn.h"
#include "..\minimize.h"

double unifrand () ;
double normal () ;
void qsortd ( int first , int last , double *data ) ;
void qsortdsi ( int first , int last , double *data , int *slave ) ;

static MLFN **models ;
static double criter ( double *params ) ;
static double univar_criter ( double param ) ;

/*
--------------------------------------------------------------------------------

   GRNNgate - GRNN gating model combination

   This class lets the model be an external reference.
   This is generally regarded as sloppy coding, but there is a reason here.
   By making the models external, we avoid the need to pass a typed identifier
   in the parameter lists.  An alternative method would be to define a parent
   "Model" class, but this would mean redefining the specific model to reflect
   its parentage.  The method shown here fits in best with the other code.
   Feel free to modify it as desired.

--------------------------------------------------------------------------------
*/

class GRNNgate {

public:

   GRNNgate ( int n , int n_gates , int nmods , double *gates ,
              double *contenders , double *trueval ) ;
   ~GRNNgate () ;
   double trial ( double *gates , double *contenders ,
                  int i_exclude , int n_exclude ) ;
   void predict ( double *gates , double *contenders , double *output ) ;

private:
   int ncases ;       // Number of cases
   int ngates ;       // Number of gate variables
   int nmodels ;      // Number of models (nmods in constructor call)
   double *tset ;     // A copy of the trainng set is perpetually needed
   double *sigma ;    // Learned sigma weights
   double *errvals ;  // Used in trial()

   friend double criter ( double *params ) ;
} ;

static GRNNgate *local_model ;  // Used in criter()

GRNNgate::GRNNgate (
   int n ,            // Number of training cases
   int n_gates ,      // Number of gate variables
   int nmods ,        // Number of models in 'models' array
   double *gates ,    // Training cases: gate variables, n by n_gates
   double *contenders , // Training cases: contender variables, n by nmods
   double *trueval    // Training cases: n true values
   )
{

   int i ;
   double *case_ptr, *params, *base, *p0, *direc ;
   double x1, x2, x3, y1, y2, y3, err ;

   ncases = n ;
   ngates = n_gates ;
   nmodels = nmods ;

/*
   Allocate memory that is needed throughout the life of the GRNNgate object.
   Then copy the training data to a local area.
*/

   tset = (double *) malloc ( ncases * (ngates+nmodels+1) * sizeof(double) ) ;
   sigma = (double *) malloc ( ngates * sizeof(double) ) ;
   errvals = (double *) malloc ( nmodels * sizeof(double) ) ;

   case_ptr = tset ;
   for (i=0 ; i<ncases ; i++) {
      memcpy ( case_ptr , gates + i * ngates , ngates * sizeof(double) ) ;
      case_ptr += ngates ;
      memcpy ( case_ptr , contenders + i * nmodels , nmodels * sizeof(double) );
      case_ptr += nmodels ;
      *case_ptr++ = trueval[i] ;
      }

/*
   Allocate scratch memory that is only needed during construction:
*/

   params = (double *) malloc ( ngates * sizeof(double) ) ;

/*
   Keep a local copy of this class pointer (needed by the criterion function).
   Call the optimizer, then call criter() one last time to set sigmas.
*/

   local_model = this ;

   if (ngates == 1) {
      glob_min ( -3.0 , 3.0 , 15 , 0 , 0.0 , univar_criter ,
                 &x1 , &y1 , &x2 , &y2 , &x3 , &y3 ) ;
      brentmin ( 10 , 0.0 , 1.e-5 , 1.e-5 ,
                 univar_criter , &x1 , &x2 , &x3 , y2 ) ;
      params[0] = x2 ;
      }
   else {
      base = (double *) malloc ( (2 * ngates + ngates * ngates) * sizeof(double));
      p0 = base + ngates ;
      direc = p0 + ngates ;
      for (i=0 ; i<ngates ; i++)
         params[i] = 0.0 ;
      err = criter ( params ) ;
      if (err > 0.0) // Perfect model makes powell() unnecessary
         err = powell ( 10 , 0.0 , 1.e-4 , criter , ngates , params ,
                        err , base , p0 , direc ) ;
         free ( base ) ;
      }

   criter ( params ) ;  // One last call to set sigma vector
   free ( params ) ;
}

GRNNgate::~GRNNgate ()
{
   free ( tset ) ;
   free ( sigma ) ;
   free ( errvals ) ;
}

/*
   Local routine is the criterion for the optimizer
*/

static double criter ( double *params )
{
   int i, ngates, nmodels ;
   double *inputs, out, diff, error, penalty ;

   ngates = local_model->ngates ;
   nmodels = local_model->nmodels ;

/*
   Get the sigmas from the parameter vector
*/

   penalty = 0.0 ;
   for (i=0 ; i<ngates ; i++) {
      if (params[i] > 8.0) {
         local_model->sigma[i] = exp ( 8.0 ) ;
         penalty += 10.0 * (params[i] - 8.0) ;
         }
      else if (params[i] < -8.0) {
         local_model->sigma[i] = exp ( -8.0 ) ;
         penalty += 10.0 * (-params[i] - 8.0) ;
         }
      else
         local_model->sigma[i] = exp ( params[i] ) ;
      }

/*
   This is the main loop that does the cross validation
*/

   error = 0.0 ;

   for (i=0 ; i<local_model->ncases ; i++) {  // For each case in the training set
      inputs = local_model->tset + i * (ngates+nmodels+1) ; // This case
      out = local_model->trial ( inputs , inputs+ngates , i , 0 ) ;
      diff = inputs[ngates+nmodels] - out ; // True minus predicted
      error += diff * diff ;
      } // For each case

   return error / local_model->ncases + penalty ;
}

static double univar_criter ( double param )
{
   double tempvec[1] ;

   tempvec[0] = param ;
   return criter ( tempvec ) ;
}


double GRNNgate::trial ( double *gates , double *contenders ,
                         int i_exclude , int n_exclude )
{
   int icase, imodel, ivar, idist, size ;
   double psum, *dptr, diff, dist, *outptr, err, out ;

   for (imodel=0 ; imodel<nmodels ; imodel++)
      errvals[imodel] = 0.0 ;

   size = ngates + nmodels + 1 ;              // Size of each training case

   for (icase=0 ; icase<ncases ; icase++) {  // Do all training cases

      idist = abs ( i_exclude - icase ) ;    // How close to excluded case?
      if (ncases - idist < idist)            // Also check distance going
         idist = ncases - idist ;            // around end

      if (idist <= n_exclude)                // If we are too close to excluded
         continue ;                          // Skip this tset case

      dptr = tset + size * icase ;           // Point to this case
      dist = 0.0 ;                           // Will sum distance here

      for (ivar=0 ; ivar<ngates ; ivar++) {  // All gates in this case
         diff = gates[ivar] - dptr[ivar] ;   // Input minus case
         diff /= sigma[ivar] ;               // Scale per sigma
         dist += diff * diff ;               // Cumulate Euclidean distance
         }

      dist = exp ( -dist ) ;                 // Apply the Gaussian kernel
      dptr += ngates ;                       // Contenders stored after gates
      outptr = dptr + nmodels ;              // True output is last of all

      for (ivar=0 ; ivar<nmodels ; ivar++) { // For every contender
         err = dptr[ivar] - *outptr ;        // Error for this contender
         errvals[ivar] += dist * err * err ; // Cumulate numerators
         }

      } // For all training cases

/*
   Convert unnormalized squared errors to weights
*/

   psum = 0.0 ;
   for (ivar=0 ; ivar<nmodels ; ivar++) {
      if (errvals[ivar] > 1.e-30)
         errvals[ivar] = 1.0 / errvals[ivar] ;
      else
         errvals[ivar] = 1.e30 ;   // Cheap insurance against dividing by zero
      psum += errvals[ivar] ;
      }

   for (ivar=0 ; ivar<nmodels ; ivar++)
      errvals[ivar] /= psum ;

/*
   Compute weighted output
*/

   out = 0.0 ;
   for (ivar=0 ; ivar<nmodels ; ivar++)
      out += errvals[ivar] * contenders[ivar] ;

   return out ;
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
   int i, ntries, itry, nsamps, imodel, nmodels, ndone, max_gates ;
   double *x, *xbad, *xwild, *test, std ;
   double err, diff, out, err1, err2 ;
   double *gates, *contenders, *all_gates , *all_contenders , *all_trueval ;
   double *computed_err_raw ;
   double computed_err_afterfact ;
   double computed_err_original ;
   double computed_err_random ;
   double computed_err_ratio ;
   GRNNgate *gate ;

   int nhid = 2 ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf (
         "\nUsage: GRNNGATE  nsamples  nmodels  ntries  std" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nmodels = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   std = atof ( argv[4] ) ;
#else
   nsamps = 6 ;
   nmodels = 2 ;
   ntries = 3 ;
   std = 0.5 ;
#endif

   if ((nsamps <= 0)  ||  (nmodels <= 0)  ||  (ntries <= 0)
    || (std < 0.0)) {
      printf (
         "\nUsage: GRNNGATE  nsamples  nmodels  ntries  std" ) ;
      exit ( 1 ) ;
      }

/*
   Allocate memory and initialize
*/

   max_gates = (nmodels > 2) ? nmodels : 2 ;
   models = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;
   gates = (double *) malloc ( max_gates * sizeof(double) ) ;
   contenders = (double *) malloc ( nmodels * sizeof(double) ) ;
   all_gates = (double *) malloc ( nsamps * max_gates * sizeof(double) ) ;
   all_contenders = (double *) malloc ( nsamps * nmodels * sizeof(double) ) ;
   all_trueval = (double *) malloc ( nsamps * sizeof(double) ) ;

   for (i=0 ; i<nmodels ; i++)
      models[i] = new MLFN ( nsamps , 2 , 1 , nhid ) ;

   x = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   xbad = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   xwild = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * 3 * sizeof(double) ) ;
   computed_err_raw = (double *) malloc ( nmodels * sizeof(double) ) ;

   for (imodel=0 ; imodel<nmodels ; imodel++)
      computed_err_raw[imodel] = 0.0 ;

   computed_err_afterfact = 0.0 ;
   computed_err_original = 0.0 ;
   computed_err_random = 0.0 ;
   computed_err_ratio = 0.0 ;

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
   if any, the output is occasionally wild.
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[3*i] = normal () ;                              // First predictor
         x[3*i+1] = normal () ;                            // Second predictor
         x[3*i+2] = x[3*i] - x[3*i+1] + std * normal () ;  // Predicted
         }

      if (nmodels >= 4) {                  // Generate totally random data
         for (i=0 ; i<nsamps ; i++) {
            xbad[3*i] = x[3*i] ;
            xbad[3*i+1] = x[3*i+1] ;
            xbad[3*i+2] = normal () ;
            }
         }

      if (nmodels >= 5) {                  // Generate wild data
         for (i=0 ; i<nsamps ; i++) {
            xwild[3*i] = x[3*i] ;
            xwild[3*i+1] = x[3*i+1] ;
            xwild[3*i+2] = x[3*i+2] * 1000.0 ;
            }
         }

      for (i=0 ; i<10*nsamps ; i++) {      // Build a test dataset
         test[3*i] = normal () ;
         test[3*i+1] = normal () ;
         test[3*i+2] = test[3*i] - test[3*i+1] + std * normal () ;
         }

      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->reset () ;
         if (imodel == 3) {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( xbad + 3 * i ) ;
            }
         else if (imodel == 4) {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( xwild + 3 * i ) ;
            }
         else {
            for (i=0 ; i<nsamps ; i++)
               models[imodel]->add_case ( x + 3 * i ) ;
            }

         models[imodel]->train () ;
   
         err = 0.0 ;
         for (i=0 ; i<10*nsamps ; i++) {
            models[imodel]->predict ( test + 3 * i , &out ) ;
            diff = test[3*i+2] - out ;
            err += diff * diff ;
            }
         computed_err_raw[imodel] += err / (10 * nsamps) ;
         } // For all models

/*
   Compute and print the errors of the raw models
*/

      err = 0.0 ;
      printf ( "\n\n\nDid%5d    Raw errors:", ndone ) ;
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         printf ( "  %.4lf", computed_err_raw[imodel] / ndone ) ;
         err += computed_err_raw[imodel] / ndone ;
         }
      printf ( "\n       Mean raw error = %8.5lf", err / nmodels ) ;

/*
   Pass through the training set, evaluating each model for each case.
   Save this information in the 'all_contenders' array for use by the gates.
   Also save true values.
*/

      for (i=0 ; i<nsamps ; i++) {
         all_trueval[i] = x[3*i+2] ;
         for (imodel=0 ; imodel<nmodels ; imodel++)
            models[imodel]->predict ( x + 3 * i ,
                                      all_contenders + nmodels * i + imodel ) ;
         }

/*
   Compute and print results for after-the-fact.
   This uses the contenders (model outputs) as the gates.
*/

      gate = new GRNNgate ( nsamps , nmodels , nmodels ,
                            all_contenders , all_contenders , all_trueval ) ;
      err = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         for (imodel=0 ; imodel<nmodels ; imodel++) {
            models[imodel]->predict ( test + 3 * i , &contenders[imodel] ) ;
            gates[imodel] = contenders[imodel] ;
            }
         out = gate->trial ( gates , contenders , -1 , 0 ) ;
         diff = test[3*i+2] - out ;
         err += diff * diff ;
         }
      computed_err_afterfact += err / (10 * nsamps) ;
      delete gate ;
      printf ( "\n  AfterFact error = %8.5lf", computed_err_afterfact / ndone );

/*
   Compute and print results for original variables as gates.
*/

      for (i=0 ; i<nsamps ; i++) {
         all_gates[2*i] = x[3*i] ;
         all_gates[2*i+1] = x[3*i+1] ;
         }

      gate = new GRNNgate ( nsamps , 2 , nmodels ,
                            all_gates , all_contenders , all_trueval ) ;
      err = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         for (imodel=0 ; imodel<nmodels ; imodel++)
            models[imodel]->predict ( test + 3 * i , &contenders[imodel] ) ;
         gates[0] = test[3*i] ;
         gates[1] = test[3*i+1] ;
         out = gate->trial ( gates , contenders , -1 , 0 ) ;
         diff = test[3*i+2] - out ;
         err += diff * diff ;
         }
      computed_err_original += err / (10 * nsamps) ;
      delete gate ;
      printf ( "\n  Original error = %8.5lf", computed_err_original / ndone );

/*
   Compute and print results for a random number as a gate.
*/

      for (i=0 ; i<nsamps ; i++)
         all_gates[i] = normal () ;

      gate = new GRNNgate ( nsamps , 1 , nmodels ,
                            all_gates , all_contenders , all_trueval ) ;
      err = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         for (imodel=0 ; imodel<nmodels ; imodel++)
            models[imodel]->predict ( test + 3 * i , &contenders[imodel] ) ;
         gates[0] = normal () ;
         out = gate->trial ( gates , contenders , -1 , 0 ) ;
         diff = test[3*i+2] - out ;
         err += diff * diff ;
         }
      computed_err_random += err / (10 * nsamps) ;
      delete gate ;
      printf ( "\n  Random error = %8.5lf", computed_err_random / ndone ) ;

/*
   Compute and print results for models' error ratio as a gate.
*/

      for (i=0 ; i<nsamps ; i++) {
         err1 = all_contenders[i*nmodels] - all_trueval[i] ;
         err2 = all_contenders[i*nmodels+1] - all_trueval[i] ;
         all_gates[i] = log ( (fabs(err1)+1.e-60) / (fabs(err2)+1.e-60) ) ;
         }

      gate = new GRNNgate ( nsamps , 1 , nmodels ,
                            all_gates , all_contenders , all_trueval ) ;
      err = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         for (imodel=0 ; imodel<nmodels ; imodel++)
            models[imodel]->predict ( test + 3 * i , &contenders[imodel] ) ;
         err1 = contenders[0] - test[3*i+2] ;
         err2 = contenders[1] - test[3*i+2] ;
         gates[0] = log ( (fabs(err1)+1.e-60) / (fabs(err2)+1.e-60) ) ;
         out = gate->trial ( gates , contenders , -1 , 0 ) ;
         diff = test[3*i+2] - out ;
         err += diff * diff ;
         }
      computed_err_ratio += err / (10 * nsamps) ;
      delete gate ;
      printf ( "\n  Ratio error = %8.5lf", computed_err_ratio / ndone );

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

     } // For all tries


   free ( x ) ;
   free ( xbad ) ;
   free ( xwild ) ;
   free ( test ) ;
   free ( computed_err_raw ) ;

   for (i=0 ; i<nmodels ; i++)
      delete models[i] ;
   free ( models ) ;

   _getch () ;
   return EXIT_SUCCESS ;
}

/******************************************************************************/
/*                                                                            */
/*  AFTERFAC - Test after-the-fact oracle                                     */
/*                                                                            */
/*  This uses an arbitrary number of univariate real predictors               */
/*  to make a combined univariate prediction.                                 */
/*  The number of categories may be specified,                                */
/*  though more than two risks overfitting.                                   */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "..\mlfn.h"

double unifrand () ;
double normal () ;
void qsortd ( int first , int last , double *data ) ;
void qsortdsi ( int first , int last , double *data , int *slave ) ;

static MLFN **models ;

/*
--------------------------------------------------------------------------------

   AfterFact - After-the-fact model combination via oracle based on outputs

   This class lets the model be an external reference.
   This is generally regarded as sloppy coding, but there is a reason here.
   By making the models external, we avoid the need to pass a typed identifier
   in the parameter lists.  An alternative method would be to define a parent
   "Model" class, but this would mean redefining the specific model to reflect
   its parentage.  The method shown here fits in best with the other code.
   Feel free to modify it as desired.

--------------------------------------------------------------------------------
*/

class AfterFact {

public:

   AfterFact ( int n , int ninputs , double *tset , int nmods , int ncat ) ;
   ~AfterFact () ;
   void predict ( double *input , double *output ) ;

private:
   int ncases ;       // Number of cases
   int nin ;          // Number of original case inputs
   int nmodels ;      // Number of models (nmods in constructor call)
   int ncats ;        // Number of categories
   double *outwork ;  // Work vector nmodels long holds outputs for a case
   double *thresh ;   // Nmodels by nfrac=ncats-1 array of fractile thresholds
   int *winners ;     // Index of winning model for each bin kept here
} ;

AfterFact::AfterFact (
   int n ,            // Number of training cases
   int ninputs ,      // Number of inputs
   double *tset ,     // Training cases, n by (nin+1)
   int nmods ,        // Number of models in 'models' array
   int ncat           // Number of fractile categories
   )
{

   int i, k, icase, imodel, ibin, nbins, index, klow, khigh, nthresh, ibest ;
   int *bins ;
   double *case_ptr, diff, x, best, *outptr, frac, *work ;
   double *outputs ;  // Ncases by nmodels vector saves all model outputs

   ncases = n ;
   nin = ninputs ;
   nmodels = nmods ;
   ncats = ncat ;        // Number of categories
   nthresh = ncats - 1 ; // Number of fractile thresholds

/*
   The number of bins is ncats ** nmodels.
*/

   nbins = 1 ;
   for (i=0 ; i<nmodels ; i++)
      nbins *= ncats ;

/*
   Allocate memory that is needed throughout the life of the AfterFact object
*/

   outwork = (double *) malloc ( nmodels * sizeof(double) ) ;
   thresh = (double *) malloc ( nmodels * nthresh * sizeof(double) ) ;
   winners = (int *) malloc ( nbins * sizeof(int) ) ;

/*
   Allocate scratch memory that is only needed during construction:
      outputs saves all model outputs
      work is scratch for fractile computation
      bins is a set of nbins vectors, each nmodels long
*/

   outputs = (double *) malloc ( ncases * nmodels * sizeof(double) ) ;
   work = (double *) malloc ( ncases * sizeof(double) ) ;
   bins = (int *) malloc ( nbins * nmodels * sizeof(int) ) ;

   memset ( bins , 0 , nbins * nmodels * sizeof(int) ) ;

/*
   Pass through the training set, invoking each model and saving outputs
*/

   for (icase=0 ; icase<ncases ; icase++) {
      case_ptr = tset + icase * (nin + 1) ; // Point to this case
      outptr = outputs + icase * nmodels ;
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( case_ptr , &outptr[imodel]  ) ;
      }

/*
   For each model, compute its output's fractile threshold
*/

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      for (icase=0 ; icase<ncases ; icase++)
         work[icase] = outputs[icase*nmodels+imodel] ;
      qsortd ( 0 , ncases-1 , work ) ;
      for (i=0 ; i<nthresh ; i++) {
         frac = (double) (i+1) / (double) ncats ;
         thresh[imodel*nthresh+i] = work[(int)(frac*(ncases-1))] ;
//         printf ( "\nModel %d thresh=%lf", imodel, thresh[imodel*nthresh+i] ) ; /*!!!!!!*/
         }
      }

/*
   Cumulate bin counts:
   For each case, determine which of the ncats ** nmodels bins the
   model outputs fall in.
   This points us to the nmodels-vector for that bin.
   Increment the element of this vector corresponding to the best model.
*/

   for (icase=0 ; icase<ncases ; icase++) {
      case_ptr = tset + icase * (nin + 1) ; // Point to this case
      outptr = outputs + icase * nmodels ;  // Corresponding model outputs
      ibin = 0 ;
      index = 1 ;
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         x = outptr[imodel] ;
         if (x <= thresh[imodel*nthresh])                   // Bottom bin
            k = 0 ;
         else if (x > thresh[imodel*nthresh+nthresh-1])     // Topmost bin
            k = nthresh ;
         else {  // Keep thresh[klow] < x <= thresh[khigh]  // Interior
            klow = 0 ;
            khigh = nthresh - 1 ;
            for (;;) {
               k = (klow + khigh) / 2 ;
               if (k == klow) {
                  k = khigh ;
                  break ;
                  }
               if (x <= thresh[imodel*nthresh+k])
                  khigh = k ;
               else 
                  klow = k ;
               }
            }

         ibin += k * index ;            // Locate within this model's layer
         index *= ncats ;               // Advance to next layer
         } // For all outputs (models)

// We have the bin.  Now find the best model for this case.

      for (imodel=0 ; imodel<nmodels ; imodel++) { // Determine best model
         diff = fabs ( outptr[imodel] - case_ptr[nin] ) ; // Pred minus true
         if ((imodel == 0)  ||  (diff < best)) {
            best = diff ;
            k = imodel ;
            }
         }

      ++bins[ibin*nmodels+k] ;    // Count this model's win within this bin
      } // For all cases

/*
   The bin counts are all cumulated.
   For each bin set, find the winning model.
*/

//   printf ( "\n(Bin,winner):" ) ; /*!!!!!!*/
   for (ibin=0 ; ibin<nbins ; ibin++) {  // For each bin
      k = 0 ;                            // Will point to winning model
      ibest = 0 ;                        // Best count
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         if (bins[ibin*nmodels+imodel] > ibest) {
            ibest = bins[ibin*nmodels+imodel] ;
            k = imodel ;
            }
         }
      winners[ibin] = k ;
//      printf ( " (%d %d)", ibin, k ) ; /*!!!!!!*/
      }

   free ( outputs ) ;
   free ( work ) ;
   free ( bins ) ;
//   _getch () ; /*!!!!!!*/
}

AfterFact::~AfterFact ()
{
   free ( outwork ) ;
   free ( thresh ) ;
   free ( winners ) ;
}

void AfterFact::predict ( double *input , double *output )
{
   int k, klow, khigh, imodel, ibin, index, nthresh ;
   double out ;

   nthresh = ncats - 1 ;

/*
   Invoke all models and find output of each.
   Then determine which bin the output set is in.
*/

   ibin = 0 ;
   index = 1 ;

//   printf ( "\nIn=(%.4lf %.4lf)", input[0], input[1] ) ; /*!!!!!!*/
   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &out ) ;
      outwork[imodel] = out ;
//      printf ( "  Out=%.4lf", out ) ; /*!!!!!!*/
      if (out <= thresh[imodel*nthresh])
         k = 0 ;
      else if (out > thresh[imodel*nthresh+nthresh-1])
         k = nthresh - 1 ;
      else {  // Keep thresh[klow] < out <= thresh[khigh] 
         klow = 0 ;
         khigh = nthresh - 1 ;
         for (;;) {
            k = (klow + khigh) / 2 ;
            if (k == klow) {
               k = khigh ;
               break ;
               }
            if (out <= thresh[imodel*nthresh+k])
               khigh = k ;
            else 
               klow = k ;
            }
         }
      ibin += k * index ;
      index *= ncats ;
      } // For all outputs (models)

   *output = outwork[winners[ibin]] ;
//   printf ( "   Pred bin=%d  winner=%d  out=%lf", ibin, winners[ibin], *output ) ; /*!!!!!!*/
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
   int i, ntries, itry, nsamps, nbins, imodel, nmodels, ndone ;
   double *x, *xbad, *xwild, *test, std ;
   double err, diff, out ;
   double *computed_err_raw ;
   double computed_err_afterfact ;
   AfterFact *afterfact ;

   int nhid = 2 ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 6) {
      printf (
         "\nUsage: AFTERFAC  nsamples  nbins  nmodels  ntries std" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nbins = atoi ( argv[2] ) ;
   nmodels = atoi ( argv[3] ) ;
   ntries = atoi ( argv[4] ) ;
   std = atof ( argv[5] ) ;
#else
   nsamps = 500 ;
   nbins = 4 ;
   nmodels = 3 ;
   ntries = 3 ;
   std = 0.5 ;
#endif

   if ((nsamps <= 0)  ||  (nbins <= 1)  ||  (nmodels <= 0)  ||  (ntries <= 0)
    || (std < 0.0)) {
      printf ( "\nUsage: AFTERFAC  nsamples  nbins  nmodels  ntries  std" ) ;
      exit ( 1 ) ;
      }

/*
   Allocate memory and initialize
*/

   models = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;

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
   if any, the output is occasionally wild.)
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


      afterfact = new AfterFact ( nsamps , 2 , x , nmodels , nbins ) ;
      err = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         afterfact->predict ( test + 3 * i , &out ) ;
         diff = test[3*i+2] - out ;
         err += diff * diff ;
         }
      computed_err_afterfact += err / (10 * nsamps) ;
      delete afterfact ;

/*
   Print results so far
*/

      err = 0.0 ;
      printf ( "\n\n\nDid%5d    Raw errors:", ndone ) ;
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         printf ( "  %.4lf", computed_err_raw[imodel] / ndone ) ;
         err += computed_err_raw[imodel] / ndone ;
         }
      printf ( "\n       Mean raw error = %8.5lf", err / nmodels ) ;
      printf ( "\n      AfterFact error = %8.5lf", computed_err_afterfact / ndone ) ;

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

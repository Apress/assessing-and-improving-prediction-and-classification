/******************************************************************************/
/*                                                                            */
/*  MULTPRED - Compare methods for combining multiple numeric predictors      */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "linreg.h"
#include "mlfn.h"
#include "grnn.h"

double unifrand () ;
double normal () ;

double powell ( int maxits , double critlim , double tol ,
   double (*criter) ( double * ) , int n , double *x , double ystart ,
   double *base , double *p0 , double *direc ) ;

static MLFN **puremodels ;   // Allocated and freed in main
static MLFN **bootmodels ;   // Ditto
static MLFN **models ;       // Main sets to one of above

/*
--------------------------------------------------------------------------------

   All of these routines let the model be an external reference.
   This is generally regarded as sloppy coding, but there is a reason here.
   By making the models external, we avoid the need to pass a typed identifier
   in the parameter lists.  An alternative method would be to define a parent
   "Model" class, but this would mean redefining the specific model to reflect
   its parentage.  The method shown here fits in best with the other code.
   Feel free to modify it as desired.

--------------------------------------------------------------------------------
*/

/*
--------------------------------------------------------------------------------

   Average - Compute the simple average of the predictions

--------------------------------------------------------------------------------
*/

class Average {

public:

   Average ( int n , int nin , double *tset , int nmods ) ;
   ~Average () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
} ;

Average::Average (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   nmodels = nmods ;
}

Average::~Average ()
{
}

void Average::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double outwork ;

   *output = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &outwork ) ;
      *output += outwork ; // Average output is mean across all models
      }

   *output /= nmodels ;
}

/*
--------------------------------------------------------------------------------

   Unconstrained - Compute the linear regression of the predictions

--------------------------------------------------------------------------------
*/

class Unconstrained {

public:

   Unconstrained ( int n , int nin , double *tset , int nmods ) ;
   ~Unconstrained () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   LinReg *linreg ;   // The linear regression object
   int nmodels ;      // Number of models (nmods in constructor call)
   double *coefs ;    // Computed coefficients here, constant last
} ;

Unconstrained::Unconstrained (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   int i, imodel ;
   double *casevec, *outs, *tptr ;

   nmodels = nmods ;

   casevec = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   outs = (double *) malloc ( n * sizeof(double) ) ;
   coefs = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;

   linreg = new LinReg ( n , (nmodels + 1) ) ;

   casevec[nmodels] = 1.0 ; // This is the regression constant term
   for (i=0 ; i<n ; i++) {  // Build the design matrix
      tptr = tset + i * (nin+1) ;  // This case is here
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( tptr , &casevec[imodel] ) ;
      linreg->add_case ( casevec ) ;
      outs[i] = tptr[nin] ; // Corresponding true value
      }

   linreg->solve ( 1.e-6 , outs , coefs ) ;

   free ( casevec ) ;
   free ( outs ) ;
}

Unconstrained::~Unconstrained ()
{
   free ( coefs ) ;
   delete linreg ;
}

void Unconstrained::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double out ;

   *output = coefs[nmodels] ;  // Regression constant

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &out ) ;
      *output += coefs[imodel] * out ;
      }
}

/*
--------------------------------------------------------------------------------

   Unbiased - Compute the optimal linear combination of the predictions
              subject to the constraints that the weights are all nonnegative
              and they sum to one.  This is appropriate for unbiased predictors.

--------------------------------------------------------------------------------
*/

class Unbiased {

public:

   Unbiased ( int n , int nin , double *tset , int nmods ) ;
   ~Unbiased () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *coefs ;    // Computed coefficients here
} ;

/*
   This local routine computes the optimization criterion, which is the
   squared error.  So that this routine can be passed in an optimizer's
   parameter list, it employs static variables for what it needs.
*/

static int unbiased_ncases ;
static int unbiased_nvars ;
static double *unbiased_x ;
static double *unbiased_y ;
static double *unbiased_work ;

double unbiased_crit ( double *wts )
{
   int i, j ;
   double sum, err, pred, *xptr, diff, penalty ;

   // Normalize weights to sum to one
   sum = 0.0 ;
   for (j=0 ; j<unbiased_nvars ; j++)
      sum += wts[j] ;
   if (sum < 1.e-60)   // Should almost never happen
      sum = 1.e-60 ;   // But be prepared to avoid division by zero
   for (j=0 ; j<unbiased_nvars ; j++)
      unbiased_work[j] = wts[j] / sum ;

   // Compute criterion
   err = 0.0 ;
   for (i=0 ; i<unbiased_ncases ; i++) {
      xptr = unbiased_x + i * unbiased_nvars ;  // Point to this case
      pred = 0.0 ;                              // Will cumulate prediction
      for (j=0 ; j<unbiased_nvars ; j++)        // For all model outputs
         pred += xptr[j] * unbiased_work[j] ;   // Weight them per call
      diff = pred - unbiased_y[i] ;             // Predicted minus true
      err += diff * diff ;                      // Cumulate squared error
      }

   penalty = 0.0 ;
   for (j=0 ; j<unbiased_nvars ; j++) {
      if (wts[j] < 0.0)
         penalty -= 1.e30 * wts[j] ;
      }

   return err + penalty ;
}

Unbiased::Unbiased (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   int i, imodel ;
   double *cases, *outs, *tptr, *base , *p0 , *direc, ystart, yend, sum ;

   nmodels = nmods ;

   cases = (double *) malloc ( n * nmodels * sizeof(double) ) ;
   outs = (double *) malloc ( n * sizeof(double) ) ;
   coefs = (double *) malloc ( nmodels * sizeof(double) ) ;
   base = (double *) malloc ( nmodels * sizeof(double) ) ;
   p0 = (double *) malloc ( nmodels * sizeof(double) ) ;
   direc = (double *) malloc ( nmodels * nmodels * sizeof(double) ) ;
   unbiased_work = (double *) malloc ( nmodels * sizeof(double) ) ;

/*
   Find and save each model's prediction for each case
*/

   for (i=0 ; i<n ; i++) {
      tptr = tset + i * (nin+1) ;  // This case is here
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( tptr , &cases[i*nmodels+imodel] ) ;
      outs[i] = tptr[nin] ; // Corresponding true value
      }

/*
   Optimize
*/

   unbiased_ncases = n ;
   unbiased_nvars = nmodels ;
   unbiased_x = cases ;
   unbiased_y = outs ;

   for (i=0 ; i<nmodels ; i++)
      coefs[i] = 1.0 / nmodels ;

   ystart = unbiased_crit ( coefs ) ;

   yend = powell ( 20 , 0.0 , 1.e-6 , unbiased_crit , nmodels , coefs , ystart ,
            base , p0 , direc ) ;

   // Normalize weights to sum to one
   sum = 0.0 ;
   for (i=0 ; i<nmodels ; i++)
      sum += coefs[i] ;
   if (sum < 1.e-60)   // Should almost never happen
      sum = 1.e-60 ;   // But be prepared to avoid division by zero
   for (i=0 ; i<nmodels ; i++)
      coefs[i] /= sum ;

   free ( cases ) ;
   free ( outs ) ;
   free ( base ) ;
   free ( p0 ) ;
   free ( direc ) ;
   free ( unbiased_work ) ;
}

Unbiased::~Unbiased ()
{
   free ( coefs ) ;
}

void Unbiased::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double out ;

   *output = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &out ) ;
      *output += coefs[imodel] * out ;
      }
}

/*
--------------------------------------------------------------------------------

   Biased - Compute the optimal linear combination of the predictions
            subject to the constraints that the weights are all nonnegative.
            A constant term is also included.
            This is appropriate for biased predictors.

--------------------------------------------------------------------------------
*/

class Biased {

public:

   Biased ( int n , int nin , double *tset , int nmods ) ;
   ~Biased () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *coefs ;    // Computed coefficients here
} ;

/*
   This local routine computes the optimization criterion, which is the
   squared error.  So that this routine can be passed in an optimizer's
   parameter list, it employs static variables for what it needs.
*/

static int biased_ncases ;
static int biased_nvars ;
static double *biased_x ;
static double *biased_y ;

double biased_crit ( double *wts )
{
   int i, j ;
   double err, pred, *xptr, diff, penalty ;

   err = 0.0 ;
   for (i=0 ; i<biased_ncases ; i++) {
      xptr = biased_x + i * biased_nvars ;  // Point to this case
      pred = wts[biased_nvars] ;            // Constant term
      for (j=0 ; j<biased_nvars ; j++)      // For all model outputs
         pred += xptr[j] * wts[j] ;         // Weight them per call
      diff = pred - biased_y[i] ;           // Predicted minus true
      err += diff * diff ;                  // Cumulate squared error
      }

   penalty = 0.0 ;
   for (j=0 ; j<biased_nvars ; j++) {
      if (wts[j] < 0.0)
         penalty -= 1.e30 * wts[j] ;
      }

   return err + penalty ;
}

Biased::Biased (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   int i, imodel ;
   double *cases, *outs, *tptr, *base , *p0 , *direc, ystart, yend ;

   nmodels = nmods ;

   cases = (double *) malloc ( n * nmodels * sizeof(double) ) ;
   outs = (double *) malloc ( n * sizeof(double) ) ;
   coefs = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   base = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   p0 = (double *) malloc ( (nmodels+1) * sizeof(double) ) ;
   direc = (double *) malloc ( (nmodels+1) * (nmodels+1) * sizeof(double) ) ;

/*
   Find and save each model's prediction for each case
*/

   for (i=0 ; i<n ; i++) {
      tptr = tset + i * (nin+1) ;  // This case is here
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( tptr , &cases[i*nmodels+imodel] ) ;
      outs[i] = tptr[nin] ; // Corresponding true value
      }

/*
   Optimize
*/

   biased_ncases = n ;
   biased_nvars = nmodels ;  // Don't include constant term
   biased_x = cases ;
   biased_y = outs ;

   for (i=0 ; i<nmodels ; i++)
      coefs[i] = 1.0 / nmodels ;
   coefs[nmodels] = 0.0 ;         // Constant term

   ystart = biased_crit ( coefs ) ;

   yend = powell ( 20 , 0.0 , 1.e-6 , biased_crit , nmodels+1 , coefs , ystart ,
            base , p0 , direc ) ;

   free ( cases ) ;
   free ( outs ) ;
   free ( base ) ;
   free ( p0 ) ;
   free ( direc ) ;
}

Biased::~Biased ()
{
   free ( coefs ) ;
}

void Biased::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double out ;

   *output = coefs[nmodels] ;  // Regression constant

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &out ) ;
      *output += coefs[imodel] * out ;
      }
}

/*
--------------------------------------------------------------------------------

   Weighted - Compute the variance-weighted average of the predictions

--------------------------------------------------------------------------------
*/

class Weighted {

public:

   Weighted ( int n , int nin , double *tset , int nmods ) ;
   ~Weighted () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   int nmodels ;      // Number of models (nmods in constructor call)
   double *coefs ;    // Computed coefficients here
} ;

Weighted::Weighted (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   int i, imodel ;
   double *tptr, out, diff, sum ;

   nmodels = nmods ;

   coefs = (double *) malloc ( nmodels * sizeof(double) ) ;

   for (i=0 ; i<nmodels ; i++)
      coefs[i] = 1.e-60 ;      // Theoretically zero, but prevent division by 0

   for (i=0 ; i<n ; i++) {  // Cumulate error variance of each model
      tptr = tset + i * (nin+1) ;  // This case is here
      for (imodel=0 ; imodel<nmodels ; imodel++) {
         models[imodel]->predict ( tptr , &out ) ;
         diff = out - tptr[nin] ;   // Predicted minus true
         coefs[imodel] += diff * diff ;
         }
      }
       
   sum = 0.0 ;                                   // Will sum normalizing factor
   for (imodel=0 ; imodel<nmodels ; imodel++) {
      coefs[imodel] = 1.0 / coefs[imodel] ;      // Weight is recip of variance
      sum += coefs[imodel] ;                     // Cumulate normalzier
      }

   for (imodel=0 ; imodel<nmodels ; imodel++)    // Normalize to unit sum
      coefs[imodel] /= sum ;
}

Weighted::~Weighted ()
{
   free ( coefs ) ;
}

void Weighted::numeric_predict ( double *input , double *output )
{
   int imodel ;
   double out ;

   *output = 0.0 ;

   for (imodel=0 ; imodel<nmodels ; imodel++) {
      models[imodel]->predict ( input , &out ) ;
      *output += coefs[imodel] * out ;
      }
}

/*
--------------------------------------------------------------------------------

   GenReg - Compute the General Regression of the predictions

--------------------------------------------------------------------------------
*/

class GenReg {

public:

   GenReg ( int n , int nin , double *tset , int nmods ) ;
   ~GenReg () ;
   void numeric_predict ( double *input , double *output ) ;

private:
   GRNN *grnn ;       // The GRNN object
   int nmodels ;      // Number of models (nmods in constructor call)
   double *work ;     // Work vector nmodels long
} ;

GenReg::GenReg (
   int n ,            // Number of training cases
   int nin ,          // Number of inputs
   double *tset ,     // Training cases, ntrain by (nin+nout) where nout=1 here
   int nmods          // Number of models in 'models' array
   )
{
   int i, imodel ;
   double *casevec, *tptr ;

   nmodels = nmods ;

   casevec = (double *) malloc ( (nmodels + 1) * sizeof(double) ) ;
   work = (double *) malloc ( nmodels * sizeof(double) ) ;

   grnn = new GRNN ( n , nmodels , 1 ) ;

   for (i=0 ; i<n ; i++) {  // Build the design matrix
      tptr = tset + i * (nin+1) ;     // This case is here
      for (imodel=0 ; imodel<nmodels ; imodel++)
         models[imodel]->predict ( tptr , &casevec[imodel] ) ;
      casevec[nmodels] = tptr[nin] ;  // Desired output
      grnn->add_case ( casevec ) ;
      }

   grnn->train () ;
   free ( casevec ) ;
}

GenReg::~GenReg ()
{
   delete grnn ;
   free ( work ) ;
}

void GenReg::numeric_predict ( double *input , double *output )
{
   int imodel ;

   for (imodel=0 ; imodel<nmodels ; imodel++)
      models[imodel]->predict ( input , &work[imodel] ) ;

   grnn->predict ( work , output ) ;
}


/*
--------------------------------------------------------------------------------

   Optional main to test it

   If there are four or more models, the fourth model is deliberately worthless.
   If there are five or more models, the fifth model is deliberately biased.

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, k, ntries, itry, nsamps, imodel, nmodels, divisor, ndone ;
   double *x, *xbad, *xbiased, *test ;
   double diff, var, std, temp, out ;
   double *computed_err_raw ;
   double computed_err_average ;
   double computed_err_unconstrained ;
   double computed_err_unbiased ;
   double computed_err_biased ;
   double computed_err_weighted ;
   double computed_err_bagged ;
   double computed_err_genreg ;
   Average *average ;
   Unconstrained *unconstrained ;
   Unbiased *unbiased ;
   Biased *biased ;
   Weighted *weighted ;
   GenReg *genreg ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf (
         "\nUsage: MULTPRED  nsamples  nmodels  ntries  var" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nmodels = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   var = atof ( argv[4] ) ;
#else
   nsamps = 4 ;
   nmodels = 2 ;
   ntries = 100 ;
   var = 0.1 ;
#endif

   if ((nsamps <= 0)  ||  (nmodels <= 0)  ||  (ntries <= 0)  ||  (var < 0.0)) {
      printf ( "\nUsage: MULTPRED  nsamples  nmodels  ntries  var" ) ;
      exit ( 1 ) ;
      }

   std = sqrt ( var ) ;

   divisor = 1 ;

/*
   Allocate memory and initialize
*/

   puremodels = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;
   bootmodels = (MLFN **) malloc ( nmodels * sizeof(MLFN *) ) ;

   for (i=0 ; i<nmodels ; i++) {
      puremodels[i] = new MLFN ( nsamps , 2 , 1 , 2 ) ;
      bootmodels[i] = new MLFN ( nsamps , 2 , 1 , 2 ) ;
      }

   x = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   xbad = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   xbiased = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * 3 * sizeof(double) ) ;
   computed_err_raw = (double *) malloc ( nmodels * sizeof(double) ) ;

   for (imodel=0 ; imodel<nmodels ; imodel++)
      computed_err_raw[imodel] = 0.0 ;

   computed_err_average = 0.0 ;
   computed_err_unconstrained = 0.0 ;
   computed_err_unbiased = 0.0 ;
   computed_err_biased = 0.0 ;
   computed_err_weighted = 0.0 ;
   computed_err_bagged = 0.0 ;
   computed_err_genreg = 0.0 ;

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
   if any, the output is deliberately biased.)
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[3*i] = normal () ;
         x[3*i+1] = normal () ;
         x[3*i+2] = sin ( x[3*i] ) - x[3*i+1] * x[3*i+1] + std * normal () ;
         }

      if (nmodels >= 4) {
         for (i=0 ; i<nsamps ; i++) {
            xbad[3*i] = x[3*i] ;
            xbad[3*i+1] = x[3*i+1] ;
            xbad[3*i+2] = normal () ;         // Results in a worthless model
            }
         }

      if (nmodels >= 5) {
         for (i=0 ; i<nsamps ; i++) {
            xbiased[3*i] = x[3*i] ;
            xbiased[3*i+1] = x[3*i+1] ;
            xbiased[3*i+2] = x[3*i+2] + 1.0 ; // Results in a biased model
            }
         }

      for (i=0 ; i<10*nsamps ; i++) {
         test[3*i] = normal () ;
         test[3*i+1] = normal () ;
         test[3*i+2] = sin ( test[3*i] ) - test[3*i+1] * test[3*i+1] + std * normal () ;
         }

      for (imodel=0 ; imodel<nmodels ; imodel++) {
         puremodels[imodel]->reset () ;
         bootmodels[imodel]->reset () ;
         if (imodel == 3) {
            for (i=0 ; i<nsamps ; i++) {
               puremodels[imodel]->add_case ( xbad + 3 * i ) ;
               k = (int) (unifrand() * nsamps) ;
               if (k >= nsamps)
                  k = nsamps - 1 ;
               bootmodels[imodel]->add_case ( x + 3 * k ) ; // Bagging gets good data only
               }
            }
         else if (imodel == 4) {
            for (i=0 ; i<nsamps ; i++) {
               puremodels[imodel]->add_case ( xbiased + 3 * i ) ;
               k = (int) (unifrand() * nsamps) ;
               if (k >= nsamps)
                  k = nsamps - 1 ;
               bootmodels[imodel]->add_case ( x + 3 * k ) ; // Bagging gets good data only
               }
            }
         else {
            for (i=0 ; i<nsamps ; i++) {
               puremodels[imodel]->add_case ( x + 3 * i ) ;
               k = (int) (unifrand() * nsamps) ;
               if (k >= nsamps)
                  k = nsamps - 1 ;
               bootmodels[imodel]->add_case ( x + 3 * k ) ;
               }
            }

         puremodels[imodel]->train () ;
         bootmodels[imodel]->train () ;
   
         temp = 0.0 ;
         for (i=0 ; i<10*nsamps ; i++) {
            puremodels[imodel]->predict ( test + 3 * i , &out ) ;
            diff = out - test[3*i+2] ;
            temp += diff * diff ;
            }
         computed_err_raw[imodel] += temp / (10 * nsamps) ;
         }

/*
   Average
*/

      models = puremodels ;
      average = new Average ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         average->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_average += temp / (10 * nsamps) ;
      delete average ;

/*
   Unconstrained
*/

      models = puremodels ;
      unconstrained = new Unconstrained ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         unconstrained->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_unconstrained += temp / (10 * nsamps) ;
      delete unconstrained ;

/*
   Unbiased
*/

      models = puremodels ;
      unbiased = new Unbiased ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         unbiased->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_unbiased += temp / (10 * nsamps) ;
      delete unbiased ;

/*
   Biased
*/

      models = puremodels ;
      biased = new Biased ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         biased->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_biased += temp / (10 * nsamps) ;
      delete biased ;

/*
   Weighted
*/

      models = puremodels ;
      weighted = new Weighted ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         weighted->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_weighted += temp / (10 * nsamps) ;
      delete weighted ;

/*
   Bagged
*/

      models = bootmodels ;
      average = new Average ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         average->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_bagged += temp / (10 * nsamps) ;
      delete average ;

/*
   GenReg
*/

      models = puremodels ;
      genreg = new GenReg ( nsamps , 2 , x , nmodels ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         genreg->numeric_predict ( test + 3 * i , &out ) ;
         diff = out - test[3*i+2] ;
         temp += diff * diff ;
         }
      computed_err_genreg += temp / (10 * nsamps) ;
      delete genreg ;

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
      printf ( "\n  Unconstrained error = %8.5lf", computed_err_unconstrained / ndone ) ;
      printf ( "\n       Unbiased error = %8.5lf", computed_err_unbiased / ndone ) ;
      printf ( "\n         Biased error = %8.5lf", computed_err_biased / ndone ) ;
      printf ( "\n       Weighted error = %8.5lf", computed_err_weighted / ndone ) ;
      printf ( "\n         Bagged error = %8.5lf", computed_err_bagged / ndone ) ;
      printf ( "\n         GenReg error = %8.5lf", computed_err_genreg / ndone ) ;

      if (_kbhit ()) {
         if (_getch() == 27)
            break ;
         }

     } // For all tries


   return EXIT_SUCCESS ;
}

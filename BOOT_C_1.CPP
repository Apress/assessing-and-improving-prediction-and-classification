/******************************************************************************/
/*                                                                            */
/*  BOOT_C_1 - Compare resampling methods for estimating error variance       */
/*             This uses a numeric prediction problem.                        */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

#include "linreg.h"

double unifrand () ;
double normal () ;

static LinReg *linreg_n, *linreg_nm1 ;   // Allocated and freed in main


/*
--------------------------------------------------------------------------------

   Given training and test sets, train and test a model

   In order to keep this routine general, we do not worry about passing
   as parameters special things like the LinReg object and work areas.
   Custom versions of this routine should make them static and allocate/free
   them in an outer position such as the main program.

   This avoids the need for the library routines that call a special function
   of this name to deal with changeable parameter lists.  What a pain
   that would be!  We want to keep this routine's parameter list universal.

--------------------------------------------------------------------------------
*/

static LinReg *linreg ;   // Allocated and freed in main
double *work_npredp1 ;    // Ditto, work vector npred+1 long
double *work_ntrain ;     // Ditto, work vector ntrain long

void train_test (
   int ntrain ,           // Number of training cases
   int ntest ,            // Number of test cases
   int npred ,            // Number of predictor variables
   double *train ,        // Ntrain by npred+1 matrix of predictors followed by predicted
   double *test ,         // Above is training set;  This is test set; May be same
   double *predicted      // Output of 'ntest' test set predictions
   )
{
   int i, j ;
   double *tptr ;

   linreg->reset() ;

/*
   Build the design matrix.  Place 1.0 in the last column so that the
   regression function also has a constant term.
   Simultaneously build the right hand side: the true predicted values.
*/

   work_npredp1[npred] = 1.0 ;   // This is the regression constant term
   for (i=0 ; i<ntrain ; i++) {  // Build the design matrix
      tptr = train + i * (npred+1) ;  // This case is here
      memcpy ( work_npredp1 , tptr , npred * sizeof(double) ) ;
      linreg->add_case ( work_npredp1 ) ;
      work_ntrain[i] = tptr[npred] ; // Corresponding true value
      }

   linreg->solve ( 1.e-8 , work_ntrain , work_npredp1 ) ;

/*
   Compute test set predictions
*/

   for (i=0 ; i<ntest ; i++) {
      tptr = test + i * (npred+1) ;         // This case is here
      predicted[i] = work_npredp1[npred] ;  // Constant is here
      for (j=0 ; j<npred ; j++)             // Compute predicted value
         predicted[i] += work_npredp1[j] * tptr[j] ; // Using regression equation
      }
}

/*
--------------------------------------------------------------------------------

   This function computes the error of a prediction.
   It is called q out of deference to common convention.

--------------------------------------------------------------------------------
*/

double q ( double trueval , double predicted )
{
   double diff ;

   diff = trueval - predicted ;
   return diff * diff ;
}

/*
--------------------------------------------------------------------------------

   cross_validation - Use cross validation to estimate error variance

--------------------------------------------------------------------------------
*/

void cross_validation (
   int n ,              // Number of cases in sample
   int npred ,          // Number of predictor variables
   double *data ,       // The n by npred+1 dataset of predictors followed by predicted
   void (*tt) (         // Train and test model
      int ntrain ,          // Number of training cases
      int ntest ,           // Number of test cases
      int npred ,           // Number of predictor variables
      double *train ,       // Ntrain by npred+1 matrix of predictors followed by predicted
      double *test ,        // Above is training set;  This is test set
      double *predicted ) , // Output of test set predictions
   double *mean_err     // Output of error estimate
   )
{
   int i, j ;
   double temp, err, *excluded, *test ;

   *mean_err = 0.0 ;        // Will cumulate mean error here
   test = data + (n - 1) * (npred + 1) ;  // Last case in dataset

   for (i=0 ; i<n ; i++) {            // Exclude one case at a time
      excluded = data + i * (npred + 1) ;  // Excluded case

      for (j=0 ; j<=npred ; j++) {    // Swap excluded case to test spot at end
         temp = test[j] ;
         test[j] = excluded[j] ;
         excluded[j] = temp ;
         }

      tt ( n-1 , 1 , npred , data , test , &temp ) ;

      err = q ( test[npred] , temp ) ;
      *mean_err += err ;              // Cumulate for mean error

      for (j=0 ; j<=npred ; j++) {    // Swap to restore original order
         temp = test[j] ;
         test[j] = excluded[j] ;
         excluded[j] = temp ;
         }
      }

   *mean_err /= n ;
}

/*
--------------------------------------------------------------------------------

   bootstrap - Use ordinary bootstrap to estimate error variance

--------------------------------------------------------------------------------
*/

void bootstrap (
   int n ,              // Number of cases in sample
   int npred ,          // Number of predictor variables
   double *data ,       // The n by npred+1 dataset of predictors followed by predicted
   int nboot ,          // Number of bootstrap replications
   void (*tt) (         // Train and test model
      int ntrain ,          // Number of training cases
      int ntest ,           // Number of test cases
      int npred ,           // Number of predictor variables
      double *train ,       // Ntrain by npred+1 matrix of predictors followed by predicted
      double *test ,        // Above is training set;  This is test set
      double *predicted ) , // Output of test set predictions
   double *mean_err ,   // Output of error estimate
   double *bootsamp ,   // Work area n * (npred+1) long
   double *predicted ,  // Work area n long
   int *count           // Work area n long
   )
{
   int i, rep, k ;
   double err, apparent, excess, *tptr ;

   excess = 0.0 ;                       // Cumulates excess error

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      memset ( count , 0 , n * sizeof(int) ) ;

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         memcpy ( bootsamp + i * (npred+1) ,  // Put case in bootstrap sample
                  data + k * (npred+1) , (npred+1) * sizeof(double) ) ;
         ++count[k] ;                   // Count inclusion of this case
         }

      tt ( n , n , npred , bootsamp , data , predicted ) ; // Train and predict

      for (i=0 ; i<n ; i++) {           // Compute mean error.
         tptr = data + i * (npred+1) ;  // This case is here
         err = q ( tptr[npred] , predicted[i] ) ;
         excess += (1.0 - count[i]) * err ;
         }
      } // For all bootstrap replications

   excess /= n * nboot ;  // Computed excess error is grand mean

/*
   Compute apparent error.  Add it to excess to get population error estimate.
*/

   tt ( n , n , npred , data , data , predicted ) ; // Train and predict

   apparent = 0.0 ;
   for (i=0 ; i<n ; i++) {
      tptr = data + i * (npred+1) ;        // This case is here
      err = q ( tptr[npred] , predicted[i] ) ;
      apparent += err ;
      }

   apparent /= n ;

   *mean_err = apparent + excess ;
}

/*
--------------------------------------------------------------------------------

   E0 - Use Efron's E0 bootstrap to estimate error variance

   For clarity, this implementation lets tt() predict for the entire dataset,
   even though only predictions for cases not used in the bootstrap sample
   are needed.  If speed is critical and prediction is slow, this could
   be changed easily.

--------------------------------------------------------------------------------
*/

void E0 (
   int n ,              // Number of cases in sample
   int npred ,          // Number of predictor variables
   double *data ,       // The n by npred+1 dataset of predictors followed by predicted
   int nboot ,          // Number of bootstrap replications
   void (*tt) (         // Train and test model
      int ntrain ,          // Number of training cases
      int ntest ,           // Number of test cases
      int npred ,           // Number of predictor variables
      double *train ,       // Ntrain by npred+1 matrix of predictors followed by predicted
      double *test ,        // Above is training set;  This is test set
      double *predicted ) , // Output of test set predictions
   double *mean_err ,   // Output of error estimate
   double *bootsamp ,   // Work area n * (npred+1) long
   double *predicted ,  // Work area n long
   int *count           // Work area n long
   )
{
   int i, rep, k, ntot ;
   double err, *tptr ;

   *mean_err = 0.0 ;                    // Cumulates excess error
   ntot = 0 ;

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      memset ( count , 0 , n * sizeof(int) ) ;

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         memcpy ( bootsamp + i * (npred+1) ,  // Put case in bootstrap sample
                  data + k * (npred+1) , (npred+1) * sizeof(double) ) ;
         ++count[k] ;                   // Count inclusion of this case
         }

      tt ( n , n , npred , bootsamp , data , predicted ) ; // Train and predict

      for (i=0 ; i<n ; i++) {           // Compute mean error.
         if (count[i])
            continue ;
         tptr = data + i * (npred+1) ;  // This case is here
         err = q ( tptr[npred] , predicted[i] ) ;
         *mean_err += err ;
         ++ntot ;
         }
      } // For all bootstrap replications

   if (ntot)
      *mean_err /= ntot ;
}

/*
--------------------------------------------------------------------------------

   E623 - Use Efron's E632 bootstrap to estimate error variance

--------------------------------------------------------------------------------
*/

void E632 (
   int n ,              // Number of cases in sample
   int npred ,          // Number of predictor variables
   double *data ,       // The n by npred+1 dataset of predictors followed by predicted
   int nboot ,          // Number of bootstrap replications
   void (*tt) (         // Train and test model
      int ntrain ,          // Number of training cases
      int ntest ,           // Number of test cases
      int npred ,           // Number of predictor variables
      double *train ,       // Ntrain by npred+1 matrix of predictors followed by predicted
      double *test ,        // Above is training set;  This is test set
      double *predicted ) , // Output of test set predictions
   double *mean_err ,   // Output of error estimate
   double *bootsamp ,   // Work area n * (npred+1) long
   double *predicted ,  // Work area n long
   int *count           // Work area n long
   )
{
   int i ;
   double apparent, *tptr ;

   E0 ( n , npred , data , nboot , tt , mean_err ,
        bootsamp , predicted , count ) ;

/*
   Compute apparent error.
   E632 = .632 E0  +  .368 Apparent
*/

   tt ( n , n , npred , data , data , predicted ) ; // Train and predict

   apparent = 0.0 ;
   for (i=0 ; i<n ; i++) {
      tptr = data + i * (npred+1) ;        // This case is here
      apparent += q ( tptr[npred] , predicted[i] ) ;
      }

   apparent /= n ;

   *mean_err = 0.632 * *mean_err  +  0.368 * apparent ;
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
   int i, ntries, itry, nsamps, nboot, divisor, ndone, *count ;
   double *x, *test, *bootsamp, *predicted, err, diff, var, std, temp, *tptr ;
   double *computed_err_cv, *computed_err_boot ;
   double *computed_err_E0, *computed_err_E632 ;
   double sum_observed_error, mean_computed_err, var_computed_err ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf (
         "\nUsage: BOOT_C_1  nsamples  nboot  ntries  var" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   var = atof ( argv[4] ) ;
#else
   nsamps = 4 ;
   nboot = 1000 ;
   ntries = 1000 ;
   var = 0.0 ;
#endif

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)  ||  (var < 0.0)) {
      printf ( "\nUsage: BOOT_C_1  nsamples  nboot  ntries  var" ) ;
      exit ( 1 ) ;
      }

   std = sqrt ( var ) ;

   divisor = 1000000 / (nsamps * nboot) ;  // This is for progress reports only
   if (divisor < 2)
      divisor = 2 ;

/*
   Allocate memory and initialize
*/

   linreg_n = new LinReg ( nsamps , 3 ) ;      // train_test() will need this
   linreg_nm1 = new LinReg ( nsamps-1 , 3 ) ;  // Ditto

   x = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   test = (double *) malloc ( 10 * nsamps * 3 * sizeof(double) ) ;
   computed_err_cv = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_err_boot = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_err_E0 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_err_E632 = (double *) malloc ( ntries * sizeof(double) ) ;
   work_npredp1 = (double *) malloc ( 3 * sizeof(double) ) ;
   work_ntrain = (double *) malloc ( nsamps * sizeof(double) ) ;
   bootsamp = (double *) malloc ( nsamps * 3 * sizeof(double) ) ;
   predicted = (double *) malloc ( 10 * nsamps * sizeof(double) ) ;
   count = (int *) malloc ( nsamps * sizeof(int) ) ;

/*
   Main outer loop does all tries
*/

   sum_observed_error = 0.0 ;   // For comparison purposes

   for (itry=0 ; itry<ntries ; itry++) {

/*
   Generate the data.
   The model is Y = X1 - X2 + error
   We use x as the dataset for all resampling algorithms.
   The other dataset, test, is used only to keep track of the observed
   error of the model to give us a basis of comparison.
*/

      for (i=0 ; i<nsamps ; i++) {
         x[3*i] = normal () ;
         x[3*i+1] = normal () ;
         x[3*i+2] = x[3*i] - x[3*i+1] + std * normal () ;
         }

      for (i=0 ; i<10*nsamps ; i++) {
         test[3*i] = normal () ;
         test[3*i+1] = normal () ;
         test[3*i+2] = test[3*i] - test[3*i+1] + std * normal () ;
         }

/*
   Train a model with this data and test it on an independent test set.
   This gives us a basis of comparison for the resampling methods.
*/

      linreg = linreg_n ;
      train_test ( nsamps , 10 * nsamps , 2 , x , test , predicted ) ;
      temp = 0.0 ;
      for (i=0 ; i<10*nsamps ; i++) {
         tptr = test + 3 * i ;  // This case is here
         err = q ( tptr[2] , predicted[i] ) ;
         temp += err ;
         }

      sum_observed_error += temp / (10 * nsamps) ;

/*
   Do the resampling methods
*/

      linreg = linreg_nm1 ;
      cross_validation ( nsamps , 2 , x , train_test ,
                         &computed_err_cv[itry] ) ;

      linreg = linreg_n ;
      bootstrap ( nsamps , 2 , x , nboot , train_test ,
                  &computed_err_boot[itry] , bootsamp , predicted , count ) ;

      E0 ( nsamps , 2 , x , nboot , train_test ,
           &computed_err_E0[itry] , bootsamp , predicted , count ) ;

      E632 ( nsamps , 2 , x , nboot , train_test ,
             &computed_err_E632[itry] , bootsamp , predicted , count ) ;

/*
   Periodically stop and print results for user
*/

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.

         ndone = itry + 1 ;           // This many tries done (and in arrays)
         printf ( "\n\n\nDid %d   Observed error = %.5lf",
                  ndone, sum_observed_error / ndone ) ;

/*
   Process cross validation test
*/

         mean_computed_err = 0.0 ;
         var_computed_err = 0.0 ;
         for (i=0 ; i<ndone ; i++)
            mean_computed_err += computed_err_cv[i] ;
         mean_computed_err /= ndone ;
         for (i=0 ; i<ndone ; i++) {
            diff = computed_err_cv[i] - mean_computed_err ;
            var_computed_err += diff * diff ;
            }
         var_computed_err /= ndone ;
         printf ( "\n  CV: computed error  mean=%10.5lf      std=%10.5lf",
            mean_computed_err, sqrt ( var_computed_err ) ) ;

/*
   Process bootstrap test
*/

         mean_computed_err = 0.0 ;
         var_computed_err = 0.0 ;
         for (i=0 ; i<ndone ; i++)
            mean_computed_err += computed_err_boot[i] ;
         mean_computed_err /= ndone ;
         for (i=0 ; i<ndone ; i++) {
            diff = computed_err_boot[i] - mean_computed_err ;
            var_computed_err += diff * diff ;
            }
         var_computed_err /= ndone ;
         printf ( "\nBOOT: computed error  mean=%10.5lf      std=%10.5lf",
            mean_computed_err, sqrt ( var_computed_err ) ) ;

/*
   Process E0 test
*/

         mean_computed_err = 0.0 ;
         var_computed_err = 0.0 ;
         for (i=0 ; i<ndone ; i++)
            mean_computed_err += computed_err_E0[i] ;
         mean_computed_err /= ndone ;
         for (i=0 ; i<ndone ; i++) {
            diff = computed_err_E0[i] - mean_computed_err ;
            var_computed_err += diff * diff ;
            }
         var_computed_err /= ndone ;
         printf ( "\n  E0: computed error  mean=%10.5lf      std=%10.5lf",
            mean_computed_err, sqrt ( var_computed_err ) ) ;

/*
   Process E632 test
*/

         mean_computed_err = 0.0 ;
         var_computed_err = 0.0 ;
         for (i=0 ; i<ndone ; i++)
            mean_computed_err += computed_err_E632[i] ;
         mean_computed_err /= ndone ;
         for (i=0 ; i<ndone ; i++) {
            diff = computed_err_E632[i] - mean_computed_err ;
            var_computed_err += diff * diff ;
            }
         var_computed_err /= ndone ;
         printf ( "\nE632: computed error  mean=%10.5lf      std=%10.5lf",
            mean_computed_err, sqrt ( var_computed_err ) ) ;

         } // If pausing to update statistics

      if ((itry % 10) == 1) {
         if (_kbhit ()) {
            if (_getch() == 27)
               break ;
            }
         }

     } // For all tries


   return EXIT_SUCCESS ;
}

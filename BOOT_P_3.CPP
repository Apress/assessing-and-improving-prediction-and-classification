/******************************************************************************/
/*                                                                            */
/*  BOOT_P_3 - Bootstrap estimate of bias and variance for regression coef    */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

double unifrand () ;
double normal () ;


/*
--------------------------------------------------------------------------------

   Compute the parameter

--------------------------------------------------------------------------------
*/

double param_beta ( int n , double *x , double *y )
{
   int i ;
   double xmean, ymean, xdif, ydif, xvar, covar ;

   xmean = ymean = 0.0 ;
   for (i=0 ; i<n ; i++) {
      xmean += x[i] ;
      ymean += y[i] ;
      }

   xmean /= n ;
   ymean /= n ;

   xvar = covar = 0.0 ;
   for (i=0 ; i<n ; i++) {
      xdif = x[i] - xmean ;
      ydif = y[i] - ymean ;
      xvar += xdif * xdif ;
      covar += xdif * ydif ;
      }

   if (xvar != 0.0)
      return covar / xvar ;
   return 0.0 ;
}

/*
--------------------------------------------------------------------------------

   boot_bias_var - Compute bias and variance of parameter, bivariate data

--------------------------------------------------------------------------------
*/

void boot_bias_var (
   int n ,              // Number of cases in sample
   double *x ,          // Independent variable in sample
   double *y ,          // Dependent variable in sample
   double (*user_t) (int , double * , double *) , // Compute parameter
   int nboot ,          // Number of bootstrap replications
   double *rawstat ,    // Raw statistic of sample, theta-hat
   double *bias ,       // Output of bias estimate
   double *var ,        // Output of variance estimate
   double *xwork ,      // Work area n long
   double *ywork ,      // Work area n long
   double *work2        // Work area nboot long
   )
{
   int i, rep, k ;
   double stat, mean, variance, diff ;

   mean = 0.0 ;

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         xwork[i] = x[k] ;              // Put bootstrap sample in work
         ywork[i] = y[k] ;
         }

      stat = user_t ( n , xwork , ywork ) ; // Evaluate estimator for this rep
      work2[rep] = stat ;               // Enables more accurate variance
      mean += stat ;                    // Cumulate theta-hat star dot
      }

   mean /= nboot ;
   variance = 0.0 ;
   for (rep=0 ; rep<nboot ; rep++) {    // Cumulate variance
      diff = work2[rep] - mean ;
      variance += diff * diff ;
      }

   stat = user_t ( n , x , y ) ;        // This is the final but biased estimate
   *rawstat = stat ;
   *bias = mean - stat ;
   *var = variance / (nboot - 1) ;
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
   int i, ntries, itry, nsamps, nboot, divisor, ndone ;
   double beta, *x, *y, diff, *xwork, *ywork, *work2 ;
   double *computed_param, *computed_bias, *computed_var ;
   double mean_computed_param, var_computed_param ;
   double mean_computed_bias, var_computed_bias, mean_computed_var ;

/*
   Process command line parameters
*/

   if (argc != 5) {
      printf (
         "\nUsage: BOOT_P_3  nsamples  nboot  ntries  beta" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   beta = atof ( argv[4] ) ;

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)) {
      printf ( "\nUsage: BOOT_P_3  nsamples  nboot  ntries  beta" ) ;
      exit ( 1 ) ;
      }


   divisor = 1000000 / (nsamps * nboot) ;  // This is for progress reports only
   if (divisor < 2)
      divisor = 2 ;

/*
   Allocate memory and initialize
*/

   x = (double *) malloc ( nsamps * sizeof(double) ) ;
   y = (double *) malloc ( nsamps * sizeof(double) ) ;
   xwork = (double *) malloc ( nsamps * sizeof(double) ) ;
   ywork = (double *) malloc ( nsamps * sizeof(double) ) ;
   work2 = (double *) malloc ( nboot * sizeof(double) ) ;
   computed_param = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_bias = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_var = (double *) malloc ( ntries * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   for (itry=0 ; itry<ntries ; itry++) {

      if ((itry % divisor) == 0)
         printf ( "\n\n\nTry %d", itry ) ;

      for (i=0 ; i<nsamps ; i++) {
         x[i] = normal () ;
         y[i] = beta * x[i] + 0.2 * normal () ;
         }

      boot_bias_var ( nsamps , x , y , param_beta , nboot ,
                      &computed_param[itry] , &computed_bias[itry] ,
                      &computed_var[itry] , xwork , ywork , work2 ) ;

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.
         ndone = itry + 1 ;           // This many tries done (and in arrays)

         mean_computed_param = 0.0 ;
         mean_computed_bias = 0.0 ;
         mean_computed_var = 0.0 ;
         var_computed_param = 0.0 ;
         var_computed_bias = 0.0 ;
         for (i=0 ; i<ndone ; i++) {
            mean_computed_param += computed_param[i] ;
            mean_computed_bias += computed_bias[i] ;
            mean_computed_var += computed_var[i] ;
            }
         mean_computed_param /= ndone ;  // Mean across tries of theta-hat
         mean_computed_bias /= ndone ;   // Mean of bias estimates
         mean_computed_var /= ndone ;    // Mean of variance estimates
         for (i=0 ; i<ndone ; i++) {
            diff = computed_param[i] - mean_computed_param ;
            var_computed_param += diff * diff ;
            diff = computed_bias[i] - mean_computed_bias ;
            var_computed_bias += diff * diff ;
            }
         var_computed_param /= ndone ;   // Variance of theta-hat
         var_computed_bias /= ndone ;    // Variance of bias estimate
         printf ( "\nComputed parameter (theta-hat)  mean=%.5lf  variance=%.5lf",
            mean_computed_param, var_computed_param ) ;
         printf ( "\nMean boot bias=%.5lf (its std=%.5lf)   Mean boot var=%.5lf",
            mean_computed_bias, sqrt(var_computed_bias), mean_computed_var ) ;
         }

      if ((itry % 10) == 1) {
         if (_kbhit ()) {
            if (_getch() == 27)
               break ;
            }
         }

     } // For all tries


   return EXIT_SUCCESS ;
}

/******************************************************************************/
/*                                                                            */
/*  BOOT_P_1 - Bootstrap estimate of bias and variance when s != t            */
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
void qsortd ( int istart , int istop , double *x ) ;


/*
--------------------------------------------------------------------------------

   Compute the parameter

--------------------------------------------------------------------------------
*/

double param_mean ( int n , double *x )
{
   int i ;
   double mean ;

   mean = 0.0 ;
   for (i=0 ; i<n ; i++)
      mean += x[i] ;

   return mean / n ;
}

double param_median ( int n , double *x )
{
   qsortd ( 0 , n-1 , x ) ;

   if (n % 2)
      return x[n/2] ;
   else 
      return 0.5 * (x[n/2-1] + x[n/2]) ;
}

/*
--------------------------------------------------------------------------------

   boot_bias_var - Compute bias and variance of parameter (general s != t case)

--------------------------------------------------------------------------------
*/

void boot_bias_var (
   int n ,              // Number of cases in sample
   double *data ,       // The sample
   double (*user_s) (int , double *) , // Compute param being bootstrapped, s
   double (*user_t) (int , double *) , // Compute true param, t
   int nboot ,          // Number of bootstrap replications
   double *rawstat ,    // Raw statistic of sample, theta-hat
   double *bias ,       // Output of bias estimate
   double *var ,        // Output of variance estimate
   double *work ,       // Work area n long
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
         work[i] = data[k] ;            // Put bootstrap sample in work
         }

      stat = user_s ( n , work ) ;      // Evaluate estimator for this boot rep
      work2[rep] = stat ;               // Enables more accurate variance
      mean += stat ;                    // Cumulate theta-hat star dot
      }

   mean /= nboot ;
   variance = 0.0 ;
   for (rep=0 ; rep<nboot ; rep++) {    // Cumulate variance
      diff = work2[rep] - mean ;
      variance += diff * diff ;
      }

   *rawstat = user_s ( n , data ) ;     // This is the final but biased estimate
   *bias = mean - user_t ( n , data ) ;
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
   double *x, diff, *work, *work2 ;
   double *computed_param_1, *computed_bias_1, *computed_var_1 ;
   double *computed_param_2, *computed_bias_2, *computed_var_2 ;
   double mean_computed_param, var_computed_param ;
   double mean_computed_bias, var_computed_bias, mean_computed_var ;

/*
   Process command line parameters
*/

   if (argc != 4) {
      printf (
         "\nUsage: BOOT_P_1  nsamples  nboot  ntries" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)) {
      printf ( "\nUsage: BOOT_P_1  nsamples  nboot  ntries" ) ;
      exit ( 1 ) ;
      }


   divisor = 1000000 / (nsamps * nboot) ;  // This is for progress reports only
   if (divisor < 2)
      divisor = 2 ;

/*
   Allocate memory and initialize
*/

   x = (double *) malloc ( nsamps * sizeof(double) ) ;
   work = (double *) malloc ( nsamps * sizeof(double) ) ;
   work2 = (double *) malloc ( nboot * sizeof(double) ) ;
   computed_param_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_bias_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_var_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_param_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_bias_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   computed_var_2 = (double *) malloc ( ntries * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   for (itry=0 ; itry<ntries ; itry++) {

      if ((itry % divisor) == 0)
         printf ( "\n\n\nTry %d", itry ) ;

/*
   This is the first of two tests.

   This test uses the mean of the sample to estimate the mean of the population.
   So this is actually a plug-in statistic, and a silly one at that since
   we already know that the mean is unbiased and has an optimal variance
   estimate: the sample variance divided by n-1.
   But this test is a useful verification that all is going as planned.
   It also provides for a comparison with BOOT_P_2 later.
*/

      for (i=0 ; i<nsamps ; i++)
         x[i] = normal () ;

      boot_bias_var ( nsamps , x , param_mean , param_mean , nboot ,
                      &computed_param_1[itry] , &computed_bias_1[itry] ,
                      &computed_var_1[itry] , work , work2 ) ;

/*
   This is the second of two tests.

   It uses the sample median to estimate the population mean.
   In this example, the population is defined by exponentiating a normal
   random variable.  This distribution obviously has a long right tail and
   short left tail.
*/

      for (i=0 ; i<nsamps ; i++)
         x[i] = exp ( x[i] ) ;

      boot_bias_var ( nsamps , x , param_median , param_mean , nboot ,
                      &computed_param_2[itry] , &computed_bias_2[itry] ,
                      &computed_var_2[itry] , work , work2 ) ;

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.
         ndone = itry + 1 ;           // This many tries done (and in arrays)

/*
   Process test 1 of 2
*/

         mean_computed_param = 0.0 ;
         mean_computed_bias = 0.0 ;
         mean_computed_var = 0.0 ;
         var_computed_param = 0.0 ;
         var_computed_bias = 0.0 ;
         for (i=0 ; i<ndone ; i++) {
            mean_computed_param += computed_param_1[i] ;
            mean_computed_bias += computed_bias_1[i] ;
            mean_computed_var += computed_var_1[i] ;
            }
         mean_computed_param /= ndone ;  // Mean across tries of theta-hat
         mean_computed_bias /= ndone ;   // Mean of bias estimates
         mean_computed_var /= ndone ;    // Mean of variance estimates
         for (i=0 ; i<ndone ; i++) {
            diff = computed_param_1[i] - mean_computed_param ;
            var_computed_param += diff * diff ;
            diff = computed_bias_1[i] - mean_computed_bias ;
            var_computed_bias += diff * diff ;
            }
         var_computed_param /= ndone ;   // Variance of theta-hat
         var_computed_bias /= ndone ;    // Variance of bias estimate
         printf ( "\n1: Computed parameter (theta-hat)  mean=%.5lf  variance=%.5lf",
            mean_computed_param, var_computed_param ) ;
         printf ( "\nMean boot bias=%.5lf (its std=%.5lf)   Mean boot var=%.5lf",
            mean_computed_bias, sqrt(var_computed_bias), mean_computed_var ) ;

/*
   Process test 2 of 2
*/

         mean_computed_param = 0.0 ;
         mean_computed_bias = 0.0 ;
         mean_computed_var = 0.0 ;
         var_computed_param = 0.0 ;
         var_computed_bias = 0.0 ;
         for (i=0 ; i<ndone ; i++) {
            mean_computed_param += computed_param_2[i] ;
            mean_computed_bias += computed_bias_2[i] ;
            mean_computed_var += computed_var_2[i] ;
            }
         mean_computed_param /= ndone ;  // Mean across tries of theta-hat
         mean_computed_bias /= ndone ;   // Mean of bias estimates
         mean_computed_var /= ndone ;    // Mean of variance estimates
         for (i=0 ; i<ndone ; i++) {
            diff = computed_param_2[i] - mean_computed_param ;
            var_computed_param += diff * diff ;
            diff = computed_bias_2[i] - mean_computed_bias ;
            var_computed_bias += diff * diff ;
            }
         var_computed_param /= ndone ;   // Variance of theta-hat
         var_computed_bias /= ndone ;    // Variance of bias estimate
         printf ( "\n\n2: Computed parameter (theta-hat)  mean=%.5lf  variance=%.5lf",
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

/******************************************************************************/
/*                                                                            */
/*  BOOT_P_5 - Compare jackknife and bootstrap estimates for bias/var of PF   */
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

double param_pf ( int n , double *x , double *freq )
{
   int i ;
   double sum_win, sum_loss ;

   sum_win = sum_loss = 1.e-5 ;   // Really zero, but prevent division by 0

   if (freq == NULL) {
      for (i=0 ; i<n ; i++) {
         if (x[i] > 0.0)
            sum_win += x[i] ;
         else 
            sum_loss -= x[i] ;
         }
      return sum_win / sum_loss ;
      }

   for (i=0 ; i<n ; i++) {
      if (x[i] > 0.0)
         sum_win += freq[i] * x[i] ;
      else 
         sum_loss -= freq[i] * x[i] ;
      }
   return sum_win / sum_loss ;
}

/*
--------------------------------------------------------------------------------

   boot_bias_var - Compute bias and variance of parameter (plug-in case, s=t)

--------------------------------------------------------------------------------
*/

void boot_bias_var (
   int n ,              // Number of cases in sample
   double *data ,       // The sample
   double (*user_t) (int , double * , double * ) , // Compute param
   int nboot ,          // Number of bootstrap replications
   double *rawstat ,    // Raw statistic of sample, theta-hat
   double *bias ,       // Output of bias estimate
   double *var ,        // Output of variance estimate
   double *work ,       // Work area n long
   double *work2 ,      // Work area nboot long
   double *freq         // Work area n long
   )
{
   int i, rep, k ;
   double stat, mean, variance, diff ;

   mean = 0.0 ;

   for (i=0 ; i<n ; i++)
      freq[i] = 0.0 ;

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         work[i] = data[k] ;            // Put bootstrap sample in work
         ++freq[k] ;                    // Tally for mean frequency
         }

      stat = user_t ( n , work , NULL ) ; // Evaluate estimator for this rep
      work2[rep] = stat ;               // Enables more accurate variance
      mean += stat ;                    // Cumulate theta-hat star dot
      }

   mean /= nboot ;
   variance = 0.0 ;
   for (rep=0 ; rep<nboot ; rep++) {    // Cumulate variance
      diff = work2[rep] - mean ;
      variance += diff * diff ;
      }

   for (i=0 ; i<n ; i++)                // Convert tally of useage
      freq[i] /= nboot * n ;            // To mean frequency of use

   memcpy ( work , data , n * sizeof(double) ) ; // user_t may reorder, so preserve
   *rawstat = user_t ( n , data , NULL) ;        // Final but biased estimate
   *bias = mean - user_t ( n , work , freq ) ;
   *var = variance / (nboot - 1) ;
}

/*
--------------------------------------------------------------------------------

   jack_bias_var - Compute bias and variance of parameter via jackknife

--------------------------------------------------------------------------------
*/

void jack_bias_var (
   int n ,              // Number of cases in sample
   double *data ,       // The sample
   double (*user_t) (int , double * , double * ) , // Compute param
   double *rawstat ,    // Raw statistic of sample, theta-hat
   double *bias ,       // Output of bias estimate
   double *var ,        // Output of variance estimate
   double *work         // Work area n long
   )
{
   int i ;
   double last, temp, param, diff, theta_dot ;

   last = data[n-1] ;
   theta_dot = 0.0 ;
   for (i=0 ; i<n ; i++) {
      temp = data[i] ;              // Preserve case being temporarily removed
      data[i] = last ;              // Swap in last case
      param = user_t ( n-1 , data , NULL ) ; // Param for this jackknife
      theta_dot += param ;          // Cumulate mean across jackknife
      work[i] = param ;             // Save for computing variance later
      data[i] = temp ;              // Restore original case
      }

   theta_dot /= n ;

   *rawstat = user_t ( n , data , NULL ) ;
   *bias = (n - 1) * (theta_dot - *rawstat) ;

   *var = 0.0 ;
   for (i=0 ; i<n ; i++) {
      diff = work[i] - theta_dot ;
      *var += diff * diff ;
      }

   *var *= (n - 1.0) / (double) n ;
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
   double *x, diff, *work, *work2, *freq, mean ;
   double *computed_param_1, *computed_bias_1, *computed_var_1 ;
   double *computed_param_2, *computed_bias_2, *computed_var_2 ;
   double mean_computed_param, var_computed_param ;
   double mean_computed_bias, var_computed_bias, mean_computed_var ;
   double grand_wins, grand_losses ;

/*
   Process command line parameters
*/

   if (argc != 5) {
      printf (
         "\nUsage: BOOT_P_5  nsamples  nboot  ntries  mean" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   mean = atof ( argv[4] ) ;

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)) {
      printf ( "\nUsage: BOOT_P_5  nsamples  nboot  ntries  mean" ) ;
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
   freq = (double *) malloc ( nsamps * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   grand_wins = grand_losses = 0.0 ;

   for (itry=0 ; itry<ntries ; itry++) {

      for (i=0 ; i<nsamps ; i++) {
         x[i] = 1000.0 * normal () + mean ;
         if (x[i] > 0.0)          // Cumulate so we know the true value
            grand_wins += x[i] ;
         else 
            grand_losses -= x[i] ;
         }

      boot_bias_var ( nsamps , x , param_pf , nboot ,
                      &computed_param_1[itry] , &computed_bias_1[itry] ,
                      &computed_var_1[itry] , work , work2 , freq ) ;

      jack_bias_var ( nsamps , x , param_pf ,
                      &computed_param_2[itry] , &computed_bias_2[itry] ,
                      &computed_var_2[itry] , work ) ;

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.
         ndone = itry + 1 ;           // This many tries done (and in arrays)
         printf ( "\n\n\nTry %d   True value = %lf", itry,
                  grand_wins / grand_losses ) ;

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
         printf ( "\nComputed parameter (theta-hat)  mean=%.5lf  variance=%.5lf",
            mean_computed_param, var_computed_param ) ;
         printf ( "\nBOOT: Mean bias=%.5lf (its std=%.5lf)   Mean var=%.5lf",
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
         printf ( "\nJACK: Mean bias=%.5lf (its std=%.5lf)   Mean var=%.5lf",
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

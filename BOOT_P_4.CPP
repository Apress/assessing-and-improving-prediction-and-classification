/******************************************************************************/
/*                                                                            */
/*  BOOT_P_4 - Bootstrap confidence intervals for correlation coef            */
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
double normal_cdf ( double z ) ;
double inverse_normal_cdf ( double p ) ;


/*
--------------------------------------------------------------------------------

   Compute the parameter
   This must NOT change or reorder x or y!

--------------------------------------------------------------------------------
*/

double param_corr ( int n , double *x , double *y )
{
   int i ;
   double xmean, ymean, xdif, ydif, xvar, yvar, covar ;

   xmean = ymean = 0.0 ;
   for (i=0 ; i<n ; i++) {
      xmean += x[i] ;
      ymean += y[i] ;
      }

   xmean /= n ;
   ymean /= n ;

   xvar = yvar = covar = 0.0 ;
   for (i=0 ; i<n ; i++) {
      xdif = x[i] - xmean ;
      ydif = y[i] - ymean ;
      xvar += xdif * xdif ;
      yvar += ydif * ydif ;
      covar += xdif * ydif ;
      }

   if (xvar * yvar > 0.0)
      return covar / sqrt ( xvar * yvar ) ;
   return 0.0 ;
}

/*
--------------------------------------------------------------------------------

   boot_conf_percentile - Compute confidence intervals using percentile method

--------------------------------------------------------------------------------
*/

void boot_conf_percentile (
   int n ,              // Number of cases in sample
   double *x ,          // One variable in sample
   double *y ,          // Other variable in sample
   double (*user_t) (int , double * , double *) , // Compute parameter
   int nboot ,          // Number of bootstrap replications
   double *low5 ,       // Output of lower 5% bound
   double *high5 ,      // Output of upper 5% bound
   double *low10 ,      // Output of lower 10% bound
   double *high10 ,     // Output of upper 10% bound
   double *low25 ,      // Output of lower 25% bound
   double *high25 ,     // Output of upper 25% bound
   double *xwork ,      // Work area n long
   double *ywork ,      // Work area n long
   double *work2        // Work area nboot long
   )
{
   int i, rep, k ;

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         xwork[i] = x[k] ;              // Put bootstrap sample in work
         ywork[i] = y[k] ;
         }

      work2[rep] = user_t ( n , xwork , ywork ) ;
      }

   qsortd ( 0 , nboot-1 , work2 ) ;     // Sort ascending

   k = (int) (0.05 * (nboot + 1)) - 1 ; // Unbiased quantile estimator
   if (k < 0)
      k = 0 ;
   *low5 = work2[k] ;
   *high5 = work2[nboot-1-k] ;

   k = (int) (0.10 * (nboot + 1)) - 1 ;
   if (k < 0)
      k = 0 ;
   *low10 = work2[k] ;
   *high10 = work2[nboot-1-k] ;

   k = (int) (0.25 * (nboot + 1)) - 1 ;
   if (k < 0)
      k = 0 ;
   *low25 = work2[k] ;
   *high25 = work2[nboot-1-k] ;
}

/*
--------------------------------------------------------------------------------

   boot_conf_BCa - Compute confidence intervals using improved quantile method

--------------------------------------------------------------------------------
*/

void boot_conf_BCa (
   int n ,              // Number of cases in sample
   double *x ,          // One variable in sample
   double *y ,          // Other variable in sample
   double (*user_t) (int , double * , double *) , // Compute parameter
   int nboot ,          // Number of bootstrap replications
   double *low5 ,       // Output of lower 5% bound
   double *high5 ,      // Output of upper 5% bound
   double *low10 ,      // Output of lower 10% bound
   double *high10 ,     // Output of upper 10% bound
   double *low25 ,      // Output of lower 25% bound
   double *high25 ,     // Output of upper 25% bound
   double *xwork ,      // Work area n long
   double *ywork ,      // Work area n long
   double *work2        // Work area nboot long
   )
{
   int i, rep, k, z0_count ;
   double param, theta_hat, theta_dot, z0, zlo, zhi, alo, ahi ;
   double xtemp, ytemp, xlast, ylast, diff, numer, denom, accel ;

   theta_hat = user_t ( n , x , y ) ;   // Parameter for full set

   z0_count = 0 ;                       // Will count for computing z0 later

   for (rep=0 ; rep<nboot ; rep++) {    // Do all bootstrap reps (b from 1 to B)

      for (i=0 ; i<n ; i++) {           // Generate the bootstrap sample
         k = (int) (unifrand() * n) ;   // Select a case from the sample
         if (k >= n)                    // Should never happen, but be prepared
            k = n - 1 ;
         xwork[i] = x[k] ;              // Put bootstrap sample in work
         ywork[i] = y[k] ;
         }

      param = user_t ( n , xwork , ywork ) ;  // Param for this bootstrap rep
      work2[rep] = param ;              // Save it for CDF later
      if (param < theta_hat)            // Count how many < full set param
         ++z0_count ;                   // For computing z0 later
      }

   z0 = inverse_normal_cdf ( (double) z0_count / (double) nboot ) ;

/*
   Do the jackknife for computing accel.
   Borrow xwork for storing jackknifed parameter values.
*/

   xlast = x[n-1] ;
   ylast = y[n-1] ;
   theta_dot = 0.0 ;
   for (i=0 ; i<n ; i++) {          // Jackknife
      xtemp = x[i] ;                // Preserve case being temporarily removed
      ytemp = y[i] ;
      x[i] = xlast ;                // Swap in last case
      y[i] = ylast ;
      param = user_t ( n-1 , x , y ) ; // Param for this jackknife
      theta_dot += param ;          // Cumulate mean across lackknife
      xwork[i] = param ;            // Save for computing accel later
      x[i] = xtemp ;                // Restore original case
      y[i] = ytemp ;
      }

/*
   Compute accel
*/

   theta_dot /= n ;
   numer = denom = 0.0 ;
   for (i=0 ; i<n ; i++) {
      diff = theta_dot - xwork[i] ;
      xtemp = diff * diff ;
      denom += xtemp ;
      numer += xtemp * diff ;
      }

   denom = sqrt ( denom ) ;
   denom = denom * denom * denom ;
   accel = numer / (6.0 * denom) ;

/*
   Compute the outputs
*/

   qsortd ( 0 , nboot-1 , work2 ) ;     // Sort ascending

   zlo = inverse_normal_cdf ( 0.05 ) ;
   zhi = inverse_normal_cdf ( 0.95 ) ;
   alo = normal_cdf ( z0 + (z0 + zlo) / (1.0 - accel * (z0 + zlo)) ) ;
   ahi = normal_cdf ( z0 + (z0 + zhi) / (1.0 - accel * (z0 + zhi)) ) ;
   k = (int) (alo * (nboot + 1)) - 1 ; // Unbiased quantile estimator
   if (k < 0)
      k = 0 ;
   *low5 = work2[k] ;
   k = (int) ((1.0-ahi) * (nboot + 1)) - 1 ;
   if (k < 0)
      k = 0 ;
   *high5 = work2[nboot-1-k] ;

   zlo = inverse_normal_cdf ( 0.10 ) ;
   zhi = inverse_normal_cdf ( 0.90 ) ;
   alo = normal_cdf ( z0 + (z0 + zlo) / (1.0 - accel * (z0 + zlo)) ) ;
   ahi = normal_cdf ( z0 + (z0 + zhi) / (1.0 - accel * (z0 + zhi)) ) ;
   k = (int) (alo * (nboot + 1)) - 1 ; // Unbiased quantile estimator
   if (k < 0)
      k = 0 ;
   *low10 = work2[k] ;
   k = (int) ((1.0-ahi) * (nboot + 1)) - 1 ;
   if (k < 0)
      k = 0 ;
   *high10 = work2[nboot-1-k] ;

   zlo = inverse_normal_cdf ( 0.25 ) ;
   zhi = inverse_normal_cdf ( 0.75 ) ;
   alo = normal_cdf ( z0 + (z0 + zlo) / (1.0 - accel * (z0 + zlo)) ) ;
   ahi = normal_cdf ( z0 + (z0 + zhi) / (1.0 - accel * (z0 + zhi)) ) ;
   k = (int) (alo * (nboot + 1)) - 1 ; // Unbiased quantile estimator
   if (k < 0)
      k = 0 ;
   *low25 = work2[k] ;
   k = (int) ((1.0-ahi) * (nboot + 1)) - 1 ;
   if (k < 0)
      k = 0 ;
   *high25 = work2[nboot-1-k] ;
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
   int low5, high5, low10, high10, low25, high25 ;
   double corr, *x, *y, *xwork, *ywork, *work2, *param, x1, x2 ;
   double *low5_1, *high5_1, *low10_1, *high10_1, *low25_1, *high25_1 ;
   double *low5_2, *high5_2, *low10_2, *high10_2, *low25_2, *high25_2 ;
   double *low5_3, *high5_3, *low10_3, *high10_3, *low25_3, *high25_3 ;
   double mean_param ;

/*
   Process command line parameters
*/

   if (argc != 5) {
      printf (
         "\nUsage: BOOT_P_4  nsamples  nboot  ntries  corr" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   corr = atof ( argv[4] ) ;

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)
    || (corr < -1.0)  ||  (corr > 1.0)) {
      printf ( "\nUsage: BOOT_P_4  nsamples  nboot  ntries  corr" ) ;
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
   param = (double *) malloc ( ntries * sizeof(double) ) ;
   low5_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   high5_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   low10_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   high10_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   low25_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   high25_1 = (double *) malloc ( ntries * sizeof(double) ) ;
   low5_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   high5_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   low10_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   high10_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   low25_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   high25_2 = (double *) malloc ( ntries * sizeof(double) ) ;
   low5_3 = (double *) malloc ( ntries * sizeof(double) ) ;
   high5_3 = (double *) malloc ( ntries * sizeof(double) ) ;
   low10_3 = (double *) malloc ( ntries * sizeof(double) ) ;
   high10_3 = (double *) malloc ( ntries * sizeof(double) ) ;
   low25_3 = (double *) malloc ( ntries * sizeof(double) ) ;
   high25_3 = (double *) malloc ( ntries * sizeof(double) ) ;

/*
   Main outer loop does all tries
*/

   for (itry=0 ; itry<ntries ; itry++) {

      if ((itry % divisor) == 0)
         printf ( "\n\n\nTry %d", itry ) ;

      for (i=0 ; i<nsamps ; i++) {
         x1 = normal () ;
         x2 = normal () ;
         x[i] = x1 ;
         y[i] = corr * x1 + sqrt ( 1.0 - corr * corr ) * x2 ;
         }

      param[itry] = param_corr ( nsamps , x , y ) ;

      boot_conf_percentile ( nsamps , x , y , param_corr , nboot ,
                   &low5_1[itry] , &high5_1[itry] , &low10_1[itry] , &high10_1[itry] , 
                   &low25_1[itry] , &high25_1[itry] , xwork , ywork , work2 ) ;

      boot_conf_BCa ( nsamps , x , y , param_corr , nboot ,
           &low5_2[itry] , &high5_2[itry] , &low10_2[itry] , &high10_2[itry] , 
           &low25_2[itry] , &high25_2[itry] , xwork , ywork , work2 ) ;

      // The inverted pivot (basic) intervals are trivially obtained from the
      // percentile intervals
      low5_3[itry] = 2.0 * param[itry] - high5_1[itry] ;
      high5_3[itry] = 2.0 * param[itry] - low5_1[itry] ;
      low10_3[itry] = 2.0 * param[itry] - high10_1[itry] ;
      high10_3[itry] = 2.0 * param[itry] - low10_1[itry] ;
      low25_3[itry] = 2.0 * param[itry] - high25_1[itry] ;
      high25_3[itry] = 2.0 * param[itry] - low25_1[itry] ;

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.
         ndone = itry + 1 ;           // This many tries done (and in arrays)

         mean_param = 0.0 ;
         for (i=0 ; i<ndone ; i++)
            mean_param += param[i] ;
         mean_param /= ndone ;

         printf ( "\nMean = %.5lf   Next line is percentile method, then BCa, then pivot",
                  mean_param ) ;

/*
   Process test 1 of 3
*/

         low5 = high5 = low10 = high10 = low25 = high25 = 0 ;
         for (i=0 ; i<ndone ; i++) {
            if (low5_1[i] > corr)
               ++low5 ;
            if (high5_1[i] < corr)
               ++high5 ;
            if (low10_1[i] > corr)
               ++low10 ;
            if (high10_1[i] < corr)
               ++high10 ;
            if (low25_1[i] > corr)
               ++low25 ;
            if (high25_1[i] < corr)
               ++high25 ;
            }
         printf ( "\n5: (%.3lf %.3lf)  10: (%.3lf %.3lf)  25: (%.3lf %.3lf)",
            100.0 * low5 / ndone , 100.0 * high5 / ndone ,
            100.0 * low10 / ndone , 100.0 * high10 / ndone ,
            100.0 * low25 / ndone , 100.0 * high25 / ndone ) ;

/*
   Process test 2 of 3
*/

         low5 = high5 = low10 = high10 = low25 = high25 = 0 ;
         for (i=0 ; i<ndone ; i++) {
            if (low5_2[i] > corr)
               ++low5 ;
            if (high5_2[i] < corr)
               ++high5 ;
            if (low10_2[i] > corr)
               ++low10 ;
            if (high10_2[i] < corr)
               ++high10 ;
            if (low25_2[i] > corr)
               ++low25 ;
            if (high25_2[i] < corr)
               ++high25 ;
            }
         printf ( "\n5: (%.3lf %.3lf)  10: (%.3lf %.3lf)  25: (%.3lf %.3lf)",
            100.0 * low5 / ndone , 100.0 * high5 / ndone ,
            100.0 * low10 / ndone , 100.0 * high10 / ndone ,
            100.0 * low25 / ndone , 100.0 * high25 / ndone ) ;

/*
   Process test 3 of 3
*/

         low5 = high5 = low10 = high10 = low25 = high25 = 0 ;
         for (i=0 ; i<ndone ; i++) {
            if (low5_3[i] > corr)
               ++low5 ;
            if (high5_3[i] < corr)
               ++high5 ;
            if (low10_3[i] > corr)
               ++low10 ;
            if (high10_3[i] < corr)
               ++high10 ;
            if (low25_3[i] > corr)
               ++low25 ;
            if (high25_3[i] < corr)
               ++high25 ;
            }
         printf ( "\n5: (%.3lf %.3lf)  10: (%.3lf %.3lf)  25: (%.3lf %.3lf)",
            100.0 * low5 / ndone , 100.0 * high5 / ndone ,
            100.0 * low10 / ndone , 100.0 * high10 / ndone ,
            100.0 * low25 / ndone , 100.0 * high25 / ndone ) ;
         } // If progress report

      if ((itry % 10) == 1) {
         if (_kbhit ()) {
            if (_getch() == 27)
               break ;
            }
         }

      } // For all tries


   return EXIT_SUCCESS ;
}

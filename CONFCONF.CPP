/******************************************************************************/
/*                                                                            */
/*  CONFCONF - Demonstrate confidence in hypothesis test confidence           */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

double unifrand () ;
void qsortd ( int istart , int istop , double *x ) ;
double quantile_conf ( int n , int m , double conf ) ;
double inverse_ks ( int n , double cdf ) ;
double ks_test ( int n , double *x , double *D_plus , double *D_minus ) ;

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )
{
   int i, ncases, irep, nreps, m, n_lower, n_upper, n_ks2, n_ks_null, n_ks_alt ;
   double *x, pval, conf, pessimistic_lower, pessimistic_upper ;
   double ks_two, ks_one, D, Dp, Dm ;

   if (argc != 5) {
      printf ( "\nUsage: ConfConf  ncases  pval  conf  nreps" ) ;
      printf ( "\n  ncases - Number of cases in the sample" ) ;
      printf ( "\n  pval - Probability value (<0.5) for quantile test" ) ;
      printf ( "\n  conf - Desired confidence value (<0.5) for both tests" ) ;
      printf ( "\n  nreps - Number of replications" ) ;
      exit ( 1 ) ;
      }

   ncases = atoi ( argv[1] ) ;
   pval = atof ( argv[2] ) ;
   conf = atof ( argv[3] ) ;
   nreps = atoi ( argv[4] ) ;

   if (ncases < 10) {
      printf ( "\nERROR.. Must have at least 10 cases" ) ;
      exit ( 1 ) ;
      }

   if (pval * ncases < 1.0  ||  pval >= 0.5) {
      printf ( "\nERROR.. Pval too small or too large" ) ;
      exit ( 1 ) ;
      }

   if (conf <= 0.0  ||  conf >= 0.5) {
      printf ( "\nERROR.. Conf must be greater than 0 and less than 0.5" ) ;
      exit ( 1 ) ;
      }

   if (nreps < 1) {
      printf ( "\nERROR.. Must have at least 1 replication" ) ;
      exit ( 1 ) ;
      }


/*
   Allocate memory and initialize
*/

   x = (double *) malloc ( ncases * sizeof(double) ) ;

   m = (int) (pval * ncases) ;  // Conservative order statistic for bound
   pessimistic_lower = quantile_conf ( ncases , m , conf ) ;
   pessimistic_upper = 1.0 - pessimistic_lower ;
   ks_two = inverse_ks ( ncases , 1.0 - conf ) ;       // Two-tailed test
   ks_one = inverse_ks ( ncases , 1.0 - 2.0 * conf ) ; // One-tailed test

   printf ( "\nSuppose the model predicts values near 0 for the null hypothesis" ) ;
   printf ( "\nand values near 1 for the alternative hypothesis." ) ;

   printf ( "\n\nIf the dataset represents the null hypothesis, the threshold" ) ;
   printf ( "\nfor rejecting the null at p=%.4lf is given by the %d'th order statistic.",
            pval, ncases - m + 1 ) ;
   printf ( "\nThis is a conservative estimate of the %.4lf quantile", 1.0-pval ) ;
   printf ( "\nThere is only a %.4lf chance that it will really be the %.4lf quantile or worse.",
            conf, pessimistic_upper ) ;

   printf ( "\n\nIf the dataset represents the alternative hypothesis, the threshold" ) ;
   printf ( "\nfor rejecting the alt at p=%.4lf is given by the %d'th order statistic.",
            pval, m ) ;
   printf ( "\nThis is a conservative estimate of the %.4lf quantile", pval ) ;
   printf ( "\nThere is only a %.4lf chance that it will really be the %.4lf quantile or worse.",
            conf, pessimistic_lower) ;

   printf ( "\n\nKS thresholds: two-tailed KS = %.4lf  one-tailed KS = %.4lf",
            ks_two, ks_one ) ;

/*
   Now generate nreps samples.  Verify that our required confidence level
   is observed.  Note that the fact that this test uses a uniform distribution
   does not in any way limit its applicability to uniform distributions.
   If one were to generate cases from any other reasonable distribtion,
   the pessimistic quantile bounds would have to be transformed similarly.
   The result is that the inequalities below would pass or fail identically.
   We count the number of times 'disaster' happens.
   Disaster is when the order statistic used for the threshold is toward the
   inside (center) of the distribution, meaning that if this order statistic
   had been used as a threshold, more of the distribution would be outside
   the threshold than the user expected.  We expect disaster to happen with
   probability equal to the specified conf parameter.

   For the two-tailed Kolmogorov-Smirnov test, disaster is when the empirical
   CDF deviates (above or below) from the correct value by more than the
   conf-inspired value.  For the one-tailed test in which the dataset is from
   the NULL distribution, disaster is when the empirical CDF exceeds the true
   CDF, a situation that would encourage false rejection of the null hypothesis.
   This is measured by D+.  For the one-tailed test in which the dataset is from
   the ALT distribution, disaster is when the empirical CDF is less than the
   true CDF, a situation that would encourage false rejection of the alternative
   hypothesis.  This is measured by D-.
*/

   n_lower = n_upper = n_ks2 = n_ks_null = n_ks_alt = 0 ;

   for (irep=0 ; irep<nreps ; irep++) {

      for (i=0 ; i<ncases ; i++)
         x[i] = unifrand () ;
      qsortd ( 0 , ncases-1 , x ) ;

      if (x[m-1] > pessimistic_lower)
         ++n_lower ;

      if (x[ncases-m] < pessimistic_upper)
         ++n_upper ;

      D = ks_test ( ncases , x , &Dp , &Dm ) ;
      if (D > ks_two)
         ++n_ks2 ;
      if (Dp > ks_one)
         ++n_ks_null ;
      if (Dm > ks_one)
         ++n_ks_alt ;
      }

   printf ( "\nPoint failure (expected=%.4lf)  Lower=%.4lf  Upper=%.4lf",
            conf, (double) n_lower / nreps, (double) n_upper / nreps) ;
   printf ( "\nKS failure:  two-tailed = %.4lf  NULL = %.4lf  ALT = %.4lf",
            (double) n_ks2 / nreps, (double) n_ks_null / nreps,
            (double) n_ks_alt / nreps) ;

   free ( x ) ;
   return ( 0 ) ;
}

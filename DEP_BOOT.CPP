/******************************************************************************/
/*                                                                            */
/*  DEP_BOOT - Dependent bootstrap routines                                   */
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

#define PI 3.141592653589793

/*
--------------------------------------------------------------------------------

   Compute the parameter, which in this case is the sample mean

--------------------------------------------------------------------------------
*/

double param_mean ( int n , double *x )
{
   int i ;
   double sum ;

   sum = 0.0 ;
   for (i=0 ; i<n ; i++)
      sum += x[i] ;

   return sum / n ;
}

/*
--------------------------------------------------------------------------------

   Compute the influence function for the Tapered Block Bootstrap
   This takes the original dataset as input and outputs the set of
   influence function values.
   This particular example demonstrates the sample mean.

--------------------------------------------------------------------------------
*/

void influence_mean ( int n , double *x , double *xinf )
{
   int i ;
   double mean ;

   mean = 0.0 ;
   for (i=0 ; i<n ; i++)
      mean += x[i] ;

   mean /= n ;

   for (i=0 ; i<n ; i++)
      xinf[i] = x[i] - mean ;
}

/*
--------------------------------------------------------------------------------

   Create an AR(1) series with specified coefficient

--------------------------------------------------------------------------------
*/

void create_AR1 ( int n , double coef , double *x )
{
   int i ;

   x[0] = normal() * sqrt ( 1.0 / (1.0 - coef * coef) ) ;
   for (i=1 ; i<n ; i++)
      x[i] = coef * x[i-1] + normal() ;
}

/*
--------------------------------------------------------------------------------

   Compute the standard error of the mean of a sample from an AR1 series.
   This algorithm is the brute-force method that sums historical variances.
   I suspect that there is a simple direct formula, but I checked a vast
   number of references and couldn't find any.

--------------------------------------------------------------------------------
*/

double StdOfAR1mean ( int n , double coef )
{
   int i ;
   double total, factor, mult ;

   total = 0.0 ;
   factor = mult = 1.0 ;

   i = n ;
   while (i--) {
      total += factor * factor ;
      mult *= coef ;
      factor += mult ;
      }

   i = n ;
   factor = 0.0 ;
   mult = 1.0 ;
   while (i--) {
      mult *= coef ;
      factor += mult ;
      }

   for (;;) {
      total += factor * factor ;
      factor *= coef ;
      if (factor < 1.e-6)
         break ;
      }

   return sqrt ( total / (n * n) ) ;
}

/*
--------------------------------------------------------------------------------

   Draw a Stationary Bootstrap sample

--------------------------------------------------------------------------------
*/

void SBsample ( int n , int blocksize , double *x , double *bootsamp )
{
   int i, pos ;
   double q ;

   q = 1.0 / blocksize ;              // Parameter for geometric distribution

   pos = (int) (unifrand() * n) ;     // Pick a random starting point
   if (pos >= n)                      // Should never happen
      pos = n - 1 ;                   // But avoid disaster

   for (i=0 ; i<n ; i++) {            // Build the bootstrap sample
      bootsamp[i] = x[pos] ;          // Get a case
      if (unifrand() < q) {           // Implement the geometric distribution
         pos = (int) (unifrand() * n) ; // We may choose a new random position
         if (pos >= n)                // (protected as before)
            pos = n - 1 ;
         }
      else
         pos = (pos + 1) % n ;        // Or we may simply advance circularly
      }
}

/*
--------------------------------------------------------------------------------

   Compute the taper window for the Tapered Block Bootstrap

--------------------------------------------------------------------------------
*/

void make_taper ( int blocksize , double *window )
{
   int i, low, high ;
   double w, sum ;

   low = 0 ;
   high = blocksize - 1 ;

   for (;;) {
      w = (low + 0.5) / (double) blocksize ;  // Position in block
      if (w < 0.43) {                // If near edge
         w /= 0.43 ;                 // Taper upwards
         window[low++] = w ;         // Insert the taper
         window[high--] = w ;        // It is symmetric
         }
      else {                         // Come here when well away from edge
         while (low <= high)         // Fill in the center
            window[low++] = 1.0 ;    // With full value
         break ;                     // Done
         }
      }

   // Normalize the length of the window

   sum = 0.0 ;    // Will cumulate squared length here
   for (i=0 ; i<blocksize ; i++)
      sum += window[i] * window[i] ;

   w = sqrt ( blocksize / sum ) ;
   for (i=0 ; i<blocksize ; i++)
      window[i] *= w ;
}

/*
--------------------------------------------------------------------------------

   Draw a Tapered Block Bootstrap sample
   It's best (but not required) if n is an integer multiple of blocksize

--------------------------------------------------------------------------------
*/

void TBBsample ( int n , int blocksize , double *window ,
                 double *x , double *bootsamp )
{
   int i, j, k, pos ;

   j = 0 ;                            // Will index bootsamp
   k = (int) (n / blocksize) ;        // Number of blocks

   while (k--) {                      // Count blocks done
      pos = (int) (unifrand() * (n-blocksize+1)) ; // Pick a random starting point
      if (pos > (n - blocksize))      // Should never happen
         pos = n - blocksize ;        // But avoid disaster
      for (i=0 ; i<blocksize ; i++)   // Build the bootstrap sample
         bootsamp[j++] = x[pos+i] * window[i] ; // Get a windowed case
      }
}

/*
--------------------------------------------------------------------------------

   Estimate the standard error of the mean using the Stationary Bootstrap

--------------------------------------------------------------------------------
*/

double StdErrMeanSB (
   int n ,         // Number of cases in sample
   double *x ,     // The sample
   int blocksize , // Block size
   int nboot ,     // Number of bootstrap replications to do
   double *bs      // Work area n long for bootstrap sample
   )
{
   int i, iboot ;
   double mean, grandmean, diff, sumsq ;

   grandmean = 0.0 ;       // Compute mean of original sample
   for (i=0 ; i<n ; i++)   // We will use it instead of the mean of the
      grandmean += x[i] ;  // bootstrap samples because it is slightly
   grandmean /= n ;        // more accurate

   sumsq = 0.0 ;           // Sum squared deviations from grand mean
   for (iboot=0 ; iboot<nboot ; iboot++) {
      SBsample ( n , blocksize , x , bs ) ;
      mean = 0.0 ;         // Compute mean of bootstrap sample
      for (i=0 ; i<n ; i++)
         mean += bs[i] ;
      mean /= n ;
      diff = mean - grandmean ;  // Deviation from grand mean
      sumsq += diff * diff ;     // Cumulate squared deviations
      }

   return sqrt ( sumsq / nboot ) ;
}

/*
--------------------------------------------------------------------------------

   Estimate the standard error of the mean using the Tapered Block Bootstrap

   Note that since we know this will be used for the sample mean, we could
   easily compute the influence function here instead of using a function
   call.  But the explicit call shows how it could be generalized.

--------------------------------------------------------------------------------
*/

double StdErrMeanTBB (
   int n ,         // Number of cases in sample
   double *x ,     // The sample
   int blocksize , // Block size
   int nboot ,     // Number of bootstrap replications to do
   double *xinf ,  // Work area n long for influence function values
   double *bs ,    // Work area n long for bootstrap sample
   double *window  // Work area blocksize long for window
   )
{
   int i, k, iboot ;
   double mean, sumsq ;

   influence_mean ( n , x , xinf ) ; // Compute influence function for each case
   make_taper ( blocksize , window ) ;  // Compute the tapered window
   k = blocksize * (int) (n / blocksize) ; // Length of TBB sample (<=n)

   sumsq = 0.0 ;           // Sum squared deviations from zero (mean of xinf)
   for (iboot=0 ; iboot<nboot ; iboot++) {
      TBBsample ( n , blocksize , window , xinf , bs ) ;
      mean = 0.0 ;         // Compute mean of bootstrap sample
      for (i=0 ; i<k ; i++)
         mean += bs[i] ;
      mean /= k ;
      sumsq += mean * mean ;     // Cumulate squared deviations (from zero)
      }

   return sqrt ( (double) k / (double) n * sumsq / nboot ) ;
}

/*
--------------------------------------------------------------------------------

   Estimate a quantile of the deviation of the mean using the
   Stationary Bootstrap

--------------------------------------------------------------------------------
*/

double QuantileMeanSB (
   int n ,         // Number of cases in sample
   double *x ,     // The sample
   int blocksize , // Block size
   int nboot ,     // Number of bootstrap replications to do
   double q ,      // Desired quantile, 0-1
   double *bs ,    // Work area n long for bootstrap sample
   double *reps    // Work area nboot long for replications
   )
{
   int i, iboot, subscript ;
   double mean, grandmean ;

   if (q <= 0.5)   // Formula for unbiased subscript only works if q<=.5
      subscript = (int) (q * (nboot + 1)) - 1 ;
   else 
      subscript = nboot - (int) ((1.0 - q) * (nboot + 1)) ;

   if (subscript < 0)  // Ensure that a silly user didn't put us outside bounds
      subscript = 0 ;
   else if (subscript >= nboot)
      subscript = nboot - 1 ;

   grandmean = 0.0 ;       // Compute mean of original sample
   for (i=0 ; i<n ; i++)   // We will use it instead of the mean of the
      grandmean += x[i] ;  // bootstrap samples because it is slightly
   grandmean /= n ;        // more accurate

   for (iboot=0 ; iboot<nboot ; iboot++) {
      SBsample ( n , blocksize , x , bs ) ;
      mean = 0.0 ;         // Compute mean of bootstrap sample
      for (i=0 ; i<n ; i++)
         mean += bs[i] ;
      mean /= n ;
      reps[iboot] = mean - grandmean ;
      }

   qsortd ( 0 , nboot-1 , reps ) ;
   return reps[subscript] ;
}

/*
--------------------------------------------------------------------------------

   Estimate a quantile of the deviation of the mean using the
   Tapered Block Bootstrap

   Note that since we know this will be used for the sample mean, we could
   easily compute the influence function here instead of using a function
   call.  But the explicit call shows how it could be generalized.

--------------------------------------------------------------------------------
*/

double QuantileMeanTBB (
   int n ,         // Number of cases in sample
   double *x ,     // The sample
   int blocksize , // Block size
   int nboot ,     // Number of bootstrap replications to do
   double q ,      // Desired quantile, 0-1
   double *xinf ,  // Work area n long for influence function values
   double *bs ,    // Work area n long for bootstrap sample
   double *reps ,  // Work area nboot long for replications
   double *window  // Work area blocksize long for window
   )
{
   int i, k, iboot, subscript ;
   double mean ;

   if (q <= 0.5)   // Formula for unbiased subscript only works if q<=.5
      subscript = (int) (q * (nboot + 1)) - 1 ;
   else 
      subscript = nboot - (int) ((1.0 - q) * (nboot + 1)) ;

   if (subscript < 0)  // Ensure that a silly user didn't put us outside bounds
      subscript = 0 ;
   else if (subscript >= nboot)
      subscript = nboot - 1 ;

   influence_mean ( n , x , xinf ) ; // Compute influence function for each case
   make_taper ( blocksize , window ) ;  // Compute the tapered window
   k = blocksize * (int) (n / blocksize) ; // Length of TBB sample (<=n)

   for (iboot=0 ; iboot<nboot ; iboot++) {
      TBBsample ( n , blocksize , window , xinf , bs ) ;
      mean = 0.0 ;         // Compute mean of bootstrap sample
      for (i=0 ; i<k ; i++)
         mean += bs[i] ;
      mean /= k ;
      reps[iboot] = mean ;
      }

   qsortd ( 0 , nboot-1 , reps ) ;

   return sqrt ( (double) k / (double) n ) * reps[subscript] ;
}

/*
--------------------------------------------------------------------------------

   Find the smallest positive integer m after which the autocorrelation
   appears negligible.  This is from [Politis and White, 2003]

--------------------------------------------------------------------------------
*/

int correlation_extent (
   int n ,           // Length of series
   int maxlag ,      // Max lag to test
   double *x ,       // The series
   double *autocov   // Output of autocovariance, maxlag+1 long
   )
{
   int i, lag, kn, nsmall, m ;
   double mean, sum, thresh ;

   // Compute the autocovariance of the sample.

   mean = 0.0 ;
   for (i=0 ; i<n ; i++)
      mean += x[i] ;
   mean /= n ;

   for (lag=0 ; lag<=maxlag ; lag++) {
      sum = 0.0 ;
      for (i=lag ; i<n ; i++)
         sum += (x[i] - mean) * (x[i-lag] - mean) ;
      autocov[lag] = sum / n ;
      }

   // Use the Politis and White heuristics

   kn = (int) (sqrt ( log10 ( (double) n ) ) + 0.5) ;
   if (kn < 5)
      kn = 5 ;

   thresh = 2.0 * sqrt ( log10 ( (double) n ) / n ) ;

   nsmall = 0 ;
   m = maxlag ;

   for (i=1 ; i<=maxlag ; i++) {
      if (fabs ( autocov[i] / autocov[0] ) < thresh)
         ++nsmall ;
      else
         nsmall = 0 ;
      if (nsmall >= kn) {
         m = i - kn ;
         break ;
         }
      }

   if (m < 1)
      m = 1 ;

   return m ;
}

/*
--------------------------------------------------------------------------------

   Estimate the optimal block length for the Stationary Bootstrap

   The first little function here is g-hat(w) per
   page 7 of [Politis and White, 2003].
   The second is the integrand that will be needed later.
   Note that we use statics to pass the two nuisance parameters, which
   simplifies use of canned integration if we wish.

--------------------------------------------------------------------------------
*/


static double gw ( double w , int M , double *autocov )
{
   int k ;
   double lambda, sum ;

   sum = 0.0 ;
   for (k=1 ; k<M ; k++) {
      lambda = (double) k / M ;
      lambda = (lambda < 0.5)  ?  1.0  :  (2.0 * (1.0 - lambda)) ;
      sum += lambda * autocov[k] * cos ( k * w ) ;
      }
   // The summation is symmetric, so double one side and add center (k=0) term
   return 2.0 * sum + autocov[0] ;
}

static int this_M ;             // Facilitates future use of canned integration
static double *this_autocov ;

static double integrand ( double w )
{
   double temp ;
   temp = gw ( w , this_M , this_autocov ) ;
   return (1.0 + cos ( w )) * temp * temp ;
}

int optimal_SB_size (
   int n ,              // Number of cases in sample
   double *x ,          // Variable in sample
   double *autocov      // Work area n long (Actually, not all n used)
   )
{
   int i, k, maxlag, m, nint ;
   double sum, prior_sum, ghat, lambda, dhat, term, space ;

/*
   Compute m as the minimum integer after which correlation appears negligible.
   Also compute the autocovariance.
*/

   maxlag = n / 4 ;   // This is very conservative
   if (maxlag > 50)   // And this is fairly conservative
      maxlag = 50 ;

   m = correlation_extent ( n , maxlag , x , autocov ) ;

/*
   Compute G-hat per page 7 of Politis and White.
   Note that the sum is symmetric, so we only need to evaluate one side,
   double it, and add in the center (k=0).  But we don't even need to worry
   about the center, because this term is zero.
   Also note that lambda(1)=0 so we do not need to do k=M.
*/

   ghat = 0.0 ;
   for (k=1 ; k<2*m ; k++) {
      lambda = (double) k / (2.0 * m) ;
      lambda = (lambda < 0.5)  ?  1.0  :  (2.0 * (1.0 - lambda)) ;
      ghat += lambda * k * autocov[k] ;
      }
   ghat *= 2.0 ;

/*
   Compute the integral.  Because the integrand is not terribly slow
   to compute, and because this routine will not be used in a time-critical
   application, we use crude trapezoidal subdivision.
   Fanatics could easily use a canned package instead.
   Note that the actual integral is from -PI to PI.  But the integrand
   is even in w, so we just integrate from 0 to PI and double it.
*/

   this_M = 2 * m ;         // Lets us avoid passing nuisance parameters
   this_autocov = autocov ; // In case we want to use a canned package later
   nint = 1 ;               // Number of new integration points
   sum = 0.5 * (integrand ( 0.0 ) + integrand ( PI ) ) ; // Original interval
   space = PI ;             // Spacing for that original interval

   for (;;) {               // Endless loop waits for convergence or give up
      prior_sum = sum ;     // Convergence indicator also holds old estimate
      space *= 0.5 ;        // Spacing for the upcoming subinterval
      sum = 0.0 ;           // Will cumulate subdivision here
      for (i=0 ; i<nint ; i++)  // Sum the refinement term
         sum += integrand ( space + 2 * i * space ) ; // Subdivide
      sum /= nint ;         // Refinement term
      sum = 0.5 * (prior_sum + sum) ; // This is the refined estimate
      nint *= 2 ;           // Number of terms in next refinement subdivision
      if (nint < 64)        // Avoid convergence check until stable
         continue ;
      if (nint >= 4096)     // Unlikely, but avoid huge time
         break ;
      term = fabs ( sum - prior_sum ) ;  // This is how much we just changed
      if ((term < 1.e-10)  ||  (term / (1.e-60 + fabs ( sum )) < 1.e-10))
         break ; // This is the real (incredibly conservative) convergence test
      }
   sum *= 4 ;  // Mean height * PI * 2 * (2 / PI)

/*
   Compute D-hatSB per Equation (8) Page 7 of Politis and White
*/

   term = gw ( 0.0 , this_M , autocov ) ;
   dhat = 4.0 * term * term + sum ;

/*
   Compute the optimal block size
*/

   term = 2.0 * ghat * ghat / dhat ;
   term = exp ( log ( term * n ) / 3.0 ) ;
   return (int) (term + 0.999999999) ;
}

/*
--------------------------------------------------------------------------------

   optimal_TBB_size - Compute optimal block size for Tapered Block Bootstrap

--------------------------------------------------------------------------------
*/

static int optimal_TBB_size (
   int n ,              // Number of cases in sample
   double *x ,          // Variable in sample
   double *autocov      // Work area n long (Actually, not all n used)
   )
{
   int k, maxlag, m ;
   double gamma, lambda, delta, term ;

/*
   Compute m as the minimum integer after which correlation appears negligible.
   Also compute the autocovariance.
*/

   maxlag = n / 4 ;   // This is very conservative
   if (maxlag > 50)   // And this is fairly conservative
      maxlag = 50 ;

   m = correlation_extent ( n , maxlag , x , autocov ) ;

/*
   Compute Gamma per page 138 of Paparoditis and Politis.
   Note that the sum is symmetric, so we only need to evaluate one side,
   double it, and add in the center (k=0).  But we don't even need to worry
   about the center, because this term is zero.
   Also note that lambda(1)=0 so we do not need to do k=M.
*/

   gamma = 0.0 ;
   for (k=1 ; k<2*m ; k++) {
      lambda = (double) k / (2.0 * m) ;
      lambda = (lambda < 0.5)  ?  1.0  :  (2.0 * (1.0 - lambda)) ;
      gamma += lambda * k * k * autocov[k] ;
      }
   gamma *= -10.9 ;

/*
   Compute Delta.  This, too, is symmetric, but the center term is not zero.
*/

   delta = 0.0 ;
   for (k=1 ; k<2*m ; k++) {
      lambda = (double) k / (2.0 * m) ;
      lambda = (lambda < 0.5)  ?  1.0  :  (2.0 * (1.0 - lambda)) ;
      delta += lambda * autocov[k] ;
      }
   delta = 2.0 * delta + autocov[0] ; // 2 * side + center
   delta = 1.1 * delta * delta ;

/*
   Compute the optimal block size
*/

   term = 4.0 * gamma * gamma / delta ;
   term = exp ( log ( term * n ) / 5.0 ) ;
   return (int) (term + 0.999999999) ;
}


/*
--------------------------------------------------------------------------------

   Optional main to test routines

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, ib, lastb, maxb, ntries, itry, nsamps, nboot, divisor, ndone ;
   int OptBminSB, OptBmaxSB, OptBminTBB, OptBmaxTBB ;
   double rb, factor, coef, *x, *xinf, *bs, *reps, *window, estimate, diff ;
   double SampleMean, CorrectStdErr, CorrectQuantile ;
   double *StdErrBiasSB, *StdErrErrSB, *QuantileBiasSB ;
   double *QuantileErrSB, *QuantileRejectSB, *StdErrBiasTBB, *StdErrErrTBB ;
   double *QuantileBiasTBB, *QuantileErrTBB, *QuantileRejectTBB ;
   double *autocov ;
   double OptBmeanSB, OptBmeanTBB ;

/*
   Process command line parameters
*/

   if (argc != 5) {
      printf ( "\nUsage: DEP_BOOT  nsamples  nboot  ntries  coef" ) ;
      exit ( 1 ) ;
      }

   nsamps = atoi ( argv[1] ) ;
   nboot = atoi ( argv[2] ) ;
   ntries = atoi ( argv[3] ) ;
   coef = atof ( argv[4] ) ;

   if ((nsamps <= 0)  ||  (nboot <= 0)  ||  (ntries <= 0)
     || (coef < 0.0)  ||  (coef >= 1.0)) {
      printf ( "\nUsage: DEP_BOOT  nsamples  nboot  ntries  coef" ) ;
      exit ( 1 ) ;
      }

   CorrectStdErr = StdOfAR1mean ( nsamps , coef ) ;
   CorrectQuantile = CorrectStdErr * inverse_normal_cdf ( 0.1 ) ;

   divisor = 1000000 / (nsamps * nboot) ;  // This is for progress reports only
   if (divisor < 2)
      divisor = 2 ;

/*
   Allocate memory and initialize
*/

   maxb = nsamps / 4 ;   // Max block size to test

   x = (double *) malloc ( nsamps * sizeof(double) ) ;
   xinf = (double *) malloc ( nsamps * sizeof(double) ) ;
   bs = (double *) malloc ( nsamps * sizeof(double) ) ;
   reps = (double *) malloc ( nboot * sizeof(double) ) ;
   window = (double *) malloc ( nsamps * sizeof(double) ) ;
   autocov = (double *) malloc ( nsamps * sizeof(double) ) ;
   StdErrBiasSB = (double *) malloc ( maxb * sizeof(double) ) ;
   StdErrErrSB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileBiasSB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileErrSB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileRejectSB = (double *) malloc ( maxb * sizeof(double) ) ;
   StdErrBiasTBB = (double *) malloc ( maxb * sizeof(double) ) ;
   StdErrErrTBB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileBiasTBB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileErrTBB = (double *) malloc ( maxb * sizeof(double) ) ;
   QuantileRejectTBB = (double *) malloc ( maxb * sizeof(double) ) ;

   for (i=0 ; i<maxb ; i++) {
      StdErrBiasSB[i] = 0.0 ;
      StdErrErrSB[i] = 0.0 ;
      QuantileBiasSB[i] = 0.0 ;
      QuantileErrSB[i] = 0.0 ;
      QuantileRejectSB[i] = 0.0 ;
      StdErrBiasTBB[i] = 0.0 ;
      StdErrErrTBB[i] = 0.0 ;
      QuantileBiasTBB[i] = 0.0 ;
      QuantileErrTBB[i] = 0.0 ;
      QuantileRejectTBB[i] = 0.0 ;
      }

/*
--------------------------------------------------------------------------------

   First we compute the bias and RMS error of the estimators of
   standard error of the mean and a low quantile (0.1 here) of
   the deviation of the mean.

   Main outer loop 'itry' does all Monte-Carlo replications

--------------------------------------------------------------------------------
*/

   for (itry=0 ; itry<ntries ; itry++) {

      if ((itry % divisor) == 0)
         printf ( "\n\n\nTry %d", itry ) ;

      // Create the sample and find its mean
      SampleMean = 0.0 ;
      create_AR1 ( nsamps , coef , x ) ;
      for (i=0 ; i<nsamps ; i++)
         SampleMean += x[i] ;
      SampleMean /= nsamps ;

      rb = 1.0 ;   // This spaces block sizes intelligently for display
      factor = exp ( log ( (double) maxb ) / 20.0 ) ;
      ib = lastb = 0 ;

      for ( ; ib < maxb ; ) {

         // This spaces block sizes intelligently for display
         ib = (int) (rb + 0.5) ;
         rb *= factor ;
         if (ib > maxb)
            ib = maxb ;
         if (ib == lastb)
            continue ;
         lastb = ib ;

         estimate = StdErrMeanSB ( nsamps , x , ib , nboot , bs ) ;
         diff = estimate - CorrectStdErr ;
         StdErrBiasSB[ib-1] += diff ;
         StdErrErrSB[ib-1] += diff * diff ;

         estimate = StdErrMeanTBB ( nsamps , x , ib , nboot , xinf , bs , window ) ;
         diff = estimate - CorrectStdErr ;
         StdErrBiasTBB[ib-1] += diff ;
         StdErrErrTBB[ib-1] += diff * diff ;

         estimate = QuantileMeanSB ( nsamps , x , ib , nboot , 0.1 , bs , reps ) ;
         diff = estimate - CorrectQuantile ;
         QuantileBiasSB[ib-1] += diff ;
         QuantileErrSB[ib-1] += diff * diff ;
//       if (SampleMean + estimate >= 0.0)    // Percentile method
         if (SampleMean <= estimate)          // Basic method
            QuantileRejectSB[ib-1] += 1.0 ;

         estimate = QuantileMeanTBB ( nsamps , x , ib , nboot , 0.1 , xinf , bs , reps , window ) ;
         diff = estimate - CorrectQuantile ;
         QuantileBiasTBB[ib-1] += diff ;
         QuantileErrTBB[ib-1] += diff * diff ;
//       if (SampleMean + estimate >= 0.0)    // Percentile method
         if (SampleMean <= estimate)          // Basic method
            QuantileRejectTBB[ib-1] += 1.0 ;

         } // For ib

      if (((itry % divisor) == 1)
       || (itry == ntries-1) ) {      // Don't do this every try!  Too slow.
         ndone = itry + 1 ;           // This many tries done (and in arrays)
         printf ( "\n\n\n" ) ;
         printf (
 "  b  SEbSB SEbTBB SEerrSB SEerrTBB QbSB  QbTBB  QerrSB QerrTBB  QrejSB QrejTBB" ) ;

         rb = 1.0 ;   // This spaces block sizes intelligently for display
         factor = exp ( log ( (double) maxb ) / 20.0 ) ;
         ib = lastb = 0 ;

         for ( ; ib < maxb ; ) {

            // This spaces block sizes intelligently for display
            ib = (int) (rb + 0.5) ;
            rb *= factor ;
            if (ib > maxb)
               ib = maxb ;
            if (ib == lastb)
               continue ;
            lastb = ib ;

            printf ( "\n%3d %6.3lf %6.3lf %6.3lf %6.3lf |",
               ib, StdErrBiasSB[ib-1]/ndone, StdErrBiasTBB[ib-1]/ndone,
               sqrt ( StdErrErrSB[ib-1]/ndone ), sqrt ( StdErrErrTBB[ib-1]/ndone ) ) ;
            printf ( " %6.3lf %6.3lf %6.3lf %6.3lf | %6.3lf %6.3lf",
               QuantileBiasSB[ib-1]/ndone, QuantileBiasTBB[ib-1]/ndone,
               sqrt ( QuantileErrSB[ib-1]/ndone ), sqrt ( QuantileErrTBB[ib-1]/ndone ),
               QuantileRejectSB[ib-1]/ndone, QuantileRejectTBB[ib-1]/ndone ) ;
            }
         }

      if ((itry % 10) == 1) {
         if (_kbhit ()) {
            if (_getch() == 27)
               break ;
            }
         }

      } // For all tries

   _getch () ;

/*
--------------------------------------------------------------------------------

   Now we compute the optimal block size as a function of the AR1 weight

--------------------------------------------------------------------------------
*/

   printf ( "\ncoef  SB: min    mean   max  |  TBB: min    mean   max" ) ;
   for (coef=0.0 ; coef < 0.91 ; coef+=0.1) {

      OptBminSB = OptBminTBB = nsamps ;
      OptBmaxSB = OptBmaxTBB = 0 ;
      OptBmeanSB = OptBmeanTBB = 0.0 ;

      for (itry=0 ; itry<ntries ; itry++) {
         create_AR1 ( nsamps , coef , x ) ;

         ib = optimal_SB_size ( nsamps , x , autocov ) ;
         if (ib < OptBminSB)
            OptBminSB = ib ;
         if (ib > OptBmaxSB)
            OptBmaxSB = ib ;
         OptBmeanSB += ib ;

         ib = optimal_TBB_size ( nsamps , x , autocov ) ;
         if (ib < OptBminTBB)
            OptBminTBB = ib ;
         if (ib > OptBmaxTBB)
            OptBmaxTBB = ib ;
         OptBmeanTBB += ib ;
         }

      OptBmeanSB /= ntries ;
      OptBmeanTBB /= ntries ;
      printf ( "\n%4.1lf      %3d  %6.2lf   %3d  |       %3d  %6.2lf   %3d",
         coef, OptBminSB,  OptBmeanSB,  OptBmaxSB,
               OptBminTBB, OptBmeanTBB, OptBmaxTBB ) ;
      }

   return EXIT_SUCCESS ;
}

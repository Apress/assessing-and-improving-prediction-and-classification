/******************************************************************************/
/*                                                                            */
/*  MC_TRAIN - Demonstrate Monte-Carlo permutation training                   */
/*                                                                            */
/*  This uses randomly generated credit card fraud data to train a linear     */
/*  regression model to detect fraud.                                         */
/*                                                                            */
/******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include "..\svdcmp.h"

extern double unifrand () ;
extern void qsortds ( int first , int last , double *data , double *slave ) ;

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, j, k, ncases, irep, nreps, mcpt_count, ibest, is_fraud ;
   double power, *dptr, *data, *work, *pred, coefs[4] ;
   double gain, best_gain, original_gain ;
   double inherent_bias, mean_inherent_bias, original_inherent_bias ;
   double mean_permuted_gain, training_bias ;
   double unbiased_actual_gain, unbiased_gain_above_inherent_bias ;
   double sum, dtemp, thresh, prior_thresh ;
   double p_fraud, p_legit, c_fraud, c_legit, gain_ll, gain_lf, gain_fl, gain_ff ;
   FILE *fp ;
   SingularValueDecomp *svdptr ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 4) {
      printf ( "\nUsage: MC_TRAIN  ncases  power  nreps" ) ;
      printf ( "\n  ncases - Number of cases" ) ;
      printf ( "\n  power - True power, zero for none, one for perfect" ) ;
      printf ( "\n  nreps - Number of Monte-Carlo permutations, including unpermuted" ) ;
      exit ( 1 ) ;
      }

   ncases = atoi ( argv[1] ) ;
   power = atof ( argv[2] ) ;
   nreps = atoi ( argv[3] ) ;
#else
   ncases = 100 ;
   power = 0.0 ;
   nreps = 10 ;
#endif


/*
   Open the text file to which results will be written
*/

   fp = fopen ( "MC_TRAIN.LOG" , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open MC_TRAIN.LOG file for writing!" ) ;
      return EXIT_FAILURE ;
      }

   fprintf ( fp, "Monte-Carlo training with %d cases, power=%.4lf, %d replications",
             ncases, power, nreps ) ;

/*
   Allocate scratch memory
*/

   data = (double *) malloc ( 4 * ncases * sizeof(double) ) ;
   assert ( data != NULL ) ;

   work = (double *) malloc ( ncases * sizeof(double) ) ;
   assert ( work != NULL ) ;

   pred = (double *) malloc ( ncases * sizeof(double) ) ;
   assert ( pred != NULL ) ;

   svdptr = new SingularValueDecomp ( ncases , 4 , 0 ) ;
   assert ( svdptr != NULL ) ;


/*
   Generate the random dataset.  It has four variables:
      THIS_CHARGE - Dollar value of charge being tested
      AVG_CHARGE - Average dollar value of historic charges
      FOREIGN - 1 if transaction in a foreign country, 0 if in native country
      FRAUD - 1 if fraudulent, 0 if legitimate
*/

   for (i=0 ; i<ncases ; i++) {
      is_fraud = (unifrand() < 0.01)  ?  1 : 0 ;        // True situation
      data[4*i+0] = 1000 + 1900 * (unifrand() - 0.5) ;  // THIS_CHARGE (never predictive)
      data[4*i+1] = 1000 + 500 * (unifrand() - 0.5) ;   // AVG_CHARGE (never predictive)
      data[4*i+2] = (unifrand() < 0.01)  ?  1 : 0 ;     // FOREIGN (modified below for predictive power)
      data[4*i+3] = is_fraud ;                          // FRAUD (true situation)
      data[4*i+2] = power * is_fraud + (1.0 - power) * data[4*i+2] ; // Use FOREIGN for predictive power
      }

/*
   Set the predfined gain for each possible outcome
*/

   gain_ll = 10 ;      // Legitimate and predicted legitimate
   gain_lf = -10 ;     // Legitimate and predicted fraud
   gain_fl = -1000 ;   // Fraud and predicted legitimate
   gain_ff = 500 ;     // Fraud and predicted fraud

/*
   Save the predicted variable in 'work'.
   Simultaneously compute the probability (based on occurrence)
   of fraud and legitimate charges.
   Compute the singular value decomposition of the predictors.
*/

      k = 0 ;
      dptr = svdptr->a ;

      for (i=0 ; i<ncases ; i++) {
         *dptr++ = data[4*i+0] ;   // THIS_CHARGE
         *dptr++ = data[4*i+1] ;   // AVG_CHARGE
         *dptr++ = data[4*i+2] ;   // FOREIGN
         *dptr++ = 1.0 ;           // Constant term
         work[i] = data[4*i+3] ;   // FRAUD (predicted variable)
         if (work[i] > 0.5)        // If this case is fraud
            ++k ;                  // Count it for probability
         }

      p_fraud = (double) k / ncases ;  // Probability that a case is fraud
      p_legit  = 1.0 - p_fraud ;       // And legitimate

      svdptr->svdcmp () ;


/*
   Replication loop is here
*/

   mean_inherent_bias = 0.0 ;   // Computed from permuted only
   mean_permuted_gain = 0.0 ;   // Ditto

   for (irep=0 ; irep<nreps ; irep++) {

      // Shuffle dependent variable if in permutation run (irep>0)

      if (irep) {                   // If doing permuted runs, shuffle
         i = ncases ;               // Number remaining to be shuffled
         while (i > 1) {            // While at least 2 left to shuffle
            j = (int) (unifrand () * i) ;
            if (j >= i)
               j = i - 1 ;
            dtemp = work[--i] ;
            work[i] = work[j] ;
            work[j] = dtemp ;
            }
         }

      // Fit a linear model and compute predictions.

      memcpy ( svdptr->b , work , ncases * sizeof(double) ) ;
      svdptr->backsub ( 1.e-8 , coefs ) ;

      for (i=0 ; i<ncases ; i++) {           // Find prediction for each case
         sum = coefs[3] ;                    // Constant term
         for (j=0 ; j<3 ; j++)               // Three predictors
            sum += coefs[j] * data[4*i+j] ;
         pred[i] = sum ;
         }
            
/*
   Compute the optimal threshold.
   Begin by computing the gain if all transactions are considered fraud.
   Then raise the threshold one step at a time, finding the threshold for maximum gain.
   Each time through this loop, the case at i-1 goes from being called fraud
   to being called legitimate.  So we have to undo whatever gain resulted from it
   being called fraud, and then include the gain from it being called legitimate.
*/

      gain = 0.0 ;
      for (i=0 ; i<ncases ; i++) {
         if (work[i] < 0.5)         // If this is a legitimate transaction
            gain += gain_lf ;       // Legitimate called fraud
         else                       // This is fraud
            gain += gain_ff ;       // Fraud called fraud
         }

      // The possible threshold are unique values of the predictions
      // So we sort the predictions to get the possible thresholds in ascending order

      qsortds ( 0 , ncases-1 , pred , work ) ;  // Sort predictions ascending, simultaneously moving true

      best_gain = gain ;  // Currently, this is the gain from calling all transactions fraud

      for (i=1 ; i<=ncases ; i++) {     // Try all possible thresholds, including calling all legitimate

         if (i < ncases)                // Usual situation
            thresh = pred[i] ;         
         else                           // Must include possibility of all called legitimate
            thresh = pred[i-1] + 1.0 ;  // Actual value added makes no difference; anything to make it greater

         prior_thresh = pred[i-1] ;     // This case will now change from predicted fraud to predicted legit

         if (work[i-1] < 0.5)           // If this transaction is legitimate
            gain += gain_ll - gain_lf ; // Went from called fraud to called legitimate
         else                           // This transaction is fraud
            gain += gain_fl - gain_ff ; // Went from called fraud to called legitimate

         if (thresh > prior_thresh) {   // Only update when threshold actually changes
            if (gain > best_gain) {     // (Must not break in the middle of a block of ties)
               best_gain = gain ;       // Keep track of best
               ibest = i ;              // Lets us later compute fraction classified as fraud
               }
            }
         } // For all cases, finding optimal threshold


      // The optimal threshold is now known.
      // Handle the p-value computations.

      if (irep == 0) {              // If doing original (unpermuted), save gain
         original_gain = best_gain ;
         mcpt_count = 1 ;           // Original gain equals or exceeds itself, so count it
         }

      else {
         if (best_gain >= original_gain)  // Count for p-value
            ++mcpt_count ;
         }


      // Handle the 'gain breakdown' computations

      c_fraud = (double) (ncases - ibest) / ncases ;  // Fraction of cases classified as fraud
      c_legit = 1.0 - c_fraud ;                       // Ditto legitimate
      inherent_bias  = p_legit * c_legit * gain_ll +  // Gain expected from a similar but worthless system
                       p_legit * c_fraud * gain_lf +
                       p_fraud * c_legit * gain_fl +
                       p_fraud * c_fraud * gain_ff ;

      if (irep == 0)
         original_inherent_bias = inherent_bias ; // Needed only to display for user; not really important
      else {
         mean_inherent_bias += inherent_bias ;    // These are cumulated for permutations only, not original
         mean_permuted_gain += best_gain ;
         }


      // Print stats for original model

      if (irep == 0) {              // If doing original (unpermuted), print information
         fprintf ( fp, "\n\nCoefficients:" ) ;
         fprintf ( fp, "\n   THIS_CHARGE %12.5lf", coefs[0] ) ;
         fprintf ( fp, "\n    AVG_CHARGE %12.5lf", coefs[1] ) ;
         fprintf ( fp, "\n       FOREIGN %12.5lf", coefs[2] ) ;
         fprintf ( fp, "\n      Constant %12.5lf", coefs[3] ) ;
         fprintf ( fp, "\n\nCalled fraud %d of %d  (%.2lf percent)",
                   ncases-ibest, ncases, 100.0 * c_fraud ) ;
         fprintf ( fp, "\nActual fraud %.2lf percent", 100.0 * p_fraud ) ;
         }

      }  // For all reps

/*
   Replications are done.  Print summary.
*/

   original_gain /= ncases ;                   // Make it per case, not total
   mean_inherent_bias /= nreps-1 ;
   mean_permuted_gain /= ncases * (nreps-1) ;  // Ditto
   training_bias = mean_permuted_gain - mean_inherent_bias ;
   unbiased_actual_gain = original_gain - training_bias ;
   unbiased_gain_above_inherent_bias = unbiased_actual_gain - mean_inherent_bias ;
   
   fprintf ( fp, "\n\np = %.5lf", (double) mcpt_count / nreps ) ;
   fprintf ( fp, "\n\nOriginal gain = %.5lf  with original inherent bias = %.5lf",
             original_gain, original_inherent_bias ) ;
   fprintf ( fp, "\nMean permuted gain = %.5lf", mean_permuted_gain ) ;
   fprintf ( fp, "\nMean permuted inherent bias = %.5lf", mean_inherent_bias ) ;
   fprintf ( fp, "\nTraining bias = %.5lf   (%.5lf minus %.5lf)",
             training_bias, mean_permuted_gain, mean_inherent_bias ) ;
   fprintf ( fp, "\nUnbiased actual gain = %.5lf   (%.5lf minus %.5lf)",
             unbiased_actual_gain, original_gain, training_bias ) ;
   fprintf ( fp, "\nUnbiased gain above inherent bias = %.5lf   (%.5lf minus %.5lf)",
             unbiased_gain_above_inherent_bias, unbiased_actual_gain, mean_inherent_bias ) ;

   fclose ( fp ) ;
   free ( data ) ;
   free ( work ) ;
   free ( pred ) ;
   delete svdptr ;

   printf ( "\n\nPress any key..." ) ;
   _getch () ;
   return EXIT_SUCCESS ;
}

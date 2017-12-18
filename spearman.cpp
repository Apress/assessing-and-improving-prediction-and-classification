/******************************************************************************/
/*                                                                            */
/*  SPEARMAN - Compute and test Spearman Rho correlation                      */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>

double unifrand () ;
void qsortds ( int first , int last , double *data , double *slave ) ;

/*
--------------------------------------------------------------------------------

   This is the subroutine that computes the Spearman rho

--------------------------------------------------------------------------------
*/

double spearman (  // Returns rho in range -1 to 1
   int n ,         // Input: Number of cases
   double *var1 ,  // Input: One variable
   double *var2 ,  // Input: Other variable
   double *x ,     // Work vector n long
   double *y       // Work vector n long
   )
{
   int j, k, ntied ;
   double val, x_tie_correc, y_tie_correc ;
   double dn, ssx, ssy, rank, diff, rankerr, rho ;

   // We need to rearrange input vectors, so copy them to work vectors
   // To avoid disturbing the caller
   memcpy ( x , var1 , n * sizeof(double) ) ;
   memcpy ( y , var2 , n * sizeof(double) ) ;

   // Compute ties in x, compute correction as SUM ( ties**3 - ties )
   // The following routine sorts x ascending and simultaneously moves y
   qsortds ( 0 , n-1 , x , y ) ;
   x_tie_correc = 0.0 ;
   for (j=0 ; j<n ; ) { // Convert x to ranks, cumulate tie corec
      val = x[j] ;
      for (k=j+1 ; k<n ; k++) {  // Find all ties
         if (x[k] > val)
            break ;
         }
      ntied = k - j ;
      x_tie_correc += (double) ntied * ntied * ntied - ntied ;
      rank = 0.5 * ((double) j + (double) k + 1.0) ;
      while (j < k)
         x[j++] = rank ;
      } // For each case in sorted x array

   // Now do same for y
   qsortds ( 0 , n-1 , y , x ) ;
   y_tie_correc = 0.0 ;
   for (j=0 ; j<n ; ) { // Convert y to ranks, cumulate tie corec
      val = y[j] ;
      for (k=j+1 ; k<n ; k++) {  // Find all ties
         if (y[k] > val)
            break ;
         }
      ntied = k - j ;
      y_tie_correc += (double) ntied * ntied * ntied - ntied ;
      rank = 0.5 * ((double) j + (double) k + 1.0) ;
      while (j < k)
         y[j++] = rank ;
      } // For each case in sorted y array

   // Final computations
   dn = n ;
   ssx = (dn * dn * dn - dn - x_tie_correc) / 12.0 ;
   ssy = (dn * dn * dn - dn - y_tie_correc) / 12.0 ;
   rankerr = 0.0 ;
   for (j=0 ; j<n ; j++) { // Cumulate squared rank differences
      diff = x[j] - y[j] ;
      rankerr += diff * diff ;
      }
   rho = 0.5 * (ssx + ssy - rankerr) / sqrt (ssx * ssy + 1.e-20) ;
   return rho ;
}

/*
--------------------------------------------------------------------------------

   Main program to test Spearman rho

--------------------------------------------------------------------------------
*/

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, n, discr ;
   double *x, *y, coef, std, rho, *work1, *work2 ;

/*
   Process command line parameters
*/

   if (argc != 5) {
      printf ( "\nUsage: SPEARMAN n coefficient stddev discr" ) ;
      exit ( 1 ) ;
      }

   n = atoi ( argv[1] ) ;
   coef = atof ( argv[2] ) ;
   std = atof ( argv[3] ) ;
   discr = atoi ( argv[4] ) ;

   if ((std < 0.0)  ||  (n <= 0)  ||  (discr < 0)) {
      printf ( "\nUsage: SPEARMAN n coefficient stddev discr" ) ;
      exit ( 1 ) ;
      }

   x = (double *) malloc ( n * sizeof(double) ) ;
   y = (double *) malloc ( n * sizeof(double) ) ;
   work1 = (double *) malloc ( n * sizeof(double) ) ;
   work2 = (double *) malloc ( n * sizeof(double) ) ;


   for (i=0 ; i<n ; i++) {
      x[i] = unifrand() ;
      y[i] = 2.0 * fabs ( coef ) + coef * x[i] + std * unifrand() ;
      if (discr > 0) {
         x[i] = (int) (discr * x[i]) ;
         y[i] = (int) (discr * y[i]) ;
         }
      }

   rho = spearman ( n , x , y , work1 , work2 ) ;

   printf ( "\nRho = %.3lf", rho ) ;

   free ( x ) ;
   free ( y ) ;
   free ( work1 ) ;
   free ( work2 ) ;
   return EXIT_SUCCESS ;
}

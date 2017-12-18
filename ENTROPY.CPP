/******************************************************************************/
/*                                                                            */
/*  ENTROPY - Compute the entropy of each of a set of variables               */
/*                                                                            */
/******************************************************************************/

#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <conio.h>
#include <ctype.h>
#include <stdlib.h>
#include "..\info.h"

/*
   These are defined in MEM.CPP
*/

extern int mem_keep_log ;      // Keep a log file?
extern char mem_file_name[] ;  // Log file name
extern int mem_max_used ;      // Maximum memory ever in use

int main (
   int argc ,    // Number of command line arguments (includes prog name)
   char *argv[]  // Arguments (prog name is argv[0])
   )

{
   int i, k, nbins, itype, nvars, ncases, ivar, *counts, ilow, ihigh, nb ;
   int istart, istop, ibest, *sortwork, n_indep_vars ;
   double *data, *work, *entropies, *proportional, p, max_entropy, low, high ;
   double dist, best_dist, factor, entropy ;
   short int *bins ;
   char filename[256], **names ;
   FILE *fp ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf ( "\nUsage: ENTROPY  datafile  nvars  nbins  type" ) ;
      printf ( "\n  datafile - name of the text file containing the data" ) ;
      printf ( "\n             The first line is variable names" ) ;
      printf ( "\n             Subsequent lines are the data" ) ;
      printf ( "\n             Delimiters can be space, comma, or tab" ) ;
      printf ( "\n  nvars - This many variables, starting with the first, will be tested" ) ;
      printf ( "\n  nbins - If the data is discrete, this must be at least the" ) ;
      printf ( "\n          number of bins.  It will be automatically reduced" ) ;
      printf ( "\n          to the exact number of bins." ) ;
      printf ( "\n          If the data is continuous, it specifies the number of bins" ) ;
      printf ( "\n          computed by linearly scaling the interior/exterior range." ) ;
      printf ( "\n          About 10 bins is generally good for continuous data," ) ;
      printf ( "\n          although values as high as 15 or, rarely, even 20, may be" ) ;
      printf ( "\n          appropriate if there are tens of thousands of cases." ) ;
      printf ( "\n  type - Type of data processing:" ) ;
      printf ( "\n    1 - The data is discrete" ) ;
      printf ( "\n    2 - The data is continuous, and the entire range is to be tested" ) ;
      printf ( "\n    3 - The data is continuous, and the extremes are to be truncated" ) ;
      exit ( 1 ) ;
      }

   strcpy ( filename , argv[1] ) ;
   n_indep_vars = atoi ( argv[2] ) ;
   nbins = atoi ( argv[3] ) ;
   itype = atoi ( argv[4] ) ;
#else
   strcpy ( filename , "..\\VARS.TXT" ) ;
   n_indep_vars = 8 ;
   nbins = 10 ;
   itype = 2 ;
#endif

   if (itype < 1  ||  itype > 3) {
      printf ( "\nERROR... type illegal" ) ;
      return EXIT_FAILURE ;
      }

   if (nbins < 2  ||  (itype > 1  &&  nbins < 3)) {
      printf ( "\nERROR... nbins illegal" ) ;
      return EXIT_FAILURE ;
      }

/*
   These are used by MEM.CPP for runtime memory validation
*/

   _fullpath ( mem_file_name , "MEM.LOG" , 256 ) ;
   fp = fopen ( mem_file_name , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open MEM.LOG file for writing!" ) ;
      return EXIT_FAILURE ;
      }
   fclose ( fp ) ;
   mem_keep_log = 1 ;
   mem_max_used = 0 ;

/*
   Open the text file to which results will be written
*/

   fp = fopen ( "ENTROPY.LOG" , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open ENTROPY.LOG file for writing!" ) ;
      return EXIT_FAILURE ;
      }

/*
   Read the file
*/

   if (readfile ( filename , &nvars , &names , &ncases , &data ))
      return EXIT_FAILURE ;

/*
   Allocate scratch memory

   bins - Bin ids for all variables
   counts - Count of cases in each bin
   entropies - Entropy of each variable
   proportional - Proportional entropy of each variable
   work - Temporary use for extracting a variable from the dataset
   sortwork - Temporary use for printing variable's information sorted
*/

   MEMTEXT ( "ENTROPY 6 allocs" ) ;
   bins = (short int *) MALLOC ( ncases * sizeof(short int) ) ;
   assert ( bins != NULL ) ;
   counts = (int *) MALLOC ( nbins * sizeof(int) ) ;
   assert ( counts != NULL ) ;
   entropies = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( entropies != NULL ) ;
   proportional = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( proportional != NULL ) ;
   work = (double *) MALLOC ( ncases * sizeof(double) ) ;
   assert ( work != NULL ) ;
   sortwork = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( sortwork != NULL ) ;

/*
   If splitting a continuous variable across interior range,
   compute things that will be needed.
*/

   if (itype > 1) {
      nb = nbins ;                    // Always needed
      ilow = (ncases + 1) / nb - 1 ;  // Needed only if itype==3
      if (ilow < 0)
         ilow = 0 ;
      ihigh = ncases - 1 - ilow ;
      }

/*
   If splitting a discrete variable, warn the user if the variable is continuous
*/

   else {
      for (ivar=0 ; ivar<n_indep_vars ; ivar++) {
         for (i=0 ; i<ncases ; i++)
            work[i] = data[i*nvars+ivar] ;
         qsortd ( 0 , ncases-1 , work ) ;
         k = 1 ;
         for (i=1 ; i<ncases ; i++) {
            if (work[i] > work[i-1])
               ++k ;
            }
         if (k > nbins)
            fprintf ( fp,
               "\nWARNING... %s has %d distinct values, not %d.  Results will be incorrect.",
               names[ivar], k, nbins ) ;
         }
      }

/*
   Compute and save the entropy for each variable.
   Print the results, sort them, and print them again, this time sorted.
*/

   fprintf ( fp , "\n                       Variable     Entropy   Proportional" ) ;

   for (ivar=0 ; ivar<n_indep_vars ; ivar++) {

      for (i=0 ; i<ncases ; i++)
         work[i] = data[i*nvars+ivar] ;

      if (itype == 1) {   // Discrete?
         nb = nbins ;
         partition ( ncases , work , &nb , NULL , bins ) ;
         for (i=0 ; i<nb ; i++)
            counts[i] = 0 ;
         for (i=0 ; i<ncases ; i++)
            ++counts[bins[i]] ;
         }

      else if (itype == 2) {   // Continuous, split across full range
         low = high = work[0] ;
         for (i=1 ; i<ncases ; i++) {
            if (work[i] > high)
               high = work[i] ;
            if (work[i] < low)
               low = work[i] ;
            }
         for (i=0 ; i<nb ; i++)
            counts[i] = 0 ;
         factor = (nb - 0.00000000001) / (high - low + 1.e-60) ;
         for (i=0 ; i<ncases ; i++) {
            k = (int) (factor * (work[i] - low)) ;
            ++counts[k] ;
            }
         }

      else {              // Continuous, split across interior range
         // Find the shortest interval containing 1-2/nbins of the distribution
         qsortd ( 0 , ncases-1 , work ) ;
         istart = 0 ;
         istop = istart + ihigh - ilow - 2 ;
         best_dist = 1.e60 ;
         while (istop < ncases) {  // Try bounds containing the same n of cases
            dist = work[istop] - work[istart] ;
            if (dist < best_dist) { // We're looking for the shortest
               best_dist = dist ;
               ibest = istart ;
               }
            ++istart ;  // Advance to the next interval
            ++istop ;
            }
         istart = ibest ;  // This is the shortest interval
         istop = istart + ihigh - ilow - 2 ;
         counts[0] = istart ;                 // The count of the leftmost bin
         counts[nb-1] = ncases - istop - 1 ;  // and rightmost are implicit
         for (i=1 ; i<nb-1 ; i++)
            counts[i] = 0 ;
         low = work[istart] ;
         high = work[istop] ;
         factor = (nb - 2.00000000001) / (high - low + 1.e-60) ;
         for (i=istart ; i<=istop ; i++) {
            k = (int) (factor * (work[i] - low)) ;
            ++counts[k+1] ;
            }
         } // Else continuous
//!!!!!!!!!!!!!!!!!!!!!!
      printf ( "\n\n%s bins:", names[ivar] ) ;
      for (i=0 ; i<nb ; i++)
         printf ( " %d", counts[i] ) ;
//!!!!!!!!!!!!!!!!!!!!!!

      // Compute the entropy
      max_entropy = log ( (double) (nb) ) ;
      entropy = 0.0 ;
      for (i=0 ; i<nb ; i++) {
         if (counts[i] > 0) {
            p = (double) counts[i] / (double) ncases ;
            entropy -= p * log(p) ;
            }
         }

      sortwork[ivar] = ivar ;
      entropies[ivar] = entropy ;
      proportional[ivar] = entropy / max_entropy ;

      fprintf ( fp , "\n%31s %11.5lf %11.5lf",
                names[ivar], entropy, proportional[ivar] ) ;
      }

   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\nEntropies, in decreasing order" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n                       Variable     Entropy    Relative" ) ;

   qsortdsi ( 0 , n_indep_vars-1 , entropies , sortwork ) ;
   for (ivar=0 ; ivar<n_indep_vars ; ivar++) { // Do all candidates
      k = sortwork[n_indep_vars-1-ivar] ;      // Index of sorted candidate
      fprintf ( fp , "\n%31s %11.5lf %11.5lf",
                names[k], entropies[n_indep_vars-1-ivar], proportional[k] ) ;
      }

   MEMTEXT ( "ENTROPY: Finish... 6 arrays plus free_data()" ) ;
   fclose ( fp ) ;
   FREE ( bins ) ;
   FREE ( counts ) ;
   FREE ( entropies ) ;
   FREE ( proportional ) ;
   FREE ( work ) ;
   FREE ( sortwork ) ;
   free_data ( nvars , names , data ) ;
   MEMCLOSE () ;
   printf ( "\n\nPress any key..." ) ;
   _getch () ;
   return EXIT_SUCCESS ;
}

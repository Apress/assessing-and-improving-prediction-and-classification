/******************************************************************************/
/*                                                                            */
/*  MI_ONLY - Mutual information ONLY for continuous predicted and predictors */
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
   int i, j, k, nvars, ncases, irep, nreps, ivar, nties, ties ;
   int n_indep_vars, idep, icand, *index, *mcpt_max_counts, *mcpt_same_counts, *mcpt_solo_counts ;
   double *data, *work, dtemp, *save_info, criterion, *crits ;
   char filename[256], **names, depname[256] ;
   FILE *fp ;
   MutualInformationAdaptive *mi_adapt ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 5) {
      printf ( "\nUsage: MI_ONLY  datafile  n_indep  depname  nreps" ) ;
      printf ( "\n  datafile - name of the text file containing the data" ) ;
      printf ( "\n             The first line is variable names" ) ;
      printf ( "\n             Subsequent lines are the data." ) ;
      printf ( "\n             Delimiters can be space, comma, or tab" ) ;
      printf ( "\n  n_indep - Number of independent vars, starting with the first" ) ;
      printf ( "\n  depname - Name of the 'dependent' variable" ) ;
      printf ( "\n            It must be AFTER the first n_indep variables" ) ;
      printf ( "\n  nreps - Number of Monte-Carlo permutations, including unpermuted" ) ;
      exit ( 1 ) ;
      }

   strcpy ( filename , argv[1] ) ;
   n_indep_vars = atoi ( argv[2] ) ;
   strcpy ( depname , argv[3] ) ;
   nreps = atoi ( argv[4] ) ;
#else
   strcpy ( filename , "..\\SYNTH.TXT" ) ;
   n_indep_vars = 7 ;
   strcpy ( depname , "SUM1234" ) ;
   nreps = 100 ;
#endif

   _strupr ( depname ) ;

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
   mem_keep_log = 0 ;  // Change this to 1 to keep a memory use log (slows execution!)
   mem_max_used = 0 ;

/*
   Open the text file to which results will be written
*/

   fp = fopen ( "MI_ONLY.LOG" , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open MI_ONLY.LOG file for writing!" ) ;
      return EXIT_FAILURE ;
      }

/*
   Read the file and locate the index of the dependent variable
*/

   if (readfile ( filename , &nvars , &names , &ncases , &data ))
      return EXIT_FAILURE ;

   for (idep=0 ; idep<nvars ; idep++) {
      if (! strcmp ( depname , names[idep] ))
         break ;
      }

   if (idep == nvars) {
      printf ( "\nERROR... Dependent variable %s is not in file", depname ) ;
      return EXIT_FAILURE ;
      }

   if (idep < n_indep_vars) {
      printf ( "\nERROR... Dependent variable %s must be beyond independent vars",
               depname ) ;
      return EXIT_FAILURE ;
      }

/*
   Check each variable for ties.  This is not needed for the algorithm,
   but it is good to warn the user, because more than a very few tied values
   in any variable seriously degrades performance of the adaptive partitioning algorithm.
*/

   MEMTEXT ( "MI_ONLY: Work" ) ;
   work = (double *) MALLOC ( ncases * sizeof(double) ) ;
   assert ( work != NULL ) ;

   ties = 0 ;
   assert ( work != NULL ) ;
   for (ivar=0 ; ivar<nvars ; ivar++) {
      if (ivar > n_indep_vars  &&  ivar != idep)
         continue ; // Check only the variables selected by the user
      for (i=0 ; i<ncases ; i++)
         work[i] = data[i*nvars+ivar] ;
      qsortd ( 0 , ncases-1 , work ) ;
      nties = 0 ;
      for (i=1 ; i<ncases ; i++) {
         if (work[i] == work[i-1])
            ++nties ;
         }
      if ((double) nties / (double) ncases > 0.05) {
         ++ties ;
         fprintf ( fp , "\nWARNING... %s has %.2lf percent ties!",
                   names[ivar], 100.0 * nties / (double) ncases ) ;
         }
      } // For all variables
   if (ties) {
      fprintf ( fp , "\nThe presence of ties will seriously degrade" ) ;
      fprintf ( fp , "\nperformance of the adaptive partitioning algorithm\n\n" ) ;
      }

/*
   Allocate scratch memory and create the MutualInformation object using the
   dependent variable

   crits - Mutual information criterion
   index - Indices that sort the criterion
   save_info - Ditto, this is univariate information, to be sorted
   mi_adapt - The MutualInformation object, constructed with the 'dependent' variable
*/

   MEMTEXT ( "MI_ONLY work allocs plus MutualInformation" ) ;
   crits = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( crits != NULL ) ;
   index = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( index != NULL ) ;
   mcpt_max_counts = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( mcpt_max_counts != NULL ) ;
   mcpt_same_counts = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( mcpt_same_counts != NULL ) ;
   mcpt_solo_counts = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( mcpt_solo_counts != NULL ) ;
   save_info = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( save_info != NULL ) ;

   for (irep=0 ; irep<nreps ; irep++) {

      for (i=0 ; i<ncases ; i++)            // Get the 'dependent' variable
         work[i] = data[i*nvars+idep] ;

//    Shuffle dependent variable if in permutation run (irep>0)

      if (irep) {                   // If doing permuted runs, shuffle
         i = ncases ;              // Number remaining to be shuffled
         while (i > 1) {            // While at least 2 left to shuffle
            j = (int) (unifrand () * i) ;
            if (j >= i)
               j = i - 1 ;
            dtemp = work[--i] ;
            work[i] = work[j] ;
            work[j] = dtemp ;
            }
         }

      // Here we use a tiny split theshold (instead of the usual 6.0) so that it picks up
      // small amounts of mutual information (perhaps including noise).
      // If we used 6.0, nearly all permutations of any reasonably sized dataset
      // would have a computed mutual information of zero.  It's safe picking up
      // some noise because the permutation test will account for this.

      mi_adapt = new MutualInformationAdaptive ( ncases , work , 1 , 0.1 ) ; // Deliberately tiny for low information
      assert ( mi_adapt != NULL ) ;

/*
   Compute and save the mutual information for the dependent variable
   with each individual independent variable candidate.
*/

      for (icand=0 ; icand<n_indep_vars ; icand++) { // Try all candidates
         for (i=0 ; i<ncases ; i++)
            work[i] = data[i*nvars+icand] ;

         criterion = mi_adapt->mut_inf ( work , 1 ) ;

         save_info[icand] = criterion ; // We will sort this when all candidates are done
                                        
         if (irep == 0) {               // If doing original (unpermuted), save criterion
            index[icand] = icand ;      // Will need original indices when criteria are sorted
            crits[icand] = criterion ;
            mcpt_max_counts[icand] = mcpt_same_counts[icand] = mcpt_solo_counts[icand] = 1 ;  // This is >= itself so count it now
            }
         else {
            if (criterion >= crits[icand])
               ++mcpt_solo_counts[icand] ;
            }
         } // Initial list of all candidates

      delete mi_adapt ;
      mi_adapt = NULL ;

      if (irep == 0)  // Find the indices that sort the candidates per criterion
         qsortdsi ( 0 , n_indep_vars-1 , save_info , index ) ;

      else {
         qsortd ( 0 , n_indep_vars-1 , save_info ) ;
         for (icand=0 ; icand<n_indep_vars ; icand++) {
            if (save_info[icand] >= crits[index[icand]])
               ++mcpt_same_counts[index[icand]] ;
            if (save_info[n_indep_vars-1] >= crits[index[icand]]) // Valid only for largest
               ++mcpt_max_counts[index[icand]] ;
            }
         }

      }  // For all reps

   fprintf ( fp , "\nAdaptive partitioning mutual information of %s", depname);

   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\nPredictors, in order of decreasing mutual information" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n                       Variable   Information   Solo pval   Min pval   Max pval" ) ;

   for (icand=0 ; icand<n_indep_vars ; icand++) { // Do all candidates
      k = index[n_indep_vars-1-icand] ;           // Index of sorted candidate
      fprintf ( fp , "\n%31s %11.5lf %12.4lf %10.4lf %10.4lf", names[k], crits[k],
                (double) mcpt_solo_counts[k] / nreps,
                (double) mcpt_same_counts[k] / nreps,
                (double) mcpt_max_counts[k] / nreps ) ;
      }

   MEMTEXT ( "MI_ONLY: Finish" ) ;
   fclose ( fp ) ;
   FREE ( work ) ;
   FREE ( crits ) ;
   FREE ( index ) ;
   FREE ( mcpt_max_counts ) ;
   FREE ( mcpt_same_counts ) ;
   FREE ( mcpt_solo_counts ) ;
   FREE ( save_info ) ;
   free_data ( nvars , names , data ) ;

   MEMCLOSE () ;
   printf ( "\n\nPress any key..." ) ;
   _getch () ;
   return EXIT_SUCCESS ;
}

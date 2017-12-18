/******************************************************************************/
/*                                                                            */
/*  MI_CONT - Mutual information for continuous predicted and predictors      */
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
   int i, j, k, nvars, ncases, ndiv, maxkept, ivar, nties, ties ;
   int n_indep_vars, idep, icand, iother, ibest, *sortwork, nkept, *kept ;
   double *data, *work ;
   double *save_info, *univar_info, *pair_info, bestredun, redun, bestcrit ;
   double criterion, relevance, redundancy, *crits, *reduns ;
   char filename[256], **names, depname[256] ;
   char trial_name[256], *pair_found ;
   FILE *fp ;
   MutualInformationParzen *mi_parzen ;
   MutualInformationAdaptive *mi_adapt ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 6) {
      printf ( "\nUsage: MI_CONT  datafile  n_indep  depname  ndiv  maxkept" ) ;
      printf ( "\n  datafile - name of the text file containing the data" ) ;
      printf ( "\n             The first line is variable names" ) ;
      printf ( "\n             Subsequent lines are the data." ) ;
      printf ( "\n             Delimiters can be space, comma, or tab" ) ;
      printf ( "\n  n_indep - Number of independent vars, starting with the first" ) ;
      printf ( "\n  depname - Name of the 'dependent' variable" ) ;
      printf ( "\n            It must be AFTER the first n_indep variables" ) ;
      printf ( "\n  ndiv - Normally zero, to employ adaptive partitioning" ) ;
      printf ( "\n         Specify 5 (for very few cases) to 15 (for an" ) ;
      printf ( "\n         enormous number of cases) to use Parzen windows" ) ;
      printf ( "\n  maxkept - Stepwise will allow at most this many predictors" ) ;
      exit ( 1 ) ;
      }

   strcpy ( filename , argv[1] ) ;
   n_indep_vars = atoi ( argv[2] ) ;
   strcpy ( depname , argv[3] ) ;
   ndiv = atoi ( argv[4] ) ;
   maxkept = atoi ( argv[5] ) ;
#else
   strcpy ( filename , "..\\VARS.TXT" ) ;
   n_indep_vars = 8 ;
   strcpy ( depname , "DAY_RETURN" ) ;
   ndiv = 0 ;
   maxkept = 5 ;
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
   mem_keep_log = 1 ;
   mem_max_used = 0 ;

/*
   Open the text file to which results will be written
*/

   fp = fopen ( "MI_CONT.LOG" , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open MI_CONT.LOG file for writing!" ) ;
      return EXIT_FAILURE ;
      }

/*
   Read the file and locate the index of the 'dependent' variable
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
   If adaptive partitioning is specified, check each variable for ties.
   This is not needed for the algorithm, but it is good to warn the
   user, because more than a very few tied values in any variable seriously
   degrades performance of the adaptive partitioning algorithm.
*/

   MEMTEXT ( "MI_CONT: Work" ) ;
   work = (double *) MALLOC ( ncases * sizeof(double) ) ;
   assert ( work != NULL ) ;

   if (ndiv == 0) {  // If adaptive partitioning, check for ties
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
      } // If adaptive partitioning, so testing for ties in the data

/*
   Allocate scratch memory and create the MutualInformation object using the
   dependent variable

   kept - Array of indices of variables kept so far
   crits - Ditto, criterion
   reduns - Ditto, redundancy
   sortwork - Temporary use for printing variable's information sorted
   save_info - Ditto, this is univariate information, to be sorted
   univar_info - Also univariate information, but not sorted, for use in stepwise
   pair_found - Flag: is there valid info in the corresponding element of the next array
   pair_info - Preserve pairwise information of indeps to avoid expensive recalculation
   mi_parzen - The MutualInformation object, constructed with the 'dependent' variable
   mi_adapt - Ditto, but used if adaptive partitioning
*/

   MEMTEXT ( "MI_CONT 6 allocs plus MutualInformation" ) ;
   kept = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( kept != NULL ) ;
   crits = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( crits != NULL ) ;
   reduns = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( reduns != NULL ) ;
   sortwork = (int *) MALLOC ( n_indep_vars * sizeof(int) ) ;
   assert ( sortwork != NULL ) ;
   save_info = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( save_info != NULL ) ;
   univar_info = (double *) MALLOC ( n_indep_vars * sizeof(double) ) ;
   assert ( univar_info != NULL ) ;
   pair_found = (char *) MALLOC ( (n_indep_vars * (n_indep_vars+1) / 2) * sizeof(char) ) ;
   assert ( pair_found != NULL ) ;
   pair_info = (double *) MALLOC ( (n_indep_vars * (n_indep_vars+1) / 2) * sizeof(double) ) ;
   assert ( pair_info != NULL ) ;

   for (i=0 ; i<ncases ; i++)            // Get the 'dependent' variable
      work[i] = data[i*nvars+idep] ;

   if (ndiv > 0) {
      mi_parzen = new MutualInformationParzen ( ncases , work , ndiv ) ;
      mi_adapt = NULL ;
      assert ( mi_parzen != NULL ) ;
      }
   else {
      mi_adapt = new MutualInformationAdaptive ( ncases , work , 0 , 6.0 ) ;
      mi_parzen = NULL ;
      assert ( mi_adapt != NULL ) ;
      }

   memset ( pair_found , 0 , (n_indep_vars * (n_indep_vars+1) / 2) * sizeof(char) ) ;

   if (ndiv > 0)
      fprintf ( fp , "\nParzen mutual information of %s (ndiv=%d)", depname, ndiv);
   else
      fprintf ( fp , "\nAdaptive partitioning mutual information of %s", depname);

   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n---------------------------------------------------------------" ) ;
   fprintf ( fp , "\n" ) ;

/*
   Compute and save the mutual information for the dependent variable with
   each individual independent variable candidate.  Print the results,
   sort them, and print them again, this time sorted.
*/

   fprintf ( fp , "\nInitial candidates, in order of appearance in data file" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n                       Variable   Information" ) ;

   for (icand=0 ; icand<n_indep_vars ; icand++) { // Try all candidates
      for (i=0 ; i<ncases ; i++)
         work[i] = data[i*nvars+icand] ;

      if (ndiv > 0)
         criterion = mi_parzen->mut_inf ( work ) ;
      else
         criterion = mi_adapt->mut_inf ( work , 0 ) ;

      printf ( "\n%s = %.5lf", names[icand], criterion ) ;
      fprintf ( fp , "\n%31s   %.5lf", names[icand], criterion ) ;

      sortwork[icand] = icand ;
      save_info[icand] = univar_info[icand] = criterion ;
      } // Initial list of all candidates

   if (mi_parzen != NULL) {
      delete mi_parzen ;
      mi_parzen = NULL ;
      }
   if (mi_adapt != NULL) {
      delete mi_adapt ;
      mi_adapt = NULL ;
      }

   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\nInitial candidates, in order of decreasing mutual information" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n                       Variable   Information" ) ;

   qsortdsi ( 0 , n_indep_vars-1 , save_info , sortwork ) ;
   for (icand=0 ; icand<n_indep_vars ; icand++) { // Do all candidates
      k = sortwork[n_indep_vars-1-icand] ;        // Index of sorted candidate
      fprintf ( fp , "\n%31s   %.5lf", names[k], save_info[n_indep_vars-1-icand] ) ;
      }

/*
   Initialize the 'kept' set to be the best variable, and then begin the
   main outer loop that adds variables one at a time
*/

   kept[0] = sortwork[n_indep_vars-1] ;  // Index of best single candidate
   crits[0] = save_info[n_indep_vars-1] ;
   reduns[0] = 0.0 ;
   nkept = 1 ;

   if (maxkept > n_indep_vars)  // Guard against silly user
      maxkept = n_indep_vars ;

   while (nkept < maxkept) {

      fprintf ( fp , "\n" ) ;
      fprintf ( fp , "\nVariables so far                 Relevance  Redundancy  Criterion" ) ;
      for (i=0 ; i<nkept ; i++)
         fprintf ( fp , "\n%31s %10.5lf %10.5lf %10.5lf",
                   names[kept[i]], crits[i] + reduns[i], reduns[i], crits[i] ) ;
      fprintf ( fp , "\n" ) ;
      fprintf ( fp , "\nSearching for an additional candidate..." ) ;
      fprintf ( fp , "\n" ) ;
      fprintf ( fp , "\n                       Variable  Relevance  Redundancy  Criterion" ) ;

      bestcrit = -1.e60 ;
      for (icand=0 ; icand<n_indep_vars ; icand++) { // Try all candidates
         for (i=0 ; i<nkept ; i++) {  // Is this candidate already kept?
            if (kept[i] == icand)
               break ;
            }
         if (i < nkept)  // If this candidate 'icand' is already kept
            continue ;   // Skip it

         strcpy ( trial_name , names[icand] ) ;   // Its name for printing
         for (i=0 ; i<ncases ; i++)               // Get its cases
            work[i] = data[i*nvars+icand] ;

         if (ndiv > 0) {
            mi_parzen = new MutualInformationParzen ( ncases , work , ndiv ) ;
            mi_adapt = NULL ;
            assert ( mi_parzen != NULL ) ;
            }
         else {
            mi_adapt = new MutualInformationAdaptive ( ncases , work , 0 , 6.0 ) ;
            mi_parzen = NULL ;
            assert ( mi_adapt != NULL ) ;
            }

         relevance = univar_info[icand] ; // We saved it during initial printing
         printf ( "\n%s relevance = %.5lf", trial_name, relevance ) ;

         // Compute the redundancy of this candidate
         // This is the mean of its redundancy with all kept variables
         redundancy = 0.0 ;
         for (iother=0 ; iother<nkept ; iother++) {  // Process entire kept set
            j = kept[iother] ;           // Index of a variable in the kept set
            if (icand > j)               // pair_found and pair_info are
               k = icand*(icand+1)/2+j ; // symmetric, so k is the index
            else                         // into them
               k = j*(j+1)/2+icand ;
            if (pair_found[k])           // If we already computed it
               redun = pair_info[k] ;    // Don't do it again
            else {                       // First time for this pair, so compute
               for (i=0 ; i<ncases ; i++)       // Get its cases
                  work[i] = data[i*nvars+j] ;   // Variable already in kept set
               if (ndiv > 0)
                  redun = mi_parzen->mut_inf ( work ) ;
               else
                  redun = mi_adapt->mut_inf ( work , 0 ) ;
               pair_found[k] = 1 ;       // Flag that this pair has been computed
               pair_info[k] = redun ;    // And save the MI for this pair
               } // Else must compute redundancy
            redundancy += redun ;
            printf ( "\n  %s <-> %s redundancy = %.5lf", names[icand], names[j], redun ) ;
            } // For all kept variables, computing mean redundancy

         if (mi_parzen != NULL) {
            delete mi_parzen ;
            mi_parzen = NULL ;
            }

         if (mi_adapt != NULL) {
            delete mi_adapt ;
            mi_adapt = NULL ;
            }

         redundancy /= nkept ;  // It is the mean across all kept
         printf ( "\nRedundancy = %.5lf", redundancy ) ;

         criterion = relevance - redundancy ;
         fprintf ( fp , "\n%31s %10.5lf %10.5lf %10.5lf",
                   trial_name, relevance, redundancy, criterion ) ;

         if (criterion > bestcrit) { // Did we just set a new record?
            bestcrit = criterion ;   // If so, update the record
            bestredun = redundancy ; // Needed for printing results later
            ibest = icand ;          // Keep track of the winning candidate
            }

         } // For all candidates

      // We now have the best candidate
      if (bestcrit <= 0.0)
         break ;
      kept[nkept] = ibest ;
      crits[nkept] = bestcrit ;
      reduns[nkept] = bestredun ;
      ++nkept ;
      } // While adding new variables

   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\nFinal set                        Relevance  Redundancy  Criterion" ) ;
   for (i=0 ; i<nkept ; i++)
      fprintf ( fp , "\n%31s %10.5lf %10.5lf %10.5lf",
                names[kept[i]], crits[i] + reduns[i], reduns[i], crits[i] ) ;


   MEMTEXT ( "MI_CONT: Finish" ) ;
   fclose ( fp ) ;
   FREE ( work ) ;
   FREE ( kept ) ;
   FREE ( crits ) ;
   FREE ( reduns ) ;
   FREE ( sortwork ) ;
   FREE ( save_info ) ;
   FREE ( univar_info ) ;
   FREE ( pair_found ) ;
   FREE ( pair_info ) ;
   if (mi_parzen != NULL)
      delete mi_parzen ;
   if (mi_adapt != NULL)
      delete mi_adapt ;
   free_data ( nvars , names , data ) ;

   MEMCLOSE () ;
   printf ( "\n\nPress any key..." ) ;
   _getch () ;
   return EXIT_SUCCESS ;
}

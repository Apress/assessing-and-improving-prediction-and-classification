/******************************************************************************/
/*                                                                            */
/*  MI_DISC - Mutual information for discrete predicted and predictors        */
/*                                                                            */
/******************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE

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
   int i, j, k, nvars, ncases, maxkept, ivar ;
   int n_indep_vars, idep, icand, iother, ibest, *sortwork, nkept, *kept ;
   int nbins_dep, nbins_indep, maxbins ;
   short int *bins_dep, *bins_indep ;
   double *data, *work ;
   double *save_info, *univar_info, *pair_info, redun, bestcrit, bestredun ;
   double criterion, entropy, bound, relevance, redundancy, *crits, *reduns ;
   char filename[256], **names, depname[256] ;
   char trial_name[256], *pair_found ;
   FILE *fp ;
   MutualInformationDiscrete *mi ;

/*
   Process command line parameters
*/

#if 1
   if (argc != 7) {
      printf ( "\nUsage: MI_DISC  datafile  n_indep  depname  nbins_dep  nbins_indep  maxkept" ) ;
      printf ( "\n  datafile - name of the text file containing the data" ) ;
      printf ( "\n             The first line is variable names" ) ;
      printf ( "\n             Subsequent lines are the data." ) ;
      printf ( "\n             Delimiters can be space, comma, or tab" ) ;
      printf ( "\n  n_indep - Number of independent vars, starting with the first" ) ;
      printf ( "\n  depname - Name of the 'dependent' variable" ) ;
      printf ( "\n            It must be AFTER the first n_indep variables" ) ;
      printf ( "\n  nbins_dep - Maximum number of dependent variable partitions" ) ;
      printf ( "\n        If specified as zero, two bins are defined (>0 and <=0)" ) ;
      printf ( "\n  nbins_indep - Ditto, but for independent variables" ) ;
      printf ( "\n        If specified as zero, two bins are defined (>0 and <=0)" ) ;
      printf ( "\n  maxkept - Stepwise will allow at most this many predictors" ) ;
      exit ( 1 ) ;
      }

   strcpy ( filename , argv[1] ) ;
   n_indep_vars = atoi ( argv[2] ) ;
   strcpy ( depname , argv[3] ) ;
   nbins_dep = atoi ( argv[4] ) ;
   nbins_indep = atoi ( argv[5] ) ;
   maxkept = atoi ( argv[6] ) ;
#else
   strcpy ( filename , "..\\VARS.TXT" ) ;
   strcpy ( depname , "DAY_RETURN" ) ;
   n_indep_vars = 8 ;
   nbins_indep = 2 ;
   nbins_dep = 0 ;
   maxkept = 99 ;
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

   fp = fopen ( "MI_DISC.LOG" , "wt" ) ;
   if (fp == NULL) { // Should never happen
      printf ( "\nCannot open MI_DISC.LOG file for writing!" ) ;
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
   Allocate scratch memory

   bins_dep - Bin ids for the dependent variable
   bins_indep - Bin ids for the independent variables
   kept - Array of indices of variables kept so far
   crits - Ditto, criterion
   reduns - Ditto, redundancy
   sortwork - Temporary use for printing variable's information sorted
   save_info - Ditto, this is univariate information, to be sorted
   univar_info - Also univariate information, but not sorted, for use in stepwise
   pair_found - Flag: is there valid info in the corresponding element of the next array
   pair_info - Preserve pairwise information of indeps to avoid expensive recalculation
*/

   MEMTEXT ( "MI_DISC 9 allocs" ) ;
   work = (double *) MALLOC ( ncases * sizeof(double) ) ;
   assert ( work != NULL ) ;
   bins_dep = (short int *) MALLOC ( ncases * sizeof(short int) ) ;
   assert ( bins_dep != NULL ) ;
   bins_indep = (short int *) MALLOC ( ncases * n_indep_vars * sizeof(short int) ) ;
   assert ( bins_indep != NULL ) ;
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

/*
   Compute the bin membership of all variables.
   If the user specified a number of bins as zero, we treat the variable
   as binary (two bins) using <=0 and >0 as the definition of bin membership.
*/

   if (nbins_dep == 0) {   // The dependent variable is binary
      nbins_dep = 2 ;
      for (i=0 ; i<ncases ; i++) {
         if (data[i*nvars+idep] > 0.0)
            bins_dep[i] = (short int) 1 ;
         else
            bins_dep[i] = (short int) 0 ;
         }
      fprintf ( fp , "\n%s has been given a binary split", names[idep] ) ;
      }
   else {                  // The dependent variable is to be partitioned
      for (i=0 ; i<ncases ; i++)
         work[i] = data[i*nvars+idep] ;
      partition ( ncases , work , &nbins_dep , NULL , bins_dep ) ;
      fprintf ( fp , "\n%s has been partitioned into %d bins",
                names[idep], nbins_dep ) ;
      }

   if (nbins_indep == 0) {
      fprintf ( fp , "\nIndependent variables have been given a binary split");
      nbins_indep = maxbins = 2 ;
      for (ivar=0 ; ivar<n_indep_vars ; ivar++) {
         for (i=0 ; i<ncases ; i++) {
            if (data[i*nvars+ivar] > 0.0)
               bins_indep[ivar*ncases+i] = (short int) 1 ;
            else
               bins_indep[ivar*ncases+i] = (short int) 0 ;
            }
         }
      }
   else {
      maxbins = 0 ;
      for (ivar=0 ; ivar<n_indep_vars ; ivar++) {
         for (i=0 ; i<ncases ; i++)
            work[i] = data[i*nvars+ivar] ;
         k = nbins_indep ;
         partition ( ncases , work , &k , NULL , bins_indep+ivar*ncases ) ;
         fprintf( fp, "\n%s has been partitioned into %d bins", names[ivar], k);
         if (k > maxbins)
            maxbins = k ;
         }
      }

/*
   The array (actually, half of a symmetric matrix) flags whether the MI of
   a pair of variables has already been found.  This avoids expensive recalcs.
*/

   memset ( pair_found , 0 , (n_indep_vars * (n_indep_vars+1) / 2) * sizeof(char) ) ;

/*
   Compute and save the mutual information for the dependent variable with
   each individual independent variable candidate.  Print the results,
   sort them, and print them again, this time sorted.
*/

   mi = new MutualInformationDiscrete ( ncases , bins_dep ) ;
   assert ( mi != NULL ) ;

   entropy = mi->entropy() ;
   fprintf ( fp , "\n\n\nMutual information of %s  (Entropy = %.4lf)",
             depname, entropy ) ;

   fprintf ( fp , "\n\nInitial candidates, in order of appearance in data file" ) ;
   fprintf ( fp , "\n" ) ;
   fprintf ( fp , "\n                       Variable   Information   Fano's bound" ) ;

   for (icand=0 ; icand<n_indep_vars ; icand++) { // Try all candidates
      criterion = mi->mut_inf ( bins_indep + icand * ncases ) ;
      if (nbins_dep <= 2)
         bound = (entropy - criterion - log ( 2.0 )) / log ( (double) nbins_dep ) ;
      else
         bound = (entropy - criterion - log ( 2.0 )) / log ( nbins_dep - 1.0 ) ;
      if (bound < 0.0)
         bound = 0.0 ;
      printf ( "\n%s = %.5lf  (%.5lf)", names[icand], criterion, bound ) ;
      fprintf ( fp , "\n%31s %11.5lf  %13.5lf",
                names[icand], criterion, bound ) ;
      sortwork[icand] = icand ;
      save_info[icand] = univar_info[icand] = criterion ;
      } // Initial list of all candidates

   delete mi ;

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
         relevance = univar_info[icand] ; // We saved it during initial printing
         printf ( "\n%s relevance = %.5lf", trial_name, relevance ) ;

         // Compute the redundancy of this candidate
         // This is the mean of its redundancy with all kept variables
         redundancy = 0.0 ;
         mi = new MutualInformationDiscrete( ncases , bins_indep + icand * ncases);
         assert ( mi != NULL ) ;
         for (iother=0 ; iother<nkept ; iother++) {  // Process entire kept set
            j = kept[iother] ;           // Index of a variable in the kept set
            if (icand > j)               // pair_found and pair_info are
               k = icand*(icand+1)/2+j ; // symmetric, so k is the index
            else                         // into them
               k = j*(j+1)/2+icand ;
            if (pair_found[k])           // If we already computed it
               redun = pair_info[k] ;    // Don't do it again
            else {                       // First time for this pair, so compute
               redun = mi->mut_inf ( bins_indep + j * ncases ) ;
               pair_found[k] = 1 ;       // Flag that this pair has been computed
               pair_info[k] = redun ;    // And save the MI for this pair
               } // Else must compute redundancy
            redundancy += redun ;
            printf ( "\n  %s <-> %s redundancy = %.5lf", names[icand], names[j], redun ) ;
            } // For all kept variables, computing mean redundancy

         delete mi ;

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


   MEMTEXT ( "MI_DISC: Finish... 11 arrays" ) ;
   fclose ( fp ) ;
   FREE ( work ) ;
   FREE ( bins_dep ) ;
   FREE ( bins_indep ) ;
   FREE ( kept ) ;
   FREE ( crits ) ;
   FREE ( reduns ) ;
   FREE ( sortwork ) ;
   FREE ( save_info ) ;
   FREE ( univar_info ) ;
   FREE ( pair_found ) ;
   FREE ( pair_info ) ;
   free_data ( nvars , names , data ) ;
   MEMCLOSE () ;
   printf ( "\n\nPress any key..." ) ;
   _getch () ;
   return EXIT_SUCCESS ;
}

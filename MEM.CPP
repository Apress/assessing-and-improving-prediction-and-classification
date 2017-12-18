/******************************************************************************/
/*                                                                            */
/*  MEM.CPP                                                                   */
/*                                                                            */
/*  This set of routines is used to intercept all calls to malloc, realloc,   */
/*  and free, checking the legality of their actions.                         */
/*  The routine memtext is also here to allow runtime commenting of the       */
/*  log file.                                                                 */
/*  Finally, memclose should be called at program completion to verify that   */
/*  no memory is still dangling.                                              */
/*                                                                            */
/*  To bypass the code given here, go to the global header file for the       */
/*  program and change #define MALLOC memalloc to #define MALLOC malloc etc.  */
/*                                                                            */
/******************************************************************************/

#define _CRT_SECURE_NO_DEPRECATE


#include <assert.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <windows.h>
#include <winbase.h>     // These headers are for GlobalAlloc etc.
#include "info.h"

#define USE_MALLOC 0
#define MAX_ALLOCS 65536

#define DEBUG_PRE_POST 1

/*
   These three globals must be initialized in the main program
*/

int mem_keep_log = 0 ;         // Keep a log file?
char mem_file_name[256] = "" ; // Log file name
int mem_max_used=0 ;           // Maximum memory ever in use

static int nallocs=0 ;                   // Number of allocations
static int allocs[MAX_ALLOCS] ;          // Ptrs to them (guaranteed mod 8) here
static int actual[MAX_ALLOCS] ;          // Actual ptrs
static int preptr[MAX_ALLOCS] ;          // Pointer to validation area before
static int postptr[MAX_ALLOCS] ;         // and after the allocated block
static unsigned alloc_size[MAX_ALLOCS] ; // Size of those allocs
static int total_use=0 ;                 // Total bytes allocated
static FILE *fp_rec ;                    // File pointer for recording actions

void *memalloc ( unsigned n )
{       
   void *ptr, *ptr8, *pre, *post ;
   union {
      void *dptr ;
      int cptr ;
      } uptr ;

   if (n == 0) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: memalloc called with length=0" ) ;
         fclose ( fp_rec ) ;
         }
      return NULL ;
      }

   if (nallocs == MAX_ALLOCS) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: memalloc MAX_ALLOCS exceeded" ) ;
         fclose ( fp_rec ) ;
         }
      return NULL ;
      }

#if USE_MALLOC
   ptr = (void *) malloc ( n + 3 * 8 ) ;
#else
   ptr = (void *) GlobalAlloc ( 0 , n + 3 * 8 ) ;
#endif

   ptr8 = NULL ;

   if (ptr != NULL) {
      if (((int) ptr) > 0) {
         uptr.dptr = ptr ;
         uptr.cptr = (uptr.cptr + 7) / 8 * 8 ;
         pre = uptr.dptr ;
         uptr.cptr = uptr.cptr + 8 ;
         ptr8 = uptr.dptr ;
         uptr.cptr = uptr.cptr + n ;
         post = uptr.dptr ;
         // Place unique flags before and after array to find under/overrun
         * (int *) pre = (int) ptr ^ 12345 ;
         * (((int *) pre)+1) = (int) ptr ^ 13579 ;
         * (int *) post = (int) ptr ^ 67890 ;
         * (((int *) post)+1) = (int) ptr ^ 24680 ;
         }
      else {
         ptr8 = ptr ;
         pre = post = NULL ;   // FREE verification does not check for this /*!!!!!!*/
         if (mem_keep_log) {
            fp_rec = fopen ( mem_file_name , "at" ) ;
            fprintf ( fp_rec , " ERROR...negative address " ) ;
            fclose ( fp_rec ) ;
            }
         }
      actual[nallocs] = (int) ptr ;
      allocs[nallocs] = (int) ptr8 ;
      preptr[nallocs] = (int) pre ;
      postptr[nallocs] = (int) post ;
      alloc_size[nallocs++] = n ;
      total_use += n ;
      if (total_use > mem_max_used)
         mem_max_used = total_use ;
      }
   else {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "  Alloc = NULL... ERROR!" ) ;
         fclose ( fp_rec ) ;
         }
      return NULL ;
      }

   if (mem_keep_log) {
      fp_rec = fopen ( mem_file_name , "at" ) ;
#if DEBUG_PRE_POST
      fprintf ( fp_rec ,
         "\nAlloc=%u (%u) %u bytes  %d allocs total memory=%d (%d %d %d %d)" ,
         (unsigned) ptr8 , (unsigned) ptr , n, nallocs, total_use,
         * (int *) pre, * (((int *) pre)+1), * (int *) post,
         * (((int *) post)+1) ) ;
#else
      fprintf ( fp_rec ,
         "\nAlloc=%u  %u bytes  %d allocs  total memory=%d" ,
         (unsigned) ptr8 , n, nallocs, total_use ) ;
#endif
      fclose ( fp_rec ) ;
      }

   return ( ptr8 ) ;
}

void memfree ( void *ptr )
{
   int i ;
   void *ptr_to_free ;
   i = nallocs ;
   while (i--) {
      if (allocs[i] == (int) ptr)
         break ;
      }

   if (i < 0) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: illegal FREE = %u", (unsigned) ptr ) ;
         fclose ( fp_rec ) ;
         }
      exit ( 1 ) ;
      }

   if ((* (int *) preptr[i] != (actual[i] ^ 12345))
    || (* (((int *) preptr[i])+1) != (actual[i] ^ 13579))) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: FREE underrun = %u (wanted %d got %d)",
            (unsigned) ptr, actual[i] ^ 12345, * (int *) preptr[i] ) ;
         fclose ( fp_rec ) ;
         }
      exit ( 1 ) ;
      }

   if ((* (int *) postptr[i] != (actual[i] ^ 67890))
    || (* (((int *) postptr[i])+1) != (actual[i] ^ 24680))) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: FREE overrun = %u (wanted %d got %d)",
            (unsigned) ptr, actual[i] ^ 67890, * (int *) postptr[i] ) ;
         fclose ( fp_rec ) ;
         }
      exit ( 1 ) ;
      }

   --nallocs ;
   total_use -= alloc_size[i] ;
   ptr_to_free = (void *) actual[i] ;

   if (mem_keep_log) {
      fp_rec = fopen ( mem_file_name , "at" ) ;
      fprintf ( fp_rec , "\nFree=%u (%u) %d allocs  total memory=%d",
                (unsigned) ptr, (unsigned) ptr_to_free, nallocs, total_use );
      fclose ( fp_rec ) ;
      }

   while (i<nallocs) {   // Remove old ptr
      alloc_size[i] = alloc_size[i+1] ;
      allocs[i] = allocs[i+1] ;
      actual[i] = actual[i+1] ;
      preptr[i] = preptr[i+1] ;
      postptr[i] = postptr[i+1] ;
      ++i ;
      }

#if USE_MALLOC
   free ( ptr_to_free ) ;
#else
   GlobalFree ( ptr_to_free ) ;
#endif

}

void *memrealloc ( void *ptr , unsigned n )
{
   int i, old_offset, new_offset ;
   void *newptr, *ptr_to_realloc, *ptr8, *pre, *post ;
   union {
      void *dptr ;
      int cptr ;
      } uptr ;

   if (ptr == NULL)
      return memalloc ( n ) ;

   i = nallocs ;
   old_offset = 0 ;  // Not needed.  Shuts up LINT.
   while (i--) {
      if (allocs[i] == (int) ptr) {
         ptr_to_realloc = (void *) actual[i] ;
         old_offset = allocs[i] - actual[i] ;
         break ;
         }
      }

   if (i < 0) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , "\nMEM.CPP: Illegal REALLOC = %u", (unsigned) ptr );
         fclose ( fp_rec ) ;
         }
      return NULL ;
      }

   if ((* (int *) preptr[i] != (actual[i] ^ 12345))
    || (* (((int *) preptr[i])+1) != (actual[i] ^ 13579))) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf( fp_rec, "\nMEM.CPP: REALLOC underrun = %u (wanted %d got %d)",
            (unsigned) ptr, actual[i] ^ 12345, * (int *) preptr[i] ) ;
         fclose ( fp_rec ) ;
         }
      exit ( 1 ) ;
      }

   if ((* (int *) postptr[i] != (actual[i] ^ 67890))
    || (* (((int *) postptr[i])+1) != (actual[i] ^ 24680))) {
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf( fp_rec , "\nMEM.CPP: REALLOC overrun = %u (wanted %d got %d)",
            (unsigned) ptr, actual[i] ^ 67890, * (int *) postptr[i] ) ;
         fclose ( fp_rec ) ;
         }
      exit ( 1 ) ;
      }

#if USE_MALLOC
   newptr = (void *) realloc ( ptr_to_realloc , n + 3 * 8 ) ;
#else
   newptr = (void *) GlobalReAlloc ( ptr_to_realloc , n + 3 * 8 , GMEM_MOVEABLE ) ;
#endif

   if (((int) newptr) > 0) {
      uptr.dptr = newptr ;
      uptr.cptr = (uptr.cptr + 7) / 8 * 8 ;
      pre = uptr.dptr ;
      uptr.cptr = uptr.cptr + 8 ;
      ptr8 = uptr.dptr ;
      uptr.cptr = uptr.cptr + n ;
      post = uptr.dptr ;
      }
   else {
      ptr8 = newptr ;
      pre = post = NULL ;
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf ( fp_rec , " ERROR...negative address " ) ;
         fclose ( fp_rec ) ;
         }
      }

   if (mem_keep_log) {
      fp_rec = fopen ( mem_file_name , "at" ) ;
      fprintf ( fp_rec , "\nRealloc=%u (%u) %u bytes",
                (unsigned) ptr, (unsigned) newptr, n ) ;
      fprintf( fp_rec, " New=%u  total memory=%d",
         (unsigned) ptr8, total_use+n-alloc_size[i]) ;
      fclose ( fp_rec ) ;
      }

   if (newptr != NULL) {
      total_use = total_use + n - alloc_size[i] ;
      alloc_size[i] = n ;
      actual[i] = (int) newptr ;
      preptr[i] = (int) pre ;
      postptr[i] = (int) post ;
      allocs[i] = (int) ptr8 ;
      new_offset = allocs[i] - actual[i] ;
      if (new_offset != old_offset) {
         memmove ( (char *) newptr + new_offset ,        // = ptr8
                   (char *) newptr + old_offset , n ) ;
         if (mem_keep_log) {
            fp_rec = fopen ( mem_file_name , "at" ) ;
            fprintf ( fp_rec , " Unequal offset " ) ;
            fclose ( fp_rec ) ;
            }
         }
      if (((int) newptr) > 0) {
         // Place unique flags before and after array to find under/overrun
         * (int *) pre = (int) newptr ^ 12345 ;
         * (((int *) pre)+1) = (int) newptr ^ 13579 ;
         * (int *) post = (int) newptr ^ 67890 ;
         * (((int *) post)+1) = (int) newptr ^ 24680 ;
         }
      if (total_use > mem_max_used)
         mem_max_used = total_use ;
#if DEBUG_PRE_POST
      if (mem_keep_log) {
         fp_rec = fopen ( mem_file_name , "at" ) ;
         fprintf( fp_rec, " (%d %d %d %d)",
            * (int *) pre, * (((int *) pre)+1), * (int *) post,
            * (((int *) post)+1) ) ;
         fclose ( fp_rec ) ;
         }
#endif
      }

   return ptr8 ;
}

void memtext ( char *text )
{
   if (mem_keep_log) {
      fp_rec = fopen ( mem_file_name , "at" ) ;
      fprintf ( fp_rec , "\n%s", text ) ;
      fclose ( fp_rec ) ;
      }
}

void notext ( char * )
{
   return ;
}

void memclose ()
{
   int i ;

   if (mem_keep_log) {
      fp_rec = fopen ( mem_file_name , "at" ) ;
      fprintf( fp_rec , "\nMax memory use=%d  Dangling allocs=%d",
               mem_max_used , nallocs ) ;
      i = nallocs ;
      while (i--)
         fprintf ( fp_rec , "\n%d", (int) allocs[i] ) ;
      fclose (fp_rec ) ;
      }
}

void nomemclose ()
{
   return ;
}

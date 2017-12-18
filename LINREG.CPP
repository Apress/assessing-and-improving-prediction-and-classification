/******************************************************************************/
/*                                                                            */
/*  LINREG - Linear regression                                                */
/*                                                                            */
/*  To use this class:                                                        */
/*    1) Construct a new instance of the class                                */
/*    2) Call add_case() exactly ncases times, each time providing the        */
/*       npred=ncols vector of predictors.                                    */
/*    3) Call solve() as many times as desired with various right hand sides  */
/*    4) Optionally, call reset() and go to step 2                            */
/*                                                                            */
/*  This does not include any checks for insufficient memory.                 */
/*  It also assumes that the user calls add_case exactly ncases times         */
/*  and does not check for failure to do so.                                  */
/*                                                                            */
/******************************************************************************/

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "linreg.h"

LinReg::LinReg ( int ncase , int ncol )
{
   ncases = ncase ;
   ncols = ncol ;
   svd = new SingularValueDecomp ( ncase , ncol ) ;
   reset () ;
}


LinReg::~LinReg ()
{
   if (svd != NULL)
      delete svd ;
}

/*
   This discards any existing design matrix.
   It does not need to be called after construction, but it must
   be called if the user wants to reuse the LinReg object for a new
   design matrix.
*/

void LinReg::reset ()
{
   nrows = 0 ;      // No rows (via add_case()) yet present
   decomp = 0 ;     // Decomposition not done yet
}

/*
   Build the design matrix one case at a time.
   The user must call this member EXACTLY ncases times after construction
   or a call to reset(), and before a call to solve().
*/

void LinReg::add_case ( double *newcase )
{
   if (nrows >= ncases)  // Careful user never lets this happen
      return ;           // But cheap insurance

   memcpy ( svd->a + nrows * ncols , newcase , ncols * sizeof(double) ) ;
   ++nrows ;
}

/*
   After add_case has been called exactly ncases times, this may be called
   as many times as desired to solve a system.
*/

void LinReg::solve (
   double eps ,       // Singularity limit, typically 1.e-8 or so
   double *rhs ,      // Right hand side, ncases long
   double *b          // Output of solution, npred=ncols long
   )
{
   int i ;

   if (nrows != ncases) { // Careful user never lets this happen
      for (i=0 ; i<ncols ; i++)  // This might be a warning flag
         b[i] = 0.0 ;     // For a watchful but careless user
      return ;            // But cheap insurance
      }

   if (! decomp) {        // If the decomposition has not been done yet (1'st call)
      svd->svdcmp () ;    // Do it now
      decomp = 1 ;        // And flag that it has been done
      }

   memcpy ( svd->b , rhs , ncases * sizeof(double) ) ;
   svd->backsub ( eps , b ) ;
}

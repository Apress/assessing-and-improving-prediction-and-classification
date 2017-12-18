/******************************************************************************/
/*                                                                            */
/*  LOGISTIC - Logistic regression                                            */
/*                                                                            */
/*  Unlike most implementations (which use iterative gradient ascent)         */
/*  this uses simulated annealing.  This is considerably slower than          */
/*  the traditional method, but more likely to find the global optimum.       */
/*  It is also more numerically stable.                                       */
/*                                                                            */
/*  To use this class:                                                        */
/*    1) Construct a new instance of the class                                */
/*    2) Call add_case() exactly ncases times, each time providing the        */
/*       nin+1 vector of inputs and output.                                   */
/*    3) Call train()                                                         */
/*    4) Call predict() as many times as desired                              */
/*    5) Optionally, call reset() and go to step 2                            */
/*                                                                            */
/*  This does not include any checks for insufficient memory.                 */
/*  It also assumes that the user calls add_case exactly ncases times         */
/*  and does not check for failure to do so.                                  */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "logistic.h"
#include "minimize.h"

double normal () ;

static double max_exp = log ( 1.e190 ) ;
inline double safe_exp ( double x )
{
   if (x <= max_exp)
      return exp ( x ) ;

   return 1.e190 ;
}

/*
--------------------------------------------------------------------------------

   Constructor, destructor, reset(), add_case()

--------------------------------------------------------------------------------
*/

Logistic::Logistic ( int ncase , int nin )
{
   ncases = ncase ;
   ninputs = nin ;
   coefs = (double *) malloc ( (ninputs +1) * sizeof(double) ) ;
   tset = (double *) malloc ( ncases * (ninputs + 1) * sizeof(double) ) ;
   reset () ;
}


Logistic::~Logistic ()
{
   if (tset != NULL)
      free ( tset ) ;
   if (coefs != NULL)
      free ( coefs ) ;
}

/*
   This discards any existing training data.
   It does not need to be called after construction, but it must
   be called if the user wants to reuse the Logistic object for a new dataset.
*/

void Logistic::reset ()
{
   nrows = 0 ;      // No rows (via add_case()) yet present
   trained = 0 ;    // Training not done yet
}

/*
   Build the training set one case at a time.
   The user must call this member EXACTLY ncases times after construction
   or a call to reset(), and before a call to train().
*/

void Logistic::add_case ( double *newcase )
{
   if (nrows >= ncases)  // Careful user never lets this happen
      return ;           // But cheap insurance

#if 0
   printf ( "\n---> " ) ;
   for (int i=0 ; i<=ninputs ; i++)
      printf ( " %.2lf", newcase[i] ) ;
#endif
   memcpy ( tset + nrows * (ninputs + 1) , newcase ,
            (ninputs + 1) * sizeof(double) ) ;
   ++nrows ;
}


/*
--------------------------------------------------------------------------------

   predict() - Given an input vector, compute output using trained model
               The output is the linear combination, the log odds ratio

--------------------------------------------------------------------------------
*/

void Logistic::predict (
   double *input ,     // Input vector
   double *output      // Returned output
   )
{
   int i ;

   *output = coefs[ninputs] ;    // Constant term
   for (i=0 ; i<ninputs ; i++)
      *output += input[i] * coefs[i] ;
}

/*
--------------------------------------------------------------------------------

   execute() - Given coefficients, pass through the training set,
               return log likelihood.

--------------------------------------------------------------------------------
*/

double Logistic::execute ()
{
   int icase ;
   double *tptr, term, sum1, sum2 ;

   sum1 = sum2 = 0.0 ;

   for (icase=0 ; icase<ncases ; icase++) {
      tptr = tset + (ninputs + 1) * icase ; // This case
      predict ( tptr , &term ) ;            // Log odds ratio
      sum1 += term * tptr[ninputs] ;        // Output stored after inputs
      sum2 += log ( 1.0 + safe_exp ( term )) ;
      } // For all training cases

   return sum1 - sum2 ;
}

/*
--------------------------------------------------------------------------------

   train() - Train the model

   After add_case has been called exactly ncases times, this must be called
   to train the model.

--------------------------------------------------------------------------------
*/

double logit_crit ( double *x ) ;   // Local criterion function for optimization
double logit_unicrit ( double x ) ; // Local criterion function for optimization
static Logistic *local_logistic ;   // Needed by above

void Logistic::train ()
{
   int i, inner, outer, first ;
   double y, best_y, std, *test_wts, *best_wts, *center ;

   local_logistic = this ;

/*
   Best_wts keeps track of the (log) best coefs.
   Center is the center around which perturbation is done.
   It starts at zero.  After completion of each pass through the inner loop
   it is changed to best_wts.
*/

   test_wts = (double *) malloc ( ninputs * sizeof(double) ) ;
   best_wts = (double *) malloc ( ninputs * sizeof(double) ) ;
   center = (double *) malloc ( ninputs * sizeof(double) ) ;

   std = 1.0 ;   // Reasonable when predictors are mean ranks

   for (i=0 ; i<ninputs ; i++)
      center[i] = 0.0 ;

   first = 1 ;
   for (outer=0 ; outer<10 ; outer++) {
      for (inner=0 ; inner<10 + 5 * ninputs * ninputs ; inner++) {

         for (i=0 ; i<ninputs ; i++)
            test_wts[i] = center[i] + std * normal() ;

         y = logit_crit ( test_wts ) ;
         if (first  ||  (y > best_y)) {
            first = 0 ;
            best_y = y ;
            memcpy ( best_wts , test_wts , ninputs * sizeof(double) ) ;
            }

         } // For inner loop iterations
      memcpy ( center , best_wts , ninputs * sizeof(double) ) ;
      std *= 0.7 ;
      } // For outer loop iterations

   logit_crit ( best_wts ) ;  // Needed to set coefs correctly
   trained = 1 ;    // Training complete

#if 0
   printf ( "\n" ) ;
   for (i=0 ; i<=ninputs ; i++)
      printf ( "  %.3lf", coefs[i] ) ;
//   getch () ;
#endif


   free ( test_wts ) ;
   free ( best_wts ) ;
   free ( center ) ;
}

static double logit_crit ( double *x )
{
   int i ;
   double x1, y1, x2, y2, x3, y3 ;

   for (i=0 ; i<local_logistic->ninputs ; i++)
       local_logistic->coefs[i] = safe_exp ( x[i] ) ;

   glob_min ( -20.0 , 20.0 , 5 , 0 , -1.e160 , logit_unicrit ,
              &x1 , &y1 , &x2 , &y2 , &x3 , &y3 ) ;

   y2 = brentmin ( 50 , -1.e160 , 1.e-10 , 1.e-10 , logit_unicrit ,
                   &x1 , &x2 , &x3 , y2 ) ;

   local_logistic->coefs[local_logistic->ninputs] = x2 ;
   return -y2 ;
}

static double logit_unicrit ( double t )
{
   double penalty ;

   penalty = 0.0 ;
   if (fabs ( t ) > 20.0)                // Rare pathological event
      penalty = 1.e10 * (fabs ( t ) - 20.0) ;

   local_logistic->coefs[local_logistic->ninputs] = t ;
   return penalty - local_logistic->execute () ;
}

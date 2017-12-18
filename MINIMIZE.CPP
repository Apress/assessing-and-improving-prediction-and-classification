/******************************************************************************/
/*                                                                            */
/*  MINIMIZE - Minimization routines                                          */
/*                                                                            */
/******************************************************************************/

#include <math.h>

/*
--------------------------------------------------------------------------------

  GLOB_MIN - Check equispaced intervals to find rough global minimum
             of a univariate function

  This is called with a lower and upper bound for the domain to be
  searched.  If the function is still decreasing at one of these endpoints
  the search will continue beyond the specified endpoint.
  The total interval is equally divided into npts-1 subintervals.
  These subintervals will be spaced arithmetically or logarithmically
  according to log_space.
  If the criterion function drops to critlim or smaller, execution will
  terminate as soon as a local minimum is found.  (Global search stops.)
  In other words, this lets the function quit early, before the entire
  interval is searched, if the function value is satisfactorily small.
  This early exit means that the true global minimum might be missed,
  so in most cases critlim should be set impossibly small.
  Three points will be returned.  The center point, (x2,y2), will have
  smaller function value (y2) than its neighbors.  (In pathological
  cases they may be equal.)

  If npts is input negative, that means the user is inputting f(low) in *y2.
  That sometimes saves a function evaluation.

--------------------------------------------------------------------------------
*/

int glob_min (
   double low ,                // Lower limit for search
   double high ,               // Upper limit
   int npts ,                  // Number of points to try
   int log_space ,             // Space by log?
   double critlim ,            // Quit global if crit drops this low
   double (*criter) (double) , // Criterion function
   double *x1 ,
   double *y1 ,           // Lower X value and function there
   double *x2 ,
   double *y2 ,           // Middle (best)
   double *x3 ,
   double *y3             // And upper
   )
{
   int i, ibest, turned_up, know_first_point ;
   double x, y, rate, previous ;

   if (npts < 0) {
      npts = -npts ;
      know_first_point = 1 ;
      }
   else 
      know_first_point = 0 ;

   if (log_space)
      rate = exp ( log (high / low) / (npts - 1) ) ;
   else 
      rate = (high - low) / (npts - 1) ;

   x = low ;

   previous = 0.0 ; // Avoids "use before set" compiler warnings
   ibest = -1 ;     // For proper critlim escape
   turned_up = 0 ;  // Must know if function increased after min

   for (i=0 ; i<npts ; i++) {

      if (i  ||  ! know_first_point)
         y = criter ( x ) ;
      else
         y = *y2 ;

      if ((i == 0)  ||  (y < *y2)) {  // Keep track of best here
         ibest = i ;
         *x2 = x ;
         *y2 = y ;
         *y1 = previous ;  // Function value to its left
         turned_up = 0 ;   // Flag that min is not yet bounded
         }

      else if (i == (ibest+1)) { // Didn't improve so this point may
         *y3 = y ;               // be the right neighbor of the best
         turned_up = 1 ;         // Flag that min is bounded
         }

      previous = y ;             // Keep track for left neighbor of best

      if ((*y2 <= critlim)  &&  (ibest > 0)  &&  turned_up)
         break ; // Done if good enough and both neighbors found

      if (log_space)
         x *= rate ;
      else 
         x += rate ;
      }

/*
   At this point we have a minimum (within low,high) at (x2,y2).
   Compute x1 and x3, its neighbors.
   We already know y1 and y3 (unless the minimum is at an endpoint!).
*/

   if (log_space) {
      *x1 = *x2 / rate ;
      *x3 = *x2 * rate ;
      }
   else {
      *x1 = *x2 - rate ;
      *x3 = *x2 + rate ;
      }

/*
   Normally we would now be done.  However, the careless user may have
   given us a bad x range (low,high) for the global search.
   If the function was still decreasing at an endpoint, bail out the
   user by continuing the search.
*/

   if (! turned_up) { // Must extend to the right (larger x)
      for (;;) {      // Endless loop goes as long as necessary

         *y3 = criter ( *x3 ) ;

         if (*y3 > *y2)  // If function increased we are done
            break ;
         if ((*y1 == *y2)  &&  (*y2 == *y3)) // Give up if flat
            break ;

         *x1 = *x2 ;      // Shift all points
         *y1 = *y2 ;
         *x2 = *x3 ;
         *y2 = *y3 ;

         rate *= 2.0 ;    // Step further each time
         if (log_space)   // And advance to new frontier
            *x3 *= rate ;
         else 
            *x3 += rate ;
         }
      }

   else if (ibest == 0) {  // Must extend to the left (smaller x)
      for (;;) {           // Endless loop goes as long as necessary

         *y1 = criter ( *x1 ) ;

         if (*y1 > *y2)   // If function increased we are done
            break ;
         if ((*y1 == *y2)  &&  (*y2 == *y3)) // Give up if flat
            break ;

         *x3 = *x2 ;      // Shift all points
         *y3 = *y2 ;
         *x2 = *x1 ;
         *y2 = *y1 ;

         rate *= 2.0 ;    // Step further each time
         if (log_space)   // And advance to new frontier
            *x1 /= rate ;
         else 
            *x1 -= rate ;
         }
      }

   return 0 ;
}

/*
--------------------------------------------------------------------------------

  BRENTMIN - Use Brent's method to find a local minimum of a
             univariate function.

  This is given three points such that the center has lesser function
  value than its neighbors.  It iteratively refines the interval.
  If the criterion function drops to critlim or smaller, execution will
  terminate.  In most cases this should be set impossibly small so that
  eps and tol define convergence.

--------------------------------------------------------------------------------
*/

double brentmin (
   int itmax ,            // Iteration limit
   double critlim ,       // Quit if crit drops this low
   double eps ,           // Function convergence tolerance
   double tol ,           // X convergence tolerance
   double (*criter) (double) , // Criterion function
   double *xa ,           // Lower X value, input and output
   double *xb ,           // Middle (best), input and output
   double *xc ,           // And upper, input and output
   double y               // Function value at xb
   )
{
   int iter ;
   double x0, x1, x2, y0, y1, y2, xleft, xmid, xright, movement, trial ;
   double small_step, small_dist, numer, denom, temp1, temp2 ;
   double testdist, this_x, this_y ;

/*
   Initialize
*/

   x0 = x1 = x2 = *xb ;
   xleft = *xa ;
   xright = *xc ;

   y0 = y1 = y2 = y ;

/*
  We want a golden-section search the first iteration.  Force this by setting
  movement equal to zero.
*/

   movement = trial = 0.0 ;

/*
   Main loop.
*/

   for (iter=0 ; iter<itmax ; iter++) {

      if (y0 < critlim)    // Done?
         break ;

/*
   This test is more sophisticated than it looks.  It tests the closeness
   of xright and xleft (relative to small_dist), AND makes sure that x0 is
   near the midpoint of that interval.
*/

      small_step = fabs ( x0 ) ;
      if (small_step < 1.0)
         small_step = 1.0 ;
      small_step *= tol ;
      small_dist = 2.0 * small_step ;

      xmid = 0.5 * (xleft + xright) ;

      if (fabs ( x0 - xmid )  <=  (small_dist - 0.5 * (xright - xleft)))
         break ;

/*
   Avoid refining function to limits of precision
*/

      if ((iter >= 4)  &&  (((y2 - y0) / (y0 + 1.0)) < eps))
         break ;

      if (fabs ( movement )  >  small_step) {  // Try parabolic only if moving
         temp1 = (x0 - x2) * (y0 - y1) ;
         temp2 = (x0 - x1) * (y0 - y2) ;
         numer = (x0 - x1) * temp2 - (x0 - x2) * temp1 ;
         denom = 2. * (temp1 - temp2) ;
         testdist = movement ;     // Intervals must get smaller
         movement = trial ;
         if (fabs(denom) > 1.e-40)
            trial = numer / denom ; // Parabolic estimate of minimum
         else 
            trial = 1.e40 ;

         temp1 = trial + x0 ;
         if ((2.0 * fabs ( trial ) < fabs ( testdist ))// If shrinking
          && (temp1 > xleft) && (temp1 < xright)) {    // And safely in bounds
            this_x = temp1 ;                           // Use parabolic estimate
            if ((this_x - xleft  <  small_dist)  ||    // Cannot get too close
                (xright - this_x <  small_dist))       // to the endpoints
               trial = (x0 < xmid)  ?  small_step  :  -small_step ;
            }
         else {  // Punt via golden section because cannot use parabolic
            movement = (xmid > x0)  ?  xright - x0  :  xleft - x0 ;
            trial = .3819660 * movement ;
            }
         }
      else { // Must use golden section due to insufficient movement
         movement = (xmid > x0)  ?  xright - x0  :  xleft - x0 ;
         trial = .3819660 * movement ;
         }

      if (fabs (trial)  >=  small_step)     // Make sure we move a good distance
         this_x = x0 + trial ;
      else
         this_x = (trial > 0.0)  ?  x0 + small_step  :  x0 - small_step ;

/*
   Evaluate the function here.
   Insert this new point in the correct position in the 'best' hierarchy
*/

      this_y = criter ( this_x ) ;

      if (this_y <= y0) {    // Improvement
         if (this_x < x0)
            xright = x0 ;
         else
            xleft = x0 ;
         x2 = x1 ;
         x1 = x0 ;
         x0 = this_x ;
         y2 = y1 ;
         y1 = y0 ;
         y0 = this_y ;
         }

      else {                  // No improvement
         if (this_x >= x0)
            xright = this_x ;
         else
            xleft = this_x ;

         if ((this_y <= y1)  ||  (x1 == x0)) {
            x2 = x1 ;
            x1 = this_x ;
            y2 = y1 ;
            y1 = this_y ;
            }
         else if ((this_y <= y2)  ||  (x2 == x0)  ||  (x2 == x1)) {
            x2 = this_x ;
            y2 = this_y ;
            }
         }
      }
 
   *xa = xleft ;
   *xb = x0 ;
   *xc = xright ;

   return y0 ;
}

/*
--------------------------------------------------------------------------------

  POWELL - Use Powell's method to find a local minimum of a function

--------------------------------------------------------------------------------
*/

/*
   This routine uses the general univariate minimizers 'glob_min' and
   'brentmin' to minimize along the gradient line. So we must have a local
   function for them to call, and it must have access to the relevant data.
   These statics handle that.
*/

static double univar_crit ( double t ) ; // Local univariate criterion
static double *local_x, *local_base, *local_direc ;   // It uses these
static int local_n ;
static double (*local_criter) ( double * ) ;

double powell (
   int maxits ,           // Iteration limit
   double critlim ,       // Quit if crit drops this low (Normally set impossibly small)
   double tol ,           // Convergence tolerance
   double (*criter) ( double * ) , // Criterion func
   int n ,                // Number of variables
   double *x ,            // In/out of independent variable
   double ystart ,        // Input of starting function value
   double *base ,         // Work vector n long
   double *p0 ,           // Work vector n long
   double *direc          // Work vector n*n long
   )
{
   int i, j, idir, iter, convergence_counter, idelta, replaced ;
   double fval, fbest, f0, test, t1, t2, t3, y1, y2, y3 ;
   double prev_best, toler, delta, scale, len, ftest, mult ;

/*
   Initialize for the local univariate criterion which may be called by
   'glob_min' and 'brentmin' to minimize along the search direction.
*/


   local_x = x ;
   local_base = base ;
   local_n = n ;
   local_criter = criter ;

/*
   Initialize the direction matrix to be an identity.
*/

   for (i=0 ; i<n ; i++) {
      for (j=0 ; j<n ; j++)
         direc[i*n+j] = (j == i) ? 1.0 : 0.0 ;
      }

/*
   Main loop.  For safety we impose a limit on iterations.
*/

   replaced = -1 ;          // Assume no replacement will be done
   prev_best = 1.e60 ;
   fbest = ystart ;
   scale = 0.2 ;
   iter = convergence_counter = 0 ;

   for (;;) {

      if ((iter++ >= maxits)  &&  (maxits > 0))
         break ;

      if (fbest < critlim)     // Do we satisfy user yet?
         break ;

/*
   Convergence check
*/

      if (fabs(prev_best) <= 1.0)            // If the function is small
         toler = tol ;                       // Work on absolutes
      else                                   // But if it is large
         toler = tol * fabs(prev_best) ;     // Keep things relative

      if ((prev_best - fbest)  <=  toler) {  // If little improvement
         if (++convergence_counter >= 2)     // Then count how many
            break ;                          // And quit if too many
         }
      else                                   // But a good iteration
         convergence_counter = 0 ;           // Resets this counter

      prev_best = fbest ;

/*
   Loop over all search directions, minimizing in each.
   Keep track of the direction that gave the most improvement.
   For efficiency, we keep 'replaced' as the vector that just got replaced
   (or -1 if no replacement was done).  Skip the first direction if that
   was the one just replaced!
*/

      for (i=0 ; i<n ; i++)            // Preserve the starting point
         p0[i] = x[i] ;
      f0 = fbest ;                     // And function value at that point
      delta = -1.0 ;                   // Will keep best improvement here
      idelta = 0 ;                     // Shuts up LINT

      for (idir=0 ; idir<n ; idir++) { // For each direction
         if ((n>1)  &&  ! idir  &&  ! replaced) // If we just replaced the
            continue ;                 // first vector, avoid waste
         for (i=0 ; i<n ; i++)         // Local criter steps out from here
            base[i] = x[i] ;           // So it must be current point
         y2 = fbest ;                  // Glob_min can use first f value
         local_direc = direc + idir * n ; // This is the idir direction
         for (mult=0.1 ; mult < 11.0 ; mult *= 4.0) {
            glob_min ( -mult * scale , mult * scale , 15 , 0 ,
                       critlim , univar_crit , &t1 , &y1 , &t2 ,
                       &y2 , &t3 , &y3 ) ;
            if ((y2 < y1)  &&  (y2 < y3)) // Loop until minimum is bounded
               break ;
            }
         if (y2 < critlim) {           // Good enough already?
            if (y2 < fbest) {          // If global caused improvement
               for (i=0 ; i<n ; i++)   // Implement that improvement
                  x[i] = base[i] + t2 * local_direc[i] ;
               fbest = y2 ;
               }
            else {                     // Else revert to starting point
               for (i=0 ; i<n ; i++)
                  x[i] = base[i] ;
               }
            goto FINISH ;
            }

         if (convergence_counter)  // If failing, try extra hard
            fval = brentmin ( 40 , critlim , tol , 1.e-7 ,
                              univar_crit , &t1 , &t2 , &t3 , y2 ) ;
         else                      // But normally refine only moderately
            fval = brentmin ( 20 , critlim , 10.0 * tol , 1.e-5 ,
                              univar_crit , &t1 , &t2 , &t3 , y2 ) ;
         scale = fabs(t2) / n  +  (1.0 - 1.0/n) * scale ; // Keep reasonable

         for (i=0 ; i<n ; i++)          // Get current point from parametric
            x[i] = base[i] + t2 * local_direc[i] ;
         if (fbest - fval > delta) {    // Keep track of best direction
            delta = fbest - fval ;
            idelta = idir ;
            }
         fbest = fval ;                  // This is always the best so far
         } // For all directions

/*
   Before looping through all n directions, we stood at point p0 with f=f0.
   We now stand at point x with f=fbest.
   It is quite possible that the average direction of motion points right
   along a ravine.  Thus, it behooves us to step out in that direction.
   Try it.  We might luck out.
*/

      for (i=0 ; i<n ; i++) {
         p0[i] = x[i] - p0[i] ;    // Preserve average direction here
         base[i] = x[i] + p0[i] ;  // Step out to this point (borrow base)
         }
      fval = criter ( base ) ;     // Evaluate function at this test point

/*
   If this step improved, and if a more sophisticated second derivative
   test is passed, we replace the old direction of max improvement with this
   average direction (after first minimizing in this great direction).
   Also, set 'replaced' to the vector that was replaced if that is done.
   Otherwise set it to -1.  We do that so on the next iteration we skip
   minimization along the first direction if that happened to be the one
   that was just replaced!
*/

      replaced = -1 ;                  // Assume no replacement will be done
      ftest = fbest ;                  // Save for replacement test
      if (fval < fbest) {              // If this trial step improved
         fbest = fval ;                // Might as well keep this better point
         for (i=0 ; i<n ; i++)         // Keep x current
            x[i] = base[i] ;
         }

      if (fval < f0) {                 // First of two tests for replacement
         test = f0 - ftest - delta ;
         test = 2.0 * (f0 - 2.0 * ftest + fval) * test * test ;
         if (test < delta * (f0 - fval) * (f0 - fval)) { // Use this direction
            replaced = idelta ;        // Record this upcoming replacement
            len = 0.0 ;                // Will cumulate length of new dir
            for (i=0 ; i<n ; i++)
               len += p0[i] * p0[i] ;
            len = sqrt ( len ) ;
            for (i=0 ; i<n ; i++)
               p0[i] /= len ;          // Keep direction unit length
            local_direc = p0 ;         // We put the average direction here 
            y2 = fbest ;               // Glob_min can use first f value
            for (i=0 ; i<n ; i++)      // Set out from here
               base[i] = x[i] ;
            for (mult=0.1 ; mult < 11.0 ; mult *= 4.0) {
               glob_min ( -mult * scale , mult * scale , 15 , 0 ,
                          critlim , univar_crit , &t1 , &y1 , &t2 ,
                          &y2 , &t3 , &y3 ) ;
               if ((y2 < y1)  &&  (y2 < y3)) // Loop until minimum is bounded
                  break ;
               }
            if (y2 < critlim) {          // Good enough already?
               if (y2 < fbest) {         // If global caused improvement
                  for (i=0 ; i<n ; i++)  // Implement that improvement
                     x[i] = base[i] + t2 * local_direc[i] ;
                  fbest = y2 ;
                  }
               else {                            // Else revert to starting point
                  for (i=0 ; i<n ; i++)
                     x[i] = base[i] ;
                  }
               break ;
               }
            if (convergence_counter)  // If failing, try extra hard
               fval = brentmin ( 40 , critlim , tol , 1.e-7 ,
                                 univar_crit , &t1 , &t2 , &t3 , y2 ) ;
            else                      // But normally refine only moderately
               fval = brentmin ( 20 , critlim , 10.0 * tol , 1.e-5 ,
                                 univar_crit , &t1 , &t2 , &t3 , y2 ) ;
            scale = fabs(t2) / n  +  (1.0 - 1.0/n) * scale ; // Scale reasonable
            for (i=0 ; i<n ; i++)          // Get current point from parametric
               x[i] = base[i] + t2 * local_direc[i] ;
            fbest = fval ;                 // This is always the best so far

            for (i=0 ; i<n ; i++)
               direc[idelta*n+i] = p0[i] ; // Replace best direction with this
            } // If second derivative test passed
         } // If test one passed
      } // Main loop

FINISH:
      return fbest ;
}


/*
--------------------------------------------------------------------------------

   local criterion for univariate minimization

--------------------------------------------------------------------------------
*/

static double univar_crit ( double t )
{
   int i ;

   for (i=0 ; i<local_n ; i++)
      local_x[i] = local_base[i] + t * local_direc[i] ;
   return local_criter ( local_x ) ;
}

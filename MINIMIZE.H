extern int glob_min ( double low , double high , int npts , int log_space ,
   double critlim , double (*criter) (double) ,
   double *x1, double *y1 , double *x2, double *y2 , double *x3, double *y3 ) ;

extern double brentmin ( int itmax , double critlim , double eps ,
                         double tol , double (*criter) (double) ,
                         double *x1 , double *x2 , double *x3 , double y ) ;

extern double powell ( int maxits , double critlim , double tol ,
   double (*criter) ( double * ) , int n , double *x , double ystart ,
   double *base , double *p0 , double *direc ) ;

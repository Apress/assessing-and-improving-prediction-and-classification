// Class headers, function declarations and constants for information code

/*
   These are for intercepting memory allocation for runtime checking
*/

#if 1
#define MALLOC memalloc
#define FREE memfree
#define REALLOC memrealloc
#define MEMTEXT memtext
#define MEMCLOSE memclose
#else
#define MALLOC malloc
#define FREE free
#define REALLOC realloc
#define MEMTEXT notext
#define MEMCLOSE nomemclose
#endif

#if ! defined ( PI )
#define PI 3.141592653589793
#endif


/*
--------------------------------------------------------------------------------

   Interpolation routines

--------------------------------------------------------------------------------
*/

class CubicSpline {

public:
   CubicSpline ( int n , double *xin , double *yin ) ;
   ~CubicSpline () ;
   double evaluate ( double x ) ;

private:
   int n ;
   double *x ;
   double *y ;
   double *y2 ;
} ;

class Bilinear {

public:
   Bilinear ( int nxin , double *xin , int nyin , double *yin , double *zin ,
              int extra ) ;
   ~Bilinear () ;
   double evaluate ( double x , double y ) ;

private:
   int quadratic ;
   int nx ;
   int ny ;
   double *x ;
   double *y ;
   double *z ;
} ;

/*
--------------------------------------------------------------------------------

   ParzenDensities ParzDens_? used for continuous mutual information

--------------------------------------------------------------------------------
*/

class ParzDens_1 {

public:
   ParzDens_1 ( int n_tset , double *tset , int n_div ) ;
   ~ParzDens_1 () ;
   double density ( double x ) ;
   double low ;     // Lowest value with significant density
   double high ;    // And highest

private:
   int nd ;         // Number of points in array below
   double *d ;      // The data on which the density is based
   double var ;     // Presumed variance
   double factor ;  // Normalizing factor to make it a density
   CubicSpline *spline ; // Used only if interpolation
} ;

class ParzDens_2 {

public:
   ParzDens_2 ( int n_tset , double *tset0 , double *tset1 , int n_div ) ;
   ~ParzDens_2 () ;
   double density ( double x0 , double x1 ) ;

private:
   int nd ;         // Number of points in arrays below
   double *d0 ;     // The data on which the density is based; first variable
   double *d1 ;     // And second
   double var0 ;    // Presumed variance of first variable
   double var1 ;    // And second
   double factor ;  // Normalizing factor to make it a density
   Bilinear *bilin ; // Used only for bilinear interpolation
} ;

class ParzDens_3 {

public:
   ParzDens_3 ( int n_tset , double *tset0 , double *tset1 , double *tset2 , int n_div ) ;
   ~ParzDens_3 () ;
   double density ( double x0 , double x1 , double x2 ) ;

private:
   int nd ;         // Number of points in arrays below
   double *d0 ;     // The data on which the density is based; first variable
   double *d1 ;     // And second
   double *d2 ;     // And third
   double var0 ;    // Presumed variance of first variable
   double var1 ;    // And second
   double var2 ;    // And third
   double factor ;  // Normalizing factor to make it a density
} ;

/*
--------------------------------------------------------------------------------

   MutualInformation

--------------------------------------------------------------------------------
*/

class MutualInformationParzen {  // Parzen window method

public:
   MutualInformationParzen ( int nn , double *dep_vals , int ndiv ) ;
   ~MutualInformationParzen () ;
   double mut_inf ( double *x ) ;

private:
   int n ;             // Number of cases
   int n_div ;         // Number of divisions of range, typically 5-10
   double *depvals ;   // 'Dependent' variable
   ParzDens_1 *dens_dep ;   // Marginal density of 'dependent' variable
} ;

class MutualInformationAdaptive {  // Adaptive partitioning method

public:
   MutualInformationAdaptive ( int nn , double *dep_vals ,
                               int respect_ties , double crit ) ;
   ~MutualInformationAdaptive () ;
   double mut_inf ( double *x , int respect_ties ) ;

private:
   int n ;             // Number of cases
   int *y ;            // 'Dependent' variable ranks
   int *y_tied ;       // tied[i] != 0 if case with rank i == case with rank i+1
   double chi_crit ;   // Chi-square test criterion
} ;

class MutualInformationDiscrete {

public:
   MutualInformationDiscrete ( int nc , short int *bins ) ;
   ~MutualInformationDiscrete () ;
   double entropy () ;
   double mut_inf ( short int *bins ) ;
   double conditional ( short int *bins ) ;
   double conditional_error ( short int *bins ) ;
   double HYe ( short int *bins ) ;
   double hPe ( short int *bins ) ;

private:
   int ncases ;         // Number of cases
   short int *bins_y ;  // They are here
   int nbins_y ;        // Number of bins
   int *marginal_y ;    // Marginal distribution
} ;


/*
--------------------------------------------------------------------------------

   Function declarations

--------------------------------------------------------------------------------
*/

extern void free_data ( int nvars , char **names , double *data ) ;
extern double trans_ent ( int n , int nbins_x , int nbins_y , short int *x , short int *y ,
                          int xlag , int xhist , int yhist , int *counts , double *ab ,
                          double *bc , double *b ) ;
extern double integrate ( double low , double high , double min_width ,
                          double acc , double tol , double (*criter) (double) );
extern double inverse_normal_cdf ( double p ) ;
extern void *memalloc ( unsigned int n ) ;
extern void nomemclose () ;
extern void memclose () ;
extern void memfree ( void *ptr ) ;
extern void *memrealloc ( void *ptr , unsigned int size ) ;
extern void notext ( char *text ) ;
extern void memtext ( char *text ) ;
extern double mutinf_b ( int n , short int *y , short int *x , short int *z ) ;
extern double normal () ;
extern void partition ( int n , double *data , int *npart ,
                        double *bnds , short int *bins ) ;
extern void qsortd ( int first , int last , double *data ) ;
extern void qsortds ( int first , int last , double *data , double *slave ) ;
extern void qsortdsi ( int first , int last , double *data , int *slave ) ;
extern unsigned int RAND32 () ;
extern int readfile ( char *name , int *nvars , char ***names ,
                      int *ncases , double **data ) ;
extern double unifrand () ;

class Logistic {

public:

   Logistic ( int ncase , int nin ) ;
   ~Logistic () ;
   void reset () ;
   void add_case ( double *newcase ) ;
   void train () ;
   void predict ( double *input , double *output ) ;

   double execute () ;

   int ncases ;     // Number of cases
   int ninputs  ;   // Number of inputs
   int nrows ;      // How many times has add_case() been called?
   int trained ;    // Has it been trained yet?
   double *tset ;   // Ncases by (ninputs+1) matrix of training data
   double *coefs ;  // Trained coefficient vector ninputs+1 long
} ;


class GRNN {

public:

   GRNN ( int ncase , int nin , int nout ) ;
   ~GRNN () ;
   void reset () ;
   void add_case ( double *newcase ) ;
   void train () ;
   void anneal_train ( int n_outer , int n_inner , double start_std ) ;
   void predict ( double *input , double *output ) ;


private:
   double execute () ;

   int ncases ;     // Number of cases
   int ninputs  ;   // Number of inputs
   int noutputs  ;  // Number of outputs
   int nrows ;      // How many times has add_case() been called?
   int trained ;    // Has it been trained yet?
   double *tset ;   // Ncases by (ninputs+noutputs) matrix of training data
   double *sigma ;  // Ninputs vector of sigma weights
   double *outwork ;// Noutputs work vector
} ;

class SingularValueDecomp {

public:

   SingularValueDecomp ( int nrows , int ncols , int save_a=0 ) ;
   ~SingularValueDecomp () ;
   void svdcmp () ;
   void backsub ( double limit , double *soln ) ;

   int ok ;         // Was everything legal and allocs successful?

/*
   These are made public to allow access if desired.
   Normally, only 'a' (the design matrix) and 'b' (the right-hand-side)
   are written by the user.  If 'save_a' is nonzero, 'a' is kept intact.
*/

   double *a ;      // nrows by ncols input of design, output of U
   double *u ;      // unless save_a nonzero, in which case U output in 'u'
   double *w ;      // Unsorted ncols vector of singular values
   double *v ;      // Ncols by ncols output of 'v'
   double *b ;      // Nrows right-hand-side for backsub


private:

   void bidiag ( double *matrix ) ;
   double bid1 ( int col , double *matrix , double scale ) ;
   double bid2 ( int col , double *matrix , double scale ) ;
   void right ( double *matrix ) ;
   void left ( double *matrix ) ;
   void cancel ( int low , int high , double *matrix ) ;
   void qr ( int low , int high , double *matrix ) ;
   void qr_mrot ( int col , double sine , double cosine , double *matrix ) ;
   void qr_vrot ( int col , double sine , double cosine ) ;

   int rows ;       // Nrows preserved here
   int cols ;       // And ncols
   double *work ;   // Scratch vector ncols long
   double norm ;    // Norm of 'a' matrix
} ;

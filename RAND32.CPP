/******************************************************************************/
/*                                                                            */
/*  RAND32 - Assorted high quality random numbers                             */
/*                                                                            */
/*  These routines are nonportable because they assume 32-bit ints.           */
/*  The following routines are here:                                          */
/*                                                                            */
/*    RAND32M - A very good 32-bit generator by Marsaglia that is fast and    */
/*      produces high quality numbers.  But it requires a compiler that can   */
/*      handle 64-bit ints.                                                   */
/*                                                                            */
/*    RAND_LECUYER - A very good generator by L'Ecuyer that returns values    */
/*      in the range 1 - 2147483562, which is almost (but not quite) 2^31.    */
/*                                                                            */
/*    RAND_KNUTH - An old but quite good generator by Knuth that returns      */
/*      values in the range 0 - 999,999,999, which is almost 2^30.            */
/*                                                                            */
/*    RAND16_LECUYER - A true 16-bit (0-65535) generator derived from above.  */
/*                                                                            */
/*    RAND16_KNUTH - A true 16-bit (0-65535) generator derived from above.    */
/*                                                                            */
/*    RAND32 -  A 32-bit generator that I believe is excellent.               */
/*      It puts together 32-bit ints from pairs of 16-bit ints, a practice    */
/*      that is admittedly suspect unless the 16-bit components are excellent.*/
/*      I believe this to be the case here.  The 16-bit randoms are created   */
/*      by adding the outputs of RAND16_LECUYER and RAND16_KNUTH mod 2^16,    */
/*      using a 64K (65536) element Bays-Duram shuffle on KNUTH.              */
/*                                                                            */
/*                                                                            */
/*   Summary:                                                                 */
/*                                                                            */
/*   If your compiler can handle 64-bit ints, and you want a very good and    */
/*   very fast generator, I believe that it is tough to beat RAND32M for      */
/*   the combination of speed and quality.  You can find generators that      */
/*   are faster, and some that are better.  But RAND32M does well in both.    */
/*                                                                            */
/*   If your goal is to use the outputs to generate floats, RAND_LECUYER      */
/*   will do well.  Subtract 1 so that zero is included, and divide by        */
/*   2147483562 so that 1.0 is excluded, which is the usual practice.         */
/*                                                                            */
/*   RAND_KNUTH will also serve this purpose, although it is probably not     */
/*   quite as good, and it has coarser granularity.  But the Knuth algorithm  */
/*   is very different from the L'Ecuyer algorithm, which means that if       */
/*   one fails for your application, the other may well succeed.              */
/*                                                                            */
/*   If your goal is simple selection, as in genetic algorithms, either       */
/*   RAND16_LECUYER or RAND16_KNUTH will do well.                             */
/*                                                                            */
/*   Finally, if your application demands extreme quality, and speed is not   */
/*   critical, I can say that I use RAND32.  It is difficult to prove quality */
/*   of any generator, but I have tested RAND32 for hundreds of hours using   */
/*   routines from Marsaglia's DIEHARD suite, the NIST cryptographic suite,   */
/*   and some of my own tests.  RAND32 has passed every test I used.          */
/*   But of course, I make no guarantees.  It is almost certainly excellent,  */
/*   and I have not been able to find any test that it fails.  Still, this    */
/*   does not mean that it will perform well with every application.          */
/*                                                                            */
/******************************************************************************/

#include <math.h>
#include <float.h>

/*
--------------------------------------------------------------------------------

   This is a generator suggested by Marsaglia in his DIEHARD suite.
   It provides a great combination of speed and quality, but it requires
   a compiler that can handle 64-bit ints.

--------------------------------------------------------------------------------
*/

#if 0

static unsigned int Q[256], carry=362436 ;
static int MWC256_initialized = 0 ;
static int MWC256_seed = 123456789 ;

void RAND32M_seed ( int iseed ) { // Optionally set seed
   MWC256_seed = iseed ;
   MWC256_initialized = 0 ;
   }

unsigned int RAND32M ()
{
   unsigned _int64 t ;
   unsigned _int64 a=809430660 ;
   static unsigned char i=255 ;

   if (! MWC256_initialized) {
      unsigned int k,j=MWC256_seed ;
      MWC256_initialized = 1 ;
      for (k=0 ; k<256 ; k++) {
         j = 69069 * j + 12345 ; // This overflows, doing an automatic mod 2^32
         Q[k] = j ;
         }
      }

   t = a * Q[++i] + carry ;  // This is the 64-bit op, forced by a being 64-bit
   carry = (unsigned int) (t >> 32) ;
   Q[i] = (unsigned int) (t & 0xFFFFFFFF) ;
   return Q[i] ;
}
#endif

/*
--------------------------------------------------------------------------------

   RAND_LECUYER

--------------------------------------------------------------------------------
*/

static int LECUYER_initialized = 0 ;
static int LECUYER_seed1 = 1 ;

void RAND_LECUYER_seed ( int iseed )  // Optionally set seed
{
   LECUYER_seed1 = iseed ;
   LECUYER_initialized = 0 ;
}

unsigned int RAND_LECUYER ()
{
   int i, k, index ;
   static int LECUYER_seed2, LECUYER_output, LECUYER_table[32] ;

   if (! LECUYER_initialized) {     // Initialize the shuffle table
      LECUYER_initialized = 1 ;
      LECUYER_seed2 = LECUYER_seed1 ;       // Will use this in 'normal' call
      for (i=0 ; i<32 ; i++) {
         k = LECUYER_seed1 / 53668 ;
         LECUYER_seed1 = 40014 * (LECUYER_seed1 - k * 53668) - k * 12211 ;
         if (LECUYER_seed1 < 0)
            LECUYER_seed1 += 2147483563 ;
         LECUYER_table[i] = LECUYER_seed1 ;
         }
      LECUYER_output = LECUYER_table[0] ;  // For first use below
      }

   // Generate the second of the two randoms that will be combined
   k = LECUYER_seed2 / 52774 ;
   LECUYER_seed2 = 40692 * (LECUYER_seed2 - k * 52774) - k * 3791 ;
   if (LECUYER_seed2 < 0)
      LECUYER_seed2 += 2147483399 ;

   index = LECUYER_output / 67108862 ;      // Compute the shuffle index and combine
   LECUYER_output = LECUYER_table[index] - LECUYER_seed2 ;  // the two randoms
   if (LECUYER_output < 0)
      LECUYER_output += 2147483563 ;

   // Generate a new 'first random' and replace the just used shuffle table entry
   k = LECUYER_seed1 / 53668 ;
   LECUYER_seed1 = 40014 * (LECUYER_seed1 - k * 53668) - k * 12211 ;
   if (LECUYER_seed1 < 0)
      LECUYER_seed1 += 2147483563 ;
   LECUYER_table[index] = LECUYER_seed1 ;

   return LECUYER_output ;
}

/*
--------------------------------------------------------------------------------

   RAND_KNUTH

--------------------------------------------------------------------------------
*/

static int KNUTH_initialized = 0 ;
static int KNUTH_seed1 = 1 ;

void RAND_KNUTH_seed ( int iseed )  // Optionally set seed
{
   KNUTH_seed1 = iseed ;
   KNUTH_initialized = 0 ;
}

unsigned int RAND_KNUTH ()
{
   int i, k, index ;
   static int KNUTH_output, KNUTH_table[55], next, nextp ;

   if (! KNUTH_initialized) {     // Initialize the shuffle table
      KNUTH_initialized = 1 ;
      KNUTH_output = KNUTH_seed1 % 1000000000 ; // Ensure seed not illegal
      KNUTH_table[54] = KNUTH_output ;
      k = 1 ;
      for (i=1 ; i<=54 ; i++) {
         index = (21 * i) % 55 - 1 ;
         KNUTH_table[index] = k ;
         k = KNUTH_output - k ;
         if (k < 0)
            k += 1000000000 ;
         KNUTH_output = KNUTH_table[index] ;
         }
      for (k=0 ; k<4 ; k++) {
         for (i=0 ; i<55 ; i++) {
            KNUTH_table[i] -= KNUTH_table[(i+31)%55] ;
            if (KNUTH_table[i] < 0)
               KNUTH_table[i] += 1000000000 ;
            }
         }

      next = 0 ;
      nextp = 31 ;
      } // If not initialized

   KNUTH_output = KNUTH_table[next] - KNUTH_table[nextp] ;
   if (KNUTH_output < 0)
      KNUTH_output += 1000000000 ;
   KNUTH_table[next] = KNUTH_output ;

   next = (next + 1) % 55 ;
   nextp = (nextp + 1) % 55 ;

   return KNUTH_output ;
}

/*
--------------------------------------------------------------------------------

   These are 16-bit subgenerators

--------------------------------------------------------------------------------
*/

unsigned int RAND16_LECUYER ()
{
   long k ;
   long mult = 2147483562 / 65536 ;
   long max = mult * 65536L ;

   for (;;) {
      k = RAND_LECUYER() - 1 ;
      if (k < max )
         return k / mult ;
      }
}

unsigned int RAND16_KNUTH ()
{
   long k ;
   long mult = 1000000000 / 65536 ;
   long max = mult * 65536 ;

   for (;;) {
      k = RAND_KNUTH() ;
      if (k < max )
         return k / mult ;
      }
}

/*
--------------------------------------------------------------------------------

   RAND32

--------------------------------------------------------------------------------
*/

static int RAND32_table[65536] ;    // Keep shuffle table here
static int RAND32_initialized = 0 ; // Has it been initialized?

/*
   Set the random seed
*/

void RAND32_seed ( unsigned int iseed )
{
   RAND_KNUTH_seed ( iseed & 65535 ) ;
   RAND_LECUYER_seed ( (iseed >> 16) & 65535 ) ;
   RAND32_initialized = 0 ;
}

unsigned int RAND32 ()
{
   int i, k1, k2 ;
   static int RAND32_randout ;

   if (! RAND32_initialized) {  // Initialize shuffle table before use
      RAND32_initialized = 1 ;  // Flag to avoid more inits
      for (i=0 ; i<65536 ; i++)       // Fill entire table
         RAND32_table[i] = RAND16_KNUTH() ;
      RAND32_randout = RAND16_KNUTH() ; // One more for first use
      }

   k1 = RAND32_randout = (RAND32_table[RAND32_randout] + RAND16_LECUYER()) % 65536 ;
   RAND32_table[RAND32_randout] = RAND16_KNUTH () ;

   k2 = RAND32_randout = (RAND32_table[RAND32_randout] + RAND16_LECUYER()) % 65536 ;
   RAND32_table[RAND32_randout] = RAND16_KNUTH () ;

   return (k1 << 16)  |  k2 ;
}

/*
--------------------------------------------------------------------------------

   Generate a uniform in [0, 1).

--------------------------------------------------------------------------------
*/

double unifrand ()
{
   double r1, r2 ;
   double denom = 0x7FFFFFFFL + 1.0 ;

   r1 = RAND32 () & 0x7FFFFFFFL ;
   r2 = RAND32 () & 0x7FFFFFFFL ;
   return (r1 + r2 / denom) / denom ;
}

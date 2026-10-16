#ifndef MISC_MACROS_H
#define MISC_MACROS_H

#include "random.h"

extern Random  rnd;

// ----------------------------------------------------------------
// random-number generation macros

/* Random number generator - double from [0,1]
	 rand() is the only random number generator in C++
	 Note that rand() is from 0 to RAND_MAX inclusive!
	 */
// #define frand (((double) rand() + 0.5) / ((double) RAND_MAX + 1.0))
#define frand (rnd.getRandomDoubl01())

/* Random number generator - integer from 0 to Max-1
	 rand() is the only random number generator in VC++
	 Note that rand() is from 0 to RAND_MAX inclusive!
	 */
// #define irand(Max) (rand() % (Max))
#define irand(Max) (rnd.getRandomInt(Max))

/* Random numbers from standard normal distribution
	 using Box-Muller formula
	 */
#define Pi 3.1415926535
#define BM_normal() (double)(sqrt(-2*log(frand))*cos(2*Pi*frand))


// ----------------------------------------------------------------
// gen. mathematical macros

/* These definitions need GNU C extensions!
   With the standard definitions (see below),
   min(x,f(x)) will cause f(x) to be computed twice.
   These definitions avoid the problem...
   */
#define min(X, Y)                     \
({ typeof (X) __x = (X), __y = (Y);   \
  (__x < __y) ? __x : __y; })
#define max(X, Y)                     \
({ typeof (X) __x = (X), __y = (Y);   \
  (__x > __y) ? __x : __y; })
/* Standard definitions:
#define min(X, Y)  ((X) < (Y) ? (X) : (Y))
#define max(X, Y)  ((X) > (Y) ? (X) : (Y))
*/

// use this for comparing two float numbers
#define EPS 1e-8


// ----------------------------------------------------------------
// generic C/C++ macros

#define allocate(var, type, size) if ((var = (type *) malloc(size * sizeof(type))) == NULL) { printf("\nNOT ENOUGH MEMORY!\n\n");	exit(1); }

#endif

#ifdef MSDOS		/* 80x87 math chip microsoft C 	*/
#define EPSILON 1.0e-10
#define SQREPSILON 1.0e-20
#else			/* vax 11/750 unix cc		*/
/* POL originally #define EPSILON 1.0e-8 */
/* POL originally #define SQREPSILON 1.0e-16 */
#define EPSILON 1.0e-4
#define SQREPSILON 1.0e-8
#endif

/* dimension of arrays - can be increased*/
#define N 200   /* POL: originally 20 */


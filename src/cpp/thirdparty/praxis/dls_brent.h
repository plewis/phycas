/*----------------------------------------------------------------------------------------------------------------------
|
|	@file	dls_brent.h
|
|	@brief	Header file for optimization routines translated from Richard Brent's 1973 book.
|
|	This file contains declarations and definitions for a set of one-dimensional and multidimensional optimization
|	routines.  See the corresponding C file for more information.
*/
#ifndef DLS_BRENT_H_INCLUDED
#define DLS_BRENT_H_INCLUDED

#if defined(DLS_H_INCLUDED)

#	include "dls_types.h"
#	include "dls_minfunc.h"

#else	/* for standalone version... */

#	define gAborted 0		/* see explanation in dls_brent.c */

	typedef double (*MinimizeFxn)(double *, void *);
	typedef struct praxisData PraxisData;

	/* allocateBlock is a function that dynamically allocates a block of memory; replacing it with calloc is fine
	   unless you want to do additional error handling or you need to enforce non-standard alignment such as 128-byte
	   alignment for SSE/Altivec.  Note that we assume that memory is initialized to zeros, so calloc is used rather
	   than malloc.
	*/
#	define allocateBlock(nBytes)	calloc((size_t)(nBytes), (size_t)1)
#	define deallocateBlock(p)		(p == NULL) ? 0 : free(p)

#	define dls_assert assert
#	define dls_printf printf

#endif

#ifdef __cplusplus
extern "C" {
#endif

struct praxisData
	{
	double *	x;			/**< variables */
	double **	v;			/**< n x n matrix in row-ptr form */
	double *	work;		/**< work vector of size 6*n */
	double		machEps;	/**< machine epsilon (will be set based on single vs. double precision) */
	Bool		illc;		/**< set true for ill-conditioned problems, otherwise leave at false */
	double		glimit;		/**< defaults to Brent's setting of 100.0; can be changed in needed for stability */
	};

extern PraxisData *	newPraxisData(double *x, unsigned np, Bool singlePrecision);
extern void			deletePraxisData(PraxisData *praxd);

extern DoublePair	bracketMinimum(double a, double b, MinimizeFxn fxn, void *dataPtr, double glimit);
extern double		praxis(double tol, double h, int n, MinimizeFxn fxn, void *dataPtr, PraxisData *praxd);
extern double		localMin(double, double, double, double, MinimizeFxn, void *, double *);
extern double		minimizeBrent1D(double a, double b, double eps, double t, double glimit, MinimizeFxn fxn,
					                void * dataPtr, double *px);
extern int			checkPraxis(void);

#ifdef __cplusplus
}
#endif

#endif	/* DLS_BRENT_H_INCLUDED */

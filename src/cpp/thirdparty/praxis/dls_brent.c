/*----------------------------------------------------------------------------------------------------------------------
|
|	@file	dls_brent.c
|
|	@brief	Optimization routines from Richard Brent's 1973 book
|
|	Routines for function minimization based on those of Brent (1973).  Brent's routines were translated from Algol
|	to C and cleaned up (goto elimination, provision for dynamic memory allocation, etc.).
|
|	@references		Brent, R. P. 1973. Algorithms for Minimization Without Derivatives.  Prentice-Hall, Englewood
|					Cliffs, N.J.  Republished 2002 by Dover Publications, Mineola, New York.
|
|	@references		Brent, R. P. 1975. (Erratum) Mathematics of Computation Table Errata 520:1166.
|
|	Because this code is a straightforward translation of publicly available Algol code, it is hereby placed into the
|	public domain with no other licensing restrictions.  However, since I spent a lot of time doing this translation
|	and cleaning up the code, I would appreciate your acknowledgement of my contribution in any derived works.
|
|		David L. Swofford
|		Duke University
|
|	Usage notes:
|
|		The optimization routines require a function that returns the value of the function being minimized.
|		This function needs two arguments:
|			 1.	a "double *" pointer to a vector of parameter values, or, for one-dimensional optimizations, the
|			 	address of a single parameter value.
|			 2.	a "void *" pointer to a structure containing any additional data that the function being minimized
|			 	needs to compute its value.  This pointer can be NULL if no additional data are required.
|		It should return the value of the function, evaluated for the supplied vector of parameter values.
|
|		============
|		praxis usage
|		============
|
|			 1. Allocate and initialize a praxis object using newPraxisData:
|
|					PraxisData * praxd = newPraxisData(x, n, single_precision);
|
|				  		x = array for parameter vector (must be large enough to store n doubles)
|						n = number of parameters in function being optimized
|						single_precision = true if function evaluation is only done at single-precision accuracy,
|						                   false otherwise (i.e., function is evaluated in double precision).
|
|			 2.	Initialize x to contain the starting parameter values (in most cases, the better the guess, the faster
|			    the convergence will be)
|
|			 3. Invoke the praxis routine:
|
|			 		fx_opt = praxis(tol, h, n, f, data, praxd);
|
|						tol   = tolerance limit used to test for convergence (if passed as 0, a default of 1e-5
|						        will be used)
|						h     = maximum step size used to test for convergence (if passed as 0, a default of 1.0
|							    will be used)
|						n     = the number of parameters in the function being minimized
|						f	  = the function being minimized (see above)
|						data  = an option pointer to arbitrary data that the function may need to compute its value;
|						        if the function needs no additional data, pass NULL
|						praxd = a praxis object previously allocated by newPraxisData
|
|					On return, fx_opt will be the value of the minimized function, and praxd->h will contain the
|					corresponding optimized parameter values.
|
|			 4. Destroy the praxis object (to avoid leaking memory);
|
|					deletePraxisData(praxd);
|
|	@note	'gAborted' is a variable that can be set to request cancellation without terminating the program (e.g., by
|			a signal handler invoked when ctrl-C is pressed).  By default it is just #define'd to 0 in the header; if
|			you want to provide this capability, you will need to provide a mechanism for setting gAborted externally.
|
|	Most recent revision: 9 February 2012
|
|	@todo	Need writeups like "praxis usage" above for localMin and bracketMinimum
*/

/* This file is ordinarily used in conjunction with a package that uses a shared set of standard headers as well as
   some additional support utility routines.  If you want to build it stand-alone, remove the line below that includes
   the "dls.h" header.
*/
//#include "dls.h"
#if defined(DLS_H_INCLUDED)
#	define USING_DLS
#	include "dls_math.h"
#	include "dls_matrices.h"
#	include "dls_signal.h"
#	include "dls_useractivity.h"
#	include "dls_random.h"
#	include "dls_malloc.h"
#else
#	include "dls_core.h"
#	include <stdio.h>		/* for debugging output */
#	include <stdlib.h>		/* for malloc, rand */
#	include <math.h>		/* for fabs, fmin, sqrt, pow */
#	include <float.h>		/* for FLT_EPSILON, DBL_EPSILON */
#	include <assert.h>		/* for assert */
#	define STANDALONE
	static void *	newMatrix(size_t elementSize, unsigned nrows, unsigned ncols);
	static void *	setMatrixRowPtrs(void *a_arg, void *p_arg, unsigned nrows, unsigned ncols, size_t elementSize);
	static void		deleteMatrix(void *a);
	static void		setIdentityMatrix(double **a, unsigned n);
	static double	randUni(int *);
#endif

#include "dls_brent.h"

#if defined(ALLOW_DEBUG_CODE)		/* debug code generated only for items that are #define'd below */
#	undef CHECK_PROGRESS			/*	show one-line iteration-progress output */
#	undef CHECK_PROGRESS_OFTEN		/*	do above for every function evaluation */
#	undef CHECK_LOCALMIN_PROGRESS	/*	check progress in Brent's "localMin" procedure */
#	undef DEBUG_BRACKET				/*	check bracketing routine */
#	undef DEBUG_PRAXIS				/*	debugging for praxis  */
#	undef CHECK_FINAL_PRAXIS		/*	show final optimized parameter vector */
#endif

/** This structure localizes variables that had the scope of function 'praxis' in Algol (which supports nested
	functions).  If they are just made global to the file, praxis will not be reentrant, so they are stored in
	a structure that gets passed around.
*/
typedef struct
	{
	double *	xnew;
	double		dmin;
	double		ldt;
	double		qd0;
	double		qd1;
	double		qf1;
	double		machEps;
	double		m2, m4;
	double		esmall;			/* Brent's 'small' to 'esmall' due to conflict with Win32 "#define small char" */
	double		toler;
	double		htol;
	double		fx;
	int			praxisSeed;
	int			n;
	int			nl;				/* counter for number of line minimizations performed */
	int			nf;				/* counter for number of function evaluations */
	}
	PraxisLocal;

#if defined(TEST_BRENT)
#	if !defined(CHECK_PRAXIS)
#		define CHECK_PRAXIS
#	endif
#	if !defined(CHECK_PROGRESS)
#		define CHECK_PROGRESS
#	endif
#endif

#if defined(CHECK_PRAXIS) || defined(CHECK_PROGRESS)
#	if defined(CHECK_PROGRESS)
		static void showPraxisProgress(PraxisData *praxd, PraxisLocal *pld, int flag);
#	endif
	static double g_fmin;	/* estimates the value of function f at the minimum; used for printing log(fx-f_min) */
#endif

#if defined(USING_DLS)
#	define VOID_RETURN_IF_ABORTED    if ((gAborted < 0) || ((gAborted > 0) && confirmAbortRequest())) return;
#	define RETURN_IF_ABORTED(result) if ((gAborted < 0) || ((gAborted > 0) && confirmAbortRequest())) return result;
#else
#	define VOID_RETURN_IF_ABORTED    if (gAborted != 0) return;
#	define RETURN_IF_ABORTED(result) if (gAborted != 0) return result;
#endif

typedef struct
	{
	double	x;
	double	fx;
	}
	FxnPoint;

/* Debugging for bracketMinimum... */
#if defined(DEBUG_BRACKET)
	static void debug_abc(char *s, FxnPoint a, FxnPoint b, FxnPoint c)
		{ dls_printf("%s (%g, %g, %g): (%f, %f, %f)\n", s, a.x, b.x, c.x, a.fx, b.fx, c.fx); }
	static void debugShowSuccess(FxnPoint b, FxnPoint c)
		{ if (b.fx <= c.fx) dls_printf("  this is a bracket!\n"); }
#	define debug_bracket(a) dls_printf a
#else
#	define debug_abc(s, a, b, c)
#	define debugShowSuccess(b, c)
#	define debug_bracket(a)
#endif

/* Debugging for praxis... */
#if defined(DEBUG_PRAXIS)
#	define debug_praxis(a) dls_printf a
#	define DBG_SPACES 50
#else
#	define debug_praxis(a)
#endif

/*----------------------------------------------------------------------------------------------------------------------
|
|	Allocates memory and initializes data structure for Praxis optimization.
*/
PUBLIC PraxisData * newPraxisData(
  double * x,			/**< double array of size n used for parameter vector */
  unsigned n,			/**< number of parameters */
  Bool singlePrecision)	/**< if true, function is being calculated in single-precision (affects tolerances, etc.) */
	{
	PraxisData * praxd = allocateBlock(sizeof(PraxisData));
	if (praxd != NULL)
		{
		praxd->work = allocateBlock(6*n*sizeof(double));
		if (praxd->work != NULL)
			{
			praxd->v = newMatrix(sizeof(double), n, n);
			if (praxd->v != NULL)
				{
				praxd->x = x;
				praxd->machEps = singlePrecision ? FLT_EPSILON : DBL_EPSILON;
				}
			}
		if ((praxd->work == NULL) || (praxd->v == NULL))
			{
			deletePraxisData(praxd);
			praxd = NULL;
			}
		praxd->glimit = 100.0;		/* Brent's setting, but caller can modify if needed for stability */
		}

	return praxd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Deallocates memory allocated by newPraxisData.
*/
PUBLIC void deletePraxisData(
  PraxisData *praxd)		/**< a Praxis object previously allocated by newPraxisData */
	{
	if (praxd != NULL)
		{
		deallocateBlock(praxd->work);
		deleteMatrix(praxd->v);
		deallocateBlock(praxd);
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Locates a minimizing value by inverse parabolic interpolation
|
|	Returns a value x that such that f(x) is the minimum value of a parabola that passes through f(a), f(b), and f(c).
*/
PRIVATE double extrapolateParabola(FxnPoint a, FxnPoint b, FxnPoint c)
	{
	double x;

	double b_minus_a = b.x - a.x;
	double b_minus_c = b.x - c.x;
	double t = b_minus_a*(b.fx - c.fx);
	double u = b_minus_c*(b.fx - a.fx);
	double diff = t - u;

	if (fabs(diff) > 1e-37)
		{
		/* normal case */
		x = b.x - 0.5*(b_minus_a*t - b_minus_c*u)/diff;
		}
	else
		{
		/* a, b, and c are (nearly) colinear -- just return the value for which the function value is smallest */
		double min_fx = fmin(a.fx, fmin(b.fx, c.fx));
		if (a.fx == min_fx)
			x = a.x;
		else if (b.fx == min_fx)
			x = b.x;
		else
			x = c.x;
		}

	return x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Magnifies a potential bracketing interval using the golden ratio.
|
|	The input points `a' and `b' are used to compute a new value `c', and the function is evaluated at this new point.
*/
PRIVATE FxnPoint magnifyUsingGoldenRatio(FxnPoint a, FxnPoint b, MinimizeFxn fxn, void *data)
	{
	FxnPoint c;

#	define GOLDEN_RATIO 1.618033988
	c.x = b.x + GOLDEN_RATIO*(b.x - a.x);
	c.fx = fxn(&c.x, data);
	return c;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Brackets a function minimum using the method described by Press et al. (but completely written by DLS using the
|	ideas presented in the text and code).
|
|	This function is reentrant.
*/
PUBLIC DoublePair bracketMinimum(
  double		a_guess,	/**< first initial guess */
  double		b_guess,	/**< second initial guess */
  MinimizeFxn	fxn,		/**< the function being minimized */
  void *		data,		/**< data to be passed to function (can be NULL) */
  double		glimit)		/**< maximum expansion allowed for a parabolic fit step (typically 100, but can be smaller
  							     if needed for stability).  Pass 0 to get the default value, currently 100 */
	{
	FxnPoint	a, b, c, u;
	DoublePair	bracket = STRUCT_ZERO;

	if (glimit == 0.0)
		glimit = 100.0;

	a.fx = fxn(&a_guess, data); RETURN_IF_ABORTED(bracket)
	b.fx = fxn(&b_guess, data); RETURN_IF_ABORTED(bracket)
	a.x = a_guess;
	b.x = b_guess;
	if (b.fx > a.fx)
		{
		/* Swap a and b so we can go downhill from a to b */
		FxnPoint temp = a;
		a = b;
		b = temp;
		}
	c = magnifyUsingGoldenRatio(a, b, fxn, data);		/* first guess for c */

	debug_abc("\nbracketMinimum initial guess = ", a, b, c);
	debugShowSuccess(b, c);

	while (b.fx > c.fx)
		{
		/* Get a new trial point u via inverse parabolic interpolation */
		u.x = extrapolateParabola(a, b, c);
		if ((b.x > u.x) == (u.x > c.x))
			{
			/* u is between b and c: check for a minimum between a and c  */
			u.fx = fxn(&u.x, data); RETURN_IF_ABORTED(bracket)
			if (u.fx < c.fx)
				{
				/* f(b) > f(u) < f(c): a minimum exists between b and c */
				a = b;
				debug_bracket(("  1: (b, u, c) is a bracket!\n"));
				break;		/* we're done */
				}
			else if (u.fx > b.fx)
				{
				/* f(a) > f(b) < f(u): a minimum exists between a and u */
				c = u;
				debug_bracket(("  2: (a, b, u) is a bracket!\n"));
				break;		/* we're done */
				}
			else
				{
				/* Parabolic fit failed; (a, b, c) --> (b, c, new magnified value) */
				a = b;
				b = c;
				c = magnifyUsingGoldenRatio(a, b, fxn, data);
				debug_abc("  3: (a, b, c) from magnification of (b, c) => ", a, b, c);
				}
			}
		else
			{
			double xlim = b.x + glimit*(c.x - b.x);
			if ((c.x > u.x) == (u.x > xlim))
				{
				/* u from parabolic fit is between c and allowed limit */
				u.fx = fxn(&u.x, data); RETURN_IF_ABORTED(bracket)
				if (u.fx < c.fx)
					{
					/* (a, b, c) --> (c, u, new magnified value) */
					a = c;
					b = u;
					c = magnifyUsingGoldenRatio(a, b, fxn, data);
					debug_abc("  4a: (a, b, c) from magnification of (c, u) => ", a, b, c);
					}
				else
					{
					/* Discard oldest point and continue iteration; (a, b, c) --> (b, c, u) */
					a = b;
					b = c;
					c = u;
					debug_abc("  4b: (a, b, c) from discard of old a => ", a, b, c);
					}
				}
			else
				{
				if ((u.x > xlim) == (xlim > c.x))
					{
					/* Limit parabolic fit to maximum allowed value; (a, b, c) --> (b, c, xlim) */
					u.x = xlim;
					u.fx = fxn(&u.x, data); RETURN_IF_ABORTED(bracket)
					debug_abc("  5a: (a, b, c) from pinning parabolic u to boundary => ", b, c, u);
					}
				else
					{
					/* Reject parabolic u; (a, b, c) --> (b, c, new magnified value) */
					u = magnifyUsingGoldenRatio(b, c, fxn, data);
					debug_abc("  5b: (a, b, c) from magnification of (b, c) => ", b, c, u);
					}
				a = b;
				b = c;
				c = u;
				}
			debugShowSuccess(b, c);
			}
		}
	if (a.x < c.x)
		{
		bracket.a = a.x;
		bracket.b = c.x;
		}
	else
		{
		bracket.a = c.x;
		bracket.b = a.x;
		}
	return bracket;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Finds a local minimum of a function of one variable using Brent's (1973) method.
|
|	This a translation of Brent's (1973; p 79-80) ALGOL 60 procedure localmin.
|
|	This function is reentrant.
|
|	@return		The optimized value of the function.  The value of x that optimizes f(x) is returned to *px.
*/
PUBLIC double localMin(
  double		a,		/**< on input, (a,b) must bracket a local minimum (be as pessimistic as you need to, but
  						     not more than necessary) */
  double		b,		/**< see `a' */
  double		eps,	/**< t and eps define tol = eps|x| + t; the function is never evaluated at two points closer
							 together than tol. eps should be at least as large as sqrt(machine_epsilon), where
							 machine_epsilon is typically DBL_EPSILON or FLT_EPSILON defined in float.h.  If eps is
							 passed as 0.0, it is set to sqrt(DBL_EPSILON), which would not be appropriate if the
							 function evaluation is done to single-precision accuracy, even if the results are
							 reported as double precision. */
  double		t,		/**< see `eps' */
  MinimizeFxn	f,		/**< the function to minimize */
  void *		data,	/**< pointer to data to be passed to 'f' */
  double *		px)		/**< value of x when f(x) is minimal */
	{
#	define MAX_ITER 50
#	define CGOLD    0.3819660113	/* = (3 - sqrt(5))/2 */

	FxnPoint	u, v, w, x;
	double		e, p, q, r, tol, t2,
				d = 0.0;	/* shuts up lint false-positive */
	int			iter;

#	if defined(CHECK_LOCALMIN_PROGRESS)
		dls_printf("\nEntering localMin with a=%.6g b=%.6g\n", a, b);
#	endif

	if (eps == 0.0)
		{
		/* REMEMBER: This default epsilon value assumes double precision */
		eps = sqrt(DBL_EPSILON);
		}

	x.x = a + CGOLD*(b - a);
	x.fx = f(&x.x, data);	RETURN_IF_ABORTED(0.0)
	v = w = x;
	e = 0.0;
#	if defined(CHECK_LOCALMIN_PROGRESS)
		dls_printf("x <= %.6g, f(x)=%f\n", x.x, x.fx);
#	endif

	/* Main loop */

	for (iter = 0; iter < MAX_ITER; iter++)
		{
		double m = 0.5*(a + b);
#		if defined(CHECK_LOCALMIN_PROGRESS)
			dls_printf("iter %d, a=%f b=%f\n", iter, a, b);
#		endif
		tol = eps*fabs(x.x) + t;
		t2 = 2.0*tol;

		/* Check stopping criterion */
		if (fabs(x.x - m) <= t2 - 0.5*(b - a))
			break;

		p = q = r = 0.0;
		if (fabs(e) > tol)
			{
			/* Fit parabola (trial) */
			r = (x.x - w.x)*(x.fx - v.fx);
			q = (x.x - v.x)*(x.fx - w.fx);
			p = (x.x - v.x)*q - (x.x - w.x)*r;
			q = 2.0*(q - r);
			if (q > 0.0)
				p = -p;
			else
				q = -q;
			r = e;	/* save 'e' for test below before changing it */
			e = d;	/* lint complains about possible use of 'd' before being set, but it's not a problem because e=0
			           first time through, so "fabs(e) > tol" test fails */
			}

		/* Take parabolic-interpolation or golden-section step (note that Brent's Algol procedure had an error;
		   p > q*(a-x) is correct) */
		if ((fabs(p) < fabs(0.5*q*r)) && (p > q*(a - x.x)) && (p < q*(b - x.x)))
			{
			/* Parabolic interpolation step */
			d = p/q;
			u.x = x.x + d;
			/* Don't evaluate f too close to a or b */
			if ((u.x - a < t2) || (b - u.x < t2))
				d = (x.x < m) ? tol : -tol;
			}
		else
			{
			/* "Golden section" step */
			e = ((x.x < m) ? b : a) - x.x;
			d = CGOLD*e;
			}

		/* Don't evaluate function too close to x */
		if (fabs(d) >= tol)
			u.x = x.x + d;
		else if (d > 0.0)
			u.x = x.x + tol;
		else
			u.x = x.x - tol;

		u.fx = f(&u.x, data);	RETURN_IF_ABORTED(0.0)
#		if defined(CHECK_LOCALMIN_PROGRESS)
			dls_printf("f(%f) = %f\n", u.x, u.fx);
#		endif
		/* Update a, b, v, w, and x */
		if (u.fx <= x.fx)
			{
			if (u.x < x.x)
				b = x.x;
			else
				a = x.x;
			v = w;
			w = x;
			x = u;
			}
		else
			{
			if (u.x < x.x)
				a = u.x;
			else
				b = u.x;
			if ((u.fx <= w.fx) || (w.x == x.x))
				{
				v = w;
				w = u;
				}
			else if ((u.fx <= v.fx) || (v.x == x.x) || (v.x == w.x))
				v = u;
			}
		}
	*px = x.x;
	return x.fx;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Front end to Brent's "localmin" method that establishes a bracketing interval first.
|
|	Returns the function value at the optimum.
*/
PUBLIC double minimizeBrent1D(
  double		a,			/**< on input, (a,b) must bracket a local minimum */
  double		b,			/**< see `a' */
  double		eps,		/**< t and eps define tol = eps|x| + t; see localMin for details */
  double		t,			/**< see `eps' */
  double		glimit,		/**< value passed to to 'bracketMinimum' (see that function for details); can be 0 for
  							     a default value */
  MinimizeFxn	f,			/**< function to minimize */
  void *		data,		/**< pointer to data to be passed to 'f' */
  double *		px)			/**<value of x when f(x) is minimal */
	{
	DoublePair bracket = bracketMinimum(a, b, f, data, glimit);
	if (gAborted)
		return 0.0;

	return localMin(bracket.a, bracket.b, eps, t, f, data, px);
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Evaluates a function of one variable after moving trial point a distance lambda from the initial point in the
|	direction v[j] (if j >= 0) or performs a curvilinear extrapolation (if j < 0).
*/
PRIVATE double praxis_FLin(PraxisLocal *pld, int j, double lambda, double *x, double *xnew, double *q0, double *q1,
                           double **v, MinimizeFxn f, void *data)
	{
	int i;
	int n = pld->n;

	if (j >= 0)
		{
		/* Try point in linear search */
		for (i = 0; i < n; i++)
			xnew[i] = x[i] + lambda*v[i][j];
		}
	else
		{
		/* Search along a parabolic space curve */
		double qd0 = pld->qd0;
		double qd1 = pld->qd1;
		double qa = lambda*(lambda - qd1)/(qd0*(qd0 + qd1));
		double qb = (lambda + qd0)*(qd1 - lambda)/(qd0*qd1);
		double qc = lambda*(lambda + qd0)/(qd1*(qd0 + qd1));

		/* Previous three points were stored as follows: x' in q0, x'' in x, and x''' in q1 */
		for (i = 0; i < n; i++)
			xnew[i] = qa*q0[i] + qb*x[i] + qc*q1[i];
		}
	pld->nf++;
	return (*f)(xnew, data);
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Minimizes the objective function from point x in the direction v(*,j) (if j >= 0), or performs a quadratic search
|	in the plane defined by q0, q1, and x (if j < 0).
*/
PRIVATE void praxis_Min(
  PraxisLocal * pld,	/* global data */
  int			j,		/* column of direction matrix (or <0 flag for quadratic search) */
  int			nits,	/* number of times an attempt is made to halve the interval */
  double *		x,		/* current point */
  double *		q0,		/* previous point */
  double *		q1,		/* second previous point */
  double **		v,		/* direction-set matrix */
  double *		pd2,	/* approximation to half f'' (or zero) */
  double *		px1,	/* x1 = estimate of distance to minimum, returned as actual distance found */
  double		f1,		/* if fk=true, FLin(x1), otherwise ignored */
  Bool			fk,		/* flag (see `f1') */
  MinimizeFxn	f,		/* function to minimize (see dls_brent.h) */
  void * data)			/* pointer to data to be passed to 'f' */
	{
	int		i, k;
	Bool	success, dz;
	double	x1, x2, xm, f0, f2, fm, d1, d2, t2, s, sf1, sx1;

	/* Copy args passed by reference to locals (will pass back at end) */
	d2 = *pd2;
	x1 = *px1;

	sf1 = f1;
	sx1 = x1;
	k = 0;
	xm = 0.0;
	f0 = fm = pld->fx;
	dz = (d2 < pld->machEps);	/* if true, we need f''(0) */

	/* Find step size */
	s = 0.0;
	for (i = 0; i < pld->n; i++)
		{
		double xi = x[i];
		s += xi*xi;
		}
	s = sqrt(s);
	t2 = pld->m4*sqrt(fabs(pld->fx)/(dz ? pld->dmin : d2) + s*pld->ldt) + pld->m2*pld->ldt;
	s = pld->m4*s + pld->toler;
	if (dz && (t2 > s))
		t2 = s;
	if (t2 < pld->esmall)
		t2 = pld->esmall;
	if (t2 > 0.01*pld->htol)
		t2 = 0.01*pld->htol;

	if (fk && (f1 <= fm))
		{
		xm = x1;
		fm = f1;
		}
	if (!fk || (fabs(x1) < t2))
		{
		x1 = (x1 >= 0.0) ? t2 : -t2;
		f1 = praxis_FLin(pld, j, x1, x, pld->xnew, q0, q1, v, f, data);
		VOID_RETURN_IF_ABORTED
		debug_praxis(("%*cf1 from praxis_FLin 1 => %g\n", DBG_SPACES, ' ', f1));
		}
	if (f1 <= fm)
		{
		xm = x1;
		fm = f1;
		}

	/* Find a distance x2 ("lambda*") that approximately minimizes f in the chosen direction */
	do	{
		if (dz)
			{
			/* Evaluate FLin at another point and estimate the second derivative */
			x2 = (f0 < f1) ? -x1 : 2.0*x1;

			debug_praxis(("%*ccalling praxis_FLin with x1=%g x2=%g\n", DBG_SPACES, ' ', x1, x2));
			f2 = praxis_FLin(pld, j, x2, x, pld->xnew, q0, q1, v, f, data);
			VOID_RETURN_IF_ABORTED
			debug_praxis(("%*cf2 from praxis_FLin 2 => %g\n", DBG_SPACES, ' ', f2));
			if (f2 <= fm)
				{
				xm = x2;
				fm = f2;
				}
			d2 = (x2*(f1 - f0) - x1*(f2 - f0))/(x1*x2*(x1 - x2));
			}

		/* Estimate first derivative at 0 */
		d1 = (f1 - f0)/x1 - x1*d2;
		dz = TRUE;

		/* Predict minimum */
		if (d2 <= pld->esmall)
			x2  = (d1 < 0.0) ? pld->htol : -pld->htol;
		else
			{
			x2 = -0.5*d1/d2;
			if (fabs(x2) > pld->htol)
				x2 = (x2 > 0.0) ? pld->htol : -pld->htol;
			}

		/* Evaluate f at predicted minimum */
		do	{
			f2 = praxis_FLin(pld, j, x2, x, pld->xnew, q0, q1, v, f, data);
			VOID_RETURN_IF_ABORTED
			debug_praxis(("%*cf2 from praxis_FLin 3 => %g\n", DBG_SPACES, ' ', f2));
			success = TRUE;
			if ((k < nits) && (f2 > f0))
				{
				/* No success so halve interval and try again */
				success = FALSE;
				k++;
				if ((f0 < f1) && (x1*x2 > 0.0))
					break;
				x2 *= 0.5;
				}
			}
			while (!success);
		}
		while (!success);

	pld->nl++;	/* increment counter for number of line searches */

	if (f2 > fm)
		x2 = xm;
	else
		fm = f2;

	/* Get new estimate of second derivative */
	if (fabs(x2*(x2 - x1)) > pld->esmall)
		d2 = (x2*(f1 - f0) - x1*(fm - f0))/(x1*x2*(x1 - x2));
	else if (k > 0)
		d2 = 0.0;
	if (d2 < pld->esmall)
		d2 = pld->esmall;

	x1 = x2;
	pld->fx = fm;
	if (sf1 < fm)
		{
		pld->fx = sf1;
		x1 = sx1;
		}

	/* Update x for linear search but not for parabolic search */
	if (j >= 0)
		{
		for (i = 0; i < pld->n; i++)
			x[i] += x1*v[i][j];
		}

	*px1 = x1;
	*pd2 = d2;
	debug_praxis(("%*creturning from praxis_Min with x1=%g d2=%g\n", DBG_SPACES, ' ', x1, d2));
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Looks for the function minimum along a curve defined by q0, q1, and x.
*/
PRIVATE void praxis_Quad(PraxisLocal *pld, double *x, double *q0, double *q1, double **v, MinimizeFxn f, void *data)
	{
	int		i,
			n = pld->n;
	double	lambda, s, qa, qb, qc, qd0, qd1;

	/* q0 and q1 contain previous two points */

	s = pld->fx;
	pld->fx = pld->qf1;
	pld->qf1 = s;

	qd1 = 0.0;
	for (i = 0; i < n; i++)
		{
		s = x[i];
		x[i] = q1[i];
		q1[i] = s;
		s = q1[i] - x[i];
		qd1 += s*s;
		}
	pld->qd1 = qd1 = sqrt(qd1);
	pld->qd0 = qd0 = pld->qd0;
	if ((qd0 > 0.0) && (qd1 > 0.0) && (pld->nl >= 3*n*n))
		{
		s = 0.0;
		lambda = qd1;

		debug_praxis(("%*ccalling praxis_Min from Quad 1\n", DBG_SPACES, ' '));
		praxis_Min(pld, -1, 2, x, q0, q1, v, &s, &lambda, pld->qf1, TRUE, f, data);
		VOID_RETURN_IF_ABORTED
		debug_praxis(("%*clambda from praxis_Min => %g\n", DBG_SPACES, ' ', lambda));
		qa = lambda*(lambda - qd1)/(qd0*(qd0 + qd1));
		qb = (lambda + qd0)*(qd1 - lambda)/(qd0*qd1);
		qc = lambda*(lambda + qd0)/(qd1*(qd0 + qd1));
		}
	else
		{
		debug_praxis(("%*cno praxis_Min call from Quad 2\n", DBG_SPACES, ' '));
		pld->fx = pld->qf1;
		qa = qb = 0.0;
		qc = 1.0;
		}
	qd0 = qd1;
	for (i = 0; i < n; i++)
		{
		s = q0[i];
		q0[i] = x[i];
		x[i] = qa*s + qb*x[i] + qc*q1[i];
		}
	pld->qd0 = qd0;
	pld->qd1 = qd1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Gets singular value decomposition using a modified version of Golub and Reinsch's (1969) routine, restricted to
|	m=n.  The singular values of the array 'ab' are returned in 'q', and 'ab' is overwritten with the orthogonal
|	matrix V such that U.diag(Q) = AB.V, where U is another orthogonal matrix.
*/
PRIVATE void praxis_MinFit(double eps, double tol, double **ab, double *q,
                         double *e,	/* work vector of size n */
                         int n)
	{
	int			i, j, k, l, l2, kt;
	Bool		doCancellation;
	double		c, f, g, h, s, x, y, z;

	dls_assert(n > 0);		/* to suppress static analyzer warnings about uninitialized variables */

	/* Householder's reduction to bidiagonal form */
	g = x = 0.0;
	for (i = 0; i < n; i++)
		{
		e[i] = g;
		s = 0.0;
		l = i + 1;
		for (j = i; j < n; j++)
			s += ab[j][i]*ab[j][i];
		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[i][i];
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[i][i] = f - g;
			for (j = l; j < n; j++)
				{
				f = 0.0;
				for (k = i; k < n; k++)
					f += ab[k][i]*ab[k][j];
				f /= h;
				for (k = i; k < n; k++)
					ab[k][j] += f*ab[k][i];
				}
			}
		q[i] = g;
		s = 0.0;
		for (j = l; j < n; j++)
			s += ab[i][j]*ab[i][j];

		if (s < tol)
			g = 0.0;
		else
			{
			f = ab[i][i + 1];		/* (DLS note: when i=n-1, s is always 0) */
			g = (f < 0.0) ? sqrt(s) : -sqrt(s);
			h = f*g - s;
			ab[i][i + 1] = f - g;
			for (j = l; j < n; j++)
				e[j] = ab[i][j]/h;
			for (j = l; j < n; j++)
				{
				s = 0.0;
				for (k = l; k < n; k++)
					s += ab[j][k]*ab[i][k];
				for (k = l; k < n; k++)
					ab[j][k] += s*e[k];
				}
			}
		y = fabs(q[i]) + fabs(e[i]);
		if (y > x)
			x = y;
		}

	/* Accumulation of right-hand transformations */
	for (i = n - 1; i >= 0; i--)
		{
		if (g != 0.0)
			{
			h = ab[i][i + 1]*g;
			for (j = l; j < n; j++)
				ab[j][i] = ab[i][j]/h;
			for (j = l; j < n; j++)
				{
				s = 0.0;
				for (k = l; k < n; k++)
					s += ab[i][k]*ab[k][j];
				for (k = l; k < n; k++)
					ab[k][j] += s*ab[k][i];
				}
			}
		for (j = l; j < n; j++)
			ab[i][j] = ab[j][i] = 0.0;
		ab[i][i] = 1.0;
		g = e[i];
		l = i;
		}

	/* Diagonalization of the bidiagonal form */
	eps *= x;
	for (k = n - 1; k >= 0; k--)
		{
		kt = 0;

		test_splitting:

		if (++kt > 30)
			{
			e[k] = 0.0;
#			if defined(TEST_BRENT)
				puts("QR failed");
#			endif
			}
		doCancellation = TRUE;
		for (l2 = k; l2 >= 0; l2--)
			{
			l = l2;
			if (fabs(e[l]) <= eps)
				{
				doCancellation = FALSE;
				break;
				}
			else if (fabs(q[l - 1]) <= eps)
				break;
			}

		if (doCancellation)
			{
			/* Cancellation of e[l] if l > 1 */
			c = 0.0;
			s = 1.0;
			for (i = l; i <= k; i++)
				{
				f = s*e[i];
				e[i] *= c;
				if (fabs(f) <= eps)
					break;
				g = q[i];
				if (fabs(f) < fabs(g))						/* (implies g > eps) */
					h = fabs(g)*sqrt(1.0 + (f/g)*(f/g));
				else if (f != 0.0)
					h = fabs(f)*sqrt(1.0 + (g/f)*(g/f));
				else
					h = 0.0;
				q[i] = h;
				if (h == 0.0)
					g = h = 1.0;	/* Note: this replaces q[i]=h=sqrt(g*g+f*f) which may give incorrect results if the
									         squares underflow or if f=g=0 */
				c = g/h;
				s = -f/h;
				}
			}

		z = q[k];
		if (l != k)
			{
			/* Shift from bottom 2x2 minor */
			x = q[l];
			y = q[k - 1];
			g = e[k - 1];
			h = e[k];
			f = ((y - z)*(y + z) + (g - h)*(g + h))/(2.0*h*y);
			g = sqrt(f*f + 1.0);
			s = (f < 0.0) ? f - g : f + g;
			f = ((x - z)*(x + z) + h*(y/s - h))/x;

			/* Next QR transformation */
			c = s = 1.0;
			for (i = l + 1; i <= k; i++)
				{
				g = e[i];
				y = q[i];
				h = s*g;
				g *= c;
				/* DLS note: Golob and Reinsch had "z = sqrt(f*f + h*h))" below */
				if (fabs(f) < fabs(h))
					z = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z = 0.0;
				e[i - 1] = z;
				if (z == 0.0)
					z = f = 1.0;
				c = f/z;
				s = h/z;
				f = x*c + g*s;
				g = -x*s + g*c;
				h = y*s;
				y *= c;
				for (j = 0; j < n; j++)
					{
					x = ab[j][i - 1];
					z = ab[j][i];
					ab[j][i - 1] = x*c + z*s;
					ab[j][i] = -x*s + z*c;
					}
				/* DLS note: Golob and Reinsch had "z = sqrt(f*f + h*h))" below */
				if (fabs(f) < fabs(h))
					z = fabs(h)*sqrt(1.0 + (f/h)*(f/h));
				else if (f != 0.0)
					z = fabs(f)*sqrt(1.0 + (h/f)*(h/f));
				else
					z = 0.0;
				q[i - 1] = z;
				if (z == 0.0)
					z = f = 1.0;
				c = f/z;
				s = h/z;
				f = c*g + s*y;
				x = -s*g + c*y;
				}
			e[l] = 0.0;
			e[k] = f;
			q[k] = x;
			goto test_splitting;
			}

		if (z < 0.0)
			{
			/* q[k] is made non-negative */
			q[k] = -z;
			for (j = 0; j < n; j++)
				ab[j][k] = -ab[j][k];
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Sorts the elements of 'd' and corresponding elements of 'v' into descending order.
*/
PRIVATE void praxis_SortDV(double *d, double **v, unsigned n)
	{
	unsigned	i, j, k;
	double		s;

	for (i = 0; i < n - 1; i++)
		{
		k = i;
		s = d[i];
		for (j = i + 1; j < n; j++)
			{
			if (d[j] > s)
				{
				k = j;
				s = d[j];
				}
			}
		if (k > i)
			{
			d[k] = d[i];
			d[i] = s;
			for (j = 0; j < n; j++)
				{
				s = v[j][i];
				v[j][i] = v[j][k];
				v[j][k] = s;
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Minimizes a function f of n variables using Brent's (1973) principal axis method, which is a modification of
|	Powell's (1964) method.
|
|	This function is reentrant.
*/
PUBLIC double praxis(
  double		tol,		/* tolerance used for convergence criterion */
  double		h,			/* maximum step size (e.g., 1.0) */
  int			n,			/* number of parameters to function */
  MinimizeFxn	f,			/* the function to be minimized, declared as "double fxn(double *x, void *data)" */
  void *		data,		/* pointer to data to be passed to 'f' */
  PraxisData *	praxd)		/* an object already allocated by newPraxisData */
	{
	int			i, j, k, k2, kl, kt, ktm;
	Bool		illc;
	double		vsmall, large, vlarge, scbd, ldfac, df, f1, lds, t2, sl, dn, sz, *y, *z, *q0, *q1;
	PraxisLocal	pld = STRUCT_ZERO;

	double **	v = praxd->v;
	double *	d = praxd->work;
	double *	x = praxd->x;

	if (tol <= 0.0)
		tol = 1e-5;
	if (h <= 0.0)
		h = 1.0;

	/* Partition work space provided by caller into arrays used here */
	q0 = d + n;
	q1 = q0 + n;
	pld.xnew = q1 + n;
	y = pld.xnew + n;
	z = y + n;

	/* Store args in struct for reference elsewhere */
	pld.toler = tol;
	pld.htol = h;
	pld.n = (int)n;

	/* Machine-dependent initializations */
	pld.machEps = praxd->machEps;
	pld.m2 = sqrt(pld.machEps);
	pld.m4 = sqrt(pld.m2);
	pld.esmall = pld.machEps*pld.machEps;
	vsmall = pld.esmall*pld.esmall;
	large = 1.0/pld.esmall;
	vlarge = 1.0/vsmall;

	pld.praxisSeed = 1;		/* starting seed for random number generator */

	/* Heuristic numbers:
	   - if axes may be badly scaled (which should be avoided if possible), set scbd=10, otherwise 1
	   - if the problem is known to be ill-conditioned, set illc=TRUE, otherwise FALSE
	   - ktm+1 is the number of iterations without improvement before the algorithm terminates (see section 7.6 of
	     Brent's book).  ktm=4 is very cautious; usually ktm=1 is satisfactory.
	*/
	scbd = 1.0;
	illc = praxd->illc;
	ktm = 1;

	ldfac = illc ? 0.1 : 0.01;
	kt = pld.nl = 0;
	debug_praxis(("%*ccalling f from praxis 1\n", DBG_SPACES, ' '));
	pld.qf1 = pld.fx = f(x, data);
	RETURN_IF_ABORTED(0.0)
	debug_praxis(("%*cfx => %g\n", DBG_SPACES, ' ', pld.fx));
	pld.nf = 1;
	pld.toler = t2 = pld.esmall + fabs(pld.toler);
	pld.dmin = pld.esmall;
	if (pld.htol < 100.0*pld.toler)
		pld.htol = 100.0*pld.toler;
	pld.ldt = pld.htol;
	setIdentityMatrix(v, (unsigned)n);

	d[0] = pld.qd0 = 0.0;
	for (i = 0; i < n; i++)
		{
		q0[i] = 0.0;	/* DLS: this wasn't included in Brent's code, but it's important */
		q1[i] = x[i];
		}

#	if defined(CHECK_PROGRESS)
		showPraxisProgress(praxd, &pld, 0);
#	endif

	/* ------ main loop ------ */

	for (;;)
		{
		double sf = d[0];
		double s = d[0] = 0.0;

		/* Minimize along first direction */
		debug_praxis(("%*ccalling praxis_Min at 200\n", DBG_SPACES, ' '));
		praxis_Min(&pld, 0, 2, x, q0, q1, v, &d[0], &s, pld.fx, FALSE, f, data);
		RETURN_IF_ABORTED(0.0)
		if (s < 0.0)
			{
			for (i = 0; i < n; i++)
				v[i][0] = -v[i][0];
			}
		if ((sf <= 0.9*d[0]) || (0.9*sf >= d[0]))
			{
			for (i = 1; i < n; i++)
				d[i] = 0.0;
			}

		for (k = 1; k < n; k++)
			{
			for (i = 0; i < n; i++)
				y[i] = x[i];
			sf = pld.fx;
			illc = illc || (kt > 0);
			for (;;)
				{
				kl = k;
				df = 0.0;
				if (illc)
					{
					/* Random step to get off resolution valley */
					debug_praxis(("%*crandom step\n", DBG_SPACES, ' '));
					for (i = 0; i < n; i++)
						{
						s = z[i] = (0.1*pld.ldt + t2*pow(10.0, (double)kt))*(randUni(&pld.praxisSeed) - 0.5);
						for (j = 0; j < n; j++)
							x[j] += s*v[j][i];
						}
					pld.fx = f(x, data);
					RETURN_IF_ABORTED(0.0)
					pld.nf++;
#					if defined(CHECK_PROGRESS_OFTEN)
						showPraxisProgress(praxd, &pld, 'r');
#					endif
					}
				for (k2 = k; k2 < n; k2++)
					{
					sl = pld.fx;
					s = 0.0;
					/* Minimize along "non-conjugate" directions */
					debug_praxis(("%*ccalling praxis_Min at 201\n", DBG_SPACES, ' '));
					praxis_Min(&pld, k2, 2, x, q0, q1, v, &d[k2], &s, pld.fx, FALSE, f, data);
					RETURN_IF_ABORTED(0.0)
					if (illc)
						{
						sz = s + z[k2];
						s = d[k2]*sz*sz;
						}
					else
						s = sl - pld.fx;
					if (df < s)
						{
						df = s;
						kl = k2;
						}
					}
				if (!illc && (df < fabs(100.0*pld.machEps*pld.fx)))
					illc = TRUE;	/* no success with illc=false so try once with illc=true */
				else
					break;
				}

			for (k2 = 0; k2 < k; k2++)
				{
				/* Minimize along "conjugate" directions */
				s = 0.0;
				debug_praxis(("%*ccalling praxis_Min at 202\n", DBG_SPACES, ' '));
				praxis_Min(&pld, k2, 2, x, q0, q1, v, &d[k2], &s, pld.fx, FALSE, f, data);
				RETURN_IF_ABORTED(0.0)
				}

			f1 = pld.fx;
			pld.fx = sf;
			lds = 0.0;
			for (i = 0; i < n; i++)
				{
				sl = x[i];
				x[i] = y[i];
				y[i] = (sl -= y[i]);
				lds += sl*sl;
				}
			lds = sqrt(lds);
			if (lds > pld.esmall)
				{
				/* Throw away direction kl */
				for (i = kl - 1; i >= k; i--)
					{
					for (j = 0; j < n; j++)
						v[j][i + 1] = v[j][i];
					d[i + 1] = d[i];
					}

				/* Set new "conjugate" direction ... */
				d[k] = 0.0;
				for (i = 0; i < n; i++)
					v[i][k] = y[i]/lds;

				/* ... and minimize along it */
				debug_praxis(("%*ccalling praxis_Min at 203\n", DBG_SPACES, ' '));
				praxis_Min(&pld, k, 4, x, q0, q1, v, &d[k], &lds, f1, TRUE, f, data);
				RETURN_IF_ABORTED(0.0)
				if (lds <= 0.0)
					{
					lds = -lds;
					for (i = 0; i < n; i++)
						v[i][k] = -v[i][k];
					}
				}
			pld.ldt *= ldfac;
			if (pld.ldt < lds)
				pld.ldt = lds;
#			if defined(CHECK_PROGRESS)
				showPraxisProgress(praxd, &pld, 0);
#			endif
			t2 = 0.0;
			for (i = 0; i < n; i++)
				t2 += x[i]*x[i];
			t2 = pld.m2*sqrt(t2) + pld.toler;

			/* See if step length exceeds half the tolerance (stopping criterion) */
			kt = (pld.ldt > 0.5*t2) ? 0 : kt + 1;
			debug_praxis(("%*cChecking convergence: ldt=%g t2=%g kt=%d ktm=%d\n", DBG_SPACES, ' ',
						  pld.ldt, t2, kt, ktm));
			if (kt > ktm)
				{
#				if defined(DEBUG_PRAXIS) || defined(CHECK_FINAL_PRAXIS)
					debug_praxis(("%*cConvergence criterion satisfied\n", DBG_SPACES, ' '));
					dls_printf("\n+++++ praxis returning fx=%.14g (number of function evals=%d)\n", pld.fx, pld.nf);
					dls_printf("     final x = (");
					for (i = 0; i < n; i++)
						dls_printf("%13g", x[i]);
					dls_printf("\n");
#				endif
				return pld.fx;
				}
			}

		/* Try quadratic extrapolation in case we are stuck in a curved valley */
		debug_praxis(("%*cQuad\n", DBG_SPACES, ' '));
		praxis_Quad(&pld, x, q0, q1, v, f, data);
		RETURN_IF_ABORTED(0.0)

		/* Calculate V = U.(D^(-1/2))  (note: 'v' currently contains U) */
		dn = 0.0;
		for (i = 0; i < n; i++)
			{
			d[i] = 1.0/sqrt(d[i]);
			if (d[i] > dn)
				dn = d[i];
			}
		for (j = 0; j < n; j++)
			{
			s = d[j]/dn;
			for (i = 0; i < n; i++)
				v[i][j] *= s;
			}

		if (scbd > 1.0)
			{
			/* Scale axes in attempt to reduce condition number */
			debug_praxis(("%*cscaling\n", DBG_SPACES, ' '));
			s = vlarge;
			for (i = 0; i < n; i++)
				{
				sl = 0.0;
				for (j = 0; j < n; j++)
					sl += v[i][j]*v[i][j];
				z[i] = sqrt(sl);
				if (z[i] < pld.m4)
					z[i] = pld.m4;
				if (s > z[i])
					s = z[i];
				}
			for (i = 0; i < n; i++)
				{
				sl = s/z[i];
				z[i] = 1.0/sl;
				if (z[i] > scbd)
					{
					sl = 1.0/scbd;
					z[i] = scbd;
					}

				/* DLS 11apr02: erratum published in Mathematics of Computation TE 520 (1975):1166, found when I got
				   the 2002 Dover edition of Brent's book and checked the web site given there */
				for (j = 0; j < n; j++)
					v[i][j] *= sl;
				}
			}
		debug_praxis(("%*cresetting principal vectors\n", DBG_SPACES, ' '));

		/* Transpose v for MinFit */
		for (i = 1; i < n; i++)
			{
			for (j = 0; j < i; j++)
				{
				s = v[i][j];
				v[i][j] = v[j][i];
				v[j][i] = s;
				}
			}

		/* Find the singular value decomposition of v.  This gives the eigenvalues and principal axes of the
		   approximating quadratic form without squaring the condition number. */
		praxis_MinFit(pld.machEps, vsmall, v, d, y, (int)n);		/* ('y' is just a scratch vector) */

		if (scbd > 1.0)
			{
			/* Unscaling */
			for (i = 0; i < n; i++)
				{
				s = z[i];
				for (j = 0; j < n; j++)
					v[i][j] *= s;
				}
			for (i = 0; i < n; i++)
				{
				s = 0.0;
				for (j = 0; j < n; j++)
					s += v[j][i]*v[j][i];
				s = sqrt(s);
				d[i] *= s;
				s = 1.0/s;
				for (j = 0; j < n; j++)
					v[j][i] *= s;
				}
			}

		for (i = 0; i < n; i++)
			{
			s = dn*d[i];
			if (s > large)
				d[i] = vsmall;
			else if (s < pld.esmall)
				d[i] = vlarge;
			else
				d[i] = 1.0/(s*s);
			}

		/* Sort new eigenvalues and eigenvectors */
		praxis_SortDV(d, v, (unsigned)n);

		pld.dmin = d[n - 1];
		if (pld.dmin < pld.esmall)
			pld.dmin = pld.esmall;

		illc = (pld.m2*d[0] > pld.dmin);
		}
	}

#if defined(STANDALONE)	/* These routines are ordinarily provided by the "dls" package, but are duplicated here in
						   order to make this file distributable as a standalone version. */

/*----------------------------------------------------------------------------------------------------------------------
|
|	Initializes a matrix to the identity matrix.
*/
PRIVATE void setIdentityMatrix(double **a, unsigned n)
	{
	unsigned i, j;

	for (i = 0; i < n; i++)
		{
		for (j = 0; j < n; j++)
			a[i][j] = 0.0;
		a[i][i] = 1.0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Initializes a vector 'a_arg' of row pointers into a 'nrows' x 'ncols' matrix whose storage is at 'p_arg'.  The size
|	of each element must be passed as 'elementSize'.
|
|	This storage can be defined as matrix[M][N] in the caller or can be any buffer of size at least nrows x ncols.
|	This allows the caller to use a 2-D array declared with constant sizes, and also pass that matrix to routines that
|	accommodate variable-sized matrices by using a vector of row-pointers to access the matrix elements.
*/
PRIVATE void * setMatrixRowPtrs(void *a_arg, void *p_arg, unsigned nrows, unsigned ncols, size_t elementSize)
	{
	char ** a = a_arg;
	char *  p = p_arg;
	size_t rowBytes = ncols*elementSize;

	unsigned i;
	for (i = 0; i < nrows; i++)
		{
		a[i] = p;
		p += rowBytes;
		}

	return a_arg;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Allocates and initializes a 'nrows' x 'ncols' matrix of arbitrary type.
|
|	Storage is allocated for both a vector of row pointers and the matrix itself.  The matrix is contained in a single
|	block, so that its elements may be referenced either as a[i][j] or as (*a)[ncols*i + j].
|
|	The caller should guarantee that 'nrows' is nonzero.  If 'ncols' is zero, the vector of row pointers is allocated,
|	but no storage for the matrix is allocated (the row pointers remain NULL in this case).
|
|	It is assumed that all data pointers are the same size, which is a safe assumption on modern compilers.
*/
PRIVATE void * newMatrix(size_t elementSize, unsigned nrows, unsigned ncols)
	{
	char **a;

	dls_assert(elementSize != 0);
	dls_assert(nrows != 0);

	/* Allocate array of row pointers */
	a = allocateBlock(nrows*sizeof(void *));
	if ((a != NULL) && (ncols != 0))
		{
		/* Allocate block for data and set pointers to each matrix row */
		a[0] = allocateBlock(nrows*ncols*elementSize);
		if (a[0] == NULL)
			deallocateBlock(&a);	/* sets a = NULL */
		else
			setMatrixRowPtrs(a, a[0], nrows, ncols, elementSize);
		}
	return a;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Deallocates memory for a matrix allocated by allocMatrix or allocHalfMatrix.
*/
PRIVATE void deleteMatrix(void *pA)
	{
	void ** a = (void **)pA;
	if (a != NULL)
		{
		deallocateBlock(a[0]);
		deallocateBlock(a);
		}
	}

/* globals for random number generator */
static double	ran1;
static int		ran2;
static double	ran3[127];

/*----------------------------------------------------------------------------------------------------------------------
|
|	Initialization/seed function for Brent's random number generator.
*/
PUBLIC void ranInit(int r)
	{
	if (r < 0)
		r = -r;
	r = (r % 8190) + 1;

	ran2 = 127;
	while (ran2 > 0)
		{
		int i;
		ran2 -= 1;
		ran1 = pow(-2.0, 55.0);
		for (i = 1; i <= 7; i++)
			{
			r = (1756 * r) % 8191;
			ran1 = (ran1 + (r / 32))*(1.0/256.0);
			}
		ran3[ran2] = ran1;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Returns a pseudorandom number on U[0,1] via the multiplicative congruential method.
|
|	This is the random number generator from Brent's book.  It is used for the tests in order to reproduce the
|	results published in his tables 7.2-7.8.
|
|	Because of the globals used, it should *not* be used as is in production code unless there is no possibility of
|	reentry.
*/
PRIVATE double randUni(int *pSeed)
	{
	ran2 = (ran2 == 0) ? 126 : ran2 - 1;
	ran1 += ran3[ran2];
	ran3[ran2] = ran1 = (ran1 < 0) ? ran1 + 0.5 : ran1 - 0.5;
	return ran1 + 0.5;
	}

#endif	/* STANDALONE */

#if defined(CHECK_PROGRESS)

PRIVATE void showPraxisProgress(PraxisData *praxd, PraxisLocal *pld, int flag)
	{
	int i;

	if (pld->nf == 1)
		{
		dls_printf("%7s%7s%2s%15s", "nf", "nl", " ", "f(x)");
		if (g_fmin != DBL_MIN)
			dls_printf("%15s", "delta");
		for (i = 0; i < pld->n; i++)
			{
			char buf[32];
			sprintf(buf, "x_%d", i);
			dls_printf("%15s", buf);
			}
		printf("\n");
		}

	dls_printf("%7d%7d%2c%15g", pld->nf, pld->nl, (flag == 0) ? ' ' : flag, pld->fx);
	if (g_fmin != DBL_MIN)
		dls_printf("%15g", pld->fx - g_fmin);

	for (i = 0; i < pld->n; i++)
		dls_printf("%15f", praxd->x[i]);

	dls_printf("\n");
	}

#endif	/* CHECK_PROGRESS */

#if defined(CHECK_PRAXIS)

static int g_nl;		/* number of function evaluations */

/*====================================================================================================================*/
/*	Test functions                                                                                                    */
/*====================================================================================================================*/

PRIVATE double brentTestFunction(double *px, void *data)
	{
	int i;

	double x = *px;
	double f = 0.0;
	for (i = 1; i <= 20; i++)
		{
		double z = (2*i - 5) / (x - i*i);
		f += z*z;
		}
	g_nl++;
	return f;
	}

PRIVATE double ros(double *x, void * data)
	{
	double y = 1.0 - x[0];
	double z = x[1]-x[0]*x[0];
	return 100.0*(z*z) + y*y;
	}

PRIVATE double helix(double *x, void * data)
	{
	double t, y;

	double r = sqrt(x[0]*x[0] + x[1]*x[1]);
	if (x[0] == 0.0)
		t = 0.25;
	else
		t = atan(x[1]/x[0])/6.283185307179586;		/* 6.28... = 2 x pi */
	if (x[0] < 0.0)
		t += 0.5;

	y = x[2] - 10.0*t;
	return 100.0*(y*y + (r - 1.0)*(r - 1.0)) + x[2]*x[2];
	}

PRIVATE double sing(double *x, void * data)
	{
	double a = x[0] + 10.0*x[1];
	double b = x[2] - x[3];
	double c = x[1] - 2.0*x[2];
	double d = x[0] - x[3];
	return a*a + 5.0*b*b + c*c*c*c + 10.0*d*d*d*d;
	}

typedef struct chebyData
	{
	int      n;
	double * y;
	double * ti;
	double * tminus;
	}
	ChebyData;

PRIVATE double chebyquad(double *x, void *data)
	{
	int		i, j;
	double	f;
	Bool	even;

	ChebyData * cd = data;

	double delta = 0.0;
	for (j = 0; j < cd->n; j++)
		{
		cd->y[j] = 2.0*x[j] - 1.0;
		delta += cd->y[j];
		cd->ti[j] = cd->y[j];
		cd->tminus[j] = 1.0;
		}
	f = delta*delta;
	even = FALSE;
	for (i = 1; i < cd->n; i++)
		{
		double ii = (double)(i + 1);
		even = !even;
		delta = 0.0;
		for (j = 0; j < cd->n; j++)
			{
			double tplus = 2.0*cd->y[j]*cd->ti[j] - cd->tminus[j];
			delta += tplus;
			cd->tminus[j] = cd->ti[j];
			cd->ti[j] = tplus;
			}
		delta /= cd->n;
		if (even)
			delta -= 1.0/(1.0 - ii*ii);
		f += delta*delta;
		}
	return f;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Runs one test of the Praxis optimizer.
*/
PRIVATE int runPraxisTest(
  char *		sfunc,			/* function label for output */
  MinimizeFxn	f,				/* function to optimize */
  void *		data,			/* arbitrary data to be passed to f (can be NULL) */
  unsigned		n,				/* number of parameters to function */
  double *		x,				/* parameter vector of size n; input: starting values; output: optimized values */
  double		h,				/* initial step size (e.g., 1.0) */
  double		f_min,			/* theoretical function minimum */
  int			ranInitSeed,	/* seed to Brent's random number generator */
  Bool			illc,			/* true for ill-conditioned problem, false otherwise */
  double *		x_opt)			/* theoretical optimal parameter values */
	{
	unsigned i;
	double z;

	PraxisData * praxd = newPraxisData(x, n, FALSE);
	if (praxd == NULL)
		return RC_Error;

	ranInit(ranInitSeed);
	g_fmin = f_min;				/* for progress output */

	praxd->illc = illc;
	dls_printf("\n-------------------------\n\n%s (n=%d h=%g):\n\n", sfunc, n, h);
	dls_printf("f(x) evaluated at theoretical optimum = %.13g (should be about %g)\n\n", f(x_opt, data), f_min);

	z = praxis(1e-5, h, n, f, data, praxd);

	dls_printf("\nminimal f(x) = %.13g, final x = (", z);
	for (i = 0; i < n; i++)
		{
		if (i > 0)
			dls_printf(", ");
		dls_printf("%.13g", x[i]);
		}
	dls_printf(")\n\n");

	deletePraxisData(praxd);
	return RC_OK;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Conducts a series of tests based on examples from Brent's book.  Compare the output from this test to the
|	results summarized in his tables 7.2 - 7.8.
*/
PUBLIC int checkPraxis(void)
	{
	int i, n;
	ChebyData chebyData;

#	define N_MAX 8
	double x[N_MAX], x_opt[N_MAX], y[N_MAX], ti[N_MAX], tminus[N_MAX];

	/* Sanity check on random-number generator (doesn't evaluate quality; just makes sure nothing is blatantly wrong) */
	double sum = 0.0;
	ranInit(4);
	n = 100000;
	for (i = 0; i < n; i++)
		sum += randUni(NULL);
	dls_printf("\n-------------------------\n");
	printf("\nmean of %d random numbers = %f (should be about 0.5)\n", n, sum/n);

	/* Rosenbrock test */
	x[0] = -1.2; x_opt[0] = 1.0;
	x[1] =  1.0; x_opt[1] = 1.0;
	runPraxisTest("Rosenbrock", ros, NULL, 2, x, 1.0, 0.0, 4, FALSE, x_opt);

	/* Helix test */
	x[0] = -1.0; x_opt[0] = 1.0;
	x[1] =  0.0; x_opt[1] = 0.0;
	x[2] =  0.0; x_opt[2] = 0.0;
	runPraxisTest("Helix", helix, NULL, 3, x, 1.0, 0.0, 4, FALSE, x_opt);

	/* Singular test */
	x[0] =  3.0; x_opt[0] = 0.0;
	x[1] = -1.0; x_opt[1] = 0.0;
	x[2] =  0.0; x_opt[2] = 0.0;
	x[3] =  1.0; x_opt[3] = 0.0;
	runPraxisTest("Singular", sing, NULL, 4, x, 1.0, 0.0, 2, TRUE, x_opt);

	/* Chebyquad tests... */
	chebyData.y = y;
	chebyData.ti = ti;
	chebyData.tminus = tminus;

	/* Chebyquad N = 2 */
	x_opt[0] = 0.2113248654057;
	x_opt[1] = 0.7886751345943;
	n = chebyData.n = 2;
	for (i = 0; i < n; i++)
		x[i] = (double)(i + 1)/(n + 1);
	runPraxisTest("Chebyquad_2", chebyquad, &chebyData, n, x, 0.1, 5.611343057188e-24, 4, FALSE, x_opt);

	/* Chebyquad N = 4 */
	x_opt[0] = 0.1026727662050;
	x_opt[1] = 0.4062037685183;
	x_opt[2] = 0.5937962317961;
	x_opt[3] = 0.8973272336247;
	n = chebyData.n = 4;
	for (i = 0; i < n; i++)
		x[i] = (double)(i + 1)/(n + 1);
	runPraxisTest("Chebyquad_4", chebyquad, &chebyData, n, x, 0.1, 1.879737941753e-16, 4, FALSE, x_opt);

	/* Chebyquad N = 6 */
	x_opt[0] = 0.0668765861341;
	x_opt[1] = 0.2887406551984;
	x_opt[2] = 0.3666822994979;
	x_opt[3] = 0.6333176827606;
	x_opt[4] = 0.7112593410710;
	x_opt[5] = 0.9331234076183;
	n = chebyData.n = 6;
	for (i = 0; i < n; i++)
		x[i] = (double)(i + 1)/(n + 1);
	runPraxisTest("Chebyquad_6", chebyquad, &chebyData, n, x, 0.1, 3.767185586412e-15, 4, FALSE, x_opt);

	/* Chebyquad N = 8 */
	x_opt[0] = 0.04315264838751;
	x_opt[1] = 0.1930908158017;
	x_opt[2] = 0.2663286132896;
	x_opt[3] = 0.5;
	x_opt[4] = 0.5;
	x_opt[5] = 0.7336707838024;
	x_opt[6] = 0.806909849489;
	x_opt[7] = 0.956847623612;
	n = chebyData.n = 8;
	for (i = 0; i < n; i++)
		x[i] = (double)(i + 1)/(n + 1);
	runPraxisTest("Chebyquad_6", chebyquad, &chebyData, n, x, 0.1, 0.003516873728788, 4, FALSE, x_opt);

	return RC_OK;
	}
/*----------------------------------------------------------------------------------------------------------------------
|
|	Conducts a test corresponding to the example on page 77 of Brent's book.
*/
PRIVATE void checkLocalmin(void)
	{
	int i;
	double eps = pow(16.0, -7.0);

	dls_printf("\n-------------------------\n\n%s:\n", "Localmin test");

	dls_printf("\n%5s%15s%15s%8s\n", "i", "xmin", "f(xmin)", "n_l");
	for (i = 1; i <= 19; i++)
		{
		double fx, x;
		double a = i*i;
		double b = (i + 1)*(i + 1);
		g_nl = 0;
		fx = localMin(a, b, eps, 1e-10, brentTestFunction, NULL, &x);
		printf("%5d%15.7f%15.10f%8d\n", i, x, fx, g_nl);
		}
	}

#endif	/* defined(CHECK_PRAXIS) */

#if defined(TEST_BRENT)

/*----------------------------------------------------------------------------------------------------------------------
|
|	Runs tests corresponding to examples from Brent's book.
|
|	The tests can be be performed as follows:
|
|		cc -o testpraxis -O -DTEST_BRENT dls_brent.c
|		./testpraxis
*/
int main(void)
	{
	checkLocalmin();
	checkPraxis();
	return 0;
	}

#endif	/* defined(STANDALONE) */

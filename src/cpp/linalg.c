/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas and the PhyPy library: Python software for phylogenetic analysis    |
|  Copyright (C) 2006 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
|                                                                             |
|  This program is free software; you can redistribute it and/or modify       |
|  it under the terms of the GNU General Public License as published by       |
|  the Free Software Foundation; either version 2 of the License, or          |
|  (at your option) any later version.                                        |
|                                                                             |
|  This program is distributed in the hope that it will be useful,            |
|  but WITHOUT ANY WARRANTY; without even the implied warranty of             |
|  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              |
|  GNU General Public License for more details.                               |
|                                                                             |
|  You should have received a copy of the GNU General Public License along    |
|  with this program; if not, write to the Free Software Foundation, Inc.,    |
|  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.                |
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*	linalg.c
|
|	This file contains code for various linear algebra operations used in distance and maximum likelihood analyses.
|
|	Copyright (c) 2003 by David L. Swofford, Florida State University.
|	All rights reserved.
*/
#include "dls.h"

#if 1 /* POL modifications needed for using linalg.h and linalg.c with PHYCAS */

/* PAUP is defined in dls.h */ 
#	undef PAUP

/* In PAUP*, these are part of an enum defined in paup.h */ 
#	define RC_OK	0
#	define RC_Error	1

#endif

#if defined(PAUP)
#	define DOES_MATH
#	include "paup.h"
#else
#	include <stdio.h>
#	include <math.h>
#	include <float.h>
#endif
#include "linalg.h"

#if defined(PAUP)
#	undef REORDER_EIG					/* don't bother reordering eigenvalues/vectors */
#	if DEVELOPMENT
#		undef DEBUG_SHOW_EIGVALVECT		/* show computed eigenvalues and eigenvectors */
#	endif
#endif

#define TINY 1e-20

#if MAC_VERSION && defined(__POWERPC__)
#	define Sign(a, b)	copysign((a), (b))		/* use fast intrinsic */
#endif

#if !defined(MAX)
#	define MAX(a,b) (((a) > (b)) ? (a) : (b))
#endif
#if !defined(MIN)
#	define MIN(a,b) (((a) < (b)) ? (a) : (b))
#endif

#if !defined(ANIMATE_CURSOR)
#	define ANIMATE_CURSOR
#endif

#if defined(USE_LAPACK)
#	include "lapack.h"
	int dsyevr_(char *jobz, char *range, char *uplo, int *n, double *a, int *lda, double *vl, double *vu, int *il,
	            int *iu, double *abstol, int *m, double *w, double *z, int *ldz, int *isuppz, double *work, int *lwork,
	            int *iwork, int *liwork, integer *info);
#endif

/*----------------------------------------------------------------------------------------------------------------------
|
|	LUDecompose
|
|	Replaces matrix 'a' with its LU-decomposition, obtained using Crout's algorithm as described in Press et al.
|	(1992).  Returns RC_Error if matrix is singular, RC_OK otherwise.
*/

ENTRY int LUDecompose(
  double	**a,			/* the matrix whose LU-decomposition is wanted */
  int		n,				/* order of a */
  double	*scaling,		/* work vector of size n (stores implicit scaling of each row) */
  int		*permutation,	/* => row permutation according to partial pivoting sequence */
  double	*pd)			/* => 1 if number of row interchanges was even, -1 if odd (NULL OK) */
	{
	int			i, imax = 0, j, k;
	double		largest, sum, temp, pivotVal,
				d = 1.0;
	
	/* Get implicit scaling information */
	for (i = 0; i < n; i++)
		{
		ANIMATE_CURSOR
		largest = 0.0;
		for (j = 0; j < n; j++)
			{
			if ((temp = fabs(a[i][j])) > largest)
				largest = temp;
			}
		if (largest == 0.0)	/* matrix is all-zero */
			return RC_Error;

		scaling[i] = 1.0 / largest;
		}
		
	/* Loop over columns in Crout's method */
	for (j = 0; j < n; j++)
		{
		ANIMATE_CURSOR
		
		/* Solve for beta(i,j)'s */
		for (i = 0; i < j; i++)
			{
			sum = a[i][j];
			for (k = 0; k < i; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			}
		
		/* Solve for alpha(i,j)'s and find pivot */
		largest = 0.0;
		for (i = j; i < n; i++)
			{
			sum = a[i][j];
			for (k = 0; k < j; k++)
				sum -= a[i][k] * a[k][j];
			a[i][j] = sum;
			temp = scaling[i] * fabs(sum);
			if (temp >= largest)
				{
				/* Keep track of largest pivot element */
				largest = temp;
				imax = i;
				}
			}
		if (j != imax)
			{
			/* Exchange rows */
			for (k = 0; k < n; k++)
				{
				temp = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = temp;
				}	
			d = -d;
			scaling[imax] = scaling[j];
			}

		/* Divide through by pivot element */
		permutation[j] = imax;
		pivotVal = a[j][j];
		if (a[j][j] == 0.0)
			pivotVal = a[j][j] = TINY;			/* avoid division by zero */
		if (j != n - 1)
			{
			for (i = j + 1; i < n; i++)
				a[i][j] /= pivotVal;
			}
		}

	if (pd != NULL)
		*pd = d;

	return RC_OK;
	}

#if defined(NEED_MATINV)	/* PAUP doesn't currently need to do matrix inversion */

/*----------------------------------------------------------------------------------------------------------------------
|
|	LUBackSubst
|
|	Performs back-substition into LU-decomposed matrix in order to obtain inverse.
*/

LOCAL void LUBackSubst(double **a, int n, int *permutation, double *b)
	{
	int			i, ip, j,
				ii = -1;
	double		sum;

	for (i = 0; i < n; i++)
		{
		ip = permutation[i];
		sum = b[ip];
		b[ip] = b[i];
		if (ii >= 0)
			{
			for (j = ii; j <= i - 1; j++)
				sum -= a[i][j] * b[j];
			}
		else if (sum != 0.0)
			ii = i;
		b[i] = sum;
		}
	for (i = n - 1; i >= 0; i--)
		{
		sum = b[i];
		for (j = i + 1; j < n; j++)
			sum -= a[i][j] * b[j];
		b[i] = sum / a[i][i];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	InvertMatrix
|
|	Inverts matrix 'a' using LU-decomposition technique, storing inverse in 'a_inv'.  Matrix 'a' is destroyed.  Returns
|	RC_Error if matrix is singular, RC_OK otherwise.
*/

ENTRY int InvertMatrix(
  double	**a,			/* matrix represented as vector of row pointers */
  int		n,				/* order of matrix */
  double	*col,			/* work vector of size n */
  int		*permutation,	/* work vector of size n */
  double	**a_inv)		/* inverse of input matrix a (matrix a is destroyed) */
	{
	int			rc, i, j;
	
	rc = LUDecompose(a, n, col, permutation, NULL);
	if (rc == FALSE)
		{
		for (j = 0; j < n; j++)
			{
			for (i = 0; i < n; i++)
				col[i] = 0.0;
			col[j] = 1.0;
			LUBackSubst(a, n, permutation, col);
			for (i = 0; i < n; i++)
				a_inv[i][j] = col[i];
			}
		}
	return rc;
	}

#endif

#if defined(USE_LAPACK)

/*----------------------------------------------------------------------------------------------------------------------
|
|	EigenRealSymmetric
|
|	This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK) to
|	find the eigenvalues and eigenvectors of a real symmetric matrix.  It was converted from Fortran to C by David
|	Swofford (was originally rs.f in EISPACK).
|
|	ON INPUT:
|
|	   n  is the order of the matrix  a.
|
|	   a  contains the real symmetric matrix.
|
|	   fv  is a temporary storage array of at least n elements.
|
|	ON OUTPUT:
|
|	   w  contains the eigenvalues in ascending order.
|
|	   z  contains the eigenvectors.
|
|	   Return value is an integer output variable set equal to an error completion code described
|	   in the documentation for tql2.  The normal completion code is zero.
*/

ENTRY int EigenRealSymmetric(int n, double **a, double *w, double **z, double *fv)
	{
#	if defined(HAVE_PRAGMA_UNUSED)
#		pragma unused(fv)
#	endif

#	if !defined(LAPACK_V3)		/* This is conventional LAPACK */

		int			info, lwork;
		unsigned	i, j;
		double		work[12];	//3*n for now; need to do better

		lwork = 12;		
		(void)dsyev_("V", "L", &n, (double *)a[0], &n, w, work, &lwork, &info);
		//temp?  (could deal with this via subscript order in the caller)
		/* Copy eigenvectors before returning (have to transpose them with LAPACK */
		for (i = 0; i < n; i++)
			{
			for (j = 0; j < n; j++)
				z[i][j] = a[j][i];
			}

#	else	/* This is new dqds algorithm in LAPACK v3 ("relatively robust representations") */

		int			info, il, iu, m, lwork, liwork, isuppz[2*MAX_S], iwork[10*MAX_S];
		unsigned	i, j;
		double		vl, vu, work[26*MAX_S],
					zero = 0.0;

		lwork = 26*MAX_S;
		liwork = 10*MAX_S;
		(void)dsyevr_("V", "A", "L", &n, (double *)a[0], &n, &vl, &vu, &il, &iu, &zero, &m, w, (double *)z[0], &n, 
		              isuppz, work, &lwork, iwork, &liwork, &info);

		//transpose (TO DO (LAPACK_V3): Use BLAS?)
		for (i = 0; i < n; i++)
			{
			for (j = 0; j < i; j++)
				{
				double temp = z[i][j];
				z[i][j] = z[j][i];
				z[j][i] = temp;
				}
			}
#	endif
#	if defined(DEBUG_SHOW_EIGVALVECT)
		{
		int		i;

		Pprintf("eigenvalues=%g %g %g %g\n", w[0], w[1], w[2], w[3]);
		Pprintf("eigenvectors =\n");
		for (i = 0; i < n; i++)
			Pprintf("%g %g %g %g\n", z[i][0], z[i][1], z[i][2], z[i][3]);
		}
#	endif

	return info;
	}

#else	/* !defined(USE_LAPACK) */

#if !defined(Sign)

/*----------------------------------------------------------------------------------------------------------------------
|
|	Sign
|
|	Standard "sign" function--copies sign of y to x.
*/

LOCAL double Sign(double a, double b)
	{
	double		x;

	x = (a >= 0 ? a : -a);
	return (b >= 0 ? x : -x);
	}

#endif	/* !defined(Sign) */

/*----------------------------------------------------------------------------------------------------------------------
|
|	Pythag
|
|	Finds sqrt(a**2+b**2) without overflow or destructive underflow.
*/

LOCAL double Pythag(double a, double b)
	{
	double	p, r, s, t, u;
	
	p = fabs(a);
	if (fabs(b) > p)
		p = fabs(b);
	if (p == 0.0)
		return 0.0;

	r = fabs(a);
	if (fabs(b) < r)
		r = fabs(b);
	r /= p;
	r *= r;
	
	for (FOREVER)
		{
		t = 4.0 + r;
		if (t == 4.0)
			break;
		s = r/t;
		u = 1.0 + 2.0*s;
		p = u*p;
		s /= u;
		r = s*s * r;
		}

	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Tred2
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure TRED2, Num. Math. 11, 181-195(1968) by Martin, Reinsch, and
|	Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 212-226 (1971).
|
|	This subroutine reduces a real symmetric matrix to a symmetric tridiagonal matrix using and accumulating orthogonal
|	similarity transformations.
|
|	ON INPUT:
|
|		n is the order of the matrix.
|
|		a contains the real symmetric input matrix.  only the lower triangle of the matrix need be supplied.
|
|	ON OUTPUT:
|
|		d contains the diagonal elements of the tridiagonal matrix.
|
|		e contains the subdiagonal elements of the tridiagonal matrix in its last n-1 positions.  e(1) is set to zero.
|
|		z contains the orthogonal transformation matrix produced in the reduction.
|
|		a and z may coincide.  if distinct, a is unaltered.
*/

LOCAL void Tred2(int n, double **a, double *d, double *e, double **z)
	{
	int			i, j, k, l;
	double		f, g, h, hh, scale;

	for (i = 0; i < n; i++)
		{
		for (j = i; j < n; j++)
			z[j][i] = a[j][i];
		d[i] = a[n-1][i];
		}

	for (i = n - 1; i > 0; i--)
		{
		l = i - 1;
		h = scale = 0.0;
		if (l > 0)
			{
			/* Scale row (Algol tol then not needed) */
			for (k = 0; k < i; k++)
				scale += fabs(d[k]);
			}
		if (scale == 0.0)
			{
			e[i] = d[l];
			for (j = 0; j < i; j++)
				{
				d[j] = z[l][j];
				z[i][j] = z[j][i] = 0.0;
				}
			}
		else
			{
  			for (k = 0; k < i; k++)
  				{
				d[k] = d[k] / scale;
				h += d[k] * d[k];
				}

			f = d[l];
			g = -Sign(sqrt(h), f);
			e[i] = scale * g;
			h -= f * g;
			d[l] = f - g;

			/* Form a*u */
			for (j = 0; j < i; j++)
				e[j] = 0.0;
			for (j = 0; j < i; j++)
				{
				f = d[j];
				z[j][i] = f;
				g = e[j] + z[j][j] * f;
				if (l > j)
					{
					for (k = j + 1; k < i; k++)
						{
						g += z[k][j] * d[k];
						e[k] += z[k][j] * f;
						}
					}
				e[j] = g;
				}
				
			/* Form p */
			f = 0.0;
			for (j = 0; j < i; j++)
				{
				e[j] /= h;
				f += e[j] * d[j];
				}

			/* Form q */
			hh = f / (h + h);
			for (j = 0; j < i; j++)
				e[j] -= hh * d[j];

			/* Form reduced a */
			for (j = 0; j < i; j++)
				{
				f = d[j];
				g = e[j];
				for (k = j; k < i; k++)
					z[k][j] -= (f * e[k] + g * d[k]);
				d[j] = z[l][j];
				z[i][j] = 0.0;
				}
			}
		d[i] = h;
		}

	/* Accumulation of transformation matrices */
	for (i = 1; i < n; i++)
		{
		l = i - 1;
		z[n-1][l] = z[l][l];
		z[l][l] = 1.0;
		h = d[i];
		if (h != 0.0)
			{
			for (k = 0; k <= l; k++)
				d[k] = z[k][i] / h;
			for (j = 0; j <= l; j++)
				{
				g = 0.0;
				for (k = 0; k <= l; k++)
					g += z[k][i] * z[k][j];
				for (k = 0; k <= l; k++)
					z[k][j] -= g * d[k];
				}
			}
		for (k = 0; k <= l; k++)
			z[k][i] = 0.0;
		}

	for (i = 0; i < n; i++)
		{
		d[i] = z[n-1][i];
		z[n-1][i] = 0.0;
		}

	z[n-1][n-1] = 1.0;
	e[0] = 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Tql2
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure TQL2, Num. Math. 11, 293-306(1968) by Bbowdler, Martin,
|	Reinsch, and Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 227-240(1971).
|
|	This subroutine finds the eigenvalues and eigenvectors of a symmetric tridiagonal matrix by the QL method.  The
|	eigenvectors of a full symmetric matrix can also be found if Tred2 has been used to reduce this full matrix to
|	tridiagonal form.
|
|	ON INPUT:
|
|		n is the order of the matrix.
|
|		d contains the diagonal elements of the input matrix.
|
|		e contains the subdiagonal elements of the input matrix in its last n-1 positions.  e(1) is arbitrary.
|
|		z contains the transformation matrix produced in the reduction by  tred2, if performed. If the eigenvectors of
|		the tridiagonal matrix are desired, z must contain the identity matrix.
|
|	ON OUTPUT:
|
|		d contains the eigenvalues in ascending order.  If an error exit is made, the eigenvalues are correct but
|		unordered for indices 1,2,...,ierr-1.
|
|		e has been destroyed.
|
|		z contains orthonormal eigenvectors of the symmetric tridiagonal (or full) matrix.  If an error exit is made,
|		z contains the eigenvectors associated with the stored eigenvalues.
|
|		ierr is set to
|			zero       for normal return,
|			j          if the j-th eigenvalue has not been determined after 30 iterations.
|
|	Calls Pythag for sqrt(a*a + b*b) .
*/

LOCAL int Tql2(int n, double *d, double *e, double **z)
	{
	int		i, j, k, l, m, l1, l2;
	double	c, c2, c3=0, dl1, el1, f, g, h, p, r, s, s2=0, tst1, tst2;

	if (n == 1)
		return 0;

	for (i = 1; i < n; i++)
		e[i-1] = e[i];

	f = 0.0;
	tst1 = 0.0;
	e[n-1] = 0.0;

	for (l = 0; l < n; l++)
		{
		j = 0;
		h = fabs(d[l]) + fabs(e[l]);
		if (tst1 < h)
			tst1 = h;
		/* Look for small sub-diagonal element */
		for (m = l; m < n; m++)
			{
			tst2 = tst1 + fabs(e[m]);
			if (tst2 == tst1)
				break;
			/* e[n-1] is always zero, so there is no exit through the bottom of the loop */
			}

		if (m != l)
			{
			do	{
				if (j == 30)
					{
					/* Set error -- no convergence to an eigenvalue after 30 iterations */
					return l;
					}
				j = j + 1;

				/* Form shift */
				l1 = l + 1;
				l2 = l1 + 1;
				g = d[l];
				p = (d[l1] - g) / (2.0 * e[l]);
				r = Pythag(p, 1.0);
				d[l]  = e[l] / (p + Sign(r, p));
				d[l1] = e[l] * (p + Sign(r, p));
				dl1 = d[l1];
				h = g - d[l];
				for (i = l2; i < n; i++)
					d[i] -= h;
				f += h;

				/* QL transformation */
				p = d[m];
				c = 1.0;
				c2 = c;
				el1 = e[l1];
				s = 0.0;
				for (i = m - 1; i >= l; i--)
					{
					c3 = c2;
					c2 = c;
					s2 = s;
					g = c * e[i];
					h = c * p;
					r = Pythag(p, e[i]);
					e[i+1] = s * r;
					s = e[i] / r;
					c = p / r;
					p = c * d[i] - s * g;
					d[i+1] = h + s * (c * g + s * d[i]);
					/* Form vector */
					for (k = 0; k < n; k++)
						{
						h = z[k][i+1];
						z[k][i+1] = s * z[k][i] + c * h;
						z[k][i] = c * z[k][i] - s * h;
						}
					}

				p = -s * s2 * c3 * el1 * e[l] / dl1;
				e[l] = s * p;
				d[l] = c * p;
				tst2 = tst1 + fabs(e[l]);
				}
				while (tst2 > tst1);
			}
		d[l] += f;
		}

#	if defined(REORDER_EIG)		/* no need to do this in PAUP */

	/* Order eigenvalues and eigenvectors */
	for (l = 1; l < n; l++)
		{
		i = l - 1;
		k = i;
		p = d[i];
		for (j = l; j < n; j++)
			{
			if (d[j] < p)
				{
				k = j;
				p = d[j];
				}
			}
		if (k != i)
			{
			d[k] = d[i];
			d[i] = p;
			for (j = 0; j < n; j++)
				{
				p = z[j][i];
				z[j][i] = z[j][k];
				z[j][k] = p;
				}
			}
		}
#	endif

	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	EigenRealSymmetric
|
|	This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK) to
|	find the eigenvalues and eigenvectors of a real symmetric matrix.  It was converted from Fortran to C by David
|	Swofford (was originally rs.f in EISPACK).
|
|	ON INPUT:
|
|	   n  is the order of the matrix  a.
|
|	   a  contains the real symmetric matrix.
|
|	   fv  is a temporary storage array of at least n elements.
|
|	ON OUTPUT:
|
|	   w  contains the eigenvalues in ascending order.
|
|	   z  contains the eigenvectors.
|
|
|	   Return value is an integer output variable set equal to an error completion code described in the documentation
|	   for tql2.  The normal completion code is zero.
*/

ENTRY int EigenRealSymmetric(int n, double **a, double *w, double **z, double *fv)
	{
	int rc;

	Tred2(n, a, w, fv, z);
	rc = Tql2(n, w, fv, z);

#	if defined(DEBUG_SHOW_EIGVALVECT)
		{
		int i, j;

		Pprintf("eigenvalues =");
		for (i = 0; i < n; i++)
			Pprintf(" %g", w[i]);
		Pprintf("\n\neigenvectors =\n");
		for (i = 0; i < n; i++)
			for (j = 0; j < n; j++)
				Pprintf(" %g", z[i][j]);
		Pputc('\n');
		}
#	endif

	return rc;
	}

#endif	/* !defined(USE_LAPACK) */

#if defined(NEED_NONSYM_EIG)	/* nonsymmetric eigenvalue routines are not currently needed in PAUP */

/*----------------------------------------------------------------------------------------------------------------------
|
|	Exchange
|
|	Support function for EISPACK routine Balanc.
*/

LOCAL void Exchange(int j, int k, int l, int m, int n, double **a, double *scale)
	{
	int		i;
	double	f;

	scale[m] = (double)j;
	if (j != m)
		{
		for (i = 0; i <= l; i++)
			{
			f = a[i][j];
			a[i][j] = a[i][m];
			a[i][m] = f;
			}	
		for (i = k; i < n; i++)
			{
			f = a[j][i];
			a[j][i] = a[m][i];
			a[m][i] = f;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Balanc
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure BALANCE, Num. Math. 13, 293-304(1969) by Parlett and
|	Reinsch. Handbook for Auto. Comp., Vol. II-Linear Algebra, 315-326( 1971).
|
|	This subroutine balances a real matrix and isolates eigenvalues whenever possible.
|
|	ON INPUT:
|
|	   n is the order of the matrix.
|
|	   a contains the input matrix to be balanced.
|
|	ON OUTPUT:
|
|	   a contains the balanced matrix.
|
|	   low and high are two integers such that a(i,j) is equal to zero if
|	      (1) i is greater than j and
|	      (2) j=1,...,low-1 or i=high+1,...,n.
|
|	   scale contains information determining the permutations and scaling factors used.
|
|	Suppose that the principal submatrix in rows low through high has been balanced, that p(j) denotes the index
|	interchanged with j during the permutation step, and that the elements of the diagonal matrix used are denoted by
|	d(i,j).  Then
|	   scale(j) = p(j),    for j = 1,...,low-1
|	            = d(j,j),      j = low,...,high
|	            = p(j)         j = high+1,...,n.
|	The order in which the interchanges are made is n to high+1,  then 1 to low-1.
|
|	Note that 1 is returned for high if high is zero formally.
*/

LOCAL void Balanc(int n, double **a, int *pLow, int *pHigh, double *scale)
	{
	double		c, f, g, r, s, b2;
	int			i, j, k, l, m, noconv;

	b2 = FLT_RADIX * FLT_RADIX;
	k = 0;
	l = n - 1;
   	
	/* Search for rows isolating an eigenvalue and push them down */
	for (j = l; j >= 0; j--)
		{
		for (i = 0; i <= l; i++)
			{
			if (i != j)
				{
				if (a[j][i] != 0.0)
					goto next_j1;
				}
			}
		m = l;		
		Exchange(j, k, l, m, n, a, scale);
		if (--l < 0)
			goto leave;
		
		next_j1:
			;
		}

	/* Search for columns isolating an eigenvalue and push them left */
	for (j = k; j <= l; j++)
		{
		for (i = k; i <= l; i++)
			{
			if (i != j)
				{
				if (a[i][j] != 0.0)
					goto next_j2;
				}
			}

		m = k;
		Exchange(j, k, l, m, n, a, scale);
		k++;
		next_j2:
			;
		}

	/*------ Now balance the submatrix in rows k to l ------*/
	for (i = k; i <= l; i++)
		scale[i] = 1.0;

	/* Iterative loop for norm reduction */
	do	{
		noconv = FALSE;
	
		for (i = k; i <= l; i++)
			{
			c = 0.0;
			r = 0.0;
		
			for (j = k; j <= l; j++)
				{
				if (j != i)
					{
					c += fabs(a[j][i]);
					r += fabs(a[i][j]);
					}
				}
			/* Guard against zero c or r due to underflow */
			if ((c != 0.0) && (r != 0.0))
				{
				g = r / FLT_RADIX;
				f = 1.0;
				s = c + r;
					while (c < g)
					{
					f *= FLT_RADIX;
					c *= b2;
					}
				g = r * FLT_RADIX;
				while (c >= g)
					{
					f /= FLT_RADIX;
					c /= b2;
					}
		
				/* Now balance */
				if ((c + r) / f < s * .95)
					{
					g = 1. / f;
					scale[i] *= f;
					noconv = TRUE;				
					for (j = k; j < n; j++)
						a[i][j] *= g;
					for (j = 0; j <= l; j++)
						a[j][i] *= f;
					}
				}
			}	
		}
		while (noconv);

	leave:
		*pLow = k;
		*pHigh = l;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	ElmHes
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure ELMHES, Num. Math. 12, 349-368(1968) by Martin and
|	Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 339-358 (1971).
|
|	Given a real general matrix, this subroutine reduces a submatrix situated in rows and columns low through high to
|	upper Hessenberg form by stabilized elementary similarity transformations.
|
|	ON INPUT:
|
|		n is the order of the matrix.
|
|		low and high are integers determined by the balancing subroutine BALANC.  If BALANC has not been used, set
|		low=1, high=n.
|
|		a contains the input matrix.
|
|	ON OUTPUT:
|
|		a contains the Hessenberg matrix.  The multipliers which were used in the reduction are stored in the remaining
|		triangle under the Hessenberg matrix.
|
|		intchg contains information on the rows and columns interchanged in the reduction.  Only elements low through
|		high are used.
*/

LOCAL void ElmHes(int n, int low, int high, double **a, int *intchg)
	{
	int			i, j, m;
	double		x, y;
	int			mm1;
	
	for (m = low + 1; m < high; m++)
		{
		mm1 = m - 1;
		x = 0.0;
		i = m;
	
		for (j = m; j <= high; j++)
			{
			if (fabs(a[j][mm1]) > fabs(x))
				{
				x = a[j][mm1];
				i = j;
				}
			}
	
		intchg[m] = i;
		if (i != m)
			{
			/* Interchange rows and columns of a */
			for (j = mm1; j < n; j++)
				{
				y = a[i][j];
				a[i][j] = a[m][j];
				a[m][j] = y;
				}
			for (j = 0; j <= high; j++)
				{
				y = a[j][i];
				a[j][i] = a[j][m];
				a[j][m] = y;
				}
			}

		if (x != 0.0)
			{
			for (i = m + 1; i <= high; i++)
				{
				y = a[i][mm1];
				if (y != 0.0)
					{
					y /= x;
					a[i][mm1] = y;
					for (j = m; j < n; j++)
						a[i][j] -= y * a[m][j];
					for (j = 0; j <= high; j++)
						a[j][m] += y * a[j][i];
					}
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	ElTran
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure ELMTRANS,  Num. Math. 16, 181-204 (1970) by Peters and
|	Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 372-395 (1971).
|
|	This subroutine accumulates the stabilized elementary similarity transformations used in the reduction of a real
|	general matrix to upper Hessenberg form by  ElmHes.
|
|	ON INPUT:
|
|		n is the order of the matrix.
|
|		low and high are integers determined by the balancing subroutine Balanc.  if  Balanc has not been used, set
|		low=1, high=n.
|
|		a contains the multipliers which were used in the reduction by ElmHes in its lower triangle below the
|		subdiagonal.
|
|		intchg contains information on the rows and columns interchanged in the reduction by ElmHes.  Only elements low
|		through high are used.
|
|	ON OUTPUT:
|
|	   z contains the transformation matrix produced in the reduction by ElmHes.
*/

LOCAL void ElTran(int n, int low, int high, double **a, int *intchg, double **z)
	{
	int			i, j, mp;

	/* Initialize z to identity matrix */
	for (j = 0; j < n; j++)
		{
		for (i = 0; i < n; i++)
			z[i][j] = 0.0;
		z[j][j] = 1.0;
		}

	for (mp = high - 1; mp >= low + 1; mp--)
		{
		for (i = mp + 1; i <= high; i++)
			z[i][mp] = a[i][mp-1];
	
		i = intchg[mp];
		if (i != mp) 
			{
			for (j = mp; j <= high; j++)
				{
				z[mp][j] = z[i][j];
				z[i][j] = 0.0;
				}
			z[i][mp] = 1.0;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	CDiv
|
|	Complex division, (cr,ci) = (ar,ai)/(br,bi)
*/

LOCAL void CDiv(double ar, double ai, double br, double bi, double *cr, double *ci)
	{
	double		s, ais, bis, ars, brs;

	s = fabs(br) + fabs(bi);
	ars = ar / s;
	ais = ai / s;
	brs = br / s;
	bis = bi / s;
	s = brs*brs + bis*bis;
	*cr = (ars*brs + ais*bis) / s;
	*ci = (ais*brs - ars*bis) / s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	Hqr2
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure HQR2, Num. Math. 16, 181-204 (1970) by Peters and
|	Wilkinson.  Handbook for Auto. Comp., Vol. II-Linear Algebra, 372-395 (1971).
|
|	This subroutine finds the eigenvalues and eigenvectors of a real upper Hessenberg matrix by the QR method.  The
|	eigenvectors of a real general matrix can also be found if ElmHes and ElTran or OrtHes  and  OrTran  have been
|	used to reduce this general matrix to Hessenberg form and to accumulate the similarity transformations.
|
|	ON INPUT:
|
|		n is the order of the matrix
|
|		low and high are integers determined by the balancing subroutine Balanc.  If Balanc has not been used, set
|		low=0, high=n-1.
|
|		h contains the upper Hessenberg matrix
|
|		z contains the transformation matrix produced by ElTran after the reduction by ElmHes, or by OrTran after the
|		reduction by OrtHes, if performed.  If the eigenvectors of the Hessenberg matrix are desired, z must contain
|		the identity matrix.
|
|	ON OUTPUT:
|
|		h has been destroyed
|
|		wr and wi contain the real and imaginary parts, respectively, of the eigenvalues.  The eigenvalues are
|		unordered except that complex conjugate pairs of values appear consecutively with the eigenvalue having the
|		positive imaginary part first.  If an error exit is made, the eigenvalues should be correct for indices
|		ierr,...,n-1.
|
|		z contains the real and imaginary parts of the eigenvectors.   If the i-th eigenvalue is real, the i-th column
|		of z contains its eigenvector.  If the i-th eigenvalue is complex with positive imaginary part, the i-th and
|		(i+1)-th columns of z contain the real and imaginary parts of its eigenvector.  The eigenvectors are
|		unnormalized.  If an error exit is made, none of the eigenvectors has been found. 
|
|		Return value is set to:
|			zero	for normal return,
|			j		if the limit of 30*n iterations is exhausted while the j-th eigenvalue is being sought.
|
|	Calls CDiv for complex division.
*/

LOCAL int Hqr2(int n, int low, int high, double **h, double *wr, double *wi, double **z)
	{
	int			i, j, k, l, m, na, en, notlas, mp2, itn, its, enm2, twoRoots;
	double		norm, p, q, r, s, t, w, x, y, ra, sa, vi, vr, zz, tst1, tst2;

	/* Store roots isolated by Balanc and compute matrix norm */
	norm = 0.0;
	k = 0;
	for (i = 0; i < n; i++)
		{
		for (j = k; j < n; j++)
			norm += fabs(h[i][j]);

		k = i;
		if ((i < low) || (i > high))
			{
			wr[i] = h[i][i];
			wi[i] = 0.0;
			}
		}

	en = high;
	t = 0.0;
	itn = n * 30;

	/* Search for next eigenvalues */

	while (en >= low)
		{
		its = 0;
		na = en - 1;
		enm2 = na - 1;
		twoRoots = FALSE;

		/* Look for single small sub-diagonal element */
		for (FOREVER)
			{
			for (l = en; l > low; l--)
				{
				s = fabs(h[l-1][l-1]) + fabs(h[l][l]);
				if (s == 0.0)
					s = norm;
				tst1 = s;
				tst2 = tst1 + fabs(h[l][l-1]);
				if (tst2 == tst1)
					break;
				}
	
			/* Form shift */
			x = h[en][en];
			if (l == en)
				break;
			y = h[na][na];
			w = h[en][na] * h[na][en];
			if (l == na)
				{
				twoRoots = TRUE;
				break;
				}

			if (itn == 0)
				{
				/* Set error -- all eigenvalues have not converged after 30*n iterations */
				return en;
				}
			if ((its == 10) || (its == 20))
				{
				/* Form exceptional shift */
				t += x;
			
				for (i = low; i <= en; i++)
					h[i][i] -= x;
			
				s = fabs(h[en][na]) + fabs(h[na][enm2]);
				x = s * 0.75;
				y = x;
				w = s * -0.4375 * s;
				}
	
			its++;
			--itn;
	
			/* Look for two consecutive small sub-diagonal elements */
			for (m = enm2; m >= l; m--)
				{
				zz = h[m][m];
				r = x - zz;
				s = y - zz;
				p = (r * s - w) / h[m+1][m] + h[m][m+1];
				q = h[m+1][m+1] - zz - r - s;
				r = h[m+2][m+1];
				s = fabs(p) + fabs(q) + fabs(r);
				p /= s;
				q /= s;
				r /= s;
				if (m == l)
					break;
				tst1 = fabs(p) * (fabs(h[m-1][m-1]) + fabs(zz) + fabs(h[m+1][m+1]));
				tst2 = tst1 + fabs(h[m][m-1]) * (fabs(q) + fabs(r));
				if (tst2 == tst1)
					break;
				}
		
			mp2 = m + 2;
			for (i = mp2; i <= en; i++)
				{
				h[i][i-2] = 0.0;
				if (i != mp2)
					h[i][i-3] = 0.0;
				}
			/* Double qr step involving rows l to en and columns m to en */
			for (k = m; k <= na; k++)
				{
				notlas = (k != na);
				if (k != m)
					{
					p = h[k][k-1];
					q = h[k+1][k-1];
					r = 0.0;
					if (notlas)
						r = h[k+2][k-1];
					x = fabs(p) + fabs(q) + fabs(r);
					if (x == 0.0)
						continue;
					p /= x;
					q /= x;
					r /= x;
					}
	
				s = Sign(sqrt(p*p + q*q + r*r), p);
				if (k != m)
					h[k][k-1] = -s * x;
				else if (l != m)
					h[k][k-1] = -h[k][k-1];
				p += s;
				x = p / s;
				y = q / s;
				zz = r / s;
				q /= p;
				r /= p;
				if (!notlas)
					{
					/* Row modification */
					for (j = k; j < n; j++)
						{
						p = h[k][j] + q * h[k+1][j];
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						} 
				
					j = MIN(en, k + 3);
					/* Column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						}
					/* Accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						}
					}
				else
					{
					/* Row modification */
					for (j = k; j < n; j++)
						{
#						if defined(THINK_C) && !defined(NEED_68881)
							/* Think C 7.0.4 gives "expression too complex" */
							double temp = r * h[k+2][j];
							p = h[k][j] + q * h[k+1][j] + temp;
#						else
							p = h[k][j] + q * h[k+1][j] + r * h[k+2][j];
#						endif
						h[k][j] -= p * x;
						h[k+1][j] -= p * y;
						h[k+2][j] -= p * zz;
						}
				
					j = MIN(en, k + 3);
					/* Column modification */
					for (i = 0; i <= j; i++)
						{
						p = x * h[i][k] + y * h[i][k+1] + zz * h[i][k+2];
						h[i][k] -= p;
						h[i][k+1] -= p * q;
						h[i][k+2] -= p * r;
						}
					/* Accumulate transformations */
					for (i = low; i <= high; i++)
						{
						p = x * z[i][k] + y * z[i][k+1] + zz * z[i][k+2];
						z[i][k] -= p;
						z[i][k+1] -= p * q;
						z[i][k+2] -= p * r;
						}
					}
				}
			}

		if (twoRoots)
			{
			/* Two roots found */
			p = (y - x) / 2.0;
			q = p * p + w;
			zz = sqrt(fabs(q));
			h[en][en] = x + t;
			x = h[en][en];
			h[na][na] = y + t;
			/* DLS 28aug96: Changed "0.0" to "-1e-12" below.  Roundoff errors can cause this value to dip ever-so-
			   slightly below zero even when eigenvalue is not complex. */
			if (q >= -1e-12)
				{
				/* Real pair */
				zz = p + Sign(zz, p);
				wr[na] = x + zz;
				wr[en] = wr[na];
				if (zz != 0.0)
					wr[en] = x - w/zz;
				wi[na] = 0.0;
				wi[en] = 0.0;
				x = h[en][na];
				s = fabs(x) + fabs(zz);
				p = x / s;
				q = zz / s;
				r = sqrt(p*p + q*q);
				p /= r;
				q /= r;
				/* Row modification */
				for (j = na; j < n; j++)
					{
					zz = h[na][j];
					h[na][j] = q * zz + p * h[en][j];
					h[en][j] = q * h[en][j] - p * zz;
					}
				/* Column modification */
				for (i = 0; i <= en; i++)
					{
					zz = h[i][na];
					h[i][na] = q * zz + p * h[i][en];
					h[i][en] = q * h[i][en] - p * zz;
					}
				/* Accumulate transformations */
				for (i = low; i <= high; i++)
					{
					zz = z[i][na];
					z[i][na] = q * zz + p * z[i][en];
					z[i][en] = q * z[i][en] - p * zz;
					}
				}
			else
				{
				/* Complex pair */
				wr[na] = x + p;
				wr[en] = x + p;
				wi[na] = zz;
				wi[en] = -zz;
				}
			en = enm2;
			}
		else
			{
			/* One root found */
			h[en][en] = x + t;
			wr[en] = h[en][en];
			wi[en] = 0.0;
			en = na;
			}
		}
	
	/* All roots found.  Backsubstitute to find vectors of upper triangular form */

	if (norm == 0.0)
		return 0;

	for (en = n - 1; en >= 0; en--)
		{
		p = wr[en];
		q = wi[en];
		na = en - 1;

		/* DLS 28aug96: Changed "0.0" to -1e-12 below (see comment above) */
		if (q < -1e-12)
			{
			/* Complex vector */
			m = na;
			/* Last vector component chosen imaginary so that eigenvector matrix is triangular */
			if (fabs(h[en][na]) > fabs(h[na][en]))
				{
				h[na][na] = q / h[en][na];
				h[na][en] = -(h[en][en] - p) / h[en][na];
				}
			else
				CDiv(0.0, -h[na][en], h[na][na] - p, q, &h[na][na], &h[na][en]);

			h[en][na] = 0.0;
			h[en][en] = 1.0;
			enm2 = na - 1;
			if (enm2 >= 0)
				{
				for (i = enm2; i >= 0; i--)
					{
					w = h[i][i] - p;
					ra = 0.0;
					sa = 0.0;
			
					for (j = m; j <= en; j++)
						{
						ra += h[i][j] * h[j][na];
						sa += h[i][j] * h[j][en];
						}			
					if (wi[i] < 0.0)
						{
						zz = w;
						r = ra;
						s = sa;
						}
					else
						{
						m = i;
						if (wi[i] == 0.0)
							CDiv(-ra, -sa, w, q, &h[i][na], &h[i][en]);
						else
							{
							/* Solve complex equations */
							x = h[i][i+1];
							y = h[i+1][i];
							vr = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i] - q * q;
							vi = (wr[i] - p) * 2.0 * q;
							if ((vr == 0.0) && (vi == 0.0))
								{
								tst1 = norm * (fabs(w) + fabs(q) + fabs(x) + fabs(y) + fabs(zz));
								vr = tst1;
								do	{
									vr *= .01;
									tst2 = tst1 + vr;
									}
									while (tst2 > tst1);
								}
							CDiv(x * r - zz * ra + q * sa, x * s - zz * sa - q * ra, vr, vi, &h[i][na], &h[i][en]);
							if (fabs(x) > fabs(zz) + fabs(q))
								{
								h[i+1][na] = (-ra - w * h[i][na] + q * h[i][en]) / x;
								h[i+1][en] = (-sa - w * h[i][en] - q * h[i][na]) / x;
								}
							else
								CDiv(-r - y * h[i][na], -s - y * h[i][en], zz, q, &h[i+1][na], &h[i+1][en]);
							}
				
						/* Overflow control */
						tst1 = fabs(h[i][na]);
						tst2 = fabs(h[i][en]);
						t = MAX(tst1, tst2);
						if (t != 0.0)
							{
							tst1 = t;
							tst2 = tst1 + 1.0 / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									{
									h[j][na] /= t;
									h[j][en] /= t;
									}
								}
							}
						}
					}
				}
			/* End complex vector */
			}
		else if (q == 0.0)
			{
			/* Real vector */
			m = en;
			h[en][en] = 1.0;
			if (na >= 0)
				{
				for (i = na; i >= 0; i--)
					{
					w = h[i][i] - p;
					r = 0.0;
			
					for (j = m; j <= en; j++)
						r += h[i][j] * h[j][en];
			
					if (wi[i] < 0.0)
						{
						zz = w;
						s = r;
						continue;
						}
					else
						{
						m = i;
						if (wi[i] == 0.0)
							{
							t = w;
							if (t == 0.0)
								{
								tst1 = norm;
								t = tst1;
								do	{
									t *= .01;
									tst2 = norm + t;
									}
									while (tst2 > tst1);
								}			
							h[i][en] = -r / t;
							}
						else
							{
							/* Solve real equations */
							x = h[i][i+1];
							y = h[i+1][i];
							q = (wr[i] - p) * (wr[i] - p) + wi[i] * wi[i];
							t = (x * s - zz * r) / q;
							h[i][en] = t;
							if (fabs(x) > fabs(zz))
								h[i+1][en] = (-r - w * t) / x;
							else
								h[i+1][en] = (-s - y * t) / zz;
							}
				
						/* Overflow control */
						t = fabs(h[i][en]);
						if (t != 0.0)
							{
							tst1 = t;
							tst2 = tst1 + 1. / tst1;
							if (tst2 <= tst1)
								{
								for (j = i; j <= en; j++)
									h[j][en] /= t;
								}
							}
						}
					}
				}
			/* End real vector */
			}
		}
	/* End back substitution */
	
	/* Vectors of isolated roots */
	for (i = 0; i < n; i++)
		{
		if ((i < low) || (i > high))
			{
			for (j = i; j < n; j++)
				z[i][j] = h[i][j];
			}
		}

	/* Multiply by transformation matrix to give vectors of original full matrix */
	for (j = n - 1; j >= low; j--)
		{
		m = MIN(j, high);
		for (i = low; i <= high; i++)
			{
			zz = 0.0;
			for (k = low; k <= m; k++)
				zz += z[i][k] * h[k][j];
			z[i][j] = zz;
			}
		}

	return 0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	BalBak
|
|	EISPACK routine translated from Fortran to C by David Swofford.  Modified EISPACK comments follow.
|
|	This subroutine is a translation of the Algol procedure BALBAK, Num. Math. 13, 293-304 (1969) by Parlett and
|	Reinsch.  Handbook for Auto. Comp., vol. II-Linear Algebra, 315-326 (1971).
|
|	This subroutine forms the eigenvectors of a real general matrix by back transforming those of the corresponding
|	balanced matrix determined by  Balanc.
|
|	ON INPUT:
|
|		n is the order of the matrix.
|
|		low and high are integers determined by Balanc.
|
|		scale contains information determining the permutations and scaling factors used by Balanc.
|
|		m is the number of columns of z to be back transformed.
|
|		z contains the real and imaginary parts of the eigenvectors to be back transformed in its first m columns.
|
|	ON OUTPUT:
|
|		z contains the real and imaginary parts of the transformed eigenvectors in its first m columns.
*/

LOCAL void BalBak(int n, int low, int high, double *scale, int m, double **z)
	{
	int			i, j, k, ii;
	double		s;

	if (m != 0)
		{
		if (high != low)
			{
			for (i = low; i <= high; i++)
				{
				s = scale[i];	/* left hand eigenvectors are back transformed if this statement is replaced by
								   s = 1.0/scale[i] */
				for (j = 0; j < m; j++)
					z[i][j] *= s;
				}
			}
		for (ii = 0; ii < n; ii++)
			{
			i = ii;
			if ((i < low) || (i > high))
				{
				if (i < low)
					i = low - ii;
				k = (int)scale[i];
				if (k != i)
					{
					for (j = 0; j < m; j++)
						{
						s = z[i][j];
						z[i][j] = z[k][j];
						z[k][j] = s;
						}
					}
				}
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	EigenRG
|
|	This subroutine calls the recommended sequence of subroutines from the eigensystem subroutine package (EISPACK) to
|	find the eigenvalues of a real general matrix.  It was converted from Fortran to C by David Swofford.
|
|	ON INPUT:
|
|		n  is the order of the matrix 'a'
|
|		a  contains the real general matrix
|
|	ON OUTPUT:
|
|		wr  and  wi  contain the real and imaginary parts, respectively, of the eigenvalues. Complex conjugate pairs of
|		eigenvalues appear consecutively with the eigenvalue having the positive imaginary part first.
|
|		z  contains the real and imaginary parts of the eigenvectors.  If the j-th eigenvalue is real, the j-th column
|		of  z  contains its eigenvector.  If the j-th eigenvalue is complex with positive imaginary part, the j-th and
|		(j+1)-th columns of  z  contain the real and imaginary parts of its eigenvector.  The conjugate of this vector
|		is the eigenvector for the conjugate eigenvalue.
|
|		ierr  is an integer output variable set equal to an error completion code described in the documentation for Hqr
|		and Hqr2.  The normal completion code is zero.
|
|		iv1  and  fv1  are temporary storage vectors of size n
*/

LOCAL int EigenRG(int n, double **a, double *wr, double *wi, double **z, int *iv1, double *fv1)
	{
	int			is1, is2, ierr;

	Balanc(n, a, &is1, &is2, fv1);
	ElmHes(n, is1, is2, a, iv1);
	ElTran(n, is1, is2, a, iv1, z);
	ierr = Hqr2(n, is1, is2, a, wr, wi, z);
	if (ierr == 0)
		BalBak(n, is1, is2, fv1, n, z);

	return ierr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
|	EigenRealGeneral
|
|	Calculates eigenvalues and eigenvectors of a general real matrix assuming that all eigenvalues
|	are real, using routines from the public domain EISPACK package.
*/

ENTRY int EigenRealGeneral(
  int		n,
  double	**a,		/* input matrix in row-ptr representation; will be destroyed */
  double	*v,			/* array of size 'n' to receive eigenvalues */
  double	*vi,		/* work vector of size 'n' for imaginary components of eigenvalues */
  double	**u,		/* matrix in row-ptr representation to receive eigenvectors */
  int		*iwork,		/* work vector of size 'n' */
  double	*work)		/* work vector of size 'n' */
	{
	int		i;

	if (EigenRG(n, a, v, vi, u, iwork, work) != RC_OK)
		{
#		if defined(PAUP)
			/* This should never happen in PAUP* */
			ErrorMsg(Err_Internal_Eigen, NULL);
#		else
			puts("\nInternal error in 'EigenRealGeneral'\n");
#		endif
		return RC_Error;
		}

#	if defined(DEBUG_SHOW_EIGVALVECT)
		Pprintf("eigenvalues=%g %g %g %g\n", v[0], v[1], v[2], v[3]);
		Pprintf("eigenvectors =\n");
		for (i = 0; i < n; i++)
			Pprintf("%g %g %g %g\n", u[i][0], u[i][1], u[i][2], u[i][3]);
#	endif

#	if defined(PAUP)
	for (i = 0; i < n; i++)
		{
		if (vi[i] != 0.0)
			return RC_COMPLEX_EVAL;
		}
#	endif

	return RC_OK;
	}

#endif

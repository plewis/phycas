/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
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

/*	linalg.h
|
|	Prototypes for matrix-inversion and eigensystem functions
|
|	Copyright (c) 2002 by David L. Swofford, Florida State University.
|	All rights reserved.
*/

#ifdef __cplusplus
extern "C" {
#endif

#define NEED_MATINV		/* added by POL for Phycas */

#if defined(PAUP)
//@ below defined TEMP
#	define NEED_NONSYM_EIG		/* leave undefined for PAUP* -- doesn't need nonsymmetric eigenstuff */
#	define NEED_MATINV			/* leave undefined for PAUP* -- doesn't need matrix inversion */
#endif

#if defined(NEED_NONSYM_EIG)
#	define RC_COMPLEX_EVAL 2	/* return code that complex eigenvalue obtained */
	extern int	EigenRealGeneral(int, double **, double *, double *, double **, int *, double *);
#endif

#if defined(NEED_MATINV)
	extern int	InvertMatrix(double **, int, double *, int *, double **);
#endif

extern int	LUDecompose(double **, int, double *, int *, double *);
extern int	EigenRealSymmetric(int, double **, double *, double **, double *);

#ifdef __cplusplus
}
#endif

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

#if ! defined(Q_MATRIX_HPP)
#define Q_MATRIX_HPP

extern "C"
{
#include "linalg.h"
}

#if defined(_MSC_VER)
#	pragma warning(disable: 4267)	// warning about loss of data when converting size_t to int
#endif

#include <cmath>
#include <vector>
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/tuple.hpp>
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	The QMatrix class manages a matrix of relative instantaneous rates. The primary purpose of the class is to compute
|	transition probabilities for the GTR model (and other models in which there do not exist closed-form expressions for
|	transition probabilities). The Q matrix is a square stochastic (i.e. rows sum to 0.0) matrix that need not be
|	symmetric itself, but has the property that the matrix product Pi*Q is symmetric (where Pi is a diagonal matrix of
|	state frequencies). When setQMatrix is used to set the Q matrix (or setQMatrixElement is used to change one element
|	of the Q matrix), the QMatrix class computes the eigenvalue and eigenvectors of Q. The eigenvalues and eigenvectors
|	are used by calcPMatrix to compute transition matrices needed for likelihood calculations.
*/
class QMatrix
	{
	friend class GTR;
	friend class Codon;

	public:

										QMatrix();
										~QMatrix();

		// Accessors
		//
		unsigned						getDimension();

		// Modifiers
		//
		void							setRelativeRates(const std::vector<double> & rates);
		void							setStateFreqs(const std::vector<double> & freqs);

		// Utilities
		//
#if defined(PYTHON_ONLY)
#	if defined(USING_NUMARRAY)
		boost::python::numeric::array	getQMatrix();
		boost::python::numeric::array	getEigenVectors();
		boost::python::numeric::array	getPMatrix(double edgelen);
#	else
		std::vector<double>	getQMatrix();
		std::vector<double>	getEigenVectors();
		std::vector<double>	getPMatrix(double edgelen);
#	endif
#endif
		std::vector<double>				getEigenValues();

	private:

		void							recalcPMat(double * * pmat, double edgelen);	// used by GTR
		//void							recalcPMatrix(std::vector<double> & P, double edgelen);
		std::string						showQMatrix();
		void							clear();
		void							recalcQMatrix();
		void							recalcQMatrixImpl();
		void							clearAllExceptFreqsAndRates();
		void							redimension(unsigned new_dim);
		std::string						showMatrixImpl(const double * q) const;

		void							flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const;

	private:

		unsigned						dimension;		/**< Number of elements in any given row or column of Q matrix (e.g. 4 for a 4x4 Q matrix) */
		std::vector<int>				dim_vect;		/**< Vector of length 2 describing the shape of the Q matrix for purposes of sending matrices across to Python as NumArray objects. For a 4x4 Q matrix, dim_vect[0] would be 4 and dim_vect[1] would be 4. */

		unsigned						flat_length;	/**< Number of elements in Q matrix (e.g. 16 for a 4x4 Q matrix) */

		std::vector<double>				pi;				/**< The state frequencies (length equals dimension) */
		std::vector<double>				sqrtPi;			/**< The square roots of the state frequencies (length equals dimension) */

		std::vector<double>				rr;				/**< The relative rates (elements in the upper diagonal of the R matrix). If the R matrix is 4x4, the order of the six elements in the relrates vector should be R[0][1], R[0][2], R[0][3], R[1][2], R[1][3] and R[2][3]. The R matrix is combined with the pi vector to create the Q matrix. */

		std::vector<double>				expwv;			/**< Workspace used for storing precalculated exp(w*v), where w is an eigenvalue and v an edge length; used in recalcPMat */

		double							edgelen_scaler;	/**< factor needed */
		double * *						qmat;			/**< */

		double *						qmat_begin;		/**< */
		double *						qmat_end;		/**< */

		double *						w;				/**< matrix of eigenvalues computed by EigenRealSymmetric*/
		double *						w_begin;		/**< */
		double *						w_end;			/**< */

		double * *						z;				/**< matrix of eigenvectors computed by EigenRealSymmetric */
		double *						z_begin;		/**< */
		double *						z_end;			/**< */

		double *						fv;				/**< workspace needed by EigenRealSymmetric */

		bool							q_dirty;		/**< If true, then the Q-matrix has been changed and eigenvalues and eigenvectors need to be recalculated */
	};
}	// namespace phycas

#endif


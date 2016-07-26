/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2016 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#if ! defined(VARCOV_MATRIX_HPP)
#define VARCOV_MATRIX_HPP

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
|	The VarCovMatrix class manages a variance-covariance matrix.
*/
class VarCovMatrix
	{
	public:

										VarCovMatrix();
										~VarCovMatrix();

        void                            example();
        void                            setVarCovMatrix(std::vector<double> flatvcmat);

		// Accessors
		//
		unsigned						getDimension();
        std::vector<double>             getEigenValues();
		std::vector<double>             getEigenVectors();
        std::vector<double>             getVarCovMatrix();
		std::string						showVarCovMatrix();
		std::string						showEigenVectorsAsMatrix();

	private:

		void							redimension(unsigned new_dim);
		void							clear();
		std::string						showMatrixImpl(const double * q) const;
		void							flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const;
        void                            calcEigensystem();

	private:

		unsigned						dimension;		/**< Number of elements in any given row or column of Q matrix (e.g. 4 for a 4x4 Q matrix) */
		double * *						vcmat;			/**< the variance covariance matrix stored by this object */
		double *						w;				/**< matrix of eigenvalues computed by EigenRealSymmetric*/
		double * *						z;				/**< matrix of eigenvectors computed by EigenRealSymmetric */
		double *						fv;				/**< workspace needed by EigenRealSymmetric */
	};
}	// namespace phycas

#endif


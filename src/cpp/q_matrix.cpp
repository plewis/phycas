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

#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/format.hpp>
#include "ncl/nxsallocatematrix.h"
#include "q_matrix.hpp"
#include "xlikelihood.hpp"
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|
*/
QMatrix::QMatrix()
  : dimension(0), flat_length(0)
  , qmat(0), qmat_begin(0), qmat_end(0)
  , w(0), w_begin(0), w_end(0)
  , z(0), z_begin(0), z_end(0)
  , fv(0)
  , q_dirty(true)
	{
	// Set up Q matrix representing the JC69 model by default
	rr.assign(6, 1.0);
	pi.assign(4, 0.25);
	sqrtPi.assign(4, 0.5);
	recalcQMatrix();
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
QMatrix::~QMatrix()
	{
	clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns object to its just-constructed state, with the exception that the vectors `rr', `pi' and `sqrtPi' are not
|	cleared.
*/
void QMatrix::clearAllExceptFreqsAndRates()
	{
	dimension	= 0;
	flat_length	= 0;
	qmat_begin	= 0;
	qmat_end	= 0;
	w_begin		= 0;
	w_end		= 0;
	z_begin		= 0;
	z_end		= 0;

	dim_vect.clear();

	if (qmat)
		{
		DeleteTwoDArray<double>(qmat);
		qmat = 0;
		}

	if (w)
		{
		delete [] w;
		w = 0;
		}

	if (fv)
		{
		delete [] fv;
		fv = 0;
		}

	if (z)
		{
		DeleteTwoDArray<double>(z);
		z = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns object to its just-constructed state by calling clearAllExceptFreqsAndRates and then also clearing the
|	`rr', `pi' and `sqrtPi' vectors.
*/
void QMatrix::clear()
	{
	clearAllExceptFreqsAndRates();
	pi.clear();
	sqrtPi.clear();
	rr.clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `dimension', which is the number of rows (and
|	columns) in the Q matrix.
*/
unsigned QMatrix::getDimension()
	{
	return dimension;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies elements from the supplied vector `rates' to the data member `rr' and sets `q_dirty' to true.
*/
void QMatrix::setRelativeRates(const std::vector<double> & rates)
	{
	rr.resize(rates.size());
	std::copy(rates.begin(), rates.end(), rr.begin());
	q_dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Copies elements from the supplied vector `rates' to the data member `rr' and sets `q_dirty' to true.
*/
void QMatrix::setStateFreqs(const std::vector<double> & freqs)
	{
	pi.resize(freqs.size());
	std::copy(freqs.begin(), freqs.end(), pi.begin());

	sqrtPi.resize(freqs.size());
	std::transform(pi.begin(), pi.end(), sqrtPi.begin(), boost::lambda::bind(static_cast<double(*)(double)>(&std::sqrt), boost::lambda::_1));
	q_dirty = true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reallocates all arrays and vectors (except `rr', `pi' and `sqrtPi') according to the value of `new_dim' and sets
|	`dimension' equal to `new_dim' before returning. Assumes `new_dim' is greater than zero. Note that `qmat', `w', `z'
|	and `fv' are allocated but not filled with any values when this function returns.
*/
void QMatrix::redimension(unsigned new_dim)
	{
	PHYCAS_ASSERT(new_dim > 0);
	clearAllExceptFreqsAndRates();

	dimension = new_dim;
	flat_length = new_dim*new_dim;

	expwv.resize(new_dim);

	// Set the shape of a NumArray object that represents Q, P (transition probs), E (eigenvectors) or V (eigenvalues)
	dim_vect.push_back((int)dimension);
	dim_vect.push_back((int)dimension);

	// Create and fill qmat with default (Mk model) values
	qmat = NewTwoDArray<double>(dimension, dimension);
	PHYCAS_ASSERT(qmat);

	qmat_begin		= &qmat[0][0];
	qmat_end		= qmat_begin + flat_length;

	// Create w (array of eigenvalues)
	w = new double[dimension];
	w_begin = &w[0];
	w_end = w_begin + dimension;

	// Create z (two-dimensional array of eigenvectors)
	z = NewTwoDArray<double>(dimension, dimension);
	z_begin		= &z[0][0];
	z_end		= z_begin + flat_length;

	// Create fv (workspace used by EigenRealSymmetric)
	fv = new double[dimension];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `q_dirty' is true, calls recalcQMatrixImpl to recompute the `qmat', `z' and `w' data members. If `q_dirty' is
|	false, however, this function returns immediately and thus (because it is declared inline) is computationally
|	inexpensive.
*/
void QMatrix::recalcQMatrix()
	{
	if (!q_dirty)
		return;
	recalcQMatrixImpl();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvalues in the data member `w' as a vector. Calls recalcQMatrix first to ensure that `w' is up to
|	date.
*/
std::vector<double> QMatrix::getEigenValues()
	{
	recalcQMatrix();
	std::vector<double> v(dimension, 0.0);
	std::copy(w_begin, w_end, v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
std::string QMatrix::showQMatrix()
	{
	recalcQMatrix();
	return showMatrixImpl(qmat_begin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores a flattened version of the supplied 2-dimensional array `twoDarr', storing the result in the supplied std::vector<double>
|	reference variable `p'. The supplied `twoDarr' should be laid out so that rows occupy contiguous memory.
*/
void QMatrix::flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const
	{
	unsigned flat_length = dim*dim;
	p.resize(flat_length);
	double * twoD_begin = &twoDarr[0][0];
	double * twoD_end   = twoD_begin + flat_length;
	std::copy(twoD_begin, twoD_end, p.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Uses `pi' and `rr' vectors to create the two-dimensional array `qmat'. The data member `qmat' is reallocated if
|	changes in `pi' and `rr' imply a new dimension (or if `qmat' is NULL). Recomputes eigenvalues and eigenvectors
|	corresponding to the new Q matrix. Assumes `q_dirty' is true, and sets `q_dirty' to false when finished. Throws an
|	exception if lengths of the `pi' and `rr' vectors are incompatible.
*/
void QMatrix::recalcQMatrixImpl()
	{
	PHYCAS_ASSERT(q_dirty);

	// pi, sqrtPi and rr should all be non-empty
	PHYCAS_ASSERT(pi.size() > 0);
	PHYCAS_ASSERT(sqrtPi.size() == pi.size());
	PHYCAS_ASSERT(rr.size() > 0);

	// First check to make sure lengths of `rr' and `pi' are compatible with each other
	unsigned dim_pi = (unsigned)pi.size();
	unsigned dim_rr = (unsigned)((1.0 + std::sqrt(1.0 + 8.0*rr.size()))/2.0); // rr.size() = (n^2 - n)/2, so use quadratic formula to get n
	if (dim_pi != dim_rr)
		{
		throw XLikelihood(boost::str(boost::format("Number of relative rates (%d) and number of state frequencies (%d) specified are incompatible") % dim_rr % dim_pi));
		}

	// If qmat is NULL or dimension differs from dim_pi and dim_rr, reallocate qmat, z, w and dv
	if (!qmat || (dimension != dim_pi))
		redimension(dim_pi);

	// This vector will hold the row sums
	std::vector<double> row_sum(dimension, 0.0);

	unsigned i, j, k = 0;
	double sum_for_scaling = 0.0;
	for (i = 0; i < dimension; ++i)
		{
		double pi_i = pi[i];
		double sqrtPi_i = sqrtPi[i];
		for (j = i + 1; j < dimension; ++j, ++k)
			{
			double rr_k = rr[k];
			double pi_j = pi[j];
			double sqrtPi_j = sqrtPi[j];

			// set value in upper triangle
			qmat[i][j] = rr_k*sqrtPi_i*sqrtPi_j;

			// set value in lower triangle
			qmat[j][i] = qmat[i][j];

			// add to relevant row sums
			row_sum[i] += rr_k*pi_j;
			row_sum[j] += rr_k*pi_i;

			// add to total expected number of substitutions
			sum_for_scaling += 2.0*pi_i*rr_k*pi_j;
			}
		qmat[i][i] = -row_sum[i];
		}

	PHYCAS_ASSERT(sum_for_scaling > 0.0);
	edgelen_scaler = 1.0/sum_for_scaling;

	// Calculate eigenvalues (w) and eigenvectors (z)
	int err_code = EigenRealSymmetric(dimension, qmat, w, z, fv);

	if (err_code != 0)
		{
		clearAllExceptFreqsAndRates();
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues of the Q matrix");
		}

	q_dirty = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Recomputes a transition probability matrix for an edge length `edgelen', storing it in `pmat'. If either state
|   frequencies or relative rates have changed, the Q matrix is reconstructed and eigenvalues and eigenvectors
|   recomputed before the transition matrix is recomputed.
*/
void QMatrix::recalcPMat(
  double * * pmat,		/**< is the transition matrix to recalculate */
  double edgelen) 		/**< is the edge length */
	{
	recalcQMatrix();
    double t = edgelen;

    // The next two lines fix the "Rota" bug; see BUGS file for details
    if (t < 1.e-8)
        t = 1.e-8; //TreeNode::edgeLenEpsilon;

	// Adjust the supplied edgelen to account for the fact that the expected number of substitutions
	// implied by the Q matrix is not unity
	double v = t*edgelen_scaler;

    // Precalculate exp to avoid doing the same calculation dimension*dimension times
	for (unsigned k = 0; k < dimension; ++k)
		{
		expwv[k] = std::exp(w[k]*v);
        }

	// Exponentiate eigenvalues and put everything back together again
	// Real symmetric matrices can be diagonalized using Z*exp(D)*Z^T, where Z is the
	// orthogonal matrix of eigenvectors and D is the diagonal matrix of eigenvalues,
	// each multiplied by time (scaled to equal expected number of substitutions)
	for (unsigned i = 0; i < dimension; ++i)
		{
		double sqrtPi_i = sqrtPi[i];
		for (unsigned j = 0; j < dimension; ++j)
			{
			double factor = sqrtPi[j]/sqrtPi_i;
			double Pij = 0.0;
			for (unsigned k = 0; k < dimension; ++k)
				{
				double tmp = z[i][k]*z[j][k]*expwv[k];
				Pij +=  tmp;
				}
			pmat[i][j] = Pij*factor;
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
// void QMatrix::recalcPMatrix(
//   std::vector<double> & P,	/**< */
//   double edgelen)			/**< */
// 	{
// 	recalcQMatrix();
//
//     double t = edgelen;
//
//     // The next two lines fix the "Rota" bug; see BUGS file for details
//     if (t < 1.e-8)
//         t = 1.e-8; //TreeNode::edgeLenEpsilon;
//
//     // Adjust the supplied edgelen to account for the fact that the expected number of substitutions
// 	// implied by the Q matrix is not unity
// 	double v = t*edgelen_scaler;
//
// 	P.clear();
// 	P.reserve(flat_length);
//
// 	// Exponentiate eigenvalues and put everything back together again
// 	for (unsigned i = 0; i < dimension; ++i)
// 		{
// 		double sqrtPi_i = sqrtPi[i];
// 		for (unsigned j = 0; j < dimension; ++j)
// 			{
// 			double factor = sqrtPi[j]/sqrtPi_i;
// 			double Pij = 0.0;
// 			for (unsigned k = 0; k < dimension; ++k)
// 				{
// 				double tmp = z[i][k]*z[j][k]*std::exp(w[k]*v);
// 				Pij +=  tmp;
// 				}
// 			P.push_back(Pij*factor);
// 			}
// 		}
// 	}

/*----------------------------------------------------------------------------------------------------------------------
|
*/
std::string QMatrix::showMatrixImpl(const double * q) const
	{
	unsigned i, j;

	// Output one column for row labels and a label for every column of qmat
	std::string s = str(boost::format("%12s ") % " ");
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		}
	s += "\n";

	// Output rows of q
	double * pq = (double *)q;	//PELIGROSO
	for (i = 0; i < dimension; ++i)
		{
		s += str(boost::format("%12d ") % i);
		for (j = 0; j < dimension; ++j)
			{
			s += str(boost::format("%12.5f ") % *pq++);
			}
		s += "\n";
		}

	return s;
	}

#if defined(PYTHON_ONLY)
#if defined(USING_NUMARRAY)
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entries in the data member `qmat' as a NumArray. Calls recalcQMatrix first to ensure that `qmat' is up
|	to date.
*/
boost::python::numeric::array QMatrix::getQMatrix()
	{
	recalcQMatrix();
	return num_util::makeNum(qmat_begin, dim_vect);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `z' as a NumArray. Calls recalcQMatrix first to ensure that `z' is up
|	to date.
*/
boost::python::numeric::array QMatrix::getEigenVectors()
	{
	recalcQMatrix();
	return num_util::makeNum(z_begin, dim_vect);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the transition probability matrix as a vector of double values (where the 2-dimensional matrix is converted
|   to a 1-dimensional vector by storing rows one after the other. Intended for debugging (not fast).
*/
boost::python::numeric::array QMatrix::getPMatrix(double edgelen)
	{
	double * * pMat = NewTwoDArray<double>(dimension, dimension);
	recalcPMat(pMat, edgelen);
	return num_util::makeNum(pMat, dim_vect);
	}
#else
/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entries in the data member `qmat' as a NumArray. Calls recalcQMatrix first to ensure that `qmat' is up
|	to date.
*/
std::vector<double> QMatrix::getQMatrix()
	{
	recalcQMatrix();
	std::vector<double> p;
	flattenTwoDMatrix(p, qmat, dimension);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `z' as a NumArray. Calls recalcQMatrix first to ensure that `z' is up
|	to date.
*/
std::vector<double> QMatrix::getEigenVectors()
	{
	recalcQMatrix();
	std::vector<double> p;
	flattenTwoDMatrix(p, z, dimension);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the transition probability matrix as a vector of double values (where the 2-dimensional matrix is converted
|   to a 1-dimensional vector by storing rows one after the other. Intended for debugging (not fast).
*/
std::vector<double> QMatrix::getPMatrix(double edgelen)
	{
	double * * pMat = NewTwoDArray<double>(dimension, dimension);
	recalcPMat(pMat, edgelen);
	std::vector<double> p;
	flattenTwoDMatrix(p, pMat, dimension);
	DeleteTwoDArray<double>(pMat);
	return p;
	}
#endif
#endif


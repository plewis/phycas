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
#include <boost/foreach.hpp>
#include "ncl/nxsallocatematrix.h"
#include "pwk_marglike.hpp"
#include "xlikelihood.hpp"
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for all vectors and matrices: _V, _S, _Sinv, and _z are all p by p; _posterior_samples has length 
|   n (but each element has a _parameters vector of length p), and _w is of length p.
*/
PWKMargLike::PWKMargLike(unsigned n, unsigned p)
  : _n(0), _p(0), _V(0), _S(0), _Sinv(0), _w(0), _z(0)
	{
	PHYCAS_ASSERT(n > 0);
	PHYCAS_ASSERT(p > 0);
	_n = n;
	_p = p;

	// Create _V (variance-covariance matrix)
	_V = NewTwoDArray<double>(_p, _p);

	// Create _S (square root of variance-covariance matrix)
	_S = NewTwoDArray<double>(_p, _p);

	// Create _Sinv (inverse square root of variance-covariance matrix)
	_Sinv = NewTwoDArray<double>(_p, _p);

	// Create _w (array of eigenvalues)
	_w = new double[_p];

	// Create _z (two-dimensional array of eigenvectors)
	_z = NewTwoDArray<double>(_p, _p);

    // This vector will hold the n standardized parameter vectors for one tree topology
    _posterior_samples.resize(n);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Destroys arrays and matrices that do not automatically destroy themselves.
*/
PWKMargLike::~PWKMargLike()
	{
	if (_V)
		{
		DeleteTwoDArray<double>(_V);
		_V = 0;
		}

	if (_S)
		{
		DeleteTwoDArray<double>(_S);
		_S = 0;
		}

	if (_Sinv)
		{
		DeleteTwoDArray<double>(_Sinv);
		_Sinv = 0;
		}

	if (_w)
		{
		delete [] _w;
		_w = 0;
		}

	if (_z)
		{
		DeleteTwoDArray<double>(_z);
		_z = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `_n', which is the number of sampled parameter
|   vectors for a single tree topology.
*/
unsigned PWKMargLike::getNumSamples()
	{
	return _n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `_p', which is the number of sampled parameters
|   for a single sampled tree (includes edge lengths and model parameters).
*/
unsigned PWKMargLike::getNumParams()
	{
	return _p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies v (vector of sampled parameters) to row `row' of the `_posterior_samples' vector.
*/
void PWKMargLike::copyParamVector(unsigned row, std::vector<double> v)
    {
    assert(_posterior_samples.size() > row);
    PosteriorSample & sample = _posterior_samples[row];
    sample._parameters.resize(v.size());
    std::copy(v.begin(), v.end(), sample._parameters.begin());
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Calculates mean vector _means, variance-covariance matrix _V, standard deviation matrix _S (and also _Sinv) and uses
|   these to standardize the parameter vectors stored in `_posterior_samples'.
*/
void PWKMargLike::estimateMargLike()
    {
    // if _means is not empty, then standardization has already been done
    assert(_means.empty());

    // calculate means
    _means.resize(_p, 0.0);
    typedef PosteriorSample & sample_ref;
    BOOST_FOREACH(sample_ref sample, _posterior_samples)
        {
        unsigned i = 0;
        BOOST_FOREACH(double x, sample._parameters)
            {
            _means[i++] += x;
            }
        }
    std::transform(_means.begin(), _means.end(), _means.begin(), boost::lambda::_1/_n);


    // center observations
    BOOST_FOREACH(sample_ref sample, _posterior_samples)
        {
        unsigned i = 0;
        BOOST_FOREACH(double & x, sample._parameters)
            {
            x -= _means[i++];
            }
        }

    // construct _V = (C' C)/(n-1), where C = n by p centered observations matrix
    for (unsigned i = 0; i < _p; i++)
        {
        for (unsigned j = 0; j < _p; j++)
            {
            double xsum = 0.0;
            for (unsigned k = 0; k < _n; k++)
                {
                sample_ref sample = _posterior_samples[k];
                xsum += sample._parameters[i]*sample._parameters[j];
                }
            _V[i][j] = xsum/(_n-1);
            }
        }

    // calculate standard deviation matrix _S and its inverse, _Sinv
    calcEigensystem();
    calcSqrtVarCov();

    // standardize centered _observations by multiplying by _Sinv
    for (unsigned i = 0; i < _n; i++)
        {
        sample_ref sample = _posterior_samples[i];
        for (unsigned j = 0; j < _p; j++)
            {
            double xsum = 0.0;
            for (unsigned k = 0; k < _p; k++)
                {
                xsum += sample._parameters[k]*_Sinv[k][j];
                }
            sample._parameters[j] = xsum;
            }
        }

    }

/*----------------------------------------------------------------------------------------------------------------------
|   Computes _S and _Sinv from eigenvalues _w and eigenvectors _z, which are assumed to have already been calculated
|   from _V.
*/
void PWKMargLike::calcSqrtVarCov()
	{
    // Compute _S (square root of variance-covariance matrix)
    //
    //      | _z[0][0] _z[0][1] | | sqrt(_w[0])        0       | | _z[0][0] _z[1][0] |
    // _S = |                   | |                            | |                   |
    //      | _z[1][0] _z[1][1] | |      0         sqrt(_w[1]) | | _z[0][1] _z[1][1] |
    //
    //      | _z[0][0] sqrt(_w[0]) _z[0][0] + _z[0][1] sqrt(_w[1]) _z[0][1]     _z[0][0] sqrt(_w[0]) _z[1][0] + _z[0][1] sqrt(_w[1]) _z[1][1] |
    //    = |                                                                                                                                 |
    //      | _z[1][0] sqrt(_w[0]) _z[0][0] + _z[1][1] sqrt(_w[1]) _z[0][1]     _z[1][0] sqrt(_w[0]) _z[1][0] + _z[1][1] sqrt(_w[1]) _z[1][1] |

	double * * dbl2d = NewTwoDArray<double>(_p, _p);
    for (unsigned i = 0; i < _p; i++)
        {
        for (unsigned j = 0; j < _p; j++)
            {
            double tmp = 0.0;
            for (unsigned k = 0; k < _p; k++)
                {
                tmp += _z[i][k] * sqrt(_w[k]) * _z[j][k];
                }
            _S[i][j] = tmp;
            dbl2d[i][j] = tmp;
            }
        }

    // Compute _Sinv (inverse of square root of variance-covariance matrix)
    std::vector<double> dbl1d(_p, 0.0);
    std::vector<int> int1d(_p, 0.0);
	int err_code = InvertMatrix(dbl2d, _p, &dbl1d[0], &int1d[0], _Sinv);
    DeleteTwoDArray<double>(dbl2d);

	if (err_code != 0)
		{
		throw XLikelihood("Error encountered when inverting the square root of the variance-covariance matrix");
		}

    if (false)
        {
        // debug check: see if _S multipled by _Sinv equals I
        std::cerr << "\n_S * _Sinv:" << std::endl;
        double * * I = NewTwoDArray<double>(_p, _p);
        matrixMultiply(_S, _Sinv, I);
        std::cerr << showMatrix(I) << std::endl;

        // debug check: see if _S multipled by _Sinv equals I
        double * * V = NewTwoDArray<double>(_p, _p);
        double * * Z = NewTwoDArray<double>(_p, _p);
        matrixMultiply(_S, _S, V);
        matrixSubtract(_V, V, Z);
        std::cerr << "\n_V - _S * _S:" << std::endl;
        std::cerr << showMatrix(Z) << std::endl;

        DeleteTwoDArray<double>(I);
        DeleteTwoDArray<double>(V);
        DeleteTwoDArray<double>(Z);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Multiply u and v and return the product in u_times_v.
*/
void PWKMargLike::matrixMultiply(double * * u, double * * v, double * * u_times_v)
    {
    for (unsigned i = 0; i < _p; ++i)
        {
        for (unsigned j = 0; j < _p; ++j)
            {
            double tmp = 0.0;
            for (unsigned k = 0; k < _p; ++k)
                {
                tmp += u[i][k]*v[k][j];
                }
            u_times_v[i][j] = tmp;
            }
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Subtract matrix v from matrix u and return the product in u_minus_v.
*/
void PWKMargLike::matrixSubtract(double * * u, double * * v, double * * u_minus_v)
    {
    for (unsigned i = 0; i < _p; ++i)
        {
        for (unsigned j = 0; j < _p; ++j)
            {
            u_minus_v[i][j] = u[i][j] - v[i][j];
            }
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Add matrix u to matrix u and return the product in u_plus_v.
*/
void PWKMargLike::matrixAdd(double * * u, double * * v, double * * u_plus_v)
    {
    for (unsigned i = 0; i < _p; ++i)
        {
        for (unsigned j = 0; j < _p; ++j)
            {
            u_plus_v[i][j] = u[i][j] + v[i][j];
            }
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes eigenvalues and eigenvectors corresponding to the variance-covariance matrix.
*/
void PWKMargLike::calcEigensystem()
	{
	// Calculate eigenvalues (_w) and eigenvectors (_z)
    std::vector<double> fv(_p, 0.0);
	int err_code = EigenRealSymmetric(_p, _V, _w, _z, &fv[0]);

	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues of the variance-covariance matrix");
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvalues in the data member `_w' as a vector.
*/
std::vector<double> PWKMargLike::getEigenValues()
	{
	std::vector<double> v(_p, 0.0);
	double * w_begin = &_w[0];
	double * w_end = w_begin + _p;
	std::copy(w_begin, w_end, v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of log eigenvalues of the variance-covariance matrix.
*/
double PWKMargLike::logDet()
    {
    double logdet = 0.0;

    for (unsigned i = 0; i < _p; ++i)
        {
        logdet += log(_w[i]);
        }

    return logdet;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvectors in the data member `_z' as a NumArray. Calls recalcQMatrix first to ensure that `_z' is up
|	to date.
*/
std::vector<double> PWKMargLike::getEigenVectors()
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, _z, _p);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats `V' nicely for output.
*/
std::string PWKMargLike::showMatrix(double * * V)
	{
	double * V_begin = &V[0][0];
	return showMatrixImpl(V_begin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats the matrix of eigenvectors nicely for output.
*/
std::string PWKMargLike::showEigenVectorsAsMatrix()
	{
	double * z_begin = &_z[0][0];
	return showMatrixImpl(z_begin);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats the matrix of eigenvectors nicely for output.
*/
std::string PWKMargLike::showMatrixImpl(const double * q) const
	{
	unsigned i, j;

	// Output one column for row labels and a label for every column of _V
	std::string s = str(boost::format("%12s ") % " ");
	for (i = 0; i < _p; ++i)
		{
		s += str(boost::format("%12d ") % i);
		}
	s += "\n";

	// Output rows of q
	double * pq = (double *)q;	//PELIGROSO
	for (i = 0; i < _p; ++i)
		{
		s += str(boost::format("%12d ") % i);
		for (j = 0; j < _p; ++j)
			{
			s += str(boost::format("%12.5f ") % *pq++);
			}
		s += "\n";
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Stores a flattened version of the supplied 2-dimensional array `twoDarr', storing the result in the supplied std::vector<double>
|	reference variable `p'. The supplied `twoDarr' should be laid out so that rows occupy contiguous memory.
*/
void PWKMargLike::flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const
	{
	unsigned flat_length = dim*dim;
	p.resize(flat_length);
	double * twoD_begin = &twoDarr[0][0];
	double * twoD_end   = twoD_begin + flat_length;
	std::copy(twoD_begin, twoD_end, p.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the entries in the data member `_V' as a single dimensional vector.
*/
std::vector<double> PWKMargLike::getVarCovMatrix()
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, _V, _p);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies parameter types from `ptypes' to `_param_type' vector.
*/
void PWKMargLike::copyParamTypes(std::vector<unsigned> & ptypes)
	{
    _param_type.resize(ptypes.size());
    std::copy(ptypes.begin(), ptypes.end(), _param_type.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies supplied newick tree description to `_newick'.
*/
void PWKMargLike::copyNewick(std::string newick)
	{
    _newick = newick;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Debugging check: calculates log-likelihoods of all sampled trees for comparison with those from the original
|   parameter file.
*/
void PWKMargLike::debugCheckLogLikelihoods()
	{
	}

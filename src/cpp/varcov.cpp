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

#include <fstream>
#include <strstream>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string/join.hpp>
#include "ncl/nxsallocatematrix.h"
#include "basic_lot.hpp"
#include "basic_lot.hpp"
#include "varcov.hpp"
#include "xlikelihood.hpp"
#include <boost/math/special_functions/gamma.hpp>   // lgamma
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
| 	Allocates a two dimensional array of doubles as one contiguous block of memory
| 	the dimensions are f by s.
|
| 	The pointer should be freed by a call to DeleteTwoDArrayDouble
|
| 	The array is set up so that:
|
|	for (i = 0 ; i < f ; i++)
|		for (j = 0 ; j < s ; j++)
|			array[i][j];
|
|	would be the same order of access as:
|
|  	double * ptr = **array;
|	for (i = 0 ; i < f*s*t ; i++)
|		*ptr++;
|
|   This is a template specialization for the NewTwoDArray function defined in nxsallocatematrix.h.
*/
double * * NewTwoDArrayDouble(unsigned f, unsigned s)
	{
	assert(f > 0 && s > 0);
	double * * ptr;
	ptr = new double *[f];
	*ptr = new double [f * s];
	for (unsigned fIt = 1 ; fIt < f ; fIt ++)
		ptr[fIt] = ptr[fIt -1] +  s;
	return ptr;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Deletes a 2-dimensional array of doubles and sets ptr to NULL.
|
|   This is a template specialization for the DeleteTwoDArray function defined in nxsallocatematrix.h.
*/
void DeleteTwoDArrayDouble(double * *  & ptr)
	{
	if (ptr)
		{
		delete [] * ptr;
		delete [] ptr;
		ptr = NULL;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Allocates memory for all vectors and matrices: _V, _S, _Sinv, and _z are all p by p; _posterior_samples has length 
|   n (but each element has a _parameters vector of length p), and _w is of length p.
*/
VarCovMatrix::VarCovMatrix(std::string name, unsigned n, unsigned p)
  : _n(0), _p(0), _V(0), _S(0), _Stmp(0), _Sinv(0), _z(0), _shell_midpoint(0.0), _shell_delta(0.0)
    , _log_ratio_sum(0.0), _representative_logkernel(0.0), _representative_loglikelihood(0.0), _representative_logprior(0.0)
    , _representative_logjacobian_edgelen(0.0), _representative_logjacobian_substmodel(0.0), _representative_logjacobian_standardization(0.0)
    , _center_around_mode(false) // currently _center_around_mode cannot be modified anywhere except here
	{
	PHYCAS_ASSERT(n > 0);
	PHYCAS_ASSERT(p > 0);
	_n = n;
	_p = p;
    _name = name;

	// Create _V (variance-covariance matrix)
	_V = NewTwoDArrayDouble(_p, _p);

	// Create _S (square root of variance-covariance matrix)
	_S = NewTwoDArrayDouble(_p, _p);

	// Create _Sinv (inverse square root of variance-covariance matrix)
	_Sinv = NewTwoDArrayDouble(_p, _p);

	// Create _Stmp (inverse square root of variance-covariance matrix)
	_Stmp = NewTwoDArrayDouble(_p, _p);

	// This vector will hold eigenvalues
	_w .resize(_p);

	// Create _z (two-dimensional array of eigenvectors)
	_z = NewTwoDArrayDouble(_p, _p);

    // This vector will hold the n parameter vectors for one tree topology
    _posterior_samples.resize(n);

#if 0
    // fill _posterior_samples with multivariate standard normal vectors for debugging tests
    Lot lot = Lot();
    double debug_sigma = 5.0;
    for (unsigned i = 0; i < _n; ++i)
        {
        PosteriorSample & sample = _posterior_samples[i];
        sample._parameters.resize(_p);
        for (unsigned j = 0; j < _p; ++j)
            {
            sample._parameters[j] = debug_sigma*lot.Normal();
            }
        }
#endif
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Destroys arrays and matrices that do not automatically destroy themselves.
*/
VarCovMatrix::~VarCovMatrix()
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

	if (_Stmp)
		{
		DeleteTwoDArray<double>(_Stmp);
		_Stmp = 0;
		}

	if (_z)
		{
		DeleteTwoDArray<double>(_z);
		_z = 0;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns 0.5*logDet(), which is the log of the Jacobian for the standardization transformation (the log of the
|   determinant of the _S matrix).
*/
double VarCovMatrix::logJacobianForStandardization() const
	{
	return 0.5*logDet();
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns data member `_log_ratio_sum', which is calculated in the member function calcRepresentativeForShell.
*/
double VarCovMatrix::getLogRatioSum() const
	{
	return _log_ratio_sum;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns data member `_log_ratio_vect', which is calculated in the member function calcRepresentativeForShell.
*/
std::vector<double> VarCovMatrix::getLogRatioVect() const
	{
	return _log_ratio_vect;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_logkernel', which is calculated in the member function
|   calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogKernel() const
	{
	return _representative_logkernel;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Calculates sum of elements on the log scale, controlling for underflow, and returns sum on log scale. For example,
|   if a is the largest element of the vector {a,b,c}, then
|   sum = a + b + c
|       = a (1 + b/a + c/a)
|   log(sum) = log(a) +  log(1 + exp(log(b) - log(a)) + exp(log(c) - log(a)))
|   This approach uses only log(a), log(b) and log(c) to compute the log of the sum a + b + c.
*/
double VarCovMatrix::calcSumTermsOnLogScale(std::vector<double> & log_terms) const
    {
    assert(log_terms.size() > 0);
    double max_log_term = *std::max_element(log_terms.begin(), log_terms.end());
    double sum_of_terms = 0.0;
    BOOST_FOREACH(double logx, log_terms)
        {
        sum_of_terms += exp(logx - max_log_term);
        }
    double sum_on_log_scale = max_log_term + log(sum_of_terms);
    return sum_on_log_scale;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_loglikelihood', which is calculated in the member function
|   calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogLikelihood() const
	{
	return _representative_loglikelihood;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_logprior', which is calculated in the member function
|   calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogPrior() const
	{
	return _representative_logprior;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_logjacobian_edgelen', which is calculated in the member
|   function calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogJacobianEdgelen() const
	{
	return _representative_logjacobian_edgelen;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_logjacobian_substmodel', which is calculated in the member
|   function calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogJacobianSubstmodel() const
	{
	return _representative_logjacobian_substmodel;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_logjacobian_standardization', which is calculated in the
|   member function calcRepresentativeForShell.
*/
double VarCovMatrix::getRepresentativeLogJacobianStandardization() const
	{
	return _representative_logjacobian_standardization;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the data member `_representative_param_vector', which is calculated in the member function
|   calcRepresentativeForShell.
*/
std::vector<double> VarCovMatrix::getRepresentativeParamVect() const
	{
	return _representative_param_vector;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `_n', which is the number of sampled parameter
|   vectors for a single tree topology.
*/
unsigned VarCovMatrix::getNumSamples() const
	{
	return _n;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Accessor function that returns the current value of the data member `_p', which is the number of sampled parameters
|   for a single sampled tree (includes edge lengths and model parameters).
*/
unsigned VarCovMatrix::getNumParams() const
	{
	return _p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies v (vector of transformed sampled parameters) to row `row' of the `_posterior_samples' vector. Also copy the
|   log likelihood (lnL), log joint prior (lnP), and log jacobian for the log (or log ratio) transformation to the 
|   sample corresponding to row `row' of the `_posterior_samples' vector.
*/
void VarCovMatrix::copyParamVector(unsigned row, std::vector<double> v, double lnL, double lnP, double lnJ_edgelen, double lnJ_substmodel)
    {
    // Get object corresponding to row
    assert(_posterior_samples.size() > row);
    PosteriorSample & sample = _posterior_samples[row];

    // Copy the parameter vector
    sample._parameters.resize(v.size());
    std::copy(v.begin(), v.end(), sample._parameters.begin());

    // Copy lnL, lnP, and lnJ
    sample._log_likelihood = lnL;
    sample._log_prior = lnP;
    sample._log_jacobian_edgelen = lnJ_edgelen;
    sample._log_jacobian_substmodel = lnJ_substmodel;
    sample._log_jacobian_standardization = 0.0;  // will update this later to equal log jacobian for standardization transformation
    sample._radius = 0.0;           // cannot set this until after standardization
    }

/*----------------------------------------------------------------------------------------------------------------------
|
*/
void VarCovMatrix::setParamNames(std::vector<std::string> names)
    {
    _param_names.resize(names.size());
    std::copy(names.begin(), names.end(), _param_names.begin());
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Calculates mean vector _means, variance-covariance matrix _V, standard deviation matrix _S (and also _Sinv) and uses
|   these to standardize the parameter vectors stored in `_posterior_samples', returning Euclidian distance from the 
|   origin to the furthest sample point.
*/
double VarCovMatrix::standardizeSamples()
    {
    // if _means is not empty, then standardization has already been done
    assert(_means.empty());
    assert(_mode.empty());

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

#if 1
    // print _means vector
    std::cerr << "\n_means = c(";
    std::copy(_means.begin(), _means.end(), std::ostream_iterator<double>(std::cerr, ","));
    std::cerr << ")" << std::endl;
#endif

    // calculate median vector
    _mode.resize(_p, 0.0);
    sample_ref best = _posterior_samples[0];
    double best_log_kernel = best._log_likelihood; // + best._log_prior + best._log_jacobian_edgelen + best._log_jacobian_substmodel + best._log_jacobian_standardization;
    BOOST_FOREACH(sample_ref sample, _posterior_samples)
        {
        double log_kernel = sample._log_likelihood; // + sample._log_prior + sample._log_jacobian_edgelen + sample._log_jacobian_substmodel + sample._log_jacobian_standardization;
        if (log_kernel > best_log_kernel)
            {
            unsigned i = 0;
            BOOST_FOREACH(double x, sample._parameters)
                {
                _mode[i++] = x;
                }
            }
        }

#if 1
    // print _mode vector
    std::cerr << "\n_mode = c(";
    std::copy(_mode.begin(), _mode.end(), std::ostream_iterator<double>(std::cerr, ","));
    std::cerr << ")" << std::endl;
#endif

    if (_center_around_mode)
        {
        // center observations around median
        BOOST_FOREACH(sample_ref sample, _posterior_samples)
            {
            unsigned i = 0;
            BOOST_FOREACH(double & x, sample._parameters)
                {
                x -= _mode[i++];
                }
            }
        }
    else
        {
        // center observations around mean
        BOOST_FOREACH(sample_ref sample, _posterior_samples)
            {
            unsigned i = 0;
            BOOST_FOREACH(double & x, sample._parameters)
                {
                x -= _means[i++];
                }
            }
        }

#if 0
    // save R file for plotting density of each log-transformed and centered parameter
    assert(_p == _param_names.size());
    std::string fn = boost::str(boost::format("%s.R") % _name);

    std::ofstream sh("plots.sh", std::ios::out | std::ios::app);
    sh << "rscript " << fn << std::endl;
    sh.close();

    std::ofstream f(fn, std::ios::out);
    f << "cwd = system('cd \"$( dirname \"$0\" )\" && pwd', intern = TRUE)" << std::endl;
    f << "setwd(cwd)" << std::endl;
    f << boost::str(boost::format("pdf(\"%s.pdf\")") % _name) << std::endl;

    for (unsigned i = 0; i < _p; ++i)
        {
        f << boost::str(boost::format("p%d = c(") % i);
        unsigned j = 0;
        BOOST_FOREACH(sample_ref sample, _posterior_samples)
            {
            if (j > 0)
                f << boost::str(boost::format(",%.5f") % sample._parameters[i]);
            else
                f << boost::str(boost::format("%.5f") % sample._parameters[i]);
            ++j;
            }
        f << ")" << std::endl;
        }
    unsigned nrows = 4;
    unsigned ncols = 3;
    f << boost::str(boost::format("par(mfrow=c(%d,%d))") % nrows % ncols) << std::endl;
    for (unsigned i = 0; i < _p; ++i)
        {
        f << boost::str(boost::format("plot(density(p%d), xlab=\"\", ylab=\"\", main=\"%s\")") % i % _param_names[i]) << std::endl;
        f << "abline(v=0)" << std::endl;
        }
    f << "dev.off()" << std::endl;
    f.close();
#endif

#if 0
    // check centering by calculating mean of centered observations
    std::vector<double> centered_means;
    centered_means.resize(_p, 0.0);
    BOOST_FOREACH(sample_ref sample, _posterior_samples)
        {
        unsigned i = 0;
        BOOST_FOREACH(double x, sample._parameters)
            {
            centered_means[i++] += x;
            }
        }
    std::transform(centered_means.begin(), centered_means.end(), centered_means.begin(), boost::lambda::_1/_n);

    // print centered_means vector
    std::cerr << "centered_means = c(";
    std::copy(centered_means.begin(), centered_means.end(), std::ostream_iterator<double>(std::cerr, ","));
    std::cerr << ")" << std::endl;
#endif

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

    // ensure that _V is symmetric
    for (unsigned i = 0; i < _p - 1; i++)
        {
        for (unsigned j = i + 1; j < _p; j++)
            {
            double tmp = (_V[i][j] + _V[j][i])/2.0;
            _V[i][j] = tmp;
            _V[j][i] = tmp;
            }
        }

    // Calculate standard deviation matrix _S and its inverse, _Sinv
    //std::cerr << "Just about to call calcEigensystem()..." << std::endl;
    calcEigensystem();

    //std::cerr << "Just about to call calcSqrtVarCov()..." << std::endl;
    calcSqrtVarCov();

    //std::cerr << "Just about to standardize _observations..." << std::endl;

    // Standardize centered _observations by premultiplying by _Sinv
    // Also, update jacobian for each observation for the standardization
    double log_det_S = 0.5*logDet();
    double max_radius = 0.0;
    for (unsigned i = 0; i < _n; i++)
        {
        sample_ref sample = _posterior_samples[i];
        sample._standardized.resize(sample._parameters.size());
        sample._log_jacobian_standardization = log_det_S;
        for (unsigned j = 0; j < _p; j++)
            {
            double xsum = 0.0;
            for (unsigned k = 0; k < _p; k++)
                {
                //xsum += sample._parameters[k]*_Sinv[k][j]; // postmultiply (wrong)
                //xsum += _Sinv[j][k]*sample._parameters[k]; // premultiply
                xsum += _Stmp[j][k]*sample._parameters[k]; // premultiply using _Stmp
                }
            sample._standardized[j] = xsum;
            }

        double r = 0.0;
        BOOST_FOREACH(double v, sample._standardized)
            {
            r += v*v;
            }
        sample._radius = sqrt(r);
        if (sample._radius > max_radius)
            max_radius = sample._radius;
        }

    return max_radius;
    }


/*----------------------------------------------------------------------------------------------------------------------
|   Calculate the length of the supplied vector as the square root of the sum of squared elements.
*/
double VarCovMatrix::calcVectorNorm(std::vector<double> & v) const
    {
    //unsigned sz = (unsigned)v.size();
    double sum_of_squares = 0.0;
    BOOST_FOREACH(double x, v)
        {
        sum_of_squares += x*x;
        }
    return sqrt(sum_of_squares);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   For debugging
*/
std::vector<double> VarCovMatrix::destandardizeSample(unsigned row)
    {
    _representative_param_vector.clear();
    _representative_param_vector.resize(_p, 0.0);

    PosteriorSample & sample = _posterior_samples[row];
    std::copy(sample._standardized.begin(), sample._standardized.end(), _representative_param_vector.begin());

    // Debugging: destandardize _representative_param_vector
    // Y' = S^{-1}*(X - Xmean) (_p x _p)*(_p x 1)
    std::vector<double> Y(_p, 0.0);
    matrixTimesVector(_S, &_representative_param_vector[0], &Y[0]);
    if (_center_around_mode)
        vectorPlusVector(&Y[0], &_mode[0], &_representative_param_vector[0]);
    else
        vectorPlusVector(&Y[0], &_means[0], &_representative_param_vector[0]);

    return _representative_param_vector;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns vector of parameter vectors. The first `n' vectors are chosen, while any after that were sampled.
*/
std::vector<double> VarCovMatrix::getRepresentativesForShell(unsigned n, double r, double delta)
	{
    std::vector<double> vect_of_vects;

    Lot rng;
    std::vector<double> Y(_p, 0.0);
    std::vector<double> random_vector(_p, 0.0);
    for (unsigned i = 0; i < n; ++i)
        {
        // Choose random vector of length r by normalizing and rescaling a vector of random normal deviates
        double sumsq = 0.0;
        for (unsigned i = 0; i < _p; ++i)
            {
            double x = rng.Normal();
            random_vector[i] = x;
            sumsq += x*x;
            }
        double vector_length = sqrt(sumsq);
        std::transform(random_vector.begin(), random_vector.end(), random_vector.begin(), boost::lambda::_1*r/vector_length);

        // Destandardize vector
        matrixTimesVector(_S, &random_vector[0], &Y[0]);
        std::vector<double> v(_p, 0.0);
        if (_center_around_mode)
            vectorPlusVector(&Y[0], &_mode[0], &v[0]);
        else
            vectorPlusVector(&Y[0], &_means[0], &v[0]);
        std::copy(v.begin(), v.end(), std::back_inserter(vect_of_vects));
        }

    BOOST_FOREACH(PosteriorSample & p, _posterior_samples)
        {
        if (p._radius > r - delta && p._radius <= r + delta)
            {
            std::copy(p._standardized.begin(), p._standardized.end(), random_vector.begin());

            // Debugging: destandardize random_vector
            matrixTimesVector(_S, &random_vector[0], &Y[0]);
            std::vector<double> v(_p, 0.0);
            if (_center_around_mode)
                vectorPlusVector(&Y[0], &_mode[0], &v[0]);
            else
                vectorPlusVector(&Y[0], &_means[0], &v[0]);
            std::copy(v.begin(), v.end(), std::back_inserter(vect_of_vects));
            }
        }

    return vect_of_vects;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   If there are any samples with radius between `r' - `delta' and `r' + `delta', compute average log kernel for those
|   samples and store in `_representative_logkernel'. If not, choose random point at radius `r' and fill the vector
|   `_representative_param_vector' with the destandardized parameter values. Returns the number of samples found in the
|   shell specified by `r' and `delta' and sets `_shell_midpoint' to 'r' and `_shell_delta' to `delta'.
*/
unsigned VarCovMatrix::calcRepresentativeForShell(double r, double delta)
	{
    _shell_midpoint = r;
    _shell_delta = delta;
    _representative_param_vector.clear();
    _representative_param_vector.resize(_p, 0.0);

    // Find number of samples for this shell
    _representative_logkernel = 0.0;
    _representative_loglikelihood = 0.0;
    _representative_logprior = 0.0;
    _representative_logjacobian_edgelen = 0.0;
    _representative_logjacobian_substmodel = 0.0;
    _representative_logjacobian_standardization = 0.0;
    unsigned num_samples_this_shell = 0;
    std::vector<double> tmp_logkernelvect;
    BOOST_FOREACH(PosteriorSample & p, _posterior_samples)
        {
        if (p._radius > _shell_midpoint - _shell_delta && p._radius <= _shell_midpoint + _shell_delta)
            {
            //std::cerr << boost::str(boost::format("%6d lnL=%.5f lnP=%.5f lnJe=%.5f lnJp=%.5f lnJs=%.5f") % num_samples_this_shell % p._log_likelihood % p._log_prior % p._log_jacobian_edgelen % p._log_jacobian_substmodel % p._log_jacobian_standardization) << std::endl;
            tmp_logkernelvect.push_back(p._log_likelihood + p._log_prior + p._log_jacobian_edgelen + p._log_jacobian_substmodel + p._log_jacobian_standardization);

            _representative_logkernel += p._log_likelihood + p._log_prior + p._log_jacobian_edgelen + p._log_jacobian_substmodel + p._log_jacobian_standardization;
            num_samples_this_shell++;
            }
        }

    //std::cerr << "--> Done listing samples for this shell" << std::endl;

    if (num_samples_this_shell > 0)
        {
        // Representative log kernel value is simple average of log kernel values of samples in this shell
        _representative_logkernel /= num_samples_this_shell;

#if 0
        // Debugging: find element of _posterior_samples closest to _representative_logkernel
        PosteriorSample & closest = _posterior_samples[0];
        BOOST_FOREACH(PosteriorSample & p, _posterior_samples)
            {
            if (p._radius > _shell_midpoint - _shell_delta && p._radius <= _shell_midpoint + _shell_delta)
                {
                double plogkernel = p._log_likelihood + p._log_prior + p._log_jacobian_edgelen + p._log_jacobian_substmodel + p._log_jacobian_standardization;
                double clogkernel = closest._log_likelihood + closest._log_prior + closest._log_jacobian_edgelen + closest._log_jacobian_substmodel + closest._log_jacobian_standardization;
                if (fabs(plogkernel - _representative_logkernel) < fabs(clogkernel - _representative_logkernel))
                    {
                    closest = p;
                    }
                }
            }

        // Debugging: fill _representative_param_vector with sample closest to _representative_logkernel
        _representative_logkernel     = closest._log_likelihood + closest._log_prior + closest._log_jacobian_edgelen + closest._log_jacobian_substmodel + closest._log_jacobian_standardization;
        _representative_loglikelihood = closest._log_likelihood;
        _representative_logprior      = closest._log_prior;
        _representative_logjacobian_edgelen         = closest._log_jacobian_edgelen;
        _representative_logjacobian_substmodel      = closest._log_jacobian_substmodel;
        _representative_logjacobian_standardization = closest._log_jacobian_standardization;
        std::copy(closest._standardized.begin(), closest._standardized.end(), _representative_param_vector.begin());

        // Debugging: destandardize _representative_param_vector
        std::vector<double> Y(_p, 0.0);
        matrixTimesVector(_S, &_representative_param_vector[0], &Y[0]);
        if (_center_around_mode)
            vectorPlusVector(&Y[0], &_mode[0], &_representative_param_vector[0]);
        else
            vectorPlusVector(&Y[0], &_means[0], &_representative_param_vector[0]);

        // Spit out closest
        std::cerr << "*** comparing closest to chosen points ***" << std::endl;
        std::cerr << boost::str(boost::format("%20s %12.5f ") % "closest" % _representative_logkernel);
        BOOST_FOREACH(double v, _representative_param_vector)
            {
            std::cerr << boost::str(boost::format("%12.5f ") % v);
            }
        std::cerr << std::endl;

        // Find index of largest eigenvalue
        unsigned primary_eigenvalue_index = (unsigned)(std::max_element(_w.begin(), _w.end()) - _w.begin());

        // Copy column from eigenvector matrix representing the primary eigenvector
        std::vector<double> primary_eigenvector(_p);
        for (unsigned i = 0; i < _p; ++i)
            {
            primary_eigenvector[i] = _z[i][primary_eigenvalue_index];
            }

        // Scale primary eigenvector to have length _shell_midpoint
        std::transform(primary_eigenvector.begin(), primary_eigenvector.end(), primary_eigenvector.begin(), boost::lambda::_1*_shell_midpoint);
#pragma clang diagnostic ignored "-Wunused-variable"
        double tmp = calcVectorNorm(primary_eigenvector);
        assert(fabs(tmp - _shell_midpoint) < 1.e8);

        // First, draw point on hypersphere of radius _shell_midpoint along primary eigenvector and destandardize it
        matrixTimesVector(_S, &primary_eigenvector[0], &Y[0]);
        if (_center_around_mode)
            vectorPlusVector(&Y[0], &_mode[0], &_representative_param_vector[0]);
        else
            vectorPlusVector(&Y[0], &_means[0], &_representative_param_vector[0]);

        // Spit out point
        std::cerr << boost::str(boost::format("%20s %12s ") % "chosen" % "---");
        BOOST_FOREACH(double v, _representative_param_vector)
            {
            std::cerr << boost::str(boost::format("%12.5f ") % v);
            }
        std::cerr << std::endl;
#endif

        }
    else
        {
        // No samples available for this shell, so pick a point on the hypersphere at radius _shell_midpoint
        // and store that in _representative_param_vector so that log kernel can be computed by caller
        //double log_prior = 0.0;
        //double log_jacobian = 0.0;

        // Find index of largest eigenvalue
        unsigned primary_eigenvalue_index = (unsigned)(std::max_element(_w.begin(), _w.end()) - _w.begin());

        //std::cerr << "Eigenvalues: largest is at index " << primary_eigenvalue_index << std::endl;
        //unsigned which = 0;
        //BOOST_FOREACH(double v, _w)
        //    {
        //    std::cerr << boost::str(boost::format("%12d %12.5f") % which % v) << std::endl;
        //    which++;
        //    }
        //std::cerr << std::endl;

        // Copy column from eigenvector matrix representing the primary eigenvector
        std::vector<double> primary_eigenvector(_p);
        for (unsigned i = 0; i < _p; ++i)
            {
            primary_eigenvector[i] = _z[i][primary_eigenvalue_index];
            }

        // Scale primary eigenvector to have length _shell_midpoint
        std::transform(primary_eigenvector.begin(), primary_eigenvector.end(), primary_eigenvector.begin(), boost::lambda::_1*_shell_midpoint);
#pragma clang diagnostic ignored "-Wunused-variable"
        double tmp = calcVectorNorm(primary_eigenvector);
        assert(fabs(tmp - _shell_midpoint) < 1.e8);

        //std::cerr << "_shell_midpoint                   = " << _shell_midpoint << std::endl;
        //std::cerr << "scaled primary eigenvector length = " << tmp << std::endl;

        // First, draw point on hypersphere of radius _shell_midpoint along primary eigenvector
        // and destandardize it
        std::vector<double> Y(_p, 0.0);
        matrixTimesVector(_S, &primary_eigenvector[0], &Y[0]);
        if (_center_around_mode)
            vectorPlusVector(&Y[0], &_mode[0], &_representative_param_vector[0]);
        else
            vectorPlusVector(&Y[0], &_means[0], &_representative_param_vector[0]);
        }

    // If num_samples_this_shell > 0, then store representative kernel/sample kernel ratios
    _log_ratio_sum = 0.0;
    _log_ratio_vect.resize(num_samples_this_shell, 0.0);
    if (num_samples_this_shell > 0)
        {
        //std::vector<double> log_ratios(num_samples_this_shell, 0.0);
        unsigned i = 0;
        BOOST_FOREACH(PosteriorSample & p, _posterior_samples)
            {
            if (p._radius > _shell_midpoint - _shell_delta && p._radius <= _shell_midpoint + _shell_delta)
                {
                double log_ratio = _representative_logkernel - p._log_likelihood - p._log_prior - p._log_jacobian_edgelen - p._log_jacobian_substmodel - p._log_jacobian_standardization;
                _log_ratio_vect[i++] = log_ratio;
                }
            }
        _log_ratio_sum = calcSumTermsOnLogScale(_log_ratio_vect);

#if 0
        std::cerr << "\ncalcRepresentativeForShell: (" << (r - delta) << ", " << (r + delta) << ")" << std::endl;
        std::cerr << "log kernel values for this shell:" << std::endl;
        std::cerr << "logk <- c(";
        std::copy(tmp_logkernelvect.begin(), tmp_logkernelvect.end(), std::ostream_iterator<double>(std::cerr, ","));
        std::cerr << ")" << std::endl;
        std::cerr << "plot(density(logk))" << std::endl;
        std::cerr << "abline(v=" << _representative_logkernel << ")" << std::endl;
        std::cerr << "repr, sampled, repr-sampled, exp(repr-sampled)" << std::endl;
        std::cerr << boost::str(boost::format("%20s %20s %20s %20s") % "repr" % "sampled" % "repr-sampled" % "exp(repr-sampled)") << std::endl;
        BOOST_FOREACH(double logv, tmp_logkernelvect)
            {
            std::cerr << boost::str(boost::format("%20.5f %20.5f %20.5f %20.5f") % _representative_logkernel % logv % (_representative_logkernel - logv) % exp(_representative_logkernel - logv) ) << std::endl;
            }
        //std::cerr << "mean log kernel value for this shell:" << _representative_logkernel << std::endl;
        //std::cerr << "minimum log kernel value for this shell:" << *std::min_element(tmp_logkernelvect.begin(), tmp_logkernelvect.end()) << std::endl;
        //std::cerr << "maximum log kernel value for this shell:" << *std::max_element(tmp_logkernelvect.begin(), tmp_logkernelvect.end()) << std::endl;
#endif

        }

    return num_samples_this_shell;
}

/*----------------------------------------------------------------------------------------------------------------------
|   Computes mean of the radial error distribution (Edmundson, 1959). This is the expected distance from the origin to
|   a Normal(0, sigma^2) deviate.
*/
double VarCovMatrix::radialErrorDistMean(double sigma) const
	{

    // calculate mean of radial error distribution
    double radial_error_dist_mean = exp(0.5*log(2.0) + boost::math::lgamma<double>((_p+1.0)*0.5) - boost::math::lgamma<double>(_p*0.5));
    return radial_error_dist_mean;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Computes _S and _Sinv from eigenvalues _w and eigenvectors _z, which are assumed to have already been calculated
|   from _V.
*/
void VarCovMatrix::calcSqrtVarCov()
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

	double * * dbl2d = NewTwoDArrayDouble(_p, _p);
    for (unsigned i = 0; i < _p; i++)
        {
        for (unsigned j = 0; j < _p; j++)
            {
            double tmp = 0.0;
            double tmp2 = 0.0;
            for (unsigned k = 0; k < _p; k++)
                {
                tmp += _z[i][k] * sqrt(_w[k]) * _z[j][k];
                tmp2 += _z[i][k] * pow(_w[k],-0.5) * _z[j][k];
                }
            _S[i][j] = tmp;
            _Stmp[i][j] = tmp2;
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
        // debug check: see if _Sinv equals _Stmp
        std::cerr << "\n_Sinv - _Stmp:" << std::endl;
        double * * I = NewTwoDArrayDouble(_p, _p);
        matrixMinusMatrix(_S, _Sinv, I);
        std::cerr << showMatrix(I, _p) << std::endl;

        // debug check: see if _S multipled by _Sinv equals I
        std::cerr << "\n_S * _Sinv:" << std::endl;
        double * * J = NewTwoDArrayDouble(_p, _p);
        matrixTimesMatrix(_S, _Sinv, J);
        std::cerr << showMatrix(J, _p) << std::endl;

        // debug check: see if _S multipled by _S equals _V
        double * * V = NewTwoDArrayDouble(_p, _p);
        double * * Z = NewTwoDArrayDouble(_p, _p);
        matrixTimesMatrix(_S, _S, V);
        matrixMinusMatrix(_V, V, Z);
        std::cerr << "\n_V - _S * _S:" << std::endl;
        std::cerr << showMatrix(Z, _p) << std::endl;

        DeleteTwoDArray<double>(I);
        DeleteTwoDArray<double>(J);
        DeleteTwoDArray<double>(V);
        DeleteTwoDArray<double>(Z);
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Add vector u and vector v and return the sum in vector u_plus_v.
*/
void VarCovMatrix::vectorPlusVector(double * u, double * v, double * u_plus_v)
    {
    for (unsigned i = 0; i < _p; ++i)
        {
            u_plus_v[i] = u[i] + v[i];
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Multiply matrix u and matrix v and return the product in matrix u_times_v.
*/
void VarCovMatrix::matrixTimesMatrix(double * * u, double * * v, double * * u_times_v)
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
|	Multiply _p X _p matrix u by _p-length vector v and return the product in _p-length vector u_times_v.
*/
void VarCovMatrix::matrixTimesVector(double * * u, double * v, double * u_times_v)
    {
    for (unsigned i = 0; i < _p; ++i)
        {
        double tmp = 0.0;
        for (unsigned j = 0; j < _p; ++j)
            {
            tmp += u[i][j]*v[j];
            }
        u_times_v[i] = tmp;
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Subtract matrix v from matrix u and return the product in matrix u_minus_v.
*/
void VarCovMatrix::matrixMinusMatrix(double * * u, double * * v, double * * u_minus_v)
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
|	Add matrix u to matrix u and return the product in matrix u_plus_v.
*/
void VarCovMatrix::matrixPlusMatrix(double * * u, double * * v, double * * u_plus_v)
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
|   Raise the matrix `V' to the power `power' and place the result in A, which should have already been allocated. The
|   eigenvectors and eigenvalues will be returned in the matrix `eigenvectors' and vector `eigenvalues', respectively.
*/
void VarCovMatrix::MatrixPow(double * * V, double power, unsigned p, double * * A, double * * eigenvectors, std::vector<double> & eigenvalues)
	{
    std::cerr << "\n################## A (start MatrixPow) ##################\n\n" << showMatrix(A, p) << std::endl;
    std::cerr << "\n################## V (start MatrixPow) ##################\n\n" << showMatrix(V, p) << std::endl;
    std::cerr << "\n################## eigenvectors (start MatrixPow) ##################\n\n" << showMatrix(eigenvectors, p) << std::endl;
    std::cerr << "\n################## eigenvalues (start MatrixPow) ##################\n\n" << showVector(eigenvalues, p) << std::endl;

    // Compute eigenvalues and eigenvectors
    std::vector<double> work(p, 0.0);
	int err_code = EigenRealSymmetric(p, V, &eigenvalues[0], eigenvectors, &work[0]);
	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues for matrix V in MatrixPow");
		}

    for (unsigned i = 0; i < p; i++)
        {
        for (unsigned j = 0; j < p; j++)
            {
            double tmp = 0.0;
            for (unsigned k = 0; k < p; k++)
                {
                tmp += eigenvectors[i][k] * pow(eigenvalues[k], power) * eigenvectors[j][k];
                }
            A[i][j] = tmp;
            }
        }

    std::cerr << "\n################## A (end MatrixPow) ##################\n\n" << showMatrix(A, p) << std::endl;
    std::cerr << "\n################## V (end MatrixPow) ##################\n\n" << showMatrix(V, p) << std::endl;
    std::cerr << "\n################## eigenvectors (end MatrixPow) ##################\n\n" << showMatrix(eigenvectors, p) << std::endl;
    std::cerr << "\n################## eigenvalues (end MatrixPow) ##################\n\n" << showVector(eigenvalues, p) << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Generate random data, compute mean vector and variance-covariance matrix, compute eigenvalues and eigenvectors of
|   variance-covariance matrix, then draw random points from concentric shells and compare random points with original
|   data points to demonstrate similarity.
*/
void VarCovMatrix::debugTestStandardization()
	{
    Lot rng;

    std::cerr << "########## Entering debugTestStandardization() ##########" << std::endl;

    // Random vector length is p
    // Sample size is n
    // Each of the p variables will be iid Gamma(shape,scale) (mean = shape*scale)
    unsigned p = 2;
    unsigned n = 1000;
    unsigned shape = 2; // must be positive integer
    double   scale = 1.0;
    bool     verbose = true;

    // Output settings and dimensions
    std::cerr << "### p = " << p << std::endl;
    std::cerr << "### n = " << n << std::endl;
    std::cerr << "### shape = " << shape << std::endl;
    std::cerr << "### scale = " << scale << std::endl;
    std::cerr << "### verbose = " << (verbose ? "yes" : "no") << std::endl;

    // Generate n p-dimensional random vectors of log-transformed Gamma(shape, scale) deviates
    std::vector< std::vector<double> > D(n);

    // Allocate memory for eigenvalues and eigenvectors
    std::vector<double> eigenvalues(p, 0.0);
    double * * eigenvectors = NewTwoDArrayDouble(p, p);
    eigenvectors[0][0] = 1.0;
    eigenvectors[0][1] = 0.0;
    eigenvectors[1][0] = 0.0;
    eigenvectors[1][1] = 1.0;

#if 1
    std::vector<double> true_mean(p, 10.0);
    std::cerr << "\n################## V ##################\n\n" << showVector(true_mean, p) << std::endl;
    double * * V = NewTwoDArrayDouble(p, p);
    assert(V != NULL);
    V[0][0] = 1.0;
    V[0][1] = 2.0;
    V[1][0] = 2.0;
    V[1][1] = 4.0;
    std::cerr << "\n################## V ##################\n\n" << showMatrix(V, p) << std::endl;

    double * * A = NewTwoDArrayDouble(p, p);
    A[0][0] = 0.0;
    A[0][1] = 0.0;
    A[1][0] = 0.0;
    A[1][1] = 0.0;

    MatrixPow(V, 0.5, p, A, eigenvectors, eigenvalues);

    std::cerr << "\n################## A ##################\n\n" << showMatrix(A, p) << std::endl;

    for (unsigned i = 0; i < n; ++i)
        {
        D[i].resize(p);

        // Generate a bivariate normal random deviate
        double z0 = rng.Normal();
        double z1 = rng.Normal();
        D[i][0] = true_mean[0] + A[0][0]*z0 + A[0][1]*z1;
        D[i][1] = true_mean[1] + A[1][0]*z0 + A[1][1]*z1;
        //std::cerr << "\n################## D[" << i << "] ##################\n\n" << showVector(D[i], p) << std::endl;
        }
#endif

#if 0
    for (unsigned i = 0; i < n; ++i)
        {
        D[i].resize(p);
        for (unsigned j = 0; j < p; ++j)
            {
            // Use fact that sum of shape Exponential(1/scale) deviates ~ Gamma(shape, scale)
            double v = 0.0;
            for (unsigned k = 0; k < shape; ++k)
                {
                double u = rng.Uniform();
                v += -scale*log(1.0 - u);
                }
            D[i][j] = log(v);
            }
        }
#endif

#if 0
    for (unsigned i = 0; i < n; ++i)
        {
        D[i].resize(p);
        for (unsigned j = 0; j < p; ++j)
            {
            // use log-transformed Exponential(1.0) deviates instead of Gamma(shape, scale) deviates
            double u = rng.Uniform();
            double v = -1.0*log(1.0 - u);
            D[i][j] = log(v);
            }
        }
#endif

    if (verbose)
        {
        // Output data matrix D
        std::cerr << "### original data matrix D: " << std::endl;
        for (unsigned i = 0; i < n; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % D[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        }

    // Compute sample mean vector
    std::vector<double> mean(p, 0.0);
    for (unsigned i = 0; i < n; ++i)
        {
        for (unsigned j = 0; j < p; ++j)
            {
            mean[j] += D[i][j];
            }
        }
    std::transform(mean.begin(), mean.end(), mean.begin(), boost::lambda::_1/n);

    if (verbose)
        {
        // Output vector of sample means
        std::cerr << "### sample means: " << std::endl;
        std::copy(mean.begin(), mean.end(), std::ostream_iterator<double>(std::cerr, " "));
        std::cerr << "\n" << std::endl;
        }

    // Center observations
    std::vector< std::vector<double> > C(n);
    for (unsigned i = 0; i < n; ++i)
        {
        C[i].resize(p);
        for (unsigned j = 0; j < p; ++j)
            {
            C[i][j] = D[i][j] - mean[j];
            }
        }

    if (verbose)
        {
        // Output centered data matrix C
        std::cerr << "### centered data matrix C: " << std::endl;
        std::vector<double> centered_means(p, 0.0);
        for (unsigned i = 0; i < n; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % C[i][j]);
                centered_means[j] += C[i][j];
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        std::transform(centered_means.begin(), centered_means.end(), centered_means.begin(), boost::lambda::_1/n);
        std::cerr << boost::str(boost::format("%6s") % "mean");
        for (unsigned j = 0; j < p; ++j)
            {
            std::cerr << boost::str(boost::format("%12.5f") % centered_means[j]);
            }
        std::cerr << "\n" << std::endl;
        }

    // Compute sample variance-covariance matrix S = (C' C)/(n-1)
    double * * S = NewTwoDArrayDouble(p, p);
    for (unsigned i = 0; i < p; i++)
        {
        for (unsigned j = 0; j < p; j++)
            {
            double xsum = 0.0;
            for (unsigned k = 0; k < n; k++)
                {
                xsum += C[k][i]*C[k][j];
                }
            S[i][j] = xsum/(n-1);
            }
        }

    if (verbose)
        {
        // Output variance-covariance matrix S
        std::cerr << "### variance-covariance matrix S: " << std::endl;
        std::cerr << "      ";
        for (unsigned i = 0; i < p; ++i)
            std::cerr << boost::str(boost::format("%12d") % i);
        std::cerr << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % S[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        }

    // Compute eigenvalues and eigenvectors
    std::vector<double> work(p, 0.0);
	int err_code = EigenRealSymmetric(p, S, &eigenvalues[0], eigenvectors, &work[0]);
	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues for matrix S in debugTestStandardization");
		}

    if (verbose)
        {
        // Output vector of eigenvalues
        std::cerr << "### eigenvalues: " << std::endl;
        std::copy(eigenvalues.begin(), eigenvalues.end(), std::ostream_iterator<double>(std::cerr, " "));
        std::cerr << "\n" << std::endl;

        // Output matrix of eigenvectors
        std::cerr << "### eigenvectors: " << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % eigenvectors[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        }

    // Compute sqrtS (square root of variance-covariance matrix) and sqrtSinv (inverse square root of var-cov matrix).
    //
    //         | z[0][0] z[0][1] | | sqrt(w[0])        0       | | z[0][0] z[1][0] |
    // sqrtS = |                 | |                           | |                 |
    //         | z[1][0] z[1][1] | |      0         sqrt(w[1]) | | z[0][1] z[1][1] |
    //
    //         | z[0][0] sqrt(w[0]) z[0][0] + z[0][1] sqrt(w[1]) z[0][1]     z[0][0] sqrt(w[0]) z[1][0] + z[0][1] sqrt(w[1]) z[1][1] |
    //       = |                                                                                                                     |
    //         | z[1][0] sqrt(w[0]) z[0][0] + z[1][1] sqrt(w[1]) z[0][1]     z[1][0] sqrt(w[0]) z[1][0] + z[1][1] sqrt(w[1]) z[1][1] |
    //
    // Note: z,w substituted above for eigenvectors,eigenvalues for clarity.
	double * * sqrtS = NewTwoDArrayDouble(p, p);
	double * * sqrtSinv = NewTwoDArrayDouble(p, p);
    for (unsigned i = 0; i < p; i++)
        {
        for (unsigned j = 0; j < p; j++)
            {
            double tmp1 = 0.0;
            double tmp2 = 0.0;
            for (unsigned k = 0; k < p; k++)
                {
                tmp1 += eigenvectors[i][k] * pow(eigenvalues[k], 0.5) * eigenvectors[j][k];
                tmp2 += eigenvectors[i][k] * pow(eigenvalues[k],-0.5) * eigenvectors[j][k];
                }
            sqrtS[i][j]    = tmp1;
            sqrtSinv[i][j] = tmp2;
            }
        }

    if (verbose)
        {
        // Output sqrtS matrix
        std::cerr << "### square root of variance-covariance matrix sqrtS: " << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % sqrtS[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;

        // Output sqrtSinv matrix
        std::cerr << "### inverse of square root of variance-covariance matrix sqrtSinv: " << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % sqrtSinv[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;

        // Output sqrtS * sqrtSinv = identity matrix
        std::cerr << "### sqrtS * sqrtSinv = I (identity matrix): " << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                double xsum = 0.0;
                for (unsigned k = 0; k < p; ++k)
                    {
                    xsum += sqrtS[i][k] * sqrtSinv[k][j];
                    }
                std::cerr << boost::str(boost::format("%12.5f") % xsum);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;

        // Output sqrtS * sqrtS = S (variance covariance matrix)
        std::cerr << "### sqrtS * sqrtS = S (variance-covariance matrix): " << std::endl;
        for (unsigned i = 0; i < p; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                double xsum = 0.0;
                for (unsigned k = 0; k < p; ++k)
                    {
                    xsum += sqrtS[i][k] * sqrtS[k][j];
                    }
                std::cerr << boost::str(boost::format("%12.5f") % xsum);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        }

    // Standardize ith centered observation vector by premultiplying by sqrtSinv
    //
    // | Zi0 |   | s00 s01 s02 | | Ci0 |   | s00 Ci0 + s01 Ci1 + s02 Ci2 |
    // |     |   |             | |     |   |                             |
    // | Zi1 | = | s10 s11 s12 | | Ci1 | = | s10 Ci0 + s11 Ci1 + s12 Ci2 |
    // |     |   |             | |     |   |                             |
    // | Zi2 |   | s20 s21 s22 | | Ci2 |   | s20 Ci0 + s21 Ci1 + s22 Ci2 |
    //                                        jk  ik    jk  ik    jk  ik
    std::vector< std::vector<double> > Z(n);
    for (unsigned i = 0; i < n; i++)
        {
        Z[i].resize(p, 0.0);
        for (unsigned j = 0; j < p; j++)
            {
            double xsum = 0.0;
            for (unsigned k = 0; k < p; k++)
                {
                xsum += sqrtSinv[j][k]*C[i][k];
                }
            Z[i][j] = xsum;
            }
        }

    if (verbose)
        {
        // Output standardized data matrix Z
        std::cerr << "### standardized data matrix Z: " << std::endl;
        for (unsigned i = 0; i < n; ++i)
            {
            std::cerr << boost::str(boost::format("%6d") % i);
            for (unsigned j = 0; j < p; ++j)
                {
                std::cerr << boost::str(boost::format("%12.5f") % Z[i][j]);
                }
            std::cerr << std::endl;
            }
        std::cerr << std::endl;
        }

    // Calculate radius for each observed point
    std::vector<double> radius(n, 0.0);
    double max_radius = 0.0;
    for (unsigned i = 0; i < n; i++)
        {
        double r = 0.0;
        for (unsigned j = 0; j < p; j++)
            {
            r += pow(Z[i][j], 2.0);
            }
        radius[i] = sqrt(r);
        if (radius[i] > max_radius)
            max_radius = radius[i];
        }

    if (verbose)
        {
        // Output radius for each point
        std::cerr << "### vector of radii (maximum = " << max_radius << "): " << std::endl;
        for (unsigned i = 0; i < n; ++i)
            {
            std::cerr << boost::str(boost::format("%6d %12.5f") % i % radius[i]) << std::endl;
            }
        std::cerr << std::endl;
        }

    // For each observed radius, make up a new point having that radius
    std::vector<double> random_vector(p, 0.0);
	double * * Q = NewTwoDArrayDouble(n, p);
    for (unsigned i = 0; i < n; ++i)
        {
        // Choose random vector of length radius[i] by normalizing and rescaling a vector of random normal deviates
        double sumsq = 0.0;
        for (unsigned j = 0; j < p; ++j)
            {
            double x = rng.Normal();
            random_vector[j] = x;
            sumsq += x*x;
            }
        double vector_length = sqrt(sumsq);
        std::transform(random_vector.begin(), random_vector.end(), random_vector.begin(), boost::lambda::_1*radius[i]/vector_length);

        // Destandardize using (sqrtS * random_vector + mean) and store in Q[i]
        for (unsigned j = 0; j < p; ++j)
            {
            double tmp = 0.0;
            for (unsigned k = 0; k < p; ++k)
                {
                tmp += sqrtS[j][k]*random_vector[k];
                }
            Q[i][j] = tmp + mean[j];
            }
        }

    std::cerr << "Saving plot data..." << std::endl;

    // Compare original to chosen points at different radii
    std::ofstream outf("standardization-test.R", std::ios::out);
    outf << "# p are the original points" << std::endl;
    outf << "# q are the simulated points" << std::endl;

    for (unsigned j = 0; j < p; ++j)
        {
        // Output the original points as p0, p1, ...
        outf << "p" << j << " <- c(" << D[0][j];
        for (unsigned i = 1; i < n; ++i)
            {
            outf << "," << D[i][j];
            }
        outf << ")\n";

        // Output the simulated points as q0, q1, ...
        outf << "q" << j << " <- c(" << Q[0][j];
        for (unsigned i = 1; i < n; ++i)
            {
            outf << "," << Q[i][j];
            }
        outf << ")\n";
        }

    // Output the radius of each point
    outf << "#radius <- c(" << radius[0];
    for (unsigned i = 1; i < n; ++i)
        {
        outf << "," << radius[i];
        }
    outf << ")\n";

    outf << "#plot(density(radius),col=\"black\")\n";
    outf << "#rug(radius)\n";
    outf << std::endl;
    outf << "#plot(density(p0),col=\"black\")\n";
    outf << "#lines(density(q0),col=\"black\")\n";
    outf << "#lines(density(p1),col=\"red\")\n";
    outf << "#lines(density(q1),col=\"red\")\n";
    outf << std::endl;
    outf << "plot(p0, p1, type=\"p\",col=\"black\")\n";
    outf << "points(q0, q1, col=\"red\")\n";
    outf << std::endl;
    outf.close();

    if (verbose)
        std::cerr << "Cleaning up..." << std::endl;

    // Clean up
    //DeleteTwoDArray<double>(V);
    DeleteTwoDArray<double>(Q);
    DeleteTwoDArray<double>(S);
    DeleteTwoDArray<double>(sqrtS);
    DeleteTwoDArray<double>(sqrtSinv);
    DeleteTwoDArray<double>(eigenvectors);

    std::cerr << "########## Leaving debugTestStandardization() ##########" << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes eigenvalues and eigenvectors corresponding to the following contrived variance-covariance matrix:
|
|       | 2  0  0 |
|   V = | 0  3  4 |
|       | 0  4  9 |
|
|   Eigenvalues: 
|
|      1    2     11
|   
|   Eigenvector matrix:    Normalized eigenvectors:
|
|      0    1      0        0.000000000  1.000000000 0.000000000
|      2    0      1        0.894427191  0.000000000 0.447213595
|     -1    0      2       -0.447213595  0.000000000 0.894427191
|
|   S = V^{0.5}
|
|   Eigenvalues:
|
|      1    1.414213562     3.31662479
|   
|   Eigenvector matrix:    Normalized eigenvectors:
|
|      0    1      0        0.000000000  1.000000000 0.000000000
|      2    0      1        0.894427191  0.000000000 0.447213595
|     -1    0      2       -0.447213595  0.000000000 0.894427191
|
|   Here is the output produced:
|
|   A:
|                           0            1            2
|              0      2.00000      0.00000      0.00000
|              1      0.00000      3.00000      4.00000
|              2      0.00000      4.00000      9.00000
|
|   Eigenvalues:
|   2 1 11
|   Eigenvector matrix:
|                           0            1            2
|              0      1.00000      0.00000      0.00000
|              1      0.00000     -0.89443      0.44721
|              2      0.00000      0.44721      0.89443|
|
|   See https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors for a detailed work through of this example.
*/
void VarCovMatrix::debugEigensystem() const
	{
    double * * A = NewTwoDArrayDouble(3, 3);
    A[0][0] = 2.0;
    A[0][1] = 0.0;
    A[0][2] = 0.0;
    A[1][0] = 0.0;
    A[1][1] = 3.0;
    A[1][2] = 4.0;
    A[2][0] = 0.0;
    A[2][1] = 4.0;
    A[2][2] = 9.0;

    std::cerr << "\n~~~~~~~~~ begin debugEigensystem ~~~~~~~~~" << std::endl;
    std::cerr << "A:\n" << showMatrix(A, 3) << std::endl;

    std::vector<double> fv(3, 0.0);
    std::vector<double> v(3, 0.0);
    double * * g = NewTwoDArrayDouble(3, 3);

	// Calculate eigenvalues (v) and eigenvectors (g) (fv is just workspace)
	int err_code = EigenRealSymmetric(3, A, &v[0], g, &fv[0]);

	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues for matrix A in debugEigensystem");
		}

    std::cerr << "Eigenvalues:" << std::endl;
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cerr, " "));
    std::cerr << "\nEigenvector matrix:\n" << showMatrix(g, 3) << std::endl;

#if 0
    // Now calculate the matrix B, square root of A
	double * * B = NewTwoDArrayDouble(3, 3);
    for (unsigned i = 0; i < 3; i++)
        {
        for (unsigned j = 0; j < 3; j++)
            {
            double tmp = 0.0;
            for (unsigned k = 0; k < 3; k++)
                {
                tmp += g[i][k] * sqrt(v[k]) * g[j][k];
                }
            B[i][j] = tmp;
            }
        }

	// Calculate eigenvalues (v) and eigenvectors (g) of B
	err_code = EigenRealSymmetric(3, B, &v[0], g, &fv[0]);

	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues for matrix B in debugEigensystem");
		}

    std::cerr << "Eigenvalues:" << std::endl;
    std::copy(v.begin(), v.end(), std::ostream_iterator<double>(std::cerr, " "));
    DeleteTwoDArray<double>(B);
#endif

    std::cerr << "\n~~~~~~~~~ end debugEigensystem ~~~~~~~~~" << std::endl;

    DeleteTwoDArray<double>(A);
    DeleteTwoDArray<double>(g);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes eigenvalues and eigenvectors corresponding to the variance-covariance matrix.
*/
void VarCovMatrix::calcEigensystem()
	{
    //std::cerr << "\n~~~~~~~~~ begin ~~~~~~~~~" << _p << std::endl;
    //std::cerr << "_p = " << _p << std::endl;
    //std::cerr << "_V:\n" << showMatrix(_V, _p) << std::endl;
    //std::cerr << "\n~~~~~~~~~ end ~~~~~~~~~" << _p << std::endl;

	// Calculate eigenvalues (_w) and eigenvectors (_z)
    std::vector<double> fv(_p, 0.0);
	int err_code = EigenRealSymmetric(_p, _V, &_w[0], _z, &fv[0]);

	if (err_code != 0)
		{
		throw XLikelihood("Error in the calculation of eigenvectors and eigenvalues of the variance-covariance matrix");
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the eigenvalues in the data member `_w' as a vector.
*/
std::vector<double> VarCovMatrix::getEigenValues() const
	{
	return _w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the sum of log eigenvalues of the variance-covariance matrix, which equals the log of the determinant of the
|   variance-covariance matrix.
*/
double VarCovMatrix::logDet() const
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
std::vector<double> VarCovMatrix::getEigenVectors() const
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, _z, _p);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats `V' nicely for output.
*/
std::string VarCovMatrix::showVector(std::vector<double> & V, unsigned p) const
	{
    //std::strstream ss;
    //std::copy(V.begin(), V.end(), std::ostream_iterator<double>(ss, " "));
    //return ss.str();
    std::string s;
    for (unsigned i = 0; i < p; ++i)
        s += boost::str(boost::format("%.5f ") % V[i]);
    return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats `V' nicely for output.
*/
std::string VarCovMatrix::showMatrix(double * * V, unsigned p) const
	{
	double * V_begin = &V[0][0];
	return showMatrixImpl(V_begin, p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats the matrix of eigenvectors nicely for output.
*/
std::string VarCovMatrix::showEigenVectorsAsMatrix() const
	{
	double * z_begin = &_z[0][0];
	return showMatrixImpl(z_begin, _p);
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Creates a string that, when printed, formats the matrix of eigenvectors nicely for output.
*/
std::string VarCovMatrix::showMatrixImpl(const double * q, unsigned p) const
	{
	unsigned i, j;

	// Output one column for row labels and a label for every column of _V
	std::string s = str(boost::format("%12s ") % " ");
	for (i = 0; i < p; ++i)
		{
		s += str(boost::format("%12d ") % i);
		}
	s += "\n";

	// Output rows of q
	double * pq = (double *)q;	//PELIGROSO
	for (i = 0; i < p; ++i)
		{
		s += str(boost::format("%12d ") % i);
		for (j = 0; j < p; ++j)
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
void VarCovMatrix::flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const
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
std::vector<double> VarCovMatrix::getVarCovMatrix() const
	{
	std::vector<double> p;
	flattenTwoDMatrix(p, _V, _p);
	return p;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies parameter types from `ptypes' to `_param_type' vector.
*/
//void VarCovMatrix::copyParamTypes(std::vector<unsigned> & ptypes)
//	{
//    _param_type.resize(ptypes.size());
//    std::copy(ptypes.begin(), ptypes.end(), _param_type.begin());
//	}

/*----------------------------------------------------------------------------------------------------------------------
|   Copies supplied newick tree description to `_newick'.
*/
//void VarCovMatrix::copyNewick(std::string newick)
//	{
//    _newick = newick;
//	}

/*----------------------------------------------------------------------------------------------------------------------
|   Debugging check: calculates log-likelihoods of all sampled trees for comparison with those from the original
|   parameter file.
*/
//void VarCovMatrix::debugCheckLogLikelihoods()
//	{
//	}

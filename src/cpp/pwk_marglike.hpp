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

#if ! defined(PWKMARGLIKE_HPP)
#define PWKMARGLIKE_HPP

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
class PWKMargLike
	{
	public:

										PWKMargLike(unsigned n, unsigned p);
										~PWKMargLike();

        void                            copyNewick(std::string newick);
        void                            copyParamTypes(std::vector<unsigned> & ptypes);
        void                            copyParamVector(unsigned row, std::vector<double> v);
        void                            estimateMargLike();

		unsigned						getNumSamples();
		unsigned						getNumParams();
        void                            calcEigensystem();
        std::vector<double>             getEigenValues();
		std::vector<double>             getEigenVectors();
        std::vector<double>             getVarCovMatrix();

		std::string						showMatrix(double * * V);
		std::string						showEigenVectorsAsMatrix();

        class PosteriorSample
        {
        public:
                                            PosteriorSample() : _log_likelihood(0.0), _log_prior(0.0), _log_jacobian(0.0), _radius(0.0) {}
                                            ~PosteriorSample() {}

            bool                            operator==(const PosteriorSample & other) const {return (_radius == other._radius);}
            bool                            operator<(const PosteriorSample & other) const  {return (_radius < other._radius);}

            double                          _log_likelihood;
            double                          _log_prior;
            double                          _log_jacobian;
            double                          _radius;
            std::vector<double>             _parameters;
        };

	private:

		std::string						showMatrixImpl(const double * q) const;
		void							flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const;
        void                            calcSqrtVarCov();
        void                            matrixMultiply(double * * u, double * * v, double * * u_times_v);
        void                            matrixAdd(double * * u, double * * v, double * * u_plus_v);
        void                            matrixSubtract(double * * u, double * * v, double * * u_minus_v);
        double                          logDet();

        void                            debugCheckLogLikelihoods();

	private:

		unsigned                            _n;		/**< Number of parameter vectors in `_posterior_samples' vector */
		unsigned                            _p;		/**< Number of parameters in _posterior_samples[row]._parameters */
		double * *                          _V;     /**< the variance covariance matrix stored by this object */
		double * *                          _S;     /**< the standard deviation matrix stored by this object */
		double * *                          _Sinv;  /**< the inverse of the standard deviation matrix stored by this object */
		double *                            _w;     /**< matrix of eigenvalues computed by EigenRealSymmetric */
		double * *                          _z;     /**< matrix of eigenvectors computed by EigenRealSymmetric */

        std::string                         _newick;            /**< newick tree description for the tree topology represented by this object */
        std::vector<double>                 _means;             /**< holds means of parameters for use in centering observations */
        std::vector<unsigned>               _param_type;        /**< holds type of each parameter: 1=edgelen, 2=exchangeability, 3=frequency, 4=shape, 5=pinvar, 6=subset rel. rate */
        std::vector<PosteriorSample>        _posterior_samples; /**< holds sorted, standardized, log-transformed parameter vectors */

	};
}	// namespace phycas

#endif


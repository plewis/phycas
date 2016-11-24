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
class VarCovMatrix
	{
	public:

										VarCovMatrix(std::string name, unsigned n, unsigned p);
										~VarCovMatrix();

        void                            setParamNames(std::vector<std::string> names);

        //void                            copyNewick(std::string newick);
        //void                            copyParamTypes(std::vector<unsigned> & ptypes);
        void                            copyParamVector(unsigned row, std::vector<double> v, double lnL, double lnP, double lnJ_edgelen, double lnJ_substmodel);
        std::vector<double>             destandardizeSample(unsigned row);
        double                          standardizeSamples();

		unsigned						getNumSamples() const;
		unsigned						getNumParams() const;
        void                            debugEigensystem() const;
        void                            debugTestStandardization();
        void                            calcEigensystem();
        std::vector<double>             getEigenValues() const;
		std::vector<double>             getEigenVectors() const;
        std::vector<double>             getVarCovMatrix() const;

        std::string                     showVector(std::vector<double> & V, unsigned p) const;
		std::string						showMatrix(double * * V, unsigned p) const;
		std::string						showEigenVectorsAsMatrix() const;

        double                          logJacobianForStandardization() const;

        double                          getLogRatioSum() const;
        double                          getRepresentativeLogKernel() const;
        double                          getRepresentativeLogLikelihood() const;
        double                          getRepresentativeLogPrior() const;
        double                          getRepresentativeLogJacobianEdgelen() const;
        double                          getRepresentativeLogJacobianSubstmodel() const;
        double                          getRepresentativeLogJacobianStandardization() const;

        std::vector<double>             getLogRatioVect() const;
        std::vector<double>             getRepresentativeParamVect() const;

        unsigned                        calcRepresentativeForShell(double r, double delta);

        double                          radialErrorDistMean(double sigma) const;

        // for debugging
        std::vector<double>                 getRepresentativesForShell(unsigned n, double r, double delta);

        class PosteriorSample
        {
        public:
                                            PosteriorSample() : _log_likelihood(0.0), _log_prior(0.0), _log_jacobian_edgelen(0.0), _log_jacobian_substmodel(0.0), _log_jacobian_standardization(0.0), _radius(0.0) {}
                                            ~PosteriorSample() {}

            bool                            operator==(const PosteriorSample & other) const {return (_radius == other._radius);}
            bool                            operator<(const PosteriorSample & other) const  {return (_radius < other._radius);}

            double                          _log_likelihood;
            double                          _log_prior;
            double                          _log_jacobian_edgelen;
            double                          _log_jacobian_substmodel;
            double                          _log_jacobian_standardization;
            double                          _radius;
            std::vector<double>             _parameters;
            std::vector<double>             _standardized;
        };

	private:

		std::string						showMatrixImpl(const double * q, unsigned p) const;
		void							flattenTwoDMatrix(std::vector<double> & p, double * * twoDarr, unsigned dim) const;
        void                            calcSqrtVarCov();
        void                            vectorPlusVector(double * u, double * v, double * u_plus_v);
        void                            matrixTimesVector(double * * u, double * v, double * u_times_v);
        void                            matrixTimesMatrix(double * * u, double * * v, double * * u_times_v);
        void                            matrixPlusMatrix(double * * u, double * * v, double * * u_plus_v);
        void                            matrixMinusMatrix(double * * u, double * * v, double * * u_minus_v);
        double                          calcVectorNorm(std::vector<double> & v) const;
        double                          calcSumTermsOnLogScale(std::vector<double> & v) const;

        double                          logDet() const;

        void                            MatrixPow(double * * V, double power, unsigned p, double * * A, double * * eigenvectors, std::vector<double> & eigenvalues);

        //void                            debugCheckLogLikelihoods();

	private:

		unsigned                            _n;		/**< Number of parameter vectors in `_posterior_samples' vector */
		unsigned                            _p;		/**< Number of parameters in _posterior_samples[row]._parameters */
		double * *                          _V;     /**< the variance covariance matrix stored by this object */
		double * *                          _S;     /**< the standard deviation matrix stored by this object */
		double * *                          _Stmp;  /**< the inverse of the standard deviation matrix stored by this object (alternate calculation) */
		double * *                          _Sinv;  /**< the inverse of the standard deviation matrix stored by this object */
		double * *                          _z;     /**< matrix of eigenvectors computed by EigenRealSymmetric */

        std::string                         _name;
		std::vector<std::string>            _param_names;

        double                              _shell_midpoint;                /**< midpoint of shell last specified in call to calcRepresentativeForShell */
        double                              _shell_delta;                   /**< half the thickness of shell last specified in call to calcRepresentativeForShell */

        double                              _log_ratio_sum;                 /**< holds sum (on log scale) of ratios of representative kernel to sample kernel */
		std::vector<double>                 _log_ratio_vect;                /**< vector of ratios (on log scale) of representative kernel to sample kernel (`_log_ratio_sum' holds the sum of the elements in this vector) */

        double                              _representative_logkernel;      /**< holds last representative log kernel value determined by calcRepresentativeForShell */
        double                              _representative_loglikelihood;  /**< holds last representative log likelihood value determined by calcRepresentativeForShell */
        double                              _representative_logprior;       /**< holds last representative log prior value determined by calcRepresentativeForShell */
        double                              _representative_logjacobian_edgelen;         /**< holds last representative log jacobian value determined by calcRepresentativeForShell */
        double                              _representative_logjacobian_substmodel;      /**< holds last representative log jacobian value determined by calcRepresentativeForShell */
        double                              _representative_logjacobian_standardization; /**< holds last representative log jacobian value determined by calcRepresentativeForShell */

		std::vector<double>                 _w;                             /**< vector of eigenvalues computed by EigenRealSymmetric */
        std::vector<double>                 _means;                         /**< holds means of parameters for use in centering observations */
        std::vector<double>                 _mode;                          /**< holds mode parameter vector for use in centering observations */
        std::vector<PosteriorSample>        _posterior_samples;             /**< holds sorted, standardized, log-transformed parameter vectors */
        std::vector<double>                 _representative_param_vector;   /**< holds a representative parameter vector that has been destandardized (see calcRepresentativeForShell) */

        bool                                _center_around_mode;            /**< if true, observations are centered around median; if false, observations are centered around mean */

	};
}	// namespace phycas

#endif


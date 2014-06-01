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

#if ! defined(CODON_MODEL_HPP)
#define CODON_MODEL_HPP

#include "states_patterns.hpp"
#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/lambda/lambda.hpp>
#include <vector>
#include <numeric>
#include <algorithm>
#include "xlikelihood.hpp"
#include "mcmc_param.hpp"
#include "basic_cdf.hpp"
#include "q_matrix.hpp"
#include "model.hpp"

namespace phycas{

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents a basic codon model. The parameter kappa is the
|	transition/transversion rate ratio, and the parameter omega is the nonsynonymous/synonymous rate ratio. This model
|	also provides parameters for estimating the 61 non-stop-codon frequencies.
*/
class Codon: public Model
	{
	public:
									Codon();
									~Codon()
										{
										//std::cerr << "Codon dying..." << std::endl;
										}

		virtual std::string			getModelName() const;

		virtual void				setStateFreqUnnorm(unsigned param_index, double value);
        virtual void                setStateFreqsUnnorm(const std::vector<double> & values);
		virtual void				setAllFreqsEqual();
		virtual void				setNucleotideFreqs(double freqA, double freqC, double freqG, double freqT);

		void						fixKappa();
		void						freeKappa();

		void						fixOmega();
		void						freeOmega();

		double						getKappa();
		void						setKappa(double k);

		double						getOmega();
		void						setOmega(double w);

		void						setKappaPrior(ProbDistShPtr d);
		ProbDistShPtr				getKappaPrior();

		void						setOmegaPrior(ProbDistShPtr d);
		ProbDistShPtr				getOmegaPrior();

        void						setStateFreqPrior(MultivarProbDistShPtr d);
		MultivarProbDistShPtr		getStateFreqPrior();

		void						setStateFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr				getStateFreqParamPrior();

		virtual std::string			paramHeader() const;
		virtual std::string			paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;

        virtual unsigned            getNumFreeParameters() const;
        virtual void                appendFreeParamNames(std::vector<std::string> & names, std::string prefix = "") const;
        virtual void                appendParamNames(std::vector<std::string> & names, std::string prefix = "") const;
        virtual void                appendUntransformedParamValues(std::vector<double> & values) const;
        virtual void                appendTransformedParamValues(std::vector<double> & values) const;
        virtual bool                setParamValueFromTransformed(std::string parameter_name, double transformed_value, TreeShPtr tree);

        virtual double              calcLogDetJacobian() const;

		void						updateQMatrix() const;
		virtual void				createParameters(JointPriorManagerShPtr jpm, TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterVect & edgelen_hyperparams, MCMCUpdaterVect & parameters, int subset_pos);
		double						calcUniformizationLambda() const;
        double					    calcLMat(double * * lMat) const;
        double					    calcUMat(double * * uMat) const;
		void						calcPMat(double * * pMat, double edgeLength) const;

protected:

	ProbDistShPtr				kappa_prior;		/**< The prior distribution governing kappa */
	ProbDistShPtr				omega_prior;		/**< The prior distribution governing omega */
	ProbDistShPtr				freq_param_prior;	/**< The prior distribution governing each frequency parameter */
	MultivarProbDistShPtr		freq_prior;	        /**< The prior distribution governing the vector of frequencies (used if frequencies are updated jointly by StateFreqMove) */

	string_vect_t				freq_label;			/**< Holds names of the 60 state frequencies (e.g. freqAAC, ...) used as headers in the param file */
	string_vect_t				freq_name;			/**< deprecated: Holds names of the 61 state frequencies (e.g. freqAAA, freqAAC, ...) used as headers in the param file */

	bool						kappa_fixed;		/**< If true, the value of kappa will not change during MCMC updates */
	bool						omega_fixed;		/**< If true, the value of omega will not change during MCMC updates */

	mutable MCMCUpdaterShPtr	omega_param;		/**< Copy of the omega parameter (saved so that fixed/free status can be changed) */
	mutable MCMCUpdaterShPtr	kappa_param;		/**< Copy of the kappa parameter (saved so that fixed/free status can be changed) */

	mutable QMatrix				q_matrix;			/**< A QMatrix object used to compute transition probabilities */

	// Below here are quantities that directly affect likelihood calculations and which should increment time_stamp when modified
	double						kappa;				/**< The transition/transversion rate ratio */
	double						omega;				/**< The nonsynonymous/synonymous rate ratio */
	};

typedef boost::shared_ptr<Codon> CodonShPtr;

} // namespace phycas

#endif


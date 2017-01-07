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

#if ! defined(HKY_MODEL_HPP)
#define HKY_MODEL_HPP

#include "states_patterns.hpp"
#include "model.hpp"
#include <boost/shared_ptr.hpp>

namespace phycas{

/*----------------------------------------------------------------------------------------------------------------------
|	Specialization of the base class Model that represents the Hasegawa-Kishino-Yano (1985) model.
*/
class HKY: public Model
	{
	public:
                                        HKY();
                                        ~HKY();

		virtual std::string             getModelName() const;
        virtual void                    createParameters(JointPriorManagerShPtr jpm, TreeShPtr t, MCMCUpdaterVect & edgelens, MCMCUpdaterVect & edgelen_hyperparams, MCMCUpdaterVect & parameters, int subset_pos);
		double                          calcUniformizationLambda() const;
        double                          calcLMat(double * * lMat) const;
        double                          calcUMat(double * * uMat) const;
        void                            calcPMat(double * * pMat, double edgeLength) const;

        void                            fixKappa();
		void                            freeKappa();

        double                          getKappa();
		void                            setKappa(double k);

        void                            setKappaPrior(ProbDistShPtr d);
		ProbDistShPtr                   getKappaPrior();

        void                            setKappaFromTRatio(double tratio);
		double                          calcTRatio();

		void                            setStateFreqPrior(MultivarProbDistShPtr d);
		MultivarProbDistShPtr           getStateFreqPrior();

		void                            setStateFreqParamPrior(ProbDistShPtr d);
		ProbDistShPtr                   getStateFreqParamPrior();

		virtual std::string             paramHeader() const;
		virtual std::string             paramReport(unsigned ndecimals, bool include_edgelen_hyperparams) const;

        virtual unsigned                getNumFreeParameters() const;
        virtual void                    appendPWKParamNames(std::vector<std::string> & names, std::string prefix = "") const;
        virtual void                    appendFreeParamNames(std::vector<std::string> & names, std::string prefix = "") const;
        virtual void                    appendParamNames(std::vector<std::string> & names, std::string prefix = "") const;
        virtual void                    appendUntransformedParamValues(std::vector<double> & values) const;
        virtual void                    appendTransformedParamValues(std::vector<double> & values) const;
        virtual bool                    setParamValueFromTransformed(std::string parameter_name, double transformed_value, TreeShPtr tree);
        virtual double                  calcLogDetJacobian() const;

protected:

	ProbDistShPtr				kappa_prior;		/**< The prior distribution governing kappa */
	ProbDistShPtr				freq_param_prior;	/**< The prior distribution governing each frequency parameter (used if frequencies are updated separately by slice sampling) */
	MultivarProbDistShPtr		freq_prior;	        /**< The prior distribution governing the vector of frequencies (used if frequencies are updated jointly by StateFreqMove) */
	bool						kappa_fixed;		/**< If true, the value of kappa will not change during MCMC updates */
	mutable MCMCUpdaterShPtr	kappa_param;		/**< Copy of the kappa parameter (saved so that fixed/free status can be changed) */

	string_vect_t				freq_name;				/**< Holds names of the 4 state frequencies (e.g. freqA, freqC, freqG, freqT) used as headers in the param file */

	// Below here are quantities that directly affect likelihood calculations and which should increment time_stamp when modified
	double						kappa;				/**< The transition/transversion rate ratio */
	};

typedef boost::shared_ptr<HKY> HKYShPtr;

}

#endif

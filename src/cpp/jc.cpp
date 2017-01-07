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

#include <cmath>
#include <iostream>
#include "model.hpp"
#include "jc.hpp"
#if defined(PYTHON_ONLY) && defined(USING_NUMARRAY)
#	include <boost/python/numeric.hpp>
#	include "phycas/src/thirdparty/num_util/num_util.h"
#endif
using std::cout;
using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets `num_states' data member to 4.
*/
JC::JC()
  : Model(4)
	{
	state_repr.reserve(4);
	state_repr.push_back("A");
	state_repr.push_back("C");
	state_repr.push_back("G");
	state_repr.push_back("T");
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a string indicating the name of this model, namely "JC69", "JC69+G, "JC69+I" or "JC69+G+I".
*/
std::string JC::getModelName() const
	{
	std::string s = "JC69";
	if (is_pinvar_model)
		s += "+I";
	if (num_gamma_rates > 1)
		s += "+G";
	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function provides the names of the columns that appear in the parameter file (i.e. the "*.p" file created by
|	MrBayes). All parameter files, regardless of the model, have "Gen", "LnL" and "TL" as the first three columns (for
|	compatability with MrBayes). The JC model provides additional column labels for the gamma shape parameter (if the
|	number of rates is greater than 1) and the pinvar parameter (if an invariable sites model is being used).
*/
std::string JC::paramHeader() const	/**< is the suffix to tack onto the parameter names for this model (useful for partitioned models to show to which partition subset the parameter belongs) */
	{
	return Model::paramHeader();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	The only parameters in the JC model are edge lengths, which are not reported by this function, but we must provide
|	an override of the pure virtual base class version.
*/
std::string JC::paramReport(
  unsigned ndecimals,						/**< floating point precision to use */
  bool include_edgelen_hyperparams) const	/**< if true, include values of edge length hyperparameters */
	{
	return Model::paramReport(ndecimals, include_edgelen_hyperparams);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the transition probability matrix given an edge length. Overrides the pure virtual function inherited from
|	the base class Model. For the JC69 model, the transition probabilities are:
|>
|	Pii = 0.25 + 0.75*exp{-4*beta*t}
|	Pij = 0.25 - 0.25*exp{-4*beta*t}
|>
|	The evolutionary distance `edgeLength' = 3*beta*t, so as a function of edge length, the transition probabilities
|	are:
|>
|	Pii = 0.25 + 0.75*exp{-4*edgeLength/3}
|	Pij = 0.25 - 0.25*exp{-4*edgeLength/3}
|>
*/
void JC::calcPMat(double * * pMat, double edgeLength) const
	{
	// The next two lines fix the "Jockusch" bug (see bug 8 in the BUGS file for details)
    if (edgeLength < 1.e-8)
        edgeLength = 1.e-8; //TreeNode::edgeLenEpsilon;
	const double exp_term = exp(-(4.0/3.0)*edgeLength);
	const double prob_change = 0.25 - 0.25*exp_term;
	const double prob_nochange = 0.25 + 0.75*exp_term;
	PHYCAS_ASSERT(prob_change >= 0.0);
	pMat[0][0] = prob_nochange;
	pMat[0][1] = prob_change;
	pMat[0][2] = prob_change;
	pMat[0][3] = prob_change;
	pMat[1][0] = prob_change;
	pMat[1][1] = prob_nochange;
	pMat[1][2] = prob_change;
	pMat[1][3] = prob_change;
	pMat[2][0] = prob_change;
	pMat[2][1] = prob_change;
	pMat[2][2] = prob_nochange;
	pMat[2][3] = prob_change;
	pMat[3][0] = prob_change;
	pMat[3][1] = prob_change;
	pMat[3][2] = prob_change;
	pMat[3][3] = prob_nochange;
	//for (unsigned i = 0 ; i < 4; ++i)
	//	cout << pMat[i][0] << ' '<< pMat[i][1] << ' '<< pMat[i][2] << ' '<< pMat[i][3] << '\n';
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Needs work.
*/
double JC::calcUniformizationLambda() const
	{
    return 1.1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes the uniformized transition probability matrix given an edge length. Overrides the pure virtual function
|   inherited from the base class Model. For the JC69 model, the uniformized transition probability matrix is:
|>
|   1 - 3*beta/lambda      beta/lambda         beta/lambda        beta/lambda
|     beta/lambda       1 - 3*beta/lambda      beta/lambda        beta/lambda
|     beta/lambda          beta/lambda      1 - 3*beta/lambda     beta/lambda
|     beta/lambda          beta/lambda         beta/lambda      1 - 3*beta/lambda
|>
|	where lambda is slightly larger than the largest element in the rate matrix (i.e. slightly larger than beta).
|   See Mateiu, L., and B. Rannala (2006. Systematic Biology 55:259-269) for details.
|   For efficiency, what is actually computed is the product of lambda and the above matrix, and what is stored are
|   the logarithms of the elements (because these values are always used on the log scale): i.e.
|>
|   log(lambda - 3*beta)      log(beta)             log(beta)             log(beta)
|       log(beta)         log(lambda - 3*beta)      log(beta)             log(beta)
|       log(beta)             log(beta)         log(lambda - 3*beta)      log(beta)
|       log(beta)             log(beta)             log(beta)         log(lambda - 3*beta)
|>
|   Because rate and time are confounded, beta is arbitrarily scaled so that time is measured in
|   expected number of substitutions (i.e. set beta such that the edge length v = t):
|>
|   v = 3*beta*t
|   v = 3*(1/3)*t
|>
|   Setting beta = 1/3 results in v = t.
*/
double JC::calcLMat(double * * lMat) const
	{
    // Mateiu and Rannala:   Me:
    //  B = 4/3              beta = 1/3
    //  Qii = -1             Qii = -3*beta = -1
    //  Qij = 1/3            Qij = beta    = 1/3
    //  Pii = 3/4            Pii = 1 - 3*beta/lambda = 1 - 1/4 = 3/4
    //  Pij = 1/12           Pij = beta/lambda = 1/12
    //  lambda = 4
    double beta = 1.0/3.0;
	const double lambda = calcUniformizationLambda();
	const double diag_term = log(lambda - 3.0*beta);
	const double off_diag_term = log(beta);
	lMat[0][0] = diag_term;
	lMat[0][1] = off_diag_term;
	lMat[0][2] = off_diag_term;
	lMat[0][3] = off_diag_term;
	lMat[1][0] = off_diag_term;
	lMat[1][1] = diag_term;
	lMat[1][2] = off_diag_term;
	lMat[1][3] = off_diag_term;
	lMat[2][0] = off_diag_term;
	lMat[2][1] = off_diag_term;
	lMat[2][2] = diag_term;
	lMat[2][3] = off_diag_term;
	lMat[3][0] = off_diag_term;
	lMat[3][1] = off_diag_term;
	lMat[3][2] = off_diag_term;
	lMat[3][3] = diag_term;
    return lambda;
	}

double JC::calcUMat(double * * uMat) const
	{
    double beta = 1.0/3.0;
	const double lambda = 1.1;
	const double diag_term = 1.0 - 3.0*beta/lambda;
	const double off_diag_term = beta/lambda;
	uMat[0][0] = diag_term;
	uMat[0][1] = off_diag_term;
	uMat[0][2] = off_diag_term;
	uMat[0][3] = off_diag_term;
	uMat[1][0] = off_diag_term;
	uMat[1][1] = diag_term;
	uMat[1][2] = off_diag_term;
	uMat[1][3] = off_diag_term;
	uMat[2][0] = off_diag_term;
	uMat[2][1] = off_diag_term;
	uMat[2][2] = diag_term;
	uMat[2][3] = off_diag_term;
	uMat[3][0] = off_diag_term;
	uMat[3][1] = off_diag_term;
	uMat[3][2] = off_diag_term;
	uMat[3][3] = diag_term;
    return lambda;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the free parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void JC::appendPWKParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    // JC does not have any parameters not handled by base class Model
    Model::appendPWKParamNames(names, prefix);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the free parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendTransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void JC::appendFreeParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    // JC does not have any parameters not handled by base class Model
    Model::appendFreeParamNames(names, prefix);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the names of the parameters of this model to the supplied vector `names' in the same order used when
|   transformed parameter values are appended in the member function appendUntransformedParamValues(). Derived classes
|   should override this function, calling this version before adding the names of parameters specific to the model
|   encapsulated by the derived class.
*/
void JC::appendParamNames(
  std::vector<std::string> & names, /**< is the vector to which the parameter names should be appended */
  std::string prefix) const         /**< is the prefix (e.g. partition subset number) that should be applied to each parameter name */
	{
    // JC does not have any parameters not handled by base class Model
    Model::appendParamNames(names, prefix);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log of the determinant of the log (or log-ratio) transformation used by the
|   appendTransformedParamValues() and setParamValueFromTransformed() methods. If g(theta*) is the density function of
|   the transformed values theta*, and f(theta) is the density of the untransformed values theta, then
|   g(theta*) = f(theta) |J|, where |J| is the Jacobian computed by this function.
*/
double JC::calcLogDetJacobian() const
	{
    // JC does not have any parameters not handled by base class Model
    return Model::calcLogDetJacobian();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the untransformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void JC::appendUntransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    // JC does not have any parameters not handled by base class Model
    Model::appendUntransformedParamValues(values);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the transformed values of the parameters of this model to the supplied vector `values' in the same order
|   used when parameter names are appended in the member function appendFreeParamNames(). Derived classes should override
|   this function, calling this version before adding the values of parameters specific to the model encapsulated by the
|   derived class.
*/
void JC::appendTransformedParamValues(
  std::vector<double> & values) const   /**< is the vector to which the parameter values should be appended */
	{
    // JC does not have any parameters not handled by base class Model
    Model::appendTransformedParamValues(values);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Untransform the value `transformed_value' and use it to set the value of the parameter whose name is supplied in
|   `parameter_name'. Returns true if parameter was found, false if not.
*/
bool JC::setParamValueFromTransformed(
  std::string parameter_name,   /**< is the name of the parameter to set */
  double transformed_value,     /**< is the transformed value of the parameter to set */
  TreeShPtr tree)               /**< is the tree */
	{
    // JC does not have any parameters not handled by base class Model
    return Model::setParamValueFromTransformed(parameter_name, transformed_value, tree);
    }

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2009 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#include "joint_prior_manager.hpp"
#include "partition_model.hpp"
//#include <boost/regex.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	PartitionModel constructor.
*/
PartitionModel::PartitionModel()
    {
    _subset_proportions = SubsetProportionsShPtr(new SubsetProportions());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	PartitionModel destructor.
*/
PartitionModel::~PartitionModel()
	{
	//std::cerr << "\n>>>>> PartitionModel dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Creates a JointPriorManager object and assigned it to `_joint_prior_manager'.
*/
void PartitionModel::createJointPriorManager()
    {
    PHYCAS_ASSERT(!_joint_prior_manager);   // should not already exist
    _joint_prior_manager = JointPriorManagerShPtr(new JointPriorManager());
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current number of subsets by returning the length of the `subset_model' vector.
*/
unsigned PartitionModel::getNumSubsets() const
    {
    // In debug mode check to make sure everyone agrees about the number of subsets
    PHYCAS_ASSERT(subset_num_patterns.size() == subset_model.size());
    PHYCAS_ASSERT(subset_num_states.size() == subset_model.size());
    PHYCAS_ASSERT(subset_num_rates.size() == subset_model.size());

    return (unsigned)subset_model.size();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of patterns (sum over all partition subsets).
*/
unsigned PartitionModel::getTotalNumPatterns() const
    {
    //@POL need to make total_num_patterns data member if this function is called often!
    return (unsigned)std::accumulate(subset_num_patterns.begin(), subset_num_patterns.end(), 0);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the total number of sites (sum over all partition subsets).
*/
unsigned PartitionModel::getTotalNumSites() const
    {
    return (unsigned)std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0);
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_model' vector.
*/
const ModelVect & PartitionModel::getModelsVect() const
	{
	return subset_model;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_relrates' vector.
*/
double PartitionModel::getSubsetRelRate(unsigned i) const
	{
	PHYCAS_ASSERT(subset_relrates.size() > i);
	return subset_relrates[i];
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_relrates' vector.
*/
const std::vector<double> & PartitionModel::getSubsetRelRatesVect() const
	{
    //unsigned nsites = std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0);
    //std::vector<double> doof(subset_num_sites.size(), 0.0);
    //std::transform(subset_num_sites.begin(), subset_num_sites.end(), subset_relrates.begin(), doof.begin(), boost::lambda::_1*boost::lambda::_2/double(nsites));
    //double check = std::accumulate(doof.begin(), doof.end(), 0.0);
    //std::cerr << boost::str(boost::format("~~~~~ %.15f ~~~~~~") % check) << std::endl;
    //if (fabs(check - 1.0) > .000001) {
    //    std::cerr << "yikes!" << std::endl;
    //    }
	return subset_relrates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_relrate_prior' data member.
*/
MultivarProbDistShPtr PartitionModel::getSubsetRelRatePrior() const
	{
	return subset_relrate_prior;
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Shows current subset relative rates, proportions and transformed relative rates as well as the weighted average of the
|   relative rates, which should always be 1.0.
*/
void PartitionModel::debugShowTransformedRelativeRates(const std::string msg) const
	{
	std::cerr << "~~~~~~~~~~ " << msg << " ~~~~~~~~~~" << std::endl;
	std::cerr << "Relative rates: ";
	std::copy(subset_relrates.begin(), subset_relrates.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << boost::str(boost::format(" sum = %.15f") % std::accumulate(subset_relrates.begin(), subset_relrates.end(), 0.0)) << std::endl;

    double num_sites_total = (double)std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0.0);
    std::vector<double> relrate_proportions(subset_relrates.size(), 0.0);
    std::transform(subset_num_sites.begin(), subset_num_sites.end(), relrate_proportions.begin(), boost::lambda::_1/num_sites_total);

	std::cerr << "Subset frequencies: ";
	std::copy(subset_num_sites.begin(), subset_num_sites.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << boost::str(boost::format(" sum = %d") % std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0)) << std::endl;

	std::cerr << "Subset proportions: ";
	std::copy(relrate_proportions.begin(), relrate_proportions.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << boost::str(boost::format(" sum = %.15f") % std::accumulate(relrate_proportions.begin(), relrate_proportions.end(), 0.0)) << std::endl;

	// transform relative rates using coefficients so that elements of orig_params sum to 1
	std::cerr << "Weighted relative rates: ";
    std::vector<double> tmp(subset_relrates.size(), 0.0);
	std::transform(subset_relrates.begin(), subset_relrates.end(), relrate_proportions.begin(), tmp.begin(), boost::lambda::_1*boost::lambda::_2);
	std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<double>(std::cerr, " "));
	std::cerr << boost::str(boost::format(" sum = %.15f") %  std::accumulate(tmp.begin(), tmp.end(), 0.0)) << std::endl;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the values in the data member vector `subset_relrates' with those in the supplied `rrates' vector. Note
|   that we cannot assume that the values supplied in `rrates' are normalized. To normalize, we need to solve for the
|   value `x' in the following example involving 2 subsets, one of which contains 4/5 of the sites and an unnormalized
|   relative rate of 1, with the other subset containing 1/5 of the sites and having an unnormalized relative rate of
|   10:
|>
|   (0.8)*(1/x) + (0.2)*(10/x) = 1.0  ==> x = (0.8)*(1) + (0.2)*(10) = 2.8
|>
*/
void PartitionModel::setSubsetRelRatesVect(
  const std::vector<double> & rrates)	/**< is the vector of relative rates */
	{
	PHYCAS_ASSERT(rrates.size() == getNumSubsets());
	PHYCAS_ASSERT(subset_relrates.size() == getNumSubsets());

    double num_sites_total = (double)std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0.0);

    std::vector<double> relrate_proportion;
    relrate_proportion.resize(rrates.size());
    std::transform(rrates.begin(), rrates.end(), subset_num_sites.begin(), relrate_proportion.begin(), (boost::lambda::_1)*(boost::lambda::_2)/num_sites_total);
    double x = (double)std::accumulate(relrate_proportion.begin(), relrate_proportion.end(), 0.0);

    std::transform(rrates.begin(), rrates.end(), subset_relrates.begin(), boost::lambda::_1/x);

    //std::stringstream sbefore;
    //std::copy(rrates.begin(), rrates.end(), std::ostream_iterator<double>(sbefore, " "));
    //std::cerr << sbefore.str() << std::endl;
    //std::stringstream safter;
    //std::copy(rrates.begin(), rrates.end(), std::ostream_iterator<double>(safter, " "));
    //std::cerr << safter.str() << std::endl;
    //std::cerr << debugShowTransformedRelativeRates("PartitionModel::setSubsetRelRatesVect") << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the multivariate probability distribution (`subset_relrate_prior') used for subset relative rates.
*/
void PartitionModel::setSubsetRelRatePrior(
  MultivarProbDistShPtr rrate_prior)	/**< is the new relative rate prior */
	{
	subset_relrate_prior = rrate_prior;
    //need to set coefficients for prior here
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_patterns' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumPatternsVect() const
	{
	return subset_num_patterns;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of sites in partition subset `i' by summing pattern counts over patterns in that subset.
*/
unsigned PartitionModel::getNumSites(unsigned i) const
	{
	return std::accumulate(subset_num_sites.begin(), subset_num_sites.end(), 0);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a shared pointer to a SubsetProportions object that can provide subset proportions.
*/
SubsetProportionsShPtr PartitionModel::getSubsetProportions()
	{
    return _subset_proportions;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_sites' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumSitesVect() const
	{
	return subset_num_sites;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_states' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumStatesVect() const
	{
	return subset_num_states;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns const reference to `subset_num_rates' vector.
*/
const std::vector<unsigned> & PartitionModel::getNumRatesVect() const
	{
	return subset_num_rates;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Adds model `m' to the `subset_model' vector.
*/
void PartitionModel::addModel(ModelShPtr m)
	{
	PHYCAS_ASSERT(m);
	subset_model.push_back(m);
	subset_relrates.push_back(1.0);
	subset_num_states.push_back(m->getNumStates());
	subset_num_rates.push_back(m->getNRatesTotal());
	subset_num_patterns.push_back(0);	// can't get this from model
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_model' to a copy of the supplied vector `models'.
*/
void PartitionModel::setModelsVect(const ModelVect & models)
	{
	unsigned new_size = (unsigned)models.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_model.resize(new_size);
	std::copy(models.begin(), models.end(), subset_model.begin());

	subset_relrates.assign(new_size, 1.0);

	// Go ahead and setup subset_num_rates and subset_num_states vectors
	// according to the models
	subset_num_rates.resize(new_size);
	subset_num_states.resize(new_size);
	subset_num_patterns.resize(new_size);
	for (unsigned i = 0; i < new_size; ++i)
		{
		subset_num_states[i] = subset_model[i]->getNumStates();
		subset_num_rates[i] = subset_model[i]->getNRatesTotal();
		subset_num_patterns[i] = 0;	// can't get this from model
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_patterns' to a copy of the supplied vector `npatterns'.
*/
void PartitionModel::setNumPatternsVect(const std::vector<unsigned> & npatterns)
	{
	unsigned new_size = (unsigned)npatterns.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_patterns.resize(new_size);
	std::copy(npatterns.begin(), npatterns.end(), subset_num_patterns.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_sites' to a copy of the supplied vector `nsites'. This results in an update of the data member
|   `_subset_proportions'.
*/
void PartitionModel::setNumSitesVect(const std::vector<unsigned> & nsites)
	{
	unsigned new_size = (unsigned)nsites.size();
	//PHYCAS_ASSERT(new_size > 0);
	subset_num_sites.resize(new_size);
	std::copy(nsites.begin(), nsites.end(), subset_num_sites.begin());

    _subset_proportions->setSubsetProportionsFromNumSites(nsites);

	//std::cerr << "\nDebug output: here are the number of sites in each subset:\n**********\n";
	//std::copy(subset_num_sites.begin(), subset_num_sites.end(), std::ostream_iterator<unsigned>(std::cerr, " "));
	//std::cerr << "\n**********\n" << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_states' to a copy of the supplied vector `nstates'.
*/
void PartitionModel::setNumStatesVect(const std::vector<unsigned> & nstates)
	{
	unsigned new_size = (unsigned)nstates.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_states.resize(new_size);
	std::copy(nstates.begin(), nstates.end(), subset_num_states.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `subset_num_rates' to a copy of the supplied vector `nrates'.
*/
void PartitionModel::setNumRatesVect(const std::vector<unsigned> & nrates)
	{
	unsigned new_size = (unsigned)nrates.size();
	PHYCAS_ASSERT(new_size > 0);
	subset_num_rates.resize(new_size);
	std::copy(nrates.begin(), nrates.end(), subset_num_rates.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `site_assignments' to a copy of the supplied vector `v'.
*/
void PartitionModel::setSiteAssignments(const std::vector<unsigned> & v)
	{
	unsigned new_size = (unsigned)v.size();
	//PHYCAS_ASSERT(new_size > 0);
	site_assignments.resize(new_size);
	std::copy(v.begin(), v.end(), site_assignments.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `site_assignments' to a copy of the supplied vector `v'.
*/
const std::vector<unsigned> & PartitionModel::getSiteAssignments() const
	{
    return site_assignments;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the number of free parameters in the model. Counts 3 parameters for nucleotide frequencies, 5 parameters for
|   GTR relative rates, etc., and, importantly, does not include any edge length parameters.
*/
unsigned PartitionModel::getNumFreeParameters() const
    {
    unsigned num_subsets = (unsigned)subset_model.size();
    unsigned num_free = num_subsets - 1;
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        num_free += subset_model[i]->getNumFreeParameters();
        }
    return num_free;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of free parameter names in the same order as the transformed parameter values returned by the member
|   function getTransformedParameters().
*/
std::vector<std::string> PartitionModel::getFreeParameterNames() const
    {
    free_param_names.clear();
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets > 1)
        {
        for (unsigned i = 1; i < num_subsets; ++i)
            {
            free_param_names.push_back(boost::str(boost::format("%d_subset_rate") % (i+1)));
            }
        }
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        std::string partition_number_prefix = boost::str(boost::format("%d_") % (i+1));
        subset_model[i]->appendFreeParamNames(free_param_names, partition_number_prefix);
        }
    return free_param_names;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of all parameter names in the same order as getFreeParameterNames() except for the insertion of
|   non-free parameter names in the appropriate places.
*/
std::vector<std::string> PartitionModel::getAllParameterNames() const
    {
    param_names.clear();
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets > 1)
        {
        for (unsigned i = 0; i < num_subsets; ++i)
            {
            param_names.push_back(boost::str(boost::format("%d_subset_rate") % (i+1)));
            }
        }
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        std::string partition_number_prefix = boost::str(boost::format("%d_") % (i+1));
        subset_model[i]->appendParamNames(param_names, partition_number_prefix);
        }
    return param_names;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the untransformed values of the subset relative rate parameters of this model to the supplied vector
|   `values'.
*/
void PartitionModel::appendUntransformedSubsetRelativeRates(std::vector<double> & values) const
    {
    unsigned num_subsets = (unsigned)subset_model.size();
    PHYCAS_ASSERT(num_subsets > 1);
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        double next = subset_relrates[i];
        PHYCAS_ASSERT(next > 0.0);
        values.push_back(next);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Appends the transformed values of the subset relative rate parameters of this model to the supplied vector `values'.
*/
void PartitionModel::appendTransformedSubsetRelativeRates(std::vector<double> & values) const
    {
    unsigned num_subsets = (unsigned)subset_model.size();
    PHYCAS_ASSERT(num_subsets > 1);
    double first = subset_relrates[0];
    PHYCAS_ASSERT(first > 0.0);
    double log_first = log(first);
    for (unsigned i = 1; i < num_subsets; ++i)
        {
        double next = subset_relrates[i];
        PHYCAS_ASSERT(next > 0.0);
        double log_ratio = log(next) - log_first;
        values.push_back(log_ratio);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log of the determinant of the log-ratio transformation used for the subset relative rates. If g(theta*)
|   is the density function of the transformed values theta*, and f(theta) is the density of the untransformed values
|   theta, then g(theta*) = f(theta) |J|, where |J| is the Jacobian computed by this function. In the case of the
|   subset relative rates, f(theta) is the subset relative rate distribution.
*/
double PartitionModel::calcLogDetJacobian() const
	{
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets < 2)
        return 0.0;

    double n = (double)getTotalNumSites();
    PHYCAS_ASSERT(n > 0.0);
    double log_n = log(n);

    double log_det_jacobian = log(subset_relrates[0]);
    double ni = (double)getNumSites(0);
    log_det_jacobian += log(ni);
    log_det_jacobian -= log_n;

    for (unsigned i = 1; i < num_subsets; ++i)
        {
        log_det_jacobian += log(subset_relrates[i]);
        ni = (double)getNumSites(i);
        log_det_jacobian += log_n;
        log_det_jacobian -= log(ni);
        }

    //std::cerr << "log_det_jacobian = " << log_det_jacobian << std::endl;

    return log_det_jacobian;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of untransformed parameter values, suitable for output to a parameter sample file, for example.
|   The identity of each parameter can be ascertained from the parameter name vector returned by the member function
|   getAllParameterNames().
*/
std::vector<double> PartitionModel::getUntransformedParameters() const
    {
    param_values.clear();
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets > 1)
        appendUntransformedSubsetRelativeRates(param_values);
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        subset_model[i]->appendUntransformedParamValues(param_values);
        }
    return param_values;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns a vector of transformed parameter values. Parameters on (0, infinity) are log-transformed, while state
|   frequencies, GTR exchangeabilities and other vectors that sum to 1.0 are log-ratio transformed. The identity of each
|   parameter can be ascertained from the parameter name vector returned by the member function getFreeParameterNames().
*/
std::vector<double> PartitionModel::getTransformedParameters() const
    {
    free_param_values.clear();
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets > 1)
        appendTransformedSubsetRelativeRates(free_param_values);
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        subset_model[i]->appendTransformedParamValues(free_param_values);
        }
    return free_param_values;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns the log of the determinant of the Jacobian for the transformation used in getTransformedParameters() and
|   setTransformedParameters().
*/
double PartitionModel::getLogDetJacobian() const
    {
    double log_det_jacobian = 0.0;
    unsigned num_subsets = (unsigned)subset_model.size();
    if (num_subsets > 1)
        log_det_jacobian += calcLogDetJacobian();
    for (unsigned i = 0; i < num_subsets; ++i)
        {
        log_det_jacobian += subset_model[i]->calcLogDetJacobian();
        }
    return log_det_jacobian;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Accessor function that simply returns the data member `_joint_prior_manager'.
*/
JointPriorManagerShPtr PartitionModel::getJointPriorManager()
    {
    return _joint_prior_manager;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Expects `v' to be a vector of transformed parameter values in which parameters on (0, infinity) are log-transformed,
|   while state frequencies, GTR exchangeabilities, and other vectors that sum to 1.0 are log-ratio transformed. The
|   identity of each parameter in `v' should match the order adopted by the parameter name vector returned by the member
|   function getFreeParameterNames().
*/
void PartitionModel::setTransformedParameters(
  const std::vector<double> v,
  TreeShPtr tree)
    {
    // subset relative rate variables
    PHYCAS_ASSERT(_subset_proportions);
    std::vector<double> const & p = _subset_proportions->getSubsetProportions();
    PHYCAS_ASSERT(p.size() > 0);
    double sum_exp_terms = p[0];

    unsigned num_subset_relrates = 0;

    getFreeParameterNames();
    unsigned sz = free_param_names.size();
    PHYCAS_ASSERT(sz == v.size());
    for (unsigned i = 0; i < sz; ++i)
        {
        std::string & nm = free_param_names[i];
        //std::cerr << "nm = " << nm << std::endl;

        // The part of nm before the '_' is the subset index
        size_t pos = nm.find_first_of('_');
        PHYCAS_ASSERT(pos != std::string::npos);
        unsigned subset_index = UINT_MAX;
        try
            {
            subset_index = boost::lexical_cast<int>(nm.substr(0, pos));
            }
        catch(boost::bad_lexical_cast &)
            {
            // signal that nm could not be converted to an integer value
            subset_index = UINT_MAX;
            }
        PHYCAS_ASSERT(subset_index < UINT_MAX);
        subset_index--; // subset indices used in parameter names start at 1

        // The part of nm after the '_' is the parameter name and v[i] is the transformed parameter value
        std::string this_param_name = nm.substr(++pos,std::string::npos);
        double this_param_value = v[i];
        if (this_param_name.compare("subset_rate") == 0)
            {
            // Set the subset relative rate for subset_index. Note that
            //
            //    param_j = log(rate_j/rate_0),
            //     rate_j = rate_0*exp(param_j)
            //
            // The value of this_param_value equals param_j, where j = subset_index. The value of rate_0 can be
            // obtained from param_1, param_2, ..., param_n as follows:
            //
            //    p_0*r_0   p_1*r_1   p_2*r_2         p_n*r_n      1
            //    ------- + ------- + ------- + ... + ------- = -------
            //    p_0*r_0   p_0*r_0   p_0*r_0         p_0*r_0   p_0*r_0
            //
            //        p_1                p_2                      p_n                   1
            //    1 + --- exp(param_1) + --- exp(param_2) + ... + --- exp(param_n) = -------
            //        p_0                p_0                      p_0                p_0*r_0
            //
            //    rate_0 = 1/(p_0 [1 + p_1*exp(param_1)/p_0 + ... + p_n*exp(param_n)/p_0])
            //           = 1/(p_0*exp(param_0) + p_1*exp(param_1) + p_2*exp(param_2) + ... + p_n*exp(param_n))
            //           = 1/(p_0 + p_1*exp(param_1) + p_2*exp(param_2) + ... + p_n*exp(param_n))
            //
            // For now, add p_j*exp(this_param_value) to sum_exp_terms and set subset_relrates[j] to exp(this_param_value),
            // then later correct the values in subset_relrates after all have been encountered.
            //
            double tmp = exp(this_param_value);
            sum_exp_terms += tmp*p[subset_index];
            num_subset_relrates += 1;
            subset_relrates[subset_index] = tmp;
            }
        else
            {
            // if parameter is not a subset relative rate, let the subset model handle it
            bool found = subset_model[subset_index]->setParamValueFromTransformed(this_param_name, this_param_value, tree);
            PHYCAS_ASSERT(found);
            }
        }

        // If there are subset relative rates, correct those now
        if (num_subset_relrates > 0)
            {
            num_subset_relrates += 1;   // 0th has not yet been counted
            subset_relrates[0] = 1.0;
            double r_0 = 1.0/sum_exp_terms;
            for (unsigned j = 0; j < num_subset_relrates; ++j)
                {
                subset_relrates[j] *= r_0;
                }
            if (_joint_prior_manager)
                _joint_prior_manager->multivariateModified("subset_relrates", subset_relrates);
            }
    }

//boost::regex expression("(\\d+)_(.+)");
//boost::cmatch what;
//bool regex_ok = boost::regex_match(nm.c_str(), what, expression);
//PHYCAS_ASSERT(regex_ok);
//if (regex_ok)
//    {
//    // what[0] contains the whole string
//    // what[1] contains the subset index
//    // what[2] contains the parameter name
//    // Construct a string using characters in nm from what[1].first to what[1].second,
//    // then cast it to an integer to get the subset index
//    int sid = boost::lexical_cast<int>(std::string(what[1].first, what[1].second));
//    subset_model[sid]->setParamValueFromTransformed(std::string(what[2].first, what[2].second), v[i]);
//    }

} // namespace phycas

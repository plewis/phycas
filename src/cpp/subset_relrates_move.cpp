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
#include <boost/cast.hpp>       //temporary!
#include <sstream>
#include <boost/math/special_functions/fpclassify.hpp>
#include "probability_distribution.hpp"
#include "relative_rate_distribution.hpp"
#include "model.hpp"
#include "basic_tree_node.hpp"
#include "tree_likelihood.hpp"
#include "xlikelihood.hpp"
#include "mcmc_chain_manager.hpp"
#include "dirichlet_move.hpp"
#include "subset_relrates_move.hpp"
#include "basic_tree.hpp"
#include "tree_manip.hpp"
#include "gtr.hpp"
#include "hky.hpp"

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	The constructor simply calls the base class (DirichletMove) constructor.
*/
SubsetRelRatesMove::SubsetRelRatesMove()
  : DirichletMove()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the relative rates of the associated partition model to those in the supplied vector `v'. No transformation of
|	the relative rates is done by this function.
*/
void SubsetRelRatesMove::sendCurrValuesToModel(const double_vect_t & v)
	{
	PHYCAS_ASSERT(dim == v.size());
	partition_model->setSubsetRelRatesVect(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, storing them in the supplied vector `v'. No transformation of
|	the relative rates is done by this function.
*/
void SubsetRelRatesMove::getCurrValuesFromModel(double_vect_t & v) const
	{
	const std::vector<double> & rrates = partition_model->getSubsetRelRatesVect();

    //std::cerr << "@#@#@#@#@# SubsetRelRatesMove::getCurrValuesFromModel" << std::endl;
    //std::copy(rrates.begin(), rrates.end(), std::ostream_iterator<double>(std::cerr, " "));
    //std::cerr << std::endl;
    //std::cerr << "@#@#@#@#@# SubsetRelRatesMove::getCurrValuesFromModel" << std::endl;

	PHYCAS_ASSERT(dim == rrates.size());
	v.resize(rrates.size());
	std::copy(rrates.begin(), rrates.end(), v.begin());
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model, returning them as an anonymous vector. No transformation of
|	the relative rates is done by this function.
*/
double_vect_t SubsetRelRatesMove::listCurrValuesFromModel()
	{
	double_vect_t v(dim);
	const double_vect_t & rrates = partition_model->getSubsetRelRatesVect();
	PHYCAS_ASSERT(dim == rrates.size());
	std::copy(rrates.begin(), rrates.end(), v.begin());
	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets `partition_model' to the supplied `m'.
*/
void SubsetRelRatesMove::setPartitionModel(
  PartitionModelShPtr m)    /*< is a shared pointer to the new partition model */
	{
	partition_model = m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains the current relative rates from the model by calling SubsetRelRatesMove::getCurrValuesFromModel.
*/
void SubsetRelRatesMove::getParams()
	{
	getCurrValuesFromModel(orig_relrates);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Replaces the relative rates in the model with those supplied in the vector `v'.
*/
void SubsetRelRatesMove::setParams(
  const double_vect_t & v)    /*< is the vector of parameter values to send to the model */
	{
	PHYCAS_ASSERT(v.size() == dim);
    sendCurrValuesToModel(v);
	}

/*--------------------------------------------------------------------------------------------------------------------------
|   Returns current subset relative rates, proportions and transformed relative rates as well as the weighted average of the
|   relative rates, which should always be 1.0.
*/
std::string SubsetRelRatesMove::debugShowTransformedRelativeRates(const std::string msg) const
	{
    std::stringstream ss;
	ss << ">>>>>>>>>> " << msg << " <<<<<<<<<<" << std::endl;
	ss << "Relative rates: ";
	std::copy(orig_relrates.begin(), orig_relrates.end(), std::ostream_iterator<double>(ss, " "));
	ss << boost::str(boost::format(" sum = %.15f") % std::accumulate(orig_relrates.begin(), orig_relrates.end(), 0.0)) << std::endl;

	ss << "Subset proportions: ";
    const std::vector<unsigned> & nsites = partition_model->getNumSitesVect();
    double total_sites = (double)partition_model->getTotalNumSites();
    double_vect_t p;
    p.resize(nsites.size());
    std::transform(nsites.begin(), nsites.end(), p.begin(), boost::lambda::_1/total_sites);
	std::copy(p.begin(), p.end(), std::ostream_iterator<double>(ss, " "));
	ss << boost::str(boost::format(" sum = %.15f") % std::accumulate(p.begin(), p.end(), 0.0)) << std::endl;

	// transform relative rates using coefficients so that elements of orig_params sum to 1
	ss << "Weighted relative rates: ";
    std::vector<double> tmp(orig_relrates.size(), 0.0);
	std::transform(orig_relrates.begin(), orig_relrates.end(), p.begin(), tmp.begin(), boost::lambda::_1*boost::lambda::_2);
	std::copy(tmp.begin(), tmp.end(), std::ostream_iterator<double>(ss, " "));
	ss << boost::str(boost::format(" sum = %.15f") %  std::accumulate(tmp.begin(), tmp.end(), 0.0)) << std::endl;
    return ss.str();
    }

/*--------------------------------------------------------------------------------------------------------------------------
|   Chooses a new vector of state frequencies using a sharp Dirichlet distribution centered at the original frequencies.
*/
void SubsetRelRatesMove::proposeNewState()
	{
	// transform relative rates using coefficients so that elements of orig_params sum to 1
	orig_params.resize(orig_relrates.size());
    SubsetProportionsShPtr ptr = partition_model->getSubsetProportions();
    std::vector<double> const & p = ptr->getSubsetProportions();
    //double total_sites = (double)partition_model->getTotalNumSites();
	std::transform(orig_relrates.begin(), orig_relrates.end(), p.begin(), orig_params.begin(), boost::lambda::_1*boost::lambda::_2);

    // create vector of Dirichlet parameters for selecting new relative rate values
	// The parameter vector for this temporary distribution is obtained by multiplying
    // each of the current relative rates by the value `psi', which is usually very large (e.g. 1000)
    c_forward.resize(orig_params.size());
	std::transform(orig_params.begin(), orig_params.end(), c_forward.begin(), 1.0 + boost::lambda::_1*psi);

	// create Dirichlet distribution and sample from it
    dir_forward = DirichletShPtr(new DirichletDistribution(c_forward));
    dir_forward->SetLot(getLot().get());
	new_params.resize(orig_params.size());
    new_params = dir_forward->Sample();

    // create vector of Dirichlet parameters for selecting old frequencies (needed for Hastings ratio calculation)
    c_reverse.resize(new_params.size());
	std::transform(new_params.begin(), new_params.end(), c_reverse.begin(), 1.0 + boost::lambda::_1*psi);
    dir_reverse = DirichletShPtr(new DirichletDistribution(c_reverse));

	// transform new_params using coefficients so that elements of new_relrates are equal to proposed new relative rates
	new_relrates.resize(new_params.size());
	std::transform(new_params.begin(), new_params.end(), p.begin(), new_relrates.begin(), boost::lambda::_1/boost::lambda::_2);

#if 0
    std::cerr << "\n^^^^^ orig_relrates ^^^^^" << std::endl;
    std::copy(orig_relrates.begin(), orig_relrates.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << "\n^^^^^ orig_params ^^^^^" << std::endl;
    std::copy(orig_params.begin(), orig_params.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << "\n^^^^^ c_forward ^^^^^" << std::endl;
    std::copy(c_forward.begin(), c_forward.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << "\n^^^^^ c_reverse ^^^^^" << std::endl;
    std::copy(c_reverse.begin(), c_reverse.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << "\n^^^^^ new_params ^^^^^" << std::endl;
    std::copy(new_params.begin(), new_params.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << "\n^^^^^ new_relrates ^^^^^" << std::endl;
    std::copy(new_relrates.begin(), new_relrates.end(), std::ostream_iterator<double>(std::cerr, "|"));
    std::cerr << std::endl;
#endif
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is accepted. Simply calls the reset() function.
*/
void SubsetRelRatesMove::accept()
	{
	MCMCUpdater::accept();
	reset();
	}

/*--------------------------------------------------------------------------------------------------------------------------
|	Called if the proposed move is rejected. Calls the reset() function after reinstating the original relative rates
|   and ensuring that all conditional likelihood arrays will be recalculated when the likelihood is next calculated.
*/
void SubsetRelRatesMove::revert()
	{
	MCMCUpdater::revert();
    setParams(orig_relrates);

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();
    jpm->multivariateModified(name, orig_relrates);
    curr_ln_prior = jpm->getLogJointPrior();

    // invalidate all CLAs
    likelihood->useAsLikelihoodRoot(NULL);
    likelihood->storeAllCLAs(tree);

	reset();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls proposeNewState(), then decides whether to accept or reject the proposed new state, calling accept() or
|	revert(), whichever is appropriate.
*/
bool SubsetRelRatesMove::update()
	{
	if (is_fixed)
		return false;

	ChainManagerShPtr p = chain_mgr.lock();
	PHYCAS_ASSERT(p);
    JointPriorManagerShPtr jpm = p->getJointPriorManager();

    // copy the current subset relative rates from the model to the data member orig_params
    getParams();

	proposeNewState();

    double prev_ln_prior		= jpm->getLogJointPrior();

	double prev_ln_like			= p->getLastLnLike();
	PHYCAS_ASSERT(!use_ref_dist || mv_ref_dist);
	double prev_ln_ref_dist = (use_ref_dist ? mv_ref_dist->GetLnPDF(orig_relrates) : 0.0);

    // replace current parameter values with new ones
    setParams(new_relrates);

    jpm->multivariateModified(name, new_relrates);
    double curr_ln_prior    = jpm->getLogJointPrior();

    likelihood->useAsLikelihoodRoot(NULL);	// invalidates all CLAs
	double curr_ln_like			= (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);
	double curr_ln_ref_dist = (use_ref_dist ? mv_ref_dist->GetLnPDF(new_relrates) : 0.0);

    double prev_posterior = 0.0;
	double curr_posterior = 0.0;

	if (is_standard_heating)
		{
		prev_posterior = heating_power*(prev_ln_like + prev_ln_prior);
		curr_posterior = heating_power*(curr_ln_like + curr_ln_prior);
		if (use_ref_dist)
			{
			prev_posterior += (1.0 - heating_power)*prev_ln_ref_dist;
			curr_posterior += (1.0 - heating_power)*curr_ln_ref_dist;
			}
		}
	else
		{
		prev_posterior = heating_power*prev_ln_like + prev_ln_prior;
		curr_posterior = heating_power*curr_ln_like + curr_ln_prior;
		}

	double ln_hastings			= getLnHastingsRatio();
	double ln_accept_ratio		= curr_posterior - prev_posterior + ln_hastings;

    double lnu = std::log(rng->Uniform());

    bool accepted = false;

	if (ln_accept_ratio >= 0.0 || lnu <= ln_accept_ratio)
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("ACCEPT, prev_ln_like = %.5f, getLastLnLike() = %.5f, curr_ln_like = %.5f, prev_ln_ref_dist = %.5f, curr_ln_ref_dist = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\n") % prev_ln_like % p->getLastLnLike() % curr_ln_like % prev_ln_ref_dist % curr_ln_ref_dist % prev_ln_prior % curr_ln_prior % lnu % ln_accept_ratio);
            debug_info += debugShowTransformedRelativeRates("Details");
			//for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
			//	debug_info += boost::str(boost::format("%15.8f ") % (*it));
			//debug_info += "\nnew_params:  ";
			//for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
			//	debug_info += boost::str(boost::format("%15.8f ") % (*it));
			//debug_info += "\nc_forward:  ";
			//for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
			//	debug_info += boost::str(boost::format("%15.8f ") % (*it));
			//debug_info += "\ndir_forward:  ";
			//debug_info += dir_forward->GetDistributionDescription();
			}
		p->setLastLnLike(curr_ln_like);

		accept();
		accepted = true;
		}
	else
		{
	    if (save_debug_info)
    	    {
			debug_info = boost::str(boost::format("REJECT, prev_ln_like = %.5f, getLastLnLike() = %.5f, curr_ln_like = %.5f, prev_ln_ref_dist = %.5f, curr_ln_ref_dist = %.5f, prev_ln_prior = %.5f, curr_ln_prior = %.5f, lnu = %.5f, ln_accept_ratio = %.5f\n") % prev_ln_like % p->getLastLnLike() % curr_ln_like % prev_ln_ref_dist % curr_ln_ref_dist % prev_ln_prior % curr_ln_prior % lnu % ln_accept_ratio);
            debug_info += debugShowTransformedRelativeRates("Details");
            //for (std::vector<double>::const_iterator it = orig_params.begin(); it != orig_params.end(); ++it)
            //    debug_info += boost::str(boost::format("%15.8f ") % (*it));
            //debug_info += "\nnew_params:  ";
            //for (std::vector<double>::const_iterator it = new_params.begin(); it != new_params.end(); ++it)
            //    debug_info += boost::str(boost::format("%15.8f ") % (*it));
            //debug_info += "\nc_forward:  ";
            //for (std::vector<double>::const_iterator it = c_forward.begin(); it != c_forward.end(); ++it)
            //    debug_info += boost::str(boost::format("%15.8f ") % (*it));
            //debug_info += "\ndir_forward:  ";
            //debug_info += dir_forward->GetDistributionDescription();
			}
		curr_ln_like	= p->getLastLnLike();
		revert();
		accepted = false;
		}

    //POLTMP
    double inverse_psi = 1.0/psi;
    inverse_psi = p->adaptUpdater(inverse_psi, nattempts, accepted);
    psi = 1.0/inverse_psi;

    return accepted;
	}


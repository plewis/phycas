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

#include "model.hpp"
#include "jc.hpp"
#include "hky.hpp"
#include "gtr.hpp"
#include "codon_model.hpp"
#include "tree_likelihood.hpp"
#include "probability_distribution.hpp"
#include "mcmc_param.hpp"
#include "mcmc_chain_manager.hpp"
#include "basic_tree.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and sets `_edge_len_index' to UINT_MAX.
*/
EdgeLenParam::EdgeLenParam()
  : MCMCUpdater(), _my_node(NULL), _edge_len_index(UINT_MAX)
	{
	curr_value = 0.01;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls the base class (MCMCUpdater) constructor and sets `_edge_len_index' to the supplied value `i'.
*/
EdgeLenParam::EdgeLenParam(
  unsigned i)               /**< is the index of this edge length */
  : MCMCUpdater(), _my_node(NULL), _edge_len_index(i)
	{
	curr_value = 0.01;
	has_slice_sampler = true;
	is_move = false;
	is_master_param = false;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
EdgeLenParam::~EdgeLenParam()
	{
	//std::cerr << "\n>>>>> EdgeLenParam dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the edge managed by this EdgeLenParam is an internal edge.
*/
bool EdgeLenParam::isInternalEdge()
	{
	bool is_internal = _my_node->IsInternal();
	bool is_subroot = _my_node->IsSubroot();
	bool is_internal_edge = is_internal && !is_subroot;

	//std::cerr << "~~~~~ _my_node number = " << _my_node->GetNodeNumber() << ", is_internal = " << (is_internal ? "yes" : "no") <<  ", is_subroot = " << (is_subroot ? "yes" : "no") <<  ", is_internal_edge = " << (is_internal_edge ? "yes" : "no") << std::endl;	//@@@

	return is_internal_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns true if the edge managed by this EdgeLenParam is an external edge.
*/
bool EdgeLenParam::isExternalEdge()
	{
	bool is_tip = _my_node->IsTip();
	bool is_subroot = _my_node->IsSubroot();
	bool is_external_edge = is_tip || is_subroot;

	//std::cerr << "~~~~~ _my_node number = " << _my_node->GetNodeNumber() << ", is_tip = " << (is_tip ? "yes" : "no") <<  ", is_subroot = " << (is_subroot ? "yes" : "no") <<  ", is_external_edge = " << (is_external_edge ? "yes" : "no") << std::endl;	//@@@

	return is_external_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the sample() member function of the `slice_sampler' data member.
*/
bool EdgeLenParam::update()
	{
	if (is_fixed)
		return false;

	if (_my_node->IsInternal())
		likelihood->useAsLikelihoodRoot(_my_node);
	else
		likelihood->useAsLikelihoodRoot(_my_node->GetParent());
	likelihood->invalidateAwayFromNode(*_my_node);

	double current_brlen = getCurrValueFromModel();

	//std::cerr << "<-- EdgeLenParam::update -->" << std::endl;	//@@@

	//double last_sampled_brlen = slice_sampler->GetLastSampledXValue();
	//if (fabs(current_brlen - last_sampled_brlen) > 0.00001)
	//	{
	//	std::cerr << "BAD BAD! BAD!!" << std::endl;
	//	}
	slice_sampler->SetXValue(current_brlen);
	slice_sampler->Sample();

	ChainManagerShPtr p = chain_mgr.lock();

    if (save_debug_info)
        {
        debug_info = str(boost::format("EdgeLenParam %f") % (slice_sampler->GetLastSampledXValue()));
        }
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to set the edge length to `v'.
*/
void EdgeLenParam::sendCurrValueToModel(double v)
	{
	PHYCAS_ASSERT(_my_node != NULL);
	_my_node->SetEdgeLen(v);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Overrides base class version to return the edge length currently stored in the tree.
*/
double EdgeLenParam::getCurrValueFromModel() const
	{
	PHYCAS_ASSERT(_my_node != NULL);
	return _my_node->GetEdgeLen();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets private data member `_my_node' to point to the supplied TreeNode `nd'.
*/
void EdgeLenParam::setTreeNode(TreeNode & nd)
	{
	_my_node = &nd;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function overrides the base class version to always returns true because this derived class implements a tree
|   length prior.
*/
bool EdgeLenParam::computesTreeLengthPrior() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns string representation of split associated with `_my_node'.
*/
std::string EdgeLenParam::getSplitReprAsString() const
	{
    Split & s = _my_node->GetSplit();
    if (s.IsBitSet(0))
        s.InvertSplit();
    return s.CreatePatternRepresentation();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	EdgeLenParam is a functor whose operator() returns a value proportional to the full-conditional posterior probability
|	density for a particular value of the edge length managed by this object. If the supplied value `v' is
|	out of bounds (i.e. <= 0.0), the return value is -DBL_MAX (closest we can come to a log posterior equal to negative
|	infinity).
*/
double EdgeLenParam::operator()(
  double v)	/**< is a new value for the edge length parameter */
	{
	PHYCAS_ASSERT(_my_node != NULL);
	curr_ln_like = ln_zero;
	curr_ln_prior = 0.0;

	if (v > 0.0)
		{
		sendCurrValueToModel(v);

		curr_ln_like = (heating_power > 0.0 ? likelihood->calcLnL(tree) : 0.0);

		//likelihood->startTreeViewer(tree, boost::str(boost::format("After calcLnL: curr_ln_like = %.6f") % curr_ln_like));

		ChainManagerShPtr p = chain_mgr.lock();
		PHYCAS_ASSERT(p);
		p->setLastLnLike(curr_ln_like);

        JointPriorManagerShPtr jpm = p->getJointPriorManager();
        if (jpm->isTreeLengthPrior())
            jpm->treeLengthModified("tree_length", tree);
        else
            jpm->univariateModified(name, v);
        curr_ln_prior = jpm->getLogJointPrior();

		if (is_standard_heating)
			if (use_ref_dist)
				{
                double curr_ln_ref_dist = 0.0;
                if (likelihood->getTreeLengthRefDist())
                    {
                    curr_ln_ref_dist = likelihood->getTreeLengthRefDist()->GetLnPDF(tree);
                    }
                else
                    {
                    PHYCAS_ASSERT(ref_dist);
                    curr_ln_ref_dist = ref_dist->GetLnPDF(v);
                    }
                return heating_power*(curr_ln_like + curr_ln_prior) + (1.0 - heating_power)*curr_ln_ref_dist;
				}
			else
				return heating_power*(curr_ln_like + curr_ln_prior);
		else
			return heating_power*curr_ln_like + curr_ln_prior;
		}
    else
        return ln_zero;
	}

}	// namespace phycas

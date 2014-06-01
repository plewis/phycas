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

#include <numeric>										// for std::accumulate
#include "basic_cdf.hpp"
#include "mcmc_param.hpp"
#include "basic_tree.hpp"				// for Tree::begin() and Tree::end()
#include "mcmc_chain_manager.hpp"
#include "tree_likelihood.hpp"
#include <boost/format.hpp>

// these were at the top of basic_tree.inl
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor sets both `has_slice_sampler' and `is_move' to false. It also sets the value of `edgeLenType' data
|   member, which can be `internal', `external' or `both'. If `edgeLenType' is `internal', this EdgeLenMasterParam will
|   only calculate the prior for internal edges. If `edgeLenType' is `external', this EdgeLenMasterParam will only
|   calculate the prior for external edges. If `edgeLenType' is `both', this EdgeLenMasterParam will calculate the
|   prior for all edges in the tree.
*/
EdgeLenMasterParam::EdgeLenMasterParam(
  EdgeLenMasterParam::EdgeLenType t)    /**> is the edge length type (internal, external or both) */
  : MCMCUpdater(), edgeLenType(t), use_edge_specific_ref_dists(false)
	{
	has_slice_sampler = false;
	is_move = false;
	is_master_param = true;
	is_hyper_param = false;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor.
*/
EdgeLenMasterParam::~EdgeLenMasterParam()
	{
	//std::cerr << "\n>>>>> EdgeLenMasterParam dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	This function overrides the base class version to always returns true because this derived class implements a tree
|   length prior.
*/
bool EdgeLenMasterParam::computesTreeLengthPrior() const
	{
	return true;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	If `use_it' is true, then split-specific edge length working priors will be constructed for all edges seen during
|	an MCMC analysis. If `use_it' is false, then a single generic edge length working prior will be used for all edges.
*/
void EdgeLenMasterParam::useEdgeSpecificWorkingPriors(
  bool use_it)	/**< should be true if you wish to use split-specific edge length working priors */
	{
	use_edge_specific_ref_dists = use_it;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the probability density evaluated at the current edge length associated with `nd'.
*/
double EdgeLenMasterParam::lnPriorOneEdge(TreeNode & nd) const
	{
    bool skip = (nd.IsTipRoot())
                || ((edgeLenType == EdgeLenMasterParam::internal) && (!nd.IsInternal()))
                || ((edgeLenType == EdgeLenMasterParam::external) && (nd.IsInternal()));
	if (skip)
        {
		return 0.0;
        }
    else
        {
        double v = nd.GetEdgeLen();

	    double retval = 0.0;
	    try
		    {
		    retval = prior->GetLnPDF(v);
			//std::cerr << "OoOoOoOoO " << nd.GetNodeNumber() << " --> " << v << " --> " << retval << std::endl;	//@@@
		    }
	    catch(XProbDist &)
		    {
		    PHYCAS_ASSERT(0);
		    }
	    return retval;
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists to facilitate use of boost::lambda::bind. It returns the length of node `nd', unless the
|   node is not of the correct type (e.g. internal node when this EdgeLenMasterParam is assigned to external nodes, or
|   vice versa, or node is the tip root node, in which case it has no edge.
*/
double EdgeLenMasterParam::edgeLength(TreeNode & nd) const
	{
    bool skip = (nd.IsTipRoot())
                || ((edgeLenType == EdgeLenMasterParam::internal) && (!nd.IsInternal()))
                || ((edgeLenType == EdgeLenMasterParam::external) && (nd.IsInternal()));
	if (skip)
		return 0.0;
    else
	    return nd.GetEdgeLen();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Reports the distribution being used for both the generic `ref_dist' as well as the edge-specific working
|	priors for edges that were sampled before the working priors were finalized.
*/
std::string EdgeLenMasterParam::getWorkingPriorDescr() const
	{
	std::string s;
	std::string typestr = "all";
	if (edgeLenType == internal)
		typestr = "internal";
	else if (edgeLenType == external)
		typestr = "external";

	if (prior && ref_dist)
		{
		s += boost::str(boost::format("%s will be used for lengths of %s edges") % ref_dist->GetDistributionDescription() % typestr);
		}
	else if (mv_prior && mv_ref_dist)
		{
		s += boost::str(boost::format("%s will be used for lengths of %s edges") % mv_ref_dist->GetDistributionDescription() % typestr);
		}
	else
		{
		s += "no generic edge working prior exists at this point";
		}

	if (use_edge_specific_ref_dists)
		{
		s += "\n    Here are the working priors that will be used for edges already seen:";
		for (WorkingPriorMapConstIter it = edge_ref_dist.begin(); it != edge_ref_dist.end(); ++it)
			{
			if (it->second.wp)
				{
				s += "\n      ";
				s += it->second.wp->GetDistributionDescription();
				s += " will be used for ";
				s += it->first.CreateNewickRepresentation();
				}
			}
		}

	return s;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Member function that exists only to facilitate using boost::lambda::bind to be used in the getLnPrior() function. It
|	returns the log of the working prior probability density evaluated at the current edge length associated with `nd'.
*/
double EdgeLenMasterParam::lnWorkingPriorOneEdge(const TreeNode & nd, double v) const
	{
    bool skip = (nd.IsTipRoot())
                || ((edgeLenType == EdgeLenMasterParam::internal) && (!nd.IsInternal()))
                || ((edgeLenType == EdgeLenMasterParam::external) && (nd.IsInternal()));
	if (skip)
        {
		return 0.0;
        }
    else
        {
		double retval = 0.0;

		bool use_generic = true;
		if (use_edge_specific_ref_dists)
			{
			const Split & s = nd.GetSplitConst();
			if (edge_ref_dist.find(s) != edge_ref_dist.end())
				{
				// Found edge-specific working prior for this edge
				WorkingPriorMapConstIter ewp = edge_ref_dist.find(s);
				if ((*ewp).second.wp)
					{
					try
						{
						// edge_ref_dist[s] retrieves EdgeWorkingPrior struct for this split, of
						// which 'second' is the working prior distribution shared pointer
						retval = (*ewp).second.wp->GetLnPDF(v);
						}
					catch(XProbDist &)
						{
						PHYCAS_ASSERT(0);
						}
					use_generic = false;
					}
				}
			}
		if (use_generic)
			{
			// There is no edge-specific working prior for this edge, so use generic one
			PHYCAS_ASSERT(ref_dist);
			try
				{
				retval = ref_dist->GetLnPDF(v);
				}
			catch(XProbDist &)
				{
				PHYCAS_ASSERT(0);
				}
			}

		return retval;
		}
	}


/*----------------------------------------------------------------------------------------------------------------------
|	Loops through all nodes in the tree and computes the log of the working prior for each edge that it is responsible
|	for (according to its `edgeLenType'). Returns sum of these log working prior values.
*/
double EdgeLenMasterParam::recalcWorkingPrior() const
	{
	if (isFixed() || !isPriorSteward())
		return 0.0;

	double lnwp = 0.0;
	if (!isFixed())
		{
		for (preorder_iterator nd = tree->begin(); nd != tree->end(); ++nd)
			{
			lnwp += lnWorkingPriorOneEdge(*nd, nd->GetEdgeLen());
			}
		}
	return lnwp;
	}

}

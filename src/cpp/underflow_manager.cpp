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

#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <limits>
#include "underflow_manager.hpp"

//double tmp_max_log_f = 0.0;

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor.
*/
UnderflowManager::UnderflowManager()
: total_patterns(0), underflow_num_edges(1), underflow_max_value(10000.0)
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of the data member `underflow_max_value', which is the target value to which the largest
|	conditional likelihood for each pattern is set.
*/
double UnderflowManager::getUnderflowMaxValue() const
	{
	return underflow_max_value;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of the data member `underflow_num_edges' to `nedges'. This is the number of edges to accumulate before
|	correcting for underflow. Assumes `nedges' is greater than zero. Often several hundred taxa are required before
|	underflow becomes a problem, so a reasonable value for `underflow_num_edges' is 25.
*/
void UnderflowManager::setTriggerSensitivity(
  unsigned nedges)		/**< is the number of edges to traverse before correcting for underflow */
	{
	PHYCAS_ASSERT(nedges > 0);
	underflow_num_edges = nedges;
    std::cerr << "********** POL: setting underflow_num_edges = " << nedges << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets value of data member `underflow_max_value' to `maxval'. The largest conditional likelihood for every pattern
|	will be adjusted to a value close to this. Suppose the largest conditional likelihood for pattern i was x. A factor
|	exp(f) is found such that x*exp(f) = `underflow_max_value'. The fractional component of f is then removed, yielding
|	an unsigned int k, which is what is actually stored. The largest conditional likelihood for pattern i after the
|	correction is thus x*exp(k), which will be less than or equal to `underflow_max_value'. Assumes `maxval' is greater
|	than zero. A reasonable value for `underflow_max_value' is 10000.
*/
void UnderflowManager::setCorrectToValue(
  double maxval)	/**< is the target value to which the largest conditional likelihood for any given pattern will be scaled */
	{
	PHYCAS_ASSERT(maxval > 0.0);
	underflow_max_value = maxval;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the `num_rates', `num_patterns', and `num_states' vectors to `nr', `np' and `ns', respectively. Assumes
|	that the number of rates and states are greater than zero for every partition subset.
*/
void UnderflowManager::setDimensions(
  const uint_vect_t & np,	/**< is a vector of the number of patterns in the data for each partition subset */
  const uint_vect_t & nr, 	/**< is a vector of the number of relative rates used in modeling among-site rate variation for each partition subset */
  const uint_vect_t & ns)	/**< is a vector of the number of states for each partition subset */
	{
	// Note: num_patterns can legitimately be 0 if running with no data. In this case, no underflow
    // correction is ever needed, and all member functions are no-ops
    total_patterns = (unsigned)std::accumulate(np.begin(), np.end(), 0);
	PHYCAS_ASSERT(std::accumulate(nr.begin(), nr.end(), 0) > 0);
	PHYCAS_ASSERT(std::accumulate(ns.begin(), ns.end(), 0) > 0);

	unsigned sz = (unsigned)np.size();
	PHYCAS_ASSERT(nr.size() == (unsigned)sz);
	PHYCAS_ASSERT(ns.size() == (unsigned)sz);

	num_patterns.resize(sz);
	num_rates.resize(sz);
	num_states.resize(sz);

	for (unsigned i = 0; i < sz; ++i)
		{
		num_patterns[i]	= np[i];
		num_rates[i]	= nr[i];
		num_states[i]	= ns[i];
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of a node subtending two tips. In this case, the number of edges traversed is just two, and thus no
|	corrective action is taken other than to set the number of underflow edges traversed to 2. The zeroUF member
|	function of the supplied `cond_like' object is called, however, to ensure that the cumulative correction factor is
|	zero (this `cond_like' might have just been brought in from storage with some correction factor already in place).
*/
void UnderflowManager::twoTips(
  CondLikelihood & cond_like) /**< is the conditional likelihood array to correct */
  const
	{
    if (total_patterns > 0)
        {
	    cond_like.setUnderflowNumEdges(2);
	    cond_like.zeroUF();
        }
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then subtracts the value stored in that underflow array
|	for pattern `pat' from the site likelihood `site_like'.
*/
double UnderflowManager::correctSiteLike(
  double & site_like,						/**< is the log-site-likelihood to correct */
  unsigned pat,								/**< is the index of the pattern representing the site to correct */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	double site_uf_factor = 0.0;
    if (total_patterns > 0)
        {
	    PHYCAS_ASSERT(condlike_shptr);
	    UnderflowType const * uf = condlike_shptr->getUF();
	    PHYCAS_ASSERT(uf != NULL);
	    site_uf_factor = (double)uf[pat];
	    site_like -= site_uf_factor;
        }
	return site_uf_factor;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Obtains a pointer to the underflow array for `cond_like', then returns the value stored in that underflow array
|	for pattern `pat'.
*/
double UnderflowManager::getCorrectionFactor(
  unsigned pat,								/**< is the index of the pattern for which the correction factor is desired */
  ConstCondLikelihoodShPtr condlike_shptr)	/**< is the conditional likelihood array of the likelihood root node */
  const
	{
	double site_uf_factor = 0.0;
    if (total_patterns > 0)
        {
	    PHYCAS_ASSERT(condlike_shptr);
	    UnderflowType const * uf = condlike_shptr->getUF();
	    PHYCAS_ASSERT(uf != NULL);
	    site_uf_factor = (double)uf[pat];
        }
	return site_uf_factor;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Handles case of an internal node (`cond_like') with two internal node children (`left_cond_like' not equal to
|	`right_cond_like') or the case of an internal node (`cond_like') with one tip child and one internal node child
|	(in which case `left_cond_like' will equal `right_cond_like'). The `numEdgesSinceUnderflowProtection' data member of
|	`cond_like' is set to 2 plus the number of traversed edges stored by `left_cond_like' plus (if there are two
|	internal node children) the number of traversed edges stored by `right_cond_like'. If this new number of edges
|	tranversed is greater than or equal to `underflow_num_edges', `cond_like' is corrected for underflow.
*/
void UnderflowManager::check(
  CondLikelihood &       cond_like,			/**< the conditional likelihood array object of the focal internal node */
  const CondLikelihood & left_cond_like,	/**< the conditional likelihood array object of an internal node that is one immediate descendant of the focal node */
  const CondLikelihood & right_cond_like, 	/**< the conditional likelihood array object of an internal node that is the other immediate descendant of the focal node (if one descendant is a tip, left_cond_like and right_cond_like should refer to the same object) */
  const count_vect_t & counts,				/**< the pattern counts vector */
  bool polytomy)							/**< true if in the process of dealing with additional children (beyond first two) in a polytomy */
  const
	{
    if (total_patterns == 0)
        return;

    // Determine whether we are dealing with
    // case 1: one tip child and one internal node child
    // case 2: two internal node children (= no_tips)
    // case 3: an extra tip (can only happen in case of polytomy)
    // case 4: an extra internal node (can only happen in case of polytomy)
    unsigned which = 0;
    const unsigned char one_tip = 1;
    const unsigned char no_tips = 2;
    const unsigned char xtra_tip = 3;
    const unsigned char xtra_internal = 4;
    if (&left_cond_like == &right_cond_like)
    	which = (polytomy ? xtra_tip : one_tip);
    else
    	which = (polytomy ? xtra_internal : no_tips);

	unsigned nedges = 0;
	if (which == one_tip)
		{
		// cond_like, left_cond_like == right_cond_like
		nedges += 2 + left_cond_like.getUnderflowNumEdges();
		}
	else if (which == no_tips)
		{
		// cond_like, left_cond_like, right_cond_like
		nedges += 2 + left_cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}
	else if (which == xtra_tip)
		{
		// cond_like == left_cond_like == right_cond_like
		nedges += 1 + cond_like.getUnderflowNumEdges();
		}
	else	// xtra_internal
		{
		// cond_like == left_cond_like, right_cond_like
		nedges += 1 + cond_like.getUnderflowNumEdges() + right_cond_like.getUnderflowNumEdges();
		}

	unsigned nsubsets = (unsigned)num_patterns.size();
	bool do_correction = (nedges >= underflow_num_edges);
	if (do_correction)
		{
		// We've traversed enough edges that it is time to take another factor out for underflow control

		// Begin by finding, for each pattern, the largest conditional likelihood over all rates and states
		// (store these in underflow_work vector)
		underflow_work.resize(total_patterns);
		underflow_work.assign(total_patterns, 0.0);
		LikeFltType * cla = cond_like.getCLA();
		unsigned subset_starting_pattern = 0;
		for (unsigned k = 0; k < nsubsets; ++k)
			{
			unsigned nr = num_rates[k];
			unsigned np = num_patterns[k];
			unsigned ns = num_states[k];
			for (unsigned r = 0; r < nr; ++r)
				{
				for (unsigned pat = 0; pat < np; ++pat)
					{
					for (unsigned i = 0; i < ns; ++i)
						{
						double curr = *cla++;
						unsigned curr_pattern = subset_starting_pattern + pat;
						if (curr > underflow_work[curr_pattern])
							{
							underflow_work[curr_pattern] = curr;
							}
						}
					}
				}
			subset_starting_pattern += np;
			}

		// Assume that x is the largest of the num_rates*num_states conditional likelihoods
		// for a given pattern. Find factor g such that g*x = underflow_max_value. Suppose
		// x = 0.05 and underflow_max_value = 10000, g = 10000/0.05 = 200,000 = e^{12.206}
		// In this case, we would add the integer component of the exponent (i.e. 12) to
		// the underflow correction for this pattern and adjust all conditional likelihoods
		// by a factor f = e^{12} = 162754.791419. So, the conditional likelihood 0.05 now
		// becomes 0.05*e^{12} = 8137.739570.
		cla = cond_like.getCLA();
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_left  = left_cond_like.getUF();
		UnderflowType const * uf_right = right_cond_like.getUF();
		unsigned curr_pattern = 0;
		unsigned curr_subset_start = 0;
		for (unsigned k = 0; k < nsubsets; ++k)
			{
			unsigned nr = num_rates[k];
			unsigned np = num_patterns[k];
			unsigned ns = num_states[k];
			unsigned condlikes_per_rate = np*ns;
			for (unsigned pat = 0; pat < np; ++pat)
				{
				PHYCAS_ASSERT(underflow_work[curr_pattern] > 0.0);
				double ratio = underflow_max_value/underflow_work[curr_pattern];
				double log_ratio = std::log(ratio);
				double f = std::floor(log_ratio);
                PHYCAS_ASSERT(f < std::log(std::numeric_limits<double>::max()));
                //if (f > tmp_max_log_f)
                //    {
                //    tmp_max_log_f = f;
                //    std::cerr << "********** new tmp_max_log_f = " << tmp_max_log_f << std::endl;
                //    }
				double expf = exp(f);

				// 2 patterns, 2 rates for first (DNA) subset
				// 3 patterns, 1 rate for second (2-state morphology) subset
				// +-------------------------------+-------------------------------+--------------------------+
				// |                           subset 0                            |          subset 1        |
				// +-------------------------------+-------------------------------+--------+--------+--------+
				// |              rate 0           |             rate 1            | rate 0 | rate 0 | rate 0 |
				// +---------------+---------------+---------------+---------------+--------+--------+--------+
				// |   pattern 0   |   pattern 1   |   pattern 0   |   pattern 1   | pat. 2 | pat. 3 | pat. 4 |
				// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+----+---+----+---+----+---+
				// | A | C | G | T | A | C | G | T | A | C | G | T | A | C | G | T |  0 | 1 |  0 | 1 |  0 | 1 |
				// +---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+---+----+---+----+---+----+---+
				//                   ^ ----------------------------> ^
				// We want claPat to start out pointing to the first state for the first rate of the
				// current pattern (indicated by ^). As we loop through rates, claPat jumps np*ns element so
				// that each time through the loop it points to the first state of the current rate and pattern.
				LikeFltType * claPat = &cla[curr_subset_start + ns*pat];
				for (unsigned r = 0; r < nr; ++r, claPat += condlikes_per_rate)
					{
					for (unsigned i = 0; i < ns; ++i)
						{
						claPat[i] *= expf;
						}
					}

				if (which == one_tip)
					{
					// cond_like, left_cond_like == right_cond_like
					// uf_left is same object as uf_right in this case
					*uf = *uf_left + (UnderflowType)f;
					++uf_left;
					}
				else if (which == no_tips)
					{
					// cond_like, left_cond_like, right_cond_like
					*uf = *uf_left + *uf_right + (UnderflowType)f;
					++uf_left;
					++uf_right;
					}
				else if (which == xtra_tip)
					{
					// cond_like == left_cond_like == right_cond_like
					// nothing above this point to take account of
					*uf += (UnderflowType)f;
					}
				else	// xtra_internal
					{
					// cond_like == left_cond_like, right_cond_like
					*uf += *uf_right + (UnderflowType)f;
					++uf_right;
					}
				++uf;
				++curr_pattern;
				}
			curr_subset_start += np*ns*nr;
			}	// loop over subsets
		nedges = 0;
		}
	else
		{
		// No underflow correction needed, but need to propagate current underflow correction
		// information down the tree
		UnderflowType *	uf = cond_like.getUF();
		UnderflowType const * uf_left = left_cond_like.getUF();
		UnderflowType const * uf_right = right_cond_like.getUF();
		for (unsigned k = 0; k < nsubsets; ++k)
			{
			unsigned np = num_patterns[k];
			for (unsigned pat = 0; pat < np; ++pat, ++uf)
				{
				if (which == one_tip)
					{
					// cond_like, left_cond_like == right_cond_like
					// uf_left is same object as uf_right in this case
					*uf = *uf_left;
					++uf_left;
					}
				else if (which == no_tips)
					{
					// cond_like, left_cond_like, right_cond_like
					*uf = *uf_left + *uf_right;
					++uf_left;
					++uf_right;
					}
				else if (which == xtra_tip)
					{
					// cond_like == left_cond_like == right_cond_like
					// nothing above this point to take account of
					}
				else	// xtra_internal
					{
					// cond_like == left_cond_like, right_cond_like
					*uf += *uf_right;
					++uf_right;
					}
				}
			}
		}
	cond_like.setUnderflowNumEdges(nedges);
	}

} // namespace phycas

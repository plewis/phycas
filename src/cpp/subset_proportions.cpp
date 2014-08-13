/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2012 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#include "subset_proportions.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor does nothing.
*/
SubsetProportions::SubsetProportions()
  : _log_prod_p(0.0)
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor does nothing.
*/
SubsetProportions::~SubsetProportions()
    {
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns data member `_log_prod_p'.
*/
double SubsetProportions::getLogProdProportions() const
    {
    return _log_prod_p;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Returns data member `_subset_proportions' as a reference to a const vector of doubles.
*/
std::vector<double> const & SubsetProportions::getSubsetProportions() const
    {
    return _subset_proportions;
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Recomputes the data member `_log_prod_p' from the current subset proportions. If p1 + p2 + p3 = 1, then
|   `_log_prod_p' = log p1 + log p2 (the final proportion is not used).
*/
void SubsetProportions::recalcLogProdP()
    {
	_log_prod_p = 0.0;

    // This version generates an (incorrect) error at compile time using some older compilers, e.g.
    // "gcc version 4.0.1 (Apple Computer, Inc. build 5367)"
    //std::vector<double>::const_reverse_iterator rit = _subset_proportions.rbegin();
    //for (++rit; rit != _subset_proportions.rend(); ++rit)
    //    {
	//	_log_prod_p += log(*rit);
    //    }

    unsigned nproportions = (unsigned)_subset_proportions.size();
    std::vector<double>::const_iterator it = _subset_proportions.begin();
    for (unsigned i = 1; i < nproportions; ++it, ++i)
        {
		_log_prod_p += log(*it);
        }
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets elements of the data member `_subset_proportions' from the supplied vector.
*/
void SubsetProportions::setSubsetProportions(
  std::vector<double> const & subset_proportions)   /**< the supplied proportions */
    {
    if (subset_proportions.empty())
        {
        _subset_proportions.clear();
        return;
        }
    _subset_proportions.assign(subset_proportions.begin(), subset_proportions.end());
    recalcLogProdP();
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Sets elements of the data member `_subset_proportions' from the supplied vector, which is expected to contain the
|   number of sites for each subset.
*/
void SubsetProportions::setSubsetProportionsFromNumSites(
  std::vector<unsigned> const & subset_nsites)   /**< the supplied proportions */
    {
    if (subset_nsites.empty())
        {
        _subset_proportions.clear();
        return;
        }

    double total_sites = (double)std::accumulate(subset_nsites.begin(), subset_nsites.end(), 0);
    if (total_sites == 0)
        {
        _subset_proportions.clear();
        return;
        }

    _subset_proportions.resize(subset_nsites.size());
    std::transform(subset_nsites.begin(), subset_nsites.end(), _subset_proportions.begin(), boost::lambda::_1/total_sites);
    recalcLogProdP();
    }

} // namespace phycas



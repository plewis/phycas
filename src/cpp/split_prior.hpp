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

#ifndef PHYCAS_SPLIT_PRIOR_H
#define PHYCAS_SPLIT_PRIOR_H

#include "basic_cdf.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Provides a single member function, SplitPrior::CalcLnPriorProb, which computes the prior probability of a split
|   that is induced by a flat prior on topologies.
*/
class SplitPrior
	{
	public:

		double  CalcLnPriorProb(unsigned n, unsigned m);

    private:

        CDF     cdf;    /**< supplies LnGamma function */
    }

/*----------------------------------------------------------------------------------------------------------------------
|	Split prior probability given a flat prior on tree topologies equals:
|>
|     (no. rooted trees with n tips) * (no. rooted trees with m tips)
|     ---------------------------------------------------------------
|                  (no. unrooted trees with n + m tips)
|
|     [ (2n-3)! / 2^(n-2) (n-2)! ] [ (2m-3)! / 2^(m-2) (m-2)! ]
|     ---------------------------------------------------------
|              [ (2n + 2m - 5)! / 2^(n+m-3) (n+m-3)!]
|>
|	where n is the number of taxa on one side of split and m is the number of taxa on the other side of the split.
|	This function returns the natural logarithm of the split prior probability for this split.
*/
inline double SplitPrior::CalcLnPriorProb(
  unsigned n,  /**< is the number of taxa on one side of the split */
  unsigned m)  /**< is the number of taxa on the other side of the split */
	{
	if (n < 2 || m < 2)
		return 0.0;

	double ln_prior = log(2.0);

	CDF cdf;
	ln_prior += cdf.LnGamma((double)(2*n - 2));
	ln_prior += cdf.LnGamma((double)(2*m - 2));
	ln_prior += cdf.LnGamma((double)(n + m - 2));
	ln_prior -= cdf.LnGamma((double)(2*(n + m) - 4));
	ln_prior -= cdf.LnGamma((double)(m - 1));
	ln_prior -= cdf.LnGamma((double)(n - 1));
	return ln_prior;
	}

}   // namespace phycas

#endif


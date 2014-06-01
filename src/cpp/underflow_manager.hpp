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

#if !defined(UNDERFLOW_MANAGER_HPP)
#define UNDERFLOW_MANAGER_HPP

#include <cmath>
#include "states_patterns.hpp"
#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Underflow manager for use with TreeLikelihood class that keeps track of underflow correction factors for each data
|	pattern. It adds a vector the length of	which is the number of site patterns to each CondLikelihood object.
*/
class UnderflowManager
	{
	public:
									UnderflowManager();

		void 						setTriggerSensitivity(unsigned nedges);
		void						setCorrectToValue(double maxval);
		void 						setDimensions(const uint_vect_t & np, const uint_vect_t & nr, const uint_vect_t & ns);

		double                      getUnderflowMaxValue() const;

		void 						twoTips(CondLikelihood & cond_like) const;
		void						check(CondLikelihood & cond_like, const CondLikelihood & left_cond_like, const CondLikelihood & right_cond_like, const count_vect_t & counts, bool polytomy) const;

		double						getCorrectionFactor(unsigned pat, ConstCondLikelihoodShPtr condlike_shptr) const;
		double						correctSiteLike(double & site_like, unsigned pat, ConstCondLikelihoodShPtr condlike_shptr) const;
		void						correctLnLike(double & ln_like, ConstCondLikelihoodShPtr condlike_shptr) const;

	protected:

		uint_vect_t					num_rates;				/**< Vector of the number of among-site rate categories for each partition subset */
		uint_vect_t					num_patterns;			/**< Vector of the number of data patterns for each partition subset */
		uint_vect_t 				num_states;				/**< Vector of the number of states for each partition subset */
		unsigned					total_patterns;			/**< The total number of patterns over all partition subsets */
		unsigned					underflow_num_edges;    /**< Number of edges to traverse before underflow risk is evaluated */
		double						underflow_max_value;    /**< Maximum of the `num_states' conditional likelihoods for a given rate and pattern after underflow correction */
		mutable std::vector<double>	underflow_work;			/**< Workspace used when correcting for underflow (will have length equal to num_patterns) */
	};

} // namespace phycas

#endif

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

#ifndef NUM_SUM_HPP
#define NUM_SUM_HPP

#include <boost/shared_ptr.hpp>
#include "states_patterns.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|   Implements numerically-stable summation in which largest element of the sum is factored out and the sum returned
|   on the log scale. Used extensively in the stepping-stone and inflated density ratio methods for estimating the
|   marginal likelihood.
*/
class NumSum
	{
	public:
                    NumSum();
                    ~NumSum();

        void        clear();
        void        addValueToSum(double log_value);
        unsigned    getNSummed() const;
        double      getLogSum() const;
        double      getLogSumOfSquares() const;
        double      getMaxLogValue() const;

	private:

        double      _max_log_value;         /**< the log of the largest value in the sum */
        double      _curr_sum;              /**< the current sum of exp(log_value - `max_log_value') */
        double      _curr_sum_of_squares;   /**< the current sum of exp(log_value - `max_log_value') */
        unsigned    _num_summed;            /**< the number of elements added to `_curr_sum' */
    };

typedef boost::shared_ptr<NumSum> NumSumShPtr;

}	// namespace phycas

#endif

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
#include <cfloat>
#include "numsum.hpp"

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Constructor calls clear() to initialize.
*/
NumSum::NumSum()
	{
    clear();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Destructor has nothing to do.
*/
NumSum::~NumSum()
	{
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Returns object to just-constructed state.
*/
void NumSum::clear()
	{
    _num_summed = 0;
    _max_log_value = 0.0;
    _curr_sum = 0.0;
    _curr_sum_of_squares = 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   If `_num_summed' is greater than 0, returns `_max_log_value' plus the log of _curr_sum; otherwise, returns -DBL_MAX.
 */
double NumSum::getLogSum() const
	{
    if (_num_summed == 0)
        return -DBL_MAX;
    else
        return _max_log_value + log(_curr_sum);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   If `_num_summed' is greater than 0, returns twice `_max_log_value' plus the log of _curr_sum_of_squares; otherwise,
|   returns -DBL_MAX.
 */
double NumSum::getLogSumOfSquares() const
	{
    if (_num_summed == 0)
        return -DBL_MAX;
    else
        return 2.0*_max_log_value + log(_curr_sum_of_squares);
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of `_num_summed'.
 */
unsigned NumSum::getNSummed() const
	{
    return _num_summed;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Returns current value of the `_max_log_value' data member,
 */
double NumSum::getMaxLogValue() const
	{
    return _max_log_value;
    }

/*----------------------------------------------------------------------------------------------------------------------
|   Adds `log_value' to the current running sum. If `log_value' is larger than the current value of `_max_log_value',
|   sets `_max_log_value' to `log_value'. The basic idea is illustrated by the following sum of three values, a, b, and
|   c, where max is the largest of the three values summed:
|
|   sum = a + b + c
|   sum = max(a/max + b/max + c/max)
|   log_sum = log(max) + log( exp(log(a) - log(max)) + exp(log(b) - log(max)) + exp(log(c) - log(max)) )
|           = _max_log_value + log(_curr_sum)
|
|   Now consider updating `_max_log_value' on the fly. Suppose c = max(a, b, c, d).
|
|   Add a:
|       max = a
|       sum = a/max
|
|   Add b:
|       prevmax = a
|       max = b
|       sum = prevmax*(a/prevmax)/max + b/max
|           = a/max + b/max
|
|   Add c:
|       prevmax = b
|       max = c
|       sum = prevmax*(a/prevmax + b/prevmax)/max + c/max
|           = a/max + b/max + c/max
|
|   Add d:
|       sum = a/max + b/max + c/max
|
|   Now for what really happens:
|
|   Add a:
|       _max_log_value = log(a)
|       _curr_sum = exp[log(a) - _max_log_value]
|
|   Add b:
|       prev_max_log_value = log(a)
|       _max_log_value = log(b)
|       _curr_sum = exp[prev_max_log_value + log(_curr_sum) - _max_log_value] + exp[log(b) - _max_log_value]
|
|   Add c:
|       prev_max_log_value = log(b)
|       _max_log_valueax = log(c)
|       _curr_sum = exp[prev_max_log_value + log(_curr_sum) - _max_log_value] + exp[log(c) - _max_log_value]
|
|   Add d:
|       _curr_sum += exp[log(d) - _max_log_value]
*/
void NumSum::addValueToSum(double log_value)
	{
    if (_num_summed == 0)
        {
        _max_log_value = log_value;
        _curr_sum = 1.0;
        _curr_sum_of_squares = 1.0;
        }
    else
        {
        if (log_value > _max_log_value)
            {
            double prev_max_log_value = _max_log_value;
            _max_log_value = log_value;
            _curr_sum = exp(prev_max_log_value + log(_curr_sum) - _max_log_value);
            _curr_sum_of_squares = exp(2.0*prev_max_log_value + log(_curr_sum_of_squares) - 2.0*_max_log_value);
            }
        _curr_sum += exp(log_value - _max_log_value);
        _curr_sum_of_squares += exp(2.0*log_value - 2.0*_max_log_value);
        }
    _num_summed++;
    }

}	// namespace phycas

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
|  Phycas: Python software for phylogenetic analysis                          |
|  Copyright (C) 2007 Mark T. Holder, Paul O. Lewis and David L. Swofford     |
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

#ifndef STOP_WATCH_HPP
#define STOP_WATCH_HPP

#include <ctime>
#include <boost/shared_ptr.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Implements a simple stop watch for timing analyses. The update function should be called regularly to avoid overflow
|   problems (at 1 millions ticks/sec, it takes less than an hour to accumulate more ticks than can be stored in a
|   32-bit double!).
|
|   An Explanation of how overflow affects timing follows. When the program begins, clock ticks begin accumulating at a
|   rate of CLOCKS_PER_SEC ticks per second. CLOCKS_PER_SEC is defined in the <ctime> header. This constant can range
|   widely depending on the compiler being used. Here we will assume it is 1,000,000 ticks/sec. At this rate, if clock_t
|   is a 4-byte long, the clock() function begins returning negative numbers after about 36 minutes! Here's why. Assume
|   that a long value z holds the current number of ticks that have assumulated since the program begin.
|>
|   at start,           z =          0L = 00000000 00000000 00000000 00000000
|   after 2147 seconds, z = 2147000000L = 01111111 11111000 10011110 11000000
|   adding 1 more sec., z = 2148000000L = 10000000 00000111 11100001 00000000
|>
|   At this point (2148 seconds), the maximum size of a (signed) long value (2147483647) has been exceeded and the
|   overflow has resulted in the sign bit (most significant bit) being set. This means that the number now represents a
|   negative value. Its value can be determined by converting the sign bit back to 0, yielding 516352 (which in binary is
|   00000000 00000111 11100001 00000000), and then adding this to the most negative possible long value (-2147483648):
|>
|   min long value = 10000000 00000000 00000000 00000000
|   add 516352     = 00000000 00000111 11100001 00000000
|   = -2146967296  = 10000000 00000111 11100001 00000000
|>
|   Thus, instead of 2148000000, we get -2146967296! At this point, it is not possible to subtract the starting ticks
|   from the ending ticks to get the number of ticks that have elapsed. To get an accurate tick count, one must
|   add together the following three separate values:
|>
|   max long value minus start ticks: 2147483647 - 2147000000    = 483647
|   stop ticks minus min long value: -2146967296 - (-2147483648) = 516352
|   add 1 because min long value was not counted:                       1
|   ---------------------------------------------------------------------
|                                                                 1000000
|>
|   If the need to add that last 1 is a bit puzzling, consider an analogous situation in the decimal system. The maximum
|   value that can be stored given only one digit of storage is 9. Adding 1 more causes overflow and the value becomes
|   the minimum possible value 0 because the 1 in the second position is discarded. Now imagine the current value is 5
|   and we add 9 to this:
|>
|   start    012345
|   stop     012345678901234
|>
|   The value 9 can be obtained as follows:
|>
|   max - start = 9 - 5 =             4
|   stop - min  = 4 - 0 =             4
|   add 1 because 0 was not counted:  1
|   -----------------------------------
|                                     9
|>
*/
class StopWatch
	{
	public:
										StopWatch();
										~StopWatch();

		void							start();
		void							stop();
		double							elapsedSeconds();
		double							split();
		void							reset();
		void							normalize();

        // for debugging
		long stopTicks();
        std::string tobinary(long x);

	private:

        clock_t                         max_ticks;          /**< largest number of ticks that can accumulate before overflow */
        clock_t                         min_ticks;          /**< smallest (most negative) number of ticks that can accumulate before overflow */
        double                          ticks_per_sec;          /**< number of ticks that accumulate per second */
        double                          elapsed_seconds;    /**< keeps track of the number of seconds that have elapsed since start() was called */
        time_t                          start_time;     /**< keeps track of the starting time by storing result of calling time() */
        time_t                          stop_time;     /**< keeps track of the stopping time by storing result of calling time() */
        time_t                          last_split_time;     /**< copies elapsed_seconds whenever split() is called */
        clock_t							start_ticks;	/**< keeps track of the starting time by storing result of calling clock() */
        clock_t							stop_ticks;	    /**< keeps track of the stopping time by storing result of calling clock() */
        bool                            running;        /**< is true iff StopWatch is currently timing something */
	};

typedef boost::shared_ptr<StopWatch> StopWatchShPtr;

}	// namespace phycas

#endif

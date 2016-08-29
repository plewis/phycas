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

#if ! defined(BASIC_LOT_HPP)
#define BASIC_LOT_HPP

#include <cmath>
#include <ctime>
#include <boost/shared_ptr.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	This class was called Lot because the noun lot is defined as "an object used in deciding something by chance"
|	according to The New Merriam-Webster Dictionary.
*/
class Lot
	{
	public:
								Lot();
								Lot(unsigned);
								~Lot();

		// Accessors
		unsigned				GetSeed() const;
		unsigned 				GetInitSeed() const;

		// Modifiers
		void 					UseClockToSeed();
		void 					SetSeed(unsigned s);

		// Utilities
        unsigned                MultinomialDraw(const double * probs, unsigned n, double totalProb=1.0);
		unsigned 				SampleUInt(unsigned);
		unsigned				GetRandBits(unsigned nbits);

		double 					Uniform();
		double 					Normal();
		bool					Boolean();

	private:

		unsigned 				last_seed_setting;
		unsigned				curr_seed;
	};

typedef boost::shared_ptr<Lot> LotShPtr;

inline bool Lot::Boolean()
	{
	return (Uniform() < 0.5);
	}

inline double Lot::Normal()
	{
    double u = Uniform();
    double x = sqrt(-2.0*log(u));
    double v = Uniform();
    double y = cos(2.0*3.141592653589793238846*v);
	return x*y;
	}

} // namespace phycas

#endif


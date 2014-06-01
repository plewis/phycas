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

#if 0
#ifndef BRENT_PRAXIS_HPP
#define BRENT_PRAXIS_HPP


#include <cstdlib>
#include <climits>
#include <cfloat>
#include "basic_lot.hpp"
//#include "probability_distribution.hpp"
#include "xprobdist.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#if defined(PYTHON_ONLY)
#	include <boost/python/call_method.hpp>
#endif
//typedef std::pair<double, double> ParamAndLnProb;
//typedef std::pair<double, double> SliceInterval;

namespace phycas
{


struct FuncToMinimize 
	{
	virtual ~FuncToMinimize() 
		{
		//std::cerr << "\n>>>>> FuncToMinimize dying..." << std::endl;
		}
    virtual double operator()(double) = 0;
	};

typedef boost::shared_ptr<FuncToMinimize> FuncToMinimizeShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Implements the multivariate minimizer described by Brent, Richard P. 1973. Algorithms for Minimization Without 
|   Derivatives. Prentice-Hall. ISBN: 0-13-022335-2. This implementation adapted from Karl Gegenfurtner's C translation,
|   available from http://archives.math.utk.edu/software/msdos/numerical.analysis/praxis/
*/
class BrentPraxis
	{
	public:
								BrentPraxis();
								BrentPraxis(LotShPtr rnd, FuncToMinimizeShPtr f);
		virtual					~BrentPraxis();

		double                  Minimize(const std::vector<double> starting_point);

		void					AttachFunc(FuncToMinimizeShPtr f);
		void					AttachRandomNumberGenerator(LotShPtr rnd);
        
    protected:
		FuncToSampleShPtr		_func;          /**< is a functor representing the function to be minimized */
		LotShPtr				_r;             /**< is the random number generator */
        std::vector<double>     _minimum_point; /**< is the final point representing the minimum */
        bool                    _illc;          /**< is false by default, but can be set to true during minimization if problem is determined to be ill-conditioned */
        unsigned                _prin;          /**< controls the printed output from the routine: 0 -> no output; 1 -> print only starting and final values; 2 -> detailed map of the minimization process; 3 -> print also eigenvalues and vectors of the search directions. The default value is 1. */
        double                  _tol;           /**< is the tolerance allowed for the precision of the solution. Minimize returns if the criterion 2 * ||x[k]-x[k-1]|| <= sqrt(macheps) * ||x[k]|| + _tol is fulfilled more than ktm times. The default value depends on the machine precision. */
        unsigned                _ktm;           /**< is the number of times the tolerance criterion can be met before Minimize returns. The default value is 1, and a value of 4 leads to a very(!) cautious stopping criterion. */
        unsigned                _maxfun;        /**< is the maximum number of calls to _func allowed. Minimize will return after _maxfun calls to _func even when the minimum is not yet found. The default value of 0 indicates no limit on the number of calls. This return condition is only checked every n iterations. */
        double                  _scbd;          /**< is a scaling parameter. 1.0 is the default and indicates no scaling. If the scales for the different parameters are very different, _scbd should be set to a value of about 10.0. */
        
        unsigned                _n;             /**< number of variables (i.e. length of _minimum_point) */
        std::vector<double>     _d;             /**< eigenvalues */
        double * *              _v;             /**< matrix of eigenvectors */
	};

typedef boost::shared_ptr<BrentPraxis> BrentPraxisShPtr;


} // namespace phycas

#endif
#endif




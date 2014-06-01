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

#if ! defined(SUBSET_PROPORTIONS_HPP)
#define SUBSET_PROPORTIONS_HPP

#include <vector>
#include <numeric>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

namespace phycas
{

/*----------------------------------------------------------------------------------------------------------------------
|	Class whose objects keep subset proportions for both PartitionModel and RelativeRateDistribution.
*/
class SubsetProportions
    {
    public:
    
                                    SubsetProportions();
                                    ~SubsetProportions();
                                
        std::vector<double> const & getSubsetProportions() const;
        double                      getLogProdProportions() const;
        void                        setSubsetProportions(std::vector<double> const & subset_proportions);
        void                        setSubsetProportionsFromNumSites(std::vector<unsigned> const & subset_nsites);
    
    protected:
    
        void                        recalcLogProdP();
    
        std::vector<double>         _subset_proportions;   /**< The stored subset proportions vector */
        double                      _log_prod_p;            /**< Constant function of subset proportions needed in order to compute subset relative rates density function */
    };
    
typedef boost::shared_ptr<SubsetProportions> SubsetProportionsShPtr;
//typedef boost::shared_ptr<const SubsetProportions> SubsetProportionsConstShPtr;

} // namespace phycas
    
#endif



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

#if ! defined(SIM_DATA_HPP)
#define SIM_DATA_HPP

#include "states_patterns.hpp"

namespace phycas
{

class Tree;
typedef std::vector<std::string>				StringVect;

/*----------------------------------------------------------------------------------------------------------------------
|	Serves as a container for data simulated by the TreeLikelihood::simulate function. Has member functions for saving
|	the stored data to a nexus file.
*/
class SimData
	{
	friend class TreeLikelihood;

	public:

									SimData();


		void						clear();
		void						resizePatternVect(unsigned sz);
		void						zeroCounts();
		void						appendCountsToFile(std::string fn, bool binary);
		void						debugAppendCountsToFile(std::string row_name, std::string fn);
		std::string					createMapleTuples(unsigned row, unsigned cutoff);
		std::vector<std::string>	getPatterns(std::vector<std::string> symbols);
		void						resetPatternLength(unsigned ntaxa);
		void						wipePattern();
		void						setState(unsigned pos, int8_t state);

        void                        insertPattern(const uint_vect_t & sitelist, pattern_count_t count);
        void                        insertPatternOnly(pattern_count_t count);

        void                        buildBinVector(unsigned nstates, bool minbins);
		std::vector<double>		    getBinnedCounts();

        double						calct(unsigned nstates);
		double						calctBinned(unsigned nstates, bool minbins);
		void						includeDataFrom(SimData &);
		unsigned					getPatternLength();
		pattern_t &					getCurrPattern();
		pattern_count_t				getTotalCount();
		unsigned					getNUniquePatterns();
		void						addDataTo(SimData & other, pattern_count_t mult);
		void						setTotalCount(pattern_count_t total);
		void						multBy(pattern_count_t factor);
		void						divideBy(pattern_count_t factor);
		void						setNumAdditions(unsigned n);
		unsigned					getNumAdditions();

        std::vector<unsigned>       getPatternVectRow(unsigned i) const;

		std::string					patternTable(const StringVect & state_symbols);
		void						saveToNexusFile(const std::string filename, const StringVect & taxon_names, const std::string datatype, const StringVect & state_symbols);

		const static state_code_t	missing_state;			/**< The value that represents missing data */

		const pattern_map_t &       getSimPatternMap() const;
        pattern_to_sites_map_t &    getPatternToSitesMap();

	private:

		//void						insertPatternToRunningAverage(pattern_count_t count, pattern_count_t p);

	private:

		unsigned					pattern_length;			/**< Number of taxa, used to reserve elements in new pattern vectors */
		pattern_map_t               sim_pattern_map;		/**< Keys are patterns, values are pattern counts */
        std::vector<pattern_t>   	sim_pattern_vect;       /**< Stores patterns in the order in which they were simulated */
		std::string					outstr;					/**< Workspace for building up a tabular representation of `sim_pattern_map' (used by showPatternMap function) */

        double_vect_t		        binv;                   /**< Stores binned counts if calctBinned is called; otherwise, will be an empty vector. */
		pattern_t					tmp_pattern;			/**< Workspace for building up a pattern */
		pattern_count_t				total_count;			/**< The number of patterns inserted since sim_pattern_map was last cleared (note that this is not sim_pattern_map.size(), but instead equals the sum of all the counts) */
        pattern_to_sites_map_t      pattern_to_sites_map;   /**< stores list of 0-offset site indices (value) for each pattern (key) */
		pattern_count_t				_nchar;                 /**< This is the length of sim_pattern_vect (set by resizePatternVect) */
	};

typedef boost::shared_ptr<SimData>	SimDataShPtr;

} // namespace phycas

#endif


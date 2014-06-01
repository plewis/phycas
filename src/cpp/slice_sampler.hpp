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

#ifndef SLICE_SAMPLER_HPP
#define SLICE_SAMPLER_HPP

#include <cstdlib>
#include <climits>
#include <cfloat>
#include "basic_lot.hpp"
#include "probability_distribution.hpp"
#include "xprobdist.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#if defined(PYTHON_ONLY)
#	include <boost/python/call_method.hpp>
#endif
typedef std::pair<double, double> ParamAndLnProb;
typedef std::pair<double, double> SliceInterval;

namespace phycas
{

struct SliceStats
	{
	SliceStats() : nsamples(0), value(0.0), width(0.0), diff(0.0), failed(0.0), evals(0.0) {}

	unsigned nsamples;
	double value;
	double width;
	double diff;
	double failed;
	double evals;
	};

// Note: AdHocDensity is defined in probablity_distribution.hpp
typedef boost::shared_ptr<AdHocDensity> FuncToSampleShPtr;

#undef WEAK_FUNCTOSAMPLE
#if defined(WEAK_FUNCTOSAMPLE)
typedef boost::weak_ptr<AdHocDensity> FuncToSampleWkPtr;
#endif

/*----------------------------------------------------------------------------------------------------------------------
|	Implements the univariate slice sampler described in Neal, Radford M. 2003. Slice sampling. Annals of Statistics
|	31:705-741.
*/
class SliceSampler
	{
	public:
								SliceSampler();
#if defined(WEAK_FUNCTOSAMPLE)
								SliceSampler(LotShPtr rnd, FuncToSampleWkPtr f);
#else
								SliceSampler(LotShPtr rnd, FuncToSampleShPtr f);
#endif
		virtual					~SliceSampler();

		double					Sample();
		std::vector<double>					DebugSample();

		double					OverrelaxedSample();
		std::vector<double>					DebugOverrelaxedSample();

#if defined(WEAK_FUNCTOSAMPLE)
		void					AttachFunc(FuncToSampleWkPtr f);
#else
		void					AttachFunc(FuncToSampleShPtr f);
#endif
		void					AttachRandomNumberGenerator(LotShPtr rnd);

		void					SetXValue(double x);
		void					SetMaxUnits(unsigned umax);

		void					SetSliceUnitWidth(double uwidth);
		double					GetSliceUnitWidth() const;

		double					AdaptSimple(double multiplier);
		double					AdaptNeal(double multiplier);
		void					AdaptYConditional(double from_ends, double multiplier);
		SliceInterval			FindSliceInterval(ParamAndLnProb x0, const double ln_y0, double tol, unsigned max_steps = UINT_MAX) const;
		double					CalcW(double y0) const;

		double					GetMode() const;
		double					GetLnDensityAtMode() const;

		double					GetLastSampledXValue();
		double					GetLastSampledYValue();
		double					GetSliceYValue();

		double					GetLnZero() const {return -DBL_MAX;}
        void                    UseDoublingMethod(bool d);

		// For diagnosing inefficiency
		//
		double					GetMinX();
		double					GetMaxX();
		double					GetOrigLeftEdgeOfSlice();
		double					GetOrigRightEdgeOfSlice();
		double					GetLeftEdgeOfSlice();
		double					GetRightEdgeOfSlice();
		unsigned				GetNumFuncEvals();
		unsigned				GetNumFailedSamples();
		unsigned				GetNumUnitsRequired();
		unsigned				GetNumSamples();

		SliceStats				SummarizeDiagnostics();
		void					ResetDiagnostics();

	protected:

		void					Init();
		ParamAndLnProb			GetNextSample(const ParamAndLnProb);
		ParamAndLnProb			GetNextOverrelaxedSample(const ParamAndLnProb);
		SliceInterval			BisectionSqueeze(double left, double lnf_left, double right, double lnf_right, const double ln_y0, double tol, unsigned max_steps) const;

#if defined(WEAK_FUNCTOSAMPLE)
		FuncToSampleWkPtr		func;				/**< is a functor representing the probability distribution to be sampled */
#else
		FuncToSampleShPtr		func;				/**< is a functor representing the probability distribution to be sampled */
#endif
		LotShPtr				r;					/**< is the random number generator */
		ParamAndLnProb			lastSampled;		/**< most recent valid sample and its relative density */
		double					w;					/**< unit size for interval I */
		unsigned				maxUnits;			/**< maximum number of units of size w to use for interval I (set to UINT_MAX for unlimited) */

		// These quantities are used for diagnostic purposes in GetNextSample()
		//
		double					orig_left_edge;		/**< x coordinate of left edge of most recent slice (before cropping) */
		double					orig_right_edge;	/**< x coordinate of right edge of most recent slice (before cropping) */
		double					left_edge;			/**< x coordinate of left edge of most recent slice (after cropping) */
		double					right_edge;			/**< x coordinate of right edge of most recent slice (after cropping) */
		double					ln_y;				/**< log of f(x) representing most recent slice */

		// For diagnosing inefficiency, all reset with call to ResetDiagnostics()
		//
		double					min_x;				/**< minimum x value tried */
		double					max_x;				/**< maximum x value tried */
		double					sumValues;			/**< sum of last `num_samples' sampled values */
		double					sumWidths;			/**< sum of last `num_samples' cropped slice widths (where slice width = right_edge - left_edge) */
		double					sumDiffs;			/**< sum of last `num_samples' differences between successively sampled values */
		unsigned				func_evals;			/**< counts number of function evaluations */
		unsigned				failed_samples;		/**< counts number of samples that failed because they were not in the slice */
		unsigned				realized_m;			/**< counts number of units, each of size w, that were required to bracket the slice */
		unsigned				num_samples;		/**< counts number of times GetNextSample called */
		std::vector<ParamAndLnProb>	most_recent;	/**< vector of all (x, lnx) pairs evaluated in most recent sampling effort */

		// These are needed for y-conditional adaptation
		bool					ycond_on;			/**< if true, w chosen anew for each sample based on y-coordinate of slice */
		double					ycond_a;
		double					ycond_b;
		double					ycond_multiplier;
		ParamAndLnProb			mode;

		// These are for overrelaxed sampling
		unsigned				num_overrelaxed_samples;	/**< counts number of times GetNextOverrelaxedSample called */

        bool                    doubling;           /**< if true, doubling method will be used to increase slice interval */
	};

typedef boost::shared_ptr<SliceSampler> SliceSamplerShPtr;

} // namespace phycas

#endif

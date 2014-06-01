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

#if defined(_MSC_VER)
#	pragma warning(disable : 4267) // boost's builtin_converters.hpp casts size_t to int rather than unsigned
#endif

//#include "phycas/force_include.h"
#include <cmath>
#include "slice_sampler.hpp"
using std::ofstream;
using std::ios;
using std::vector;
#if defined(C_FUNCS_IN_STD_NAMESPACE)
	using std::sprintf;
	using std::fabs;
	using std::log;
#endif

using namespace phycas;

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Throws XPyPhy.
*/
ParamAndLnProb SliceSampler::GetNextOverrelaxedSample(const ParamAndLnProb initialPair)
	{
	PHYCAS_ASSERT(r != NULL);

#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	PHYCAS_ASSERT(fsh != NULL);
#endif

	++num_overrelaxed_samples;

	// Step 1: choose a y value for this slice uniformly from 0.0 to f(initialX)
	// where f(initialX) is the density at the original value initialX. To avoid underflow,
	// instead choose ln_y by subtracting an exponential deviate from the
	// value log(f(initialX)). If a new point x is in the slice, the density curve will
	// lie above y (and, equivalently, the log-density will be above ln_y)
	//
	const double initialLnfX = initialPair.second;
	const double initialX = initialPair.first;
	double exponential_deviate = -log(r->Uniform());
	ln_y = initialLnfX - exponential_deviate;

	// Step 2: Find slice interval at height y
	//
	SliceInterval si = FindSliceInterval(initialPair, ln_y, 1.e-6);
	left_edge = si.first;
	right_edge = si.second;

	// Step 3: Find and return new sampled point
	//
	double x = si.first + si.second - initialX;

#if defined(WEAK_FUNCTOSAMPLE)
	double ln_fx = (*fsh)(x);
#else
	double ln_fx = (*func)(x);
#endif
	++func_evals;

	ParamAndLnProb p;
	p.first  = x;
	p.second = ln_fx;
	if (ln_fx < ln_y)
		{
		// Outside slice, not a valid sample
		// Only two reasons for this: 1) bimodal distribution and 2) inadequate precision locating slice boundaries
		// XPyPhy exception thrown in both of these cases
		//
		throw XProbDist("overrelaxed sample failed: possibly not unimodal distribution, or poor precision locating slice boundaries");
		}
	else
		{
		// About to return a sample
		//
		if (x < min_x)
			min_x = x;
		if (x > max_x)
			max_x = x;
		sumValues += x;
		sumWidths += (si.second - si.first);
		sumDiffs += std::fabs(x - initialX);

		// keep track of the sample that is closest to the mode
		if (ln_fx > mode.second)
			{
			mode.first = x;
			mode.second = ln_fx;
			}

		return p;
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Takes <tt>pair<double>(x, ln(f(x)))</tt> and
|	returns <tt>pair<double>(newPoint,ln(f(newPoint)))</tt>
*/
ParamAndLnProb SliceSampler::GetNextSample(const ParamAndLnProb initialPair)
	{
	PHYCAS_ASSERT(r != NULL);

#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	PHYCAS_ASSERT(fsh != NULL);
#endif

	++num_samples;

	// Start new vector of x values tried
	//@ most_recent vector is useful for debugging, but may slow things down too much
	most_recent.clear();
	//most_recent.push_back(initialPair);

	// Step 1: choose a y value for this slice uniformly from 0.0 to f(initialX)
	// where f(initialX) is the density at the original value initialX. To avoid underflow,
	// instead choose ln_y by subtracting an exponential deviate from the
	// value log(f(initialX)). If a new point x is in the slice, the density curve will
	// lie above y (and, equivalently, the log-density will be above ln_y)
	//
	const double initialLnfX = initialPair.second;
	const double initialX = initialPair.first;
	double exponential_deviate = -log(r->Uniform());
	ln_y = initialLnfX - exponential_deviate;
	//std::cerr << "~~ ln_y = " << ln_y << std::endl;

	// Step 1.5: if using y-conditional adaptation, choose a new value for w
	//
	if (ycond_on)
		{
		w = CalcW(std::exp(ln_y));
		}

	// Step 2: randomly place interval of width w around initialX (see Fig. 3, p. 715)
	//
	double U	= r->Uniform();
	left_edge	= initialX - w*U;
	right_edge	= left_edge + w;
	++realized_m;

	// Step 3: choose maximum number of units on left (J) and right (K)
	// by dividing up the maxUnits total units randomly (see Fig. 3, p. 715)
	//
	double V	= r->Uniform();
	unsigned J	= (unsigned)(maxUnits*V);
	unsigned K	= (maxUnits - 1) - J;

	//std::cerr << "maxUnits = " << maxUnits << "\n";
	//std::cerr << "V        = " << V << "\n";
	//std::cerr << "J        = " << J << "\n";
	//std::cerr << "K        = " << K << "\n";
	//std::cerr << "-DBL_MAX = " << (-DBL_MAX) << "\n";
	//std::cerr << std::endl;
	//int ch = scanf("%c");

    if (doubling)
        {
        // Grow interval by adding the current interval width to the left or right (determined randomly)
        // (see Fig. 4, p. 715). Assumes no limit to the number of doublings (i.e. p = infinity) and that
        // the distribution is unimodal (does not perform the check in Fig. 6).
        //
        bool left_ok = false;
        bool right_ok = false;
        double new_left_edge = left_edge;
        double new_right_edge = right_edge;
#if defined(WEAK_FUNCTOSAMPLE)
      	FuncToSampleShPtr fsh = func.lock();
        double left_edge_ln_y = (*fsh)(left_edge);
        double right_edge_ln_y = (*fsh)(right_edge);
#else
        double left_edge_ln_y = (*func)(left_edge);
        double right_edge_ln_y = (*func)(right_edge);
#endif
        for (;;)
            {
		    if (left_edge_ln_y < ln_y)
                left_ok = true;
		    if (right_edge_ln_y < ln_y)
                right_ok = true;
            if (left_ok && right_ok)
                break;
        	double V = r->Uniform();
            if (V < 0.5)
                {
                left_edge -= (right_edge - left_edge);
                if (!left_ok)
                    {
#if defined(WEAK_FUNCTOSAMPLE)
                    left_edge_ln_y = (*fsh)(left_edge);
#else
                    left_edge_ln_y = (*func)(left_edge);
#endif
                    new_left_edge = left_edge;
                    }
                }
            else
                {
                right_edge += (right_edge - left_edge);
                if (!right_ok)
                    {
#if defined(WEAK_FUNCTOSAMPLE)
                    right_edge_ln_y = (*fsh)(right_edge);
#else
                    right_edge_ln_y = (*func)(right_edge);
#endif
                    new_right_edge = right_edge;
                    }
                }
            //std::cerr << "I = [" << new_left_edge << "," << new_right_edge << "]" << std::endl;
            }
        left_edge = new_left_edge;
        right_edge = new_right_edge;
        }
    else
        {
	    // Step 4: Grow interval to the left until left edge is not in the slice
	    // or J units have been added, whichever comes first
	    //
	    for (;;)
		    {
#if defined(WEAK_FUNCTOSAMPLE)
		    double left_edge_ln_y = (*fsh)(left_edge);
#else
		    double left_edge_ln_y = (*func)(left_edge);
#endif
		    //std::cerr << "~~ left edge = " << left_edge << ", left_edge_ln_y = " << left_edge_ln_y << std::endl;
		    ++func_evals;
		    if (left_edge_ln_y < ln_y || J == 0)
			    break;

		    left_edge -= w;
		    --J;
		    ++realized_m;
		    }

	    // Step 5: Grow interval to the right until right edge is not in the slice
	    // or K units have been added, whichever comes first
	    //
	    for (;;)
		    {
#if defined(WEAK_FUNCTOSAMPLE)
		    double right_edge_ln_y = (*fsh)(right_edge);
#else
		    double right_edge_ln_y = (*func)(right_edge);
#endif
		    ++func_evals;
		    //std::cerr << "~~ right edge = " << right_edge << ", right_edge_ln_y = " << right_edge_ln_y << std::endl;
		    if (right_edge_ln_y < ln_y || K == 0)
			    break;

		    right_edge += w;
		    --K;
		    ++realized_m;
		    }
        }

	orig_left_edge = left_edge;
	orig_right_edge = right_edge;

	// Step 6: Choose new x value uniformly from interval. If chosen x value
	// is in the slice, return this as the sampled value. If outside the slice,
	// adjust the interval accordingly and repeat until a sampled value is found.
	//
	ParamAndLnProb p;
	for(;;)
		{
		double x = left_edge + ((right_edge - left_edge) * r->Uniform());

#if defined(WEAK_FUNCTOSAMPLE)
		FuncToSampleShPtr fsh = func.lock();
		double ln_fx = (*fsh)(x);
#else
		double ln_fx = (*func)(x);
#endif
		//std::cerr << "~~ x = " << x << ", ln_fx = " << ln_fx << std::endl;
		++func_evals;

		p.first  = x;
		p.second = ln_fx;
		if (ln_fx < ln_y)
			{
			// Outside slice, not a valid sample
			// Adjust edge so that it is closer to slice boundary and try again
			//
			++failed_samples;
			most_recent.push_back(p);
			if (x > initialX)
				{
				right_edge = x;
				//std::cerr << "~~ failed: right_edge cropped to " << right_edge << std::endl;
				}
			else
				{
				left_edge = x;
				//std::cerr << "~~ failed: left_edge cropped to " << left_edge << std::endl;
				}
			}
		else
			{
			// About to return a sample
			//
			if (x < min_x)
				min_x = x;
			if (x > max_x)
				max_x = x;
			sumValues += x;
			sumWidths += (right_edge - left_edge);
			sumDiffs += std::fabs(x - initialX);

			// keep track of the sample that is closest to the mode
			if (ln_fx > mode.second)
				{
				mode.first = x;
				mode.second = ln_fx;
				}

			//std::cerr << "*** success: ln_y = " << ln_y << ", left = " << left_edge << ", right = " << right_edge << std::endl;
			return p;
			}
		}
	}

// beginning of functions formerly in slice_sampler.hpp

/*----------------------------------------------------------------------------------------------------------------------
|	Calls the Init() member function.
*/
SliceSampler::SliceSampler()
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the starting x-value (`lastSampled.first') to 0.1, `left_bound' to -DBL_MAX and `right_bound' to DBL_MAX.
|	In addition, sets data members `r' and `probdist' to the values passed in.
*/
SliceSampler::SliceSampler(
  LotShPtr rnd,				/**< is the random number generator object to use */
#if defined(WEAK_FUNCTOSAMPLE)
  FuncToSampleWkPtr f)		/**< is the probability distribution to sample */
#else
  FuncToSampleShPtr f)		/**< is the probability distribution to sample */
#endif
  : func(f), r(rnd)
	{
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Nothing to do.
*/
SliceSampler::~SliceSampler()
	{
	//std::cerr << "\n>>>>> SliceSampler dying..." << std::endl;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Called by constructors to initialize the object.
*/
void SliceSampler::Init()
	{
	w					= 1.0;
	ln_y				= 0.0;
	maxUnits			= UINT_MAX;

	lastSampled.first	= 0.1;
	lastSampled.second	= 0.0;

	ycond_on			= false;
	ycond_a				= 1.0;
	ycond_b				= 0.0;
	mode.first			= 0.1;
	mode.second			= -DBL_MAX;
	ycond_multiplier	= 1.0;

    doubling            = false;

	ResetDiagnostics();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Attaches a ProbabilityDistribution object representing the distribution to be sampled. Calls Init() to reset the
|	sampler to its just-constructed state.
*/
void SliceSampler::AttachFunc(
#if defined(WEAK_FUNCTOSAMPLE)
  FuncToSampleWkPtr f)
#else
  FuncToSampleShPtr f)
#endif
	{
	func = f;
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets diagnostic counters (`func_evals', `failed_samples', and `realized_m') all to 0.
*/
void SliceSampler::ResetDiagnostics()
	{
	min_x					= DBL_MAX;
	max_x					= -DBL_MAX;
	num_samples				= 0;
	num_overrelaxed_samples	= 0;
	func_evals				= 0;
	failed_samples			= 0;
	realized_m				= 0;
	sumValues				= 0.0;
	sumWidths				= 0.0;
	sumDiffs				= 0.0;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Computes means of several quantities that have accumulated since the last call to ResetDiagnostics.
*/
SliceStats SliceSampler::SummarizeDiagnostics()
	{
	SliceStats ss;
	if (num_samples > 0)
		{
		ss.nsamples	= num_samples;
		ss.value	= sumValues/(double)num_samples;
		ss.width	= sumWidths/(double)num_samples;
		ss.diff		= sumDiffs/(double)num_samples;
		ss.evals	= (double)func_evals/(double)num_samples;
		ss.failed	= (double)failed_samples/(double)num_samples;
		}

	return ss;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current value of `w'.
*/
double SliceSampler::GetSliceUnitWidth() const
	{
	return w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Resets `w' to a value equal to the mean width of a cropped slice interval times the multiplier provided.
*/
double SliceSampler::AdaptSimple(double multiplier)
	{
	PHYCAS_ASSERT(num_samples > 0);
	ycond_on = false;
	w = multiplier*sumWidths/(double)num_samples;
	return w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	`w' is set to the average distance between sampled values times the multiplier supplied. Suggested by Neal in
|	section 4.4 (p. 721).
*/
double SliceSampler::AdaptNeal(double multiplier)
	{
	PHYCAS_ASSERT(num_samples > 0);
	ycond_on = false;
	w = multiplier*sumDiffs/(double)num_samples;
	return w;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Finds edges of the slice at level y0 (the natural log of which is `ln_y0') by first stepping out in units of size
|	`w' until the edge of the slice is bracketed, then using bisection to locate the slice boundary, starting at `x0'
|	and continuing until the distance between successive points is less than `tol'. If `max_steps' is less than
|	UINT_MAX, performs exactly `max_steps' bisections and ignores `tol'. Returns outer width of slice, which is always
|	greater than or equal to the exact width of the slice by an amount that depends on `tol' or `max_steps'.
|	Throws XProbDist exception if it turns out that the mode (which is really an estimate of the mode) is not actually
|	underneath the density curve at the height y0.
*/
SliceInterval SliceSampler::FindSliceInterval(ParamAndLnProb x0, const double ln_y0, double tol, unsigned max_steps) const
	{
	// bracket left edge
	//
	double left_edge = x0.first;
	double curr_lnfx = x0.second;
	double prev_lnfx = x0.second;
	while (curr_lnfx > ln_y0)
		{
		left_edge -= w;
		prev_lnfx = curr_lnfx;
#if defined(WEAK_FUNCTOSAMPLE)
		FuncToSampleShPtr fsh = func.lock();
		curr_lnfx = (*fsh)(left_edge);
#else
		curr_lnfx = (*func)(left_edge);
#endif
		}

	// use bisection to locate left edge
	//
	SliceInterval slice_interval = BisectionSqueeze(left_edge, curr_lnfx, left_edge + w, prev_lnfx, ln_y0, tol, max_steps);
	left_edge = slice_interval.first;

	// bracket right edge
	//
	double right_edge = x0.first;
	curr_lnfx = x0.second;
	prev_lnfx = x0.second;
	while (curr_lnfx > ln_y0)
		{
		right_edge += w;
		prev_lnfx = curr_lnfx;
#if defined(WEAK_FUNCTOSAMPLE)
		FuncToSampleShPtr fsh = func.lock();
		curr_lnfx = (*fsh)(right_edge);
#else
		curr_lnfx = (*func)(right_edge);
#endif
		}

	// use bisection to locate right edge
	//
	slice_interval = BisectionSqueeze(right_edge - w, prev_lnfx, right_edge, curr_lnfx, ln_y0, tol, max_steps);
	right_edge = slice_interval.second;

	return SliceInterval(left_edge, right_edge);
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Recursive function that uses bisection to locate precisely the point at which the function crosses the value y0
|	(the natural log of which is `1n_y0'). Stops when the distance between `left' and `right' is less than `tol', or
|	when `max_steps' is 0 (each recursive call decrements `max_steps' by 1. Assumes supplied values `left' and `right'
|	bracket the crossover point.
*/
SliceInterval SliceSampler::BisectionSqueeze(double left, double lnf_left, double right, double lnf_right, const double ln_y0, double tol, unsigned max_steps) const
	{
	bool left_below = (lnf_left < ln_y0);
	bool right_below = (lnf_right < ln_y0);
	PHYCAS_ASSERT((left_below && !right_below) || (right_below && !left_below));
	double middle = (left + right)/2.0;
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	double lnf_middle = (*fsh)(middle);
#else
	double lnf_middle = (*func)(middle);
#endif
	bool middle_below = (lnf_middle < ln_y0);
	bool middle_above = !middle_below;
	double half_width = (right - left)/2.0;
	bool leave_now = (half_width < tol);
	bool go_left = false;
	if (left_below && middle_below)
		go_left = false;
	else if (left_below && middle_above)
		go_left = true;
	else if (right_below && middle_below)
		go_left = true;
	else if (right_below && middle_above)
		go_left = false;
	else
		PHYCAS_ASSERT(0);
	if (leave_now)
		{
		SliceInterval slice_interval;
		if (go_left)
			{
			slice_interval.first = left;
			slice_interval.second = middle;
			}
		else
			{
			slice_interval.first = middle;
			slice_interval.second = right;
			}
		return slice_interval;
		}
	else // half_width not yet smaller than tol
		{
		if (go_left)
			{
			return BisectionSqueeze(left, lnf_left, middle, lnf_middle, ln_y0, tol, max_steps - 1);
			}
		else
			{
			return BisectionSqueeze(middle, lnf_middle, right, lnf_right, ln_y0, tol, max_steps - 1);
			}
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns value of `w' to be used for next sample when using y-conditional adaptation. Uses the simple linear formula
|	w = m*(a + b*y), where a and b are calculated, and m is set, by calling AdaptYConditional, and y is the supplied
|	value `y0'. Note: `y0' should NOT be on the log scale.
*/
double SliceSampler::CalcW(double y0) const
	{
	double w0 = ycond_a + ycond_b*y0;
	return w0*ycond_multiplier;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns current estimate of mode based on previous sampling. Note that this is not the marginal mode! Mode here
|   means simply the value corresponding to the highest function value seen thus far.
*/
double SliceSampler::GetMode() const
	{
	return mode.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns density at current estimate of mode based on previous sampling.
*/
double SliceSampler::GetLnDensityAtMode() const
	{
	return mode.second;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	`w' is not set, but instead is recalulated for each sample based on the value of y (height of slice). The purpose
|	of this function is to parameterize the function used to compute `w' given y. The width of the slice is computed
|	using bisection at two points, one a fraction `from_ends' from the bottom and the other a fraction
|	`from_ends' from the top of the density at its mode. These two widths define a slope that can be used to
|	approximate the width of the density at any height. The actual width used for `w' is `multiplier' time the width
|	returned by the function.
*/
void SliceSampler::AdaptYConditional(double from_ends, double multiplier)
	{
	ycond_on = true;
	ycond_multiplier = multiplier;

	double ln_from_ends = std::log(from_ends);

	double ln_y0 = mode.second;
	double ln_y1 = ln_from_ends + ln_y0;
	SliceInterval s1 = FindSliceInterval(mode, ln_y1, 0.01);
	double w1 = s1.second - s1.first;

	double ln_from_ends_complement = std::log(1.0 - from_ends);

	double ln_y2 = ln_from_ends_complement + ln_y0;
	SliceInterval s2 = FindSliceInterval(mode, ln_y2, 0.01);
	double w2 = s2.second - s2.first;

	// Calculate intercept a and slope b using the two equations:
	//   w1 = a + b*y1
	//   w2 = a + b*y2
	//
	//       (w1 - w2)
	//   b = ---------
	//       (y1 - y2)
	//
	//   a = w1 - b*y1
	//
	double y1 = std::exp(ln_y1);
	double y2 = std::exp(ln_y2);
	ycond_b = (w1 - w2)/(y1 - y2);
	ycond_a = w1 - ycond_b*y1;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of times GetNextSample called. Use the value returned here to compute averages of the other
|	diagnostic counts.
*/
unsigned SliceSampler::GetNumSamples()
	{
	return num_samples;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the minimum value of x tried since last call to ResetDiagnostics function, or DBL_MAX if no sampling has
|	been attempted.
*/
double SliceSampler::GetMinX()
	{
	return min_x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the maximum value of x tried since last call to ResetDiagnostics function, or -DBL_MAX if no sampling has
|	been attempted.
*/
double SliceSampler::GetMaxX()
	{
	return max_x;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the log of the y coordinate of the most recent slice.
*/
double SliceSampler::GetSliceYValue()
	{
	return ln_y;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the left edge of the most recent slice, before cropping based on failed samples.
*/
double SliceSampler::GetOrigLeftEdgeOfSlice()
	{
	return orig_left_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the right edge of the most recent slice, before cropping based on failed samples.
*/
double SliceSampler::GetOrigRightEdgeOfSlice()
	{
	return orig_right_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the left edge of the most recent slice, after cropping based on failed samples.
*/
double SliceSampler::GetLeftEdgeOfSlice()
	{
	return left_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the x coordinate of the right edge of the most recent slice, after cropping based on failed samples.
*/
double SliceSampler::GetRightEdgeOfSlice()
	{
	return right_edge;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of function evaluations required. Divide by GetNumSamples to compute average number of function
|	evaluations required for each call to GetNextSample.
*/
unsigned SliceSampler::GetNumFuncEvals()
	{
	return func_evals;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of failed samples. Each time a sample is attempted but fails because the value is outside the slice,
|	the value of the counter `failed_samples' is incremented and the value is used to reduce the size of the sampling
|	interval. Divide by GetNumSamples to compute average number of failed samples per call to GetNextSample.
*/
unsigned SliceSampler::GetNumFailedSamples()
	{
	return failed_samples;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns number of units, each of size w, required to bracket the slice. Divide by GetNumSamples to compute average
|	number of units required per call to GetNextSample.
*/
unsigned SliceSampler::GetNumUnitsRequired()
	{
	return realized_m;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the random number generator pointer `r' to the supplied random number generator pointer `rnd'. Calls Init() to
|	reset the sampler to its just-constructed state.
*/
void SliceSampler::AttachRandomNumberGenerator(
  LotShPtr rnd)     /**< is the random number generator to attach */
	{
	r = rnd;
	Init();
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the starting x-value `lastSampled.first' to the supplied value `x'.
*/
void SliceSampler::SetXValue(
  double x)     /**< is the starting x-value to use */
	{
	if (x != lastSampled.first)
		{
		lastSampled.first = x;
#if 0
		FuncToSampleShPtr fsh = func.lock();
		PHYCAS_ASSERT(fsh);
		lastSampled.second = (*fsh)(x);
#endif
		}
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `lastSampled.first'.
*/
double SliceSampler::GetLastSampledXValue()
	{
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Returns the value of `lastSampled.second'. Note: causes function to be evaluated because, in general, we can never
|   be sure that the current value stored in lastSampled.second is correct (e.g. if func represents a posterior density,
|   then changes in other model parameters will cause the posterior density of lastSampled.first to be different than
|   it was when lastSampled.first was first evaluated).
*/
double SliceSampler::GetLastSampledYValue()
	{
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	PHYCAS_ASSERT(fsh);
	lastSampled.second = (*fsh)(lastSampled.first);
#else
	lastSampled.second = (*func)(lastSampled.first);
#endif
	return lastSampled.second;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the unit size for the interval I used for sampling from the slice S.
*/
void SliceSampler::SetSliceUnitWidth(double uwidth)
	{
	w = uwidth;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Sets the maximum number of units of size `w' to use for interval I used for sampling from the slice S. If `umax'
|	equals zero, maxUnits is set to UINT_MAX instead.
*/
void SliceSampler::SetMaxUnits(unsigned umax)
	{
	maxUnits = umax;
	if (maxUnits == 0)
		maxUnits = UINT_MAX;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Current point is `lastSampled.first'
*/
double SliceSampler::OverrelaxedSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	if (fsh == NULL)
		{
		throw XProbDist("must specify density function for slice sampler before attempting to draw samples");
		}

	// We can never be guaranteed that the full conditional density at lastSampled.first has not
	// changed since the last call during an MCMC run.
	lastSampled.second = (*fsh)(lastSampled.first);
#else
	// We can never be guaranteed that the full conditional density at lastSampled.first has not
	// changed since the last call during an MCMC run.
	lastSampled.second = (*func)(lastSampled.first);
#endif
	++func_evals;

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextOverrelaxedSample(currentPoint);

	// Let the new sampled value be the starting point for the next sample
	//
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using overrelaxed slice sampling. Same as OverrelaxedSample, but returns
|	a vector the elements of which are:
|		0: sampled x
|		1: x-coord of vertical slice
|		2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
|		3: x-coord of left edge of horizontal slice
|		4: x-coord of right edge of horizontal slice
|		5: y-coord of horizontal slice
*/
std::vector<double> SliceSampler::DebugOverrelaxedSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	if (fsh == NULL)
		{
		throw XProbDist("must attach a probability distribution to slice sampler before attempting to draw samples");
		}
#endif

	double v1 = lastSampled.first; // x-coord of vertical slice

	// We can never be guaranteed that the full conditional density at lastSampled.first has not
	// changed since the last call during an MCMC run.
#if defined(WEAK_FUNCTOSAMPLE)
	lastSampled.second = (*fsh)(lastSampled.first);
#else
	lastSampled.second = (*func)(lastSampled.first);
#endif
	++func_evals;

	double v2 = lastSampled.second; // y-coord of top of vertical slice

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextOverrelaxedSample(currentPoint);

	double v0 = lastSampled.first;	// sampled x
	double v3 = left_edge;			// x-coord of left edge of horizontal slice
	double v4 = right_edge;			// x-coord of right edge of horizontal slice
	double v5 = ln_y;				// y-coord of horizontal slice

	std::vector<double> v;
	v.push_back(v0);
	v.push_back(v1);
	v.push_back(v2);
	v.push_back(v3);
	v.push_back(v4);
	v.push_back(v5);

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Current point is `lastSampled.first'
*/
double SliceSampler::Sample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}

	// We can never be guaranteed that the full conditional density at lastSampled.first has not
	// changed since the last call during an MCMC run, so recalculate it now
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
    lastSampled.second = (*fsh)(lastSampled.first);
#else
    lastSampled.second = (*func)(lastSampled.first);
#endif
	++func_evals;

	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextSample(currentPoint);

    // Let the new sampled value be the starting point for the next sample
	//
	return lastSampled.first;
	}

/*----------------------------------------------------------------------------------------------------------------------
|	Draws a sample from the target distribution using slice sampling. Same as Sample, but returns a vector the elements
|	of which are:
|       0: sampled x
|       1: x-coord of vertical slice
|       2: y-coord of top of vertical slice (y-coord of bottom of vertical slice always 0.0)
|       3: x-coord of left edge of horizontal slice
|       4: x-coord of right edge of horizontal slice
|       5: y-coord of horizontal slice
|       6: horizontal slice interval width
|		7+: x-coord of failed sampling attempts (y-coord all equal to element 5)
*/
std::vector<double> SliceSampler::DebugSample()
	{
	if (r == NULL)
		{
		throw XProbDist("must attach random number generator to slice sampler before attempting to draw samples");
		}
#if defined(WEAK_FUNCTOSAMPLE)
	FuncToSampleShPtr fsh = func.lock();
	if (fsh == NULL)
		{
		throw XProbDist("must attach a probability distribution to slice sampler before attempting to draw samples");
		}
#endif

	// We can never be guaranteed that the full conditional density at lastSampled.first has not
	// changed since the last call during an MCMC run.
	lastSampled.second = 0.0;

#if defined(WEAK_FUNCTOSAMPLE)
	lastSampled.second = (*fsh)(lastSampled.first);
#else
	lastSampled.second = (*func)(lastSampled.first);
#endif
	++func_evals;
	//std::cerr << "~~ re-evaluating density at current point in DebugSample" << std::endl;

	double v1 = lastSampled.first; // x-coord of vertical slice
	double v2 = lastSampled.second; // y-coord of top of vertical slice

	// Let the new sampled value be the starting point for the next sample
	//
	const ParamAndLnProb currentPoint = lastSampled;
	lastSampled = GetNextSample(currentPoint);

	double v0 = lastSampled.first;	// sampled x
	double v3 = left_edge;			// x-coord of left edge of horizontal slice
	double v4 = right_edge;			// x-coord of right edge of horizontal slice
	double v5 = ln_y;				// y-coord of horizontal slice
	double v6 = w;					// horizontal slice interval width

	std::vector<double> v;
	v.push_back(v0);
	v.push_back(v1);
	v.push_back(v2);
	v.push_back(v3);
	v.push_back(v4);
	v.push_back(v5);
	v.push_back(v6);

	// let the last items in the vector be the x-coordinates of the sampling attempts (very last element
	// is the x-coord of the successful sample, which should be the same as v0)
	for (std::vector<ParamAndLnProb>::const_iterator it = most_recent.begin(); it != most_recent.end(); ++it)
		{
		v.push_back(it->first);
		}

	return v;
	}

/*----------------------------------------------------------------------------------------------------------------------
|   Sets the data member `doubling' to the value specified.
*/
void SliceSampler::UseDoublingMethod(
  bool d)   /**< is the new value for `doubling' */
    {
    doubling = d;
    }

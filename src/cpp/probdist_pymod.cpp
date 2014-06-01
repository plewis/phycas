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

#if defined(USING_NUMARRAY)
#	define PY_ARRAY_UNIQUE_SYMBOL PyArrayHandle
#endif

#include <boost/python.hpp>

#include "basic_lot.hpp"
#include "probability_distribution.hpp"
#include "mvnormal_distribution.hpp"
#include "relative_rate_distribution.hpp"
#include "lognormal.hpp"
#include "stop_watch.hpp"
#include "slice_sampler.hpp"
#include "square_matrix.hpp"
#include "rectangular_matrix.hpp"
#include "xprobdist.hpp"

using namespace boost::python;
using namespace phycas;

double getEffectiveLnZero()
	{
	return -DBL_MAX;
	}

double lnGamma(double x)
	{
	return CDF().LnGamma(x);
	}

// The following wrapper struct is needed because we will potentially be deriving Python classes from AdHocDensity
// This struct is thus only necessary if we plan to do something like this in Python:
//		MyFunctor(AdhocDensity):
//			...
// See http://wiki.python.org/moin/boost.python/InternalDataStructures
//
struct AdHocDensityWrapper : AdHocDensity
	{
	AdHocDensityWrapper(PyObject * p) : self(p) {}
	virtual ~AdHocDensityWrapper()
		{
		//std::cerr << "(ProbDist)AdHocDensityWrapper dying..." << std::endl;
		}
	double operator()(double x) { return boost::python::call_method<double>(self, "__call__", x); }
    PyObject * self;
	};

void translateXProbDist(const XProbDist &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

BOOST_PYTHON_MODULE(_ProbDistExt)
{
	def("getEffectiveLnZero", getEffectiveLnZero);
	def("lnGamma", lnGamma);

#if defined(USING_NUMARRAY)
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");
#endif

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensityWrapper> >("AdHocDensityBase")
		;

	class_<ProbabilityDistribution, bases<AdHocDensity>, boost::shared_ptr<ProbabilityDistribution>, boost::noncopyable>("ProbabilityDistribution", no_init)
        .def("lnGamma", &ProbabilityDistribution::LnGamma)
		;

	class_<MultivariateProbabilityDistribution, boost::shared_ptr<MultivariateProbabilityDistribution>, boost::noncopyable>("MultivariateProbabilityDistribution", no_init)
		;

	class_<SliceSampler, boost::shared_ptr<SliceSampler> >("SliceSamplerBase")
		.def(init<boost::shared_ptr<phycas::Lot>, boost::shared_ptr<AdHocDensity> >())
		.def("attachFunc", &SliceSampler::AttachFunc)
		.def("sample", &SliceSampler::Sample)
		.def("debugSample", &SliceSampler::DebugSample)
		.def("overrelaxedSample", &SliceSampler::OverrelaxedSample)
		.def("debugOverrelaxedSample", &SliceSampler::DebugOverrelaxedSample)
		.def("attachRandomNumberGenerator", &SliceSampler::AttachRandomNumberGenerator)
		.def("adaptSimple", &SliceSampler::AdaptSimple)
		.def("adaptNeal", &SliceSampler::AdaptNeal)
		.def("getSliceUnitWidth", &SliceSampler::GetSliceUnitWidth)
		.def("setSliceUnitWidth", &SliceSampler::SetSliceUnitWidth)
		.def("getMinX", &SliceSampler::GetMinX)
		.def("getMaxX", &SliceSampler::GetMaxX)
		.def("setMaxUnits", &SliceSampler::SetMaxUnits)
		.def("getOrigLeftEdgeOfSlice", &SliceSampler::GetOrigLeftEdgeOfSlice)
		.def("getOrigRightEdgeOfSlice", &SliceSampler::GetOrigRightEdgeOfSlice)
		.def("getLeftEdgeOfSlice", &SliceSampler::GetLeftEdgeOfSlice)
		.def("getRightEdgeOfSlice", &SliceSampler::GetRightEdgeOfSlice)
		.def("getSliceYValue", &SliceSampler::GetSliceYValue)
		.def("getNumFuncEvals", &SliceSampler::GetNumFuncEvals)
		.def("getNumFailedSamples", &SliceSampler::GetNumFailedSamples)
		.def("getNumUnitsRequired", &SliceSampler::GetNumUnitsRequired)
		.def("getNumSamples", &SliceSampler::GetNumSamples)
		.def("resetDiagnostics", &SliceSampler::ResetDiagnostics)
		.def("adaptYConditional", &SliceSampler::AdaptYConditional)
		.def("calcW", &SliceSampler::CalcW)
		.def("getMode", &SliceSampler::GetMode)
		.def("getLnDensityAtMode", &SliceSampler::GetLnDensityAtMode)
		.def("getLastSampledXValue", &SliceSampler::GetLastSampledXValue)
		.def("getLastSampledYValue", &SliceSampler::GetLastSampledYValue)
		.def("setXValue", &SliceSampler::SetXValue)
		.def("useDoublingMethod", &SliceSampler::UseDoublingMethod)
		;

	class_<phycas::StopWatch, boost::shared_ptr<phycas::StopWatch>, boost::noncopyable>("StopWatchBase")
		.def("start", &phycas::StopWatch::start)
		.def("stop", &phycas::StopWatch::stop)
		.def("reset", &phycas::StopWatch::reset)
		.def("normalize", &phycas::StopWatch::normalize)
		.def("elapsedSeconds", &phycas::StopWatch::elapsedSeconds)
		.def("split", &phycas::StopWatch::split)
		.def("stopTicks", &phycas::StopWatch::stopTicks)
		//.def("doofus", &phycas::StopWatch::doofus)
		;

	class_<phycas::SubsetProportions, boost::shared_ptr<phycas::SubsetProportions>, boost::noncopyable>("SubsetProportionsBase")
		.def("getSubsetProportions", &phycas::SubsetProportions::getSubsetProportions, return_value_policy<copy_const_reference>())
		.def("setSubsetProportions", &phycas::SubsetProportions::setSubsetProportions)
		.def("setSubsetProportionsFromNumSites", &phycas::SubsetProportions::setSubsetProportionsFromNumSites)
        .def("getLogProdProportions", &phycas::SubsetProportions::getLogProdProportions)
		;

	class_<phycas::SquareMatrix, boost::shared_ptr<phycas::SquareMatrix>, boost::noncopyable>("SquareMatrixBase", init<unsigned, double>())
		.def(init<const phycas::SquareMatrix &>())
		.def("duplicate", &phycas::SquareMatrix::Duplicate, return_value_policy<manage_new_object>())
		.def("pow", &phycas::SquareMatrix::Power, return_value_policy<manage_new_object>())
		.def("inverse", &phycas::SquareMatrix::Inverse, return_value_policy<manage_new_object>())
		.def("LUDecomposition", &phycas::SquareMatrix::LUDecomposition, return_value_policy<manage_new_object>())
		.def("CholeskyDecomposition", &phycas::SquareMatrix::CholeskyDecomposition, return_value_policy<manage_new_object>())
		.def("rightMultiplyMatrix", &phycas::SquareMatrix::RightMultiplyMatrix, return_value_policy<manage_new_object>())
		.def("leftMultiplyMatrix", &phycas::SquareMatrix::LeftMultiplyMatrix, return_value_policy<manage_new_object>())
		.def("rightMultiplyVector", &phycas::SquareMatrix::RightMultiplyVector)
		.def("leftMultiplyVector", &phycas::SquareMatrix::LeftMultiplyVector)
		.def("identity", &phycas::SquareMatrix::Identity)
		.def("trace", &phycas::SquareMatrix::Trace)
		.def("logProdMainDiag", &phycas::SquareMatrix::LogProdMainDiag)
		.def("logDeterminant", &phycas::SquareMatrix::LogDeterminant)
		.def("getDimension", &phycas::SquareMatrix::GetDimension)
		.def("__repr__", &phycas::SquareMatrix::GetStringRepresentation)
		.def("addToElement", &phycas::SquareMatrix::AddToElement)
		.def("setElement", &phycas::SquareMatrix::SetElement)
		.def("getElement", &phycas::SquareMatrix::GetElement)
		.def("setMatrix", &phycas::SquareMatrix::SetMatrix)
		.def("getMatrix", &phycas::SquareMatrix::GetMatrix)
		;

	class_<phycas::RectangularMatrix, boost::shared_ptr<phycas::RectangularMatrix>, boost::noncopyable>("RectangularMatrixBase", init<unsigned, unsigned, double>())
		.def(init<const phycas::RectangularMatrix &>())
        .def("getDimensions", &phycas::RectangularMatrix::GetDimensions)
        .def("getNRows", &phycas::RectangularMatrix::GetNRows)
        .def("getNCols", &phycas::RectangularMatrix::GetNCols)
		.def("__repr__", &phycas::RectangularMatrix::GetStringRepresentation)
        .def("addToElement", &phycas::RectangularMatrix::AddToElement)
        .def("setElement", &phycas::RectangularMatrix::SetElement)
        .def("getElement", &phycas::RectangularMatrix::GetElement)
        .def("setMatrix", &phycas::RectangularMatrix::SetMatrix)
        .def("getMatrix", &phycas::RectangularMatrix::GetMatrix)
        .def("setRow", &phycas::RectangularMatrix::SetRow)
        .def("getRow", &phycas::RectangularMatrix::GetRow)
        .def("getMean", &phycas::RectangularMatrix::GetMean)
        .def("getVarCovMatrix", &phycas::RectangularMatrix::GetVarCovMatrix, return_value_policy<manage_new_object>())
		;

//We tell boost::python the smart pointer type we're using, like this:
//class_<DrawableInterface, DrawablePtr, boost::noncopyable>
//("DrawableInterface", no_init)
//    ;
//where DrawablePtr is a typedef for the smart pointer type.

	class_<phycas::Lot, boost::shared_ptr<phycas::Lot>, boost::noncopyable>("LotBase", init<unsigned>())
		.def("getSeed", &phycas::Lot::GetSeed)
		.def("setSeed", &phycas::Lot::SetSeed)
		.def("getInitSeed", &phycas::Lot::GetInitSeed)
		.def("uniform", &phycas::Lot::Uniform)
		.def("getrandbits", &phycas::Lot::GetRandBits)
		.def("sampleUInt", &phycas::Lot::SampleUInt)
		;

	class_<MVNormalDistribution, bases<MultivariateProbabilityDistribution> >("MVNormalDistBase")
		.def(init<const std::vector<double> &, const std::vector<double> &>())
		.def(init<const MVNormalDistribution &>())
		.def("cloneAndSetLot", &MVNormalDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &MVNormalDistribution::Clone, return_value_policy<manage_new_object>())
		.def("fit", &MVNormalDistribution::Fit)
		.def("isDiscrete", &MVNormalDistribution::IsDiscrete)
		.def("getDistName", &MVNormalDistribution::GetDistributionName)
		.def("__str__", &MVNormalDistribution::GetDescriptionForPython)
		.def("__repr__", &MVNormalDistribution::GetDescriptionForPython)
		.def("setLot", &MVNormalDistribution::SetLot)
		.def("setSeed", &MVNormalDistribution::SetSeed)
		.def("resetLot", &MVNormalDistribution::ResetLot)
		.def("getMean", &MVNormalDistribution::GetMean)
		.def("getVar", &MVNormalDistribution::GetVar)
		.def("getStdDev", &MVNormalDistribution::GetStdDev)
		.def("approxCDF", &MVNormalDistribution::ApproxCDF)
		.def("sample", &MVNormalDistribution::Sample)
		.def("getLnPDF", &MVNormalDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &MVNormalDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &MVNormalDistribution::SetMeanAndVariance)
		.def("getVarCovarMatrix", &MVNormalDistribution::GetVarCovarMatrix)
		.def("getNParams", &MVNormalDistribution::GetNParams)
		.def("debugMVNorm", &MVNormalDistribution::DebugMVNorm)
		;

	class_<DirichletDistribution, bases<MultivariateProbabilityDistribution> >("DirichletDistBase")
		.def(init<const std::vector<double> &>())
		.def(init<const DirichletDistribution &>())
		.def("cloneAndSetLot", &DirichletDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &DirichletDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &DirichletDistribution::IsDiscrete)
		.def("getDistName", &DirichletDistribution::GetDistributionName)
		.def("__str__", &DirichletDistribution::GetDescriptionForPython)
		.def("__repr__", &DirichletDistribution::GetDescriptionForPython)
		.def("setLot", &DirichletDistribution::SetLot)
		.def("setSeed", &DirichletDistribution::SetSeed)
		.def("resetLot", &DirichletDistribution::ResetLot)
		.def("getMean", &DirichletDistribution::GetMean)
		.def("getVar", &DirichletDistribution::GetVar)
		.def("getStdDev", &DirichletDistribution::GetStdDev)
		.def("approxCDF", &DirichletDistribution::ApproxCDF)
		.def("sample", &DirichletDistribution::Sample)
		.def("getLnPDF", &DirichletDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &DirichletDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &DirichletDistribution::SetMeanAndVariance)
		.def("getVarCovarMatrix", &DirichletDistribution::GetVarCovarMatrix)
		.def("getNParams", &DirichletDistribution::GetNParams)
		;

	class_<RelativeRateDistribution, bases<MultivariateProbabilityDistribution> >("RelRateDistBase")
		.def(init<const std::vector<double> &, const std::vector<double> &>())
		.def(init<const RelativeRateDistribution &>())
		.def("cloneAndSetLot", &RelativeRateDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &RelativeRateDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &RelativeRateDistribution::IsDiscrete)
		.def("getDistName", &RelativeRateDistribution::GetDistributionName)
		.def("__str__", &RelativeRateDistribution::GetDescriptionForPython)
		.def("__repr__", &RelativeRateDistribution::GetDescriptionForPython)
		.def("setLot", &RelativeRateDistribution::SetLot)
		.def("setSeed", &RelativeRateDistribution::SetSeed)
		.def("resetLot", &RelativeRateDistribution::ResetLot)
		.def("getMean", &RelativeRateDistribution::GetMean)
		.def("getVar", &RelativeRateDistribution::GetVar)
		.def("getStdDev", &RelativeRateDistribution::GetStdDev)
		.def("approxCDF", &RelativeRateDistribution::ApproxCDF)
		.def("sample", &RelativeRateDistribution::Sample)
		.def("getLnPDF", &RelativeRateDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &RelativeRateDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &RelativeRateDistribution::SetMeanAndVariance)
		.def("getVarCovarMatrix", &RelativeRateDistribution::GetVarCovarMatrix)
		.def("getNParams", &RelativeRateDistribution::GetNParams)
		.def("setSubsetProportions", &RelativeRateDistribution::setSubsetProportions)
		//.def("setCoefficients", &RelativeRateDistribution::SetCoefficients)
		;

	class_<BetaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("BetaDistBase")
		.def(init<double, double>())
		.def(init<const BetaDistribution &>())
		.def("cloneAndSetLot", &BetaDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &BetaDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &BetaDistribution::IsDiscrete)
		.def("getDistName", &BetaDistribution::GetDistributionName)
		.def("__str__", &BetaDistribution::GetDistributionDescription)
		.def("__repr__", &BetaDistribution::GetDistributionDescription)
		.def("setLot", &BetaDistribution::SetLot)
		.def("setSeed", &BetaDistribution::SetSeed)
		.def("resetLot", &BetaDistribution::ResetLot)
		.def("getMean", &BetaDistribution::GetMean)
		.def("getVar", &BetaDistribution::GetVar)
		.def("getStdDev", &BetaDistribution::GetStdDev)
		.def("getCDF", &BetaDistribution::GetCDF)
		.def("getQuantile", &BetaDistribution::GetQuantile)
		.def("sample", &BetaDistribution::Sample)
		.def("getLnPDF", &BetaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BetaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BetaDistribution::SetMeanAndVariance)
		;

	class_<BetaPrimeDistribution, bases<ProbabilityDistribution, AdHocDensity> >("BetaPrimeDistBase")
		.def(init<double, double>())
		.def(init<const BetaPrimeDistribution &>())
		.def("cloneAndSetLot", &BetaPrimeDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &BetaPrimeDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &BetaPrimeDistribution::IsDiscrete)
		.def("getDistName", &BetaPrimeDistribution::GetDistributionName)
		.def("__str__", &BetaPrimeDistribution::GetDistributionDescription)
		.def("__repr__", &BetaPrimeDistribution::GetDistributionDescription)
		.def("setLot", &BetaPrimeDistribution::SetLot)
		.def("setSeed", &BetaPrimeDistribution::SetSeed)
		.def("resetLot", &BetaPrimeDistribution::ResetLot)
		.def("getMean", &BetaPrimeDistribution::GetMean)
		.def("getVar", &BetaPrimeDistribution::GetVar)
		.def("getStdDev", &BetaPrimeDistribution::GetStdDev)
		.def("getCDF", &BetaPrimeDistribution::GetCDF)
		.def("getQuantile", &BetaPrimeDistribution::GetQuantile)
		.def("sample", &BetaPrimeDistribution::Sample)
		.def("getLnPDF", &BetaPrimeDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BetaPrimeDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BetaPrimeDistribution::SetMeanAndVariance)
		;
	class_<BernoulliDistribution, bases<ProbabilityDistribution> >("BernoulliDistBase")
		.def(init<double>())
		.def(init<const BernoulliDistribution &>())
		.def("cloneAndSetLot", &BernoulliDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &BernoulliDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &BernoulliDistribution::IsDiscrete)
		.def("getDistName", &BernoulliDistribution::GetDistributionName)
		.def("__str__", &BernoulliDistribution::GetDistributionDescription)
		.def("__repr__", &BernoulliDistribution::GetDistributionDescription)
		.def("setLot", &BernoulliDistribution::SetLot)
		.def("setSeed", &BernoulliDistribution::SetSeed)
		.def("resetLot", &BernoulliDistribution::ResetLot)
		.def("getMean", &BernoulliDistribution::GetMean)
		.def("getVar", &BernoulliDistribution::GetVar)
		.def("getStdDev", &BernoulliDistribution::GetStdDev)
		.def("getCDF", &BernoulliDistribution::GetCDF)
		.def("sample", &BernoulliDistribution::Sample)
		.def("getLnPDF", &BernoulliDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BernoulliDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BernoulliDistribution::SetMeanAndVariance)
		;
	class_<BinomialDistribution, bases<ProbabilityDistribution> >("BinomialDistBase")
		.def(init<double, double>())
		.def(init<const BinomialDistribution &>())
		.def("cloneAndSetLot", &BinomialDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &BinomialDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &BinomialDistribution::IsDiscrete)
		.def("getDistName", &BinomialDistribution::GetDistributionName)
		.def("__str__", &BinomialDistribution::GetDistributionDescription)
		.def("__repr__", &BinomialDistribution::GetDistributionDescription)
		.def("setLot", &BinomialDistribution::SetLot)
		.def("setSeed", &BinomialDistribution::SetSeed)
		.def("resetLot", &BinomialDistribution::ResetLot)
		.def("getMean", &BinomialDistribution::GetMean)
		.def("getVar", &BinomialDistribution::GetVar)
		.def("getStdDev", &BinomialDistribution::GetStdDev)
		.def("getCDF", &BinomialDistribution::GetCDF)
		.def("sample", &BinomialDistribution::Sample)
		.def("getLnPDF", &BinomialDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &BinomialDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &BernoulliDistribution::SetMeanAndVariance)
		;
	class_<ImproperUniformDistribution, bases<ProbabilityDistribution, AdHocDensity> >("ImproperUniformDistBase")
		.def(init<const ImproperUniformDistribution &>())
		.def("cloneAndSetLot", &ImproperUniformDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &ImproperUniformDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &ImproperUniformDistribution::IsDiscrete)
		.def("getDistName", &ImproperUniformDistribution::GetDistributionName)
		.def("__str__", &ImproperUniformDistribution::GetDistributionDescription)
		.def("__repr__", &ImproperUniformDistribution::GetDistributionDescription)
		.def("setLot", &ImproperUniformDistribution::SetLot)
		.def("setSeed", &ImproperUniformDistribution::SetSeed)
		.def("resetLot", &ImproperUniformDistribution::ResetLot)
		.def("getMean", &ImproperUniformDistribution::GetMean)
		.def("getVar", &ImproperUniformDistribution::GetVar)
		.def("getStdDev", &ImproperUniformDistribution::GetStdDev)
		.def("getCDF", &ImproperUniformDistribution::GetCDF)
		.def("sample", &ImproperUniformDistribution::Sample)
		.def("getLnPDF", &ImproperUniformDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &ImproperUniformDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &ImproperUniformDistribution::SetMeanAndVariance)
		;
	class_<UniformDistribution, bases<ProbabilityDistribution, AdHocDensity> >("UniformDistBase")
		.def(init<double, double>())
		.def(init<const UniformDistribution &>())
		.def("cloneAndSetLot", &UniformDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &UniformDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &UniformDistribution::IsDiscrete)
		.def("getDistName", &UniformDistribution::GetDistributionName)
		.def("__str__", &UniformDistribution::GetDistributionDescription)
		.def("__repr__", &UniformDistribution::GetDistributionDescription)
		.def("setLot", &UniformDistribution::SetLot)
		.def("setSeed", &UniformDistribution::SetSeed)
		.def("resetLot", &UniformDistribution::ResetLot)
		.def("getMean", &UniformDistribution::GetMean)
		.def("getVar", &UniformDistribution::GetVar)
		.def("getStdDev", &UniformDistribution::GetStdDev)
		.def("getCDF", &UniformDistribution::GetCDF)
		.def("sample", &UniformDistribution::Sample)
		.def("getLnPDF", &UniformDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &UniformDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &UniformDistribution::SetMeanAndVariance)
		;
	class_<GammaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("GammaDistBase")
		.def(init<double, double>())
		.def(init<const GammaDistribution &>())
		.def("cloneAndSetLot", &GammaDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &GammaDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &GammaDistribution::IsDiscrete)
		.def("getDistName", &GammaDistribution::GetDistributionName)
		.def("__str__", &GammaDistribution::GetDistributionDescription)
		.def("__repr__", &GammaDistribution::GetDistributionDescription)
		.def("setLot", &GammaDistribution::SetLot)
		.def("setSeed", &GammaDistribution::SetSeed)
		.def("resetLot", &GammaDistribution::ResetLot)
		.def("getMean", &GammaDistribution::GetMean)
		.def("getVar", &GammaDistribution::GetVar)
		.def("getStdDev", &GammaDistribution::GetStdDev)
		.def("getCDF", &GammaDistribution::GetCDF)
		.def("sample", &GammaDistribution::Sample)
		.def("getLnPDF", &GammaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &GammaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &GammaDistribution::SetMeanAndVariance)
		;
	class_<ExponentialDistribution, bases<GammaDistribution, ProbabilityDistribution, AdHocDensity> >("ExponentialDistBase")
		.def(init<double>())
		.def(init<const ExponentialDistribution &>())
		.def("cloneAndSetLot", &ExponentialDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &ExponentialDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &ExponentialDistribution::IsDiscrete)
		.def("getDistName", &ExponentialDistribution::GetDistributionName)
		.def("__str__", &ExponentialDistribution::GetDistributionDescription)
		.def("__repr__", &ExponentialDistribution::GetDistributionDescription)
		.def("setLot", &ExponentialDistribution::SetLot)
		.def("setSeed", &ExponentialDistribution::SetSeed)
		.def("resetLot", &ExponentialDistribution::ResetLot)
		.def("getMean", &ExponentialDistribution::GetMean)
		.def("getVar", &ExponentialDistribution::GetVar)
		.def("getStdDev", &ExponentialDistribution::GetStdDev)
		.def("getCDF", &ExponentialDistribution::GetCDF)
		.def("sample", &ExponentialDistribution::Sample)
		.def("getLnPDF", &ExponentialDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &ExponentialDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &ExponentialDistribution::SetMeanAndVariance)
		;
	class_<InverseGammaDistribution, bases<ProbabilityDistribution, AdHocDensity> >("InverseGammaDistBase")
		.def(init<double, double>())
		.def(init<const InverseGammaDistribution &>())
		.def("cloneAndSetLot", &InverseGammaDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &InverseGammaDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &InverseGammaDistribution::IsDiscrete)
		.def("getDistName", &InverseGammaDistribution::GetDistributionName)
		.def("__str__", &InverseGammaDistribution::GetDistributionDescription)
		.def("__repr__", &InverseGammaDistribution::GetDistributionDescription)
		.def("setLot", &InverseGammaDistribution::SetLot)
		.def("setSeed", &InverseGammaDistribution::SetSeed)
		.def("resetLot", &InverseGammaDistribution::ResetLot)
		.def("getMean", &InverseGammaDistribution::GetMean)
		.def("getVar", &InverseGammaDistribution::GetVar)
		.def("getStdDev", &InverseGammaDistribution::GetStdDev)
		.def("getCDF", &InverseGammaDistribution::GetCDF)
		.def("sample", &InverseGammaDistribution::Sample)
		.def("getLnPDF", &InverseGammaDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &InverseGammaDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &InverseGammaDistribution::SetMeanAndVariance)
		;
	class_<NormalDistribution, bases<ProbabilityDistribution, AdHocDensity> >("NormalDistBase")
		.def(init<double, double>())
		.def(init<const NormalDistribution &>())
		.def("cloneAndSetLot", &NormalDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &NormalDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &NormalDistribution::IsDiscrete)
		.def("getDistName", &NormalDistribution::GetDistributionName)
		.def("__str__", &NormalDistribution::GetDistributionDescription)
		.def("__repr__", &NormalDistribution::GetDistributionDescription)
		.def("setLot", &NormalDistribution::SetLot)
		.def("setSeed", &NormalDistribution::SetSeed)
		.def("resetLot", &NormalDistribution::ResetLot)
		.def("getMean", &NormalDistribution::GetMean)
		.def("getVar", &NormalDistribution::GetVar)
		.def("getStdDev", &NormalDistribution::GetStdDev)
		.def("getCDF", &NormalDistribution::GetCDF)
		.def("sample", &NormalDistribution::Sample)
		.def("getLnPDF", &NormalDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &NormalDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &NormalDistribution::SetMeanAndVariance)
		;
	class_<LognormalDistribution, bases<ProbabilityDistribution, AdHocDensity> >("LognormalDistBase")
		.def(init<double, double>())
		.def(init<const LognormalDistribution &>())
		.def("cloneAndSetLot", &LognormalDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &LognormalDistribution::Clone, return_value_policy<manage_new_object>())
		.def("isDiscrete", &LognormalDistribution::IsDiscrete)
		.def("getDistName", &LognormalDistribution::GetDistributionName)
		.def("__str__", &LognormalDistribution::GetDistributionDescription)
		.def("__repr__", &LognormalDistribution::GetDistributionDescription)
		.def("setLot", &LognormalDistribution::SetLot)
		.def("setSeed", &LognormalDistribution::SetSeed)
		.def("resetLot", &LognormalDistribution::ResetLot)
		.def("getMean", &LognormalDistribution::GetMean)
		.def("getVar", &LognormalDistribution::GetVar)
		.def("getStdDev", &LognormalDistribution::GetStdDev)
		.def("getCDF", &LognormalDistribution::GetCDF)
		.def("sample", &LognormalDistribution::Sample)
		.def("getLnPDF", &LognormalDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &LognormalDistribution::GetRelativeLnPDF)
		.def("setMeanAndVariance", &LognormalDistribution::SetMeanAndVariance)
		;
	register_exception_translator<XProbDist>(&translateXProbDist);
}

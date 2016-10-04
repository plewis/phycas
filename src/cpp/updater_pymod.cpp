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

#include <boost/python.hpp>

#include "basic_tree_node.hpp"
#include "tree_likelihood.hpp"
#include "larget_simon_move.hpp"
#include "tree_scaler_move.hpp"
#include "dirichlet_move.hpp"
#include "state_freq_move.hpp"
#include "rel_rates_move.hpp"
#include "subset_relrates_move.hpp"
#include "bush_move.hpp"
#include "edge_move.hpp"
#include "topo_prior_calculator.hpp"

using namespace boost::python;
using namespace phycas;

void updater_pymod()
	{
   	class_<phycas::MCMCUpdater, bases<AdHocDensity>, boost::noncopyable, boost::shared_ptr<phycas::MCMCUpdater> >("MCMCUpdaterBase", no_init)
		.def("isPriorSteward", &MCMCUpdater::isPriorSteward)
		.def("getLnPrior", &MCMCUpdater::getLnPrior)
		.def("getLnLike", &MCMCUpdater::getLnLike)
		.def("computesUnivariatePrior", &MCMCUpdater::computesUnivariatePrior)
		.def("computesMultivariatePrior", &MCMCUpdater::computesMultivariatePrior)
		.def("computesTreeLengthPrior", &MCMCUpdater::computesTreeLengthPrior)
		.def("computesTopologyPrior", &MCMCUpdater::computesTopologyPrior)
		.def("getCurrValueFromModel", &MCMCUpdater::getCurrValueFromModel)
		.def("getCurrValuesFromModel", &MCMCUpdater::getCurrValuesFromModel)
        .def("sendCurrValueToModel", &MCMCUpdater::sendCurrValueToModel)
		.def("listCurrValuesFromModel", &MCMCUpdater::listCurrValuesFromModel)
		.def("useWorkingPrior", &MCMCUpdater::useWorkingPrior)
		.def("setUseWorkingPrior", &MCMCUpdater::setUseWorkingPrior)
		.def("recalcWorkingPrior", &MCMCUpdater::recalcWorkingPrior)
		.def("getWorkingPriorDescr", &MCMCUpdater::getWorkingPriorDescr)
		.def("sampleWorkingPrior", &MCMCUpdater::sampleWorkingPrior)
		.def("sampleMultivariateWorkingPrior", &MCMCUpdater::sampleMultivariateWorkingPrior)
		.def("setWorkingPrior", &MCMCUpdater::setWorkingPrior)
		.def("setMultivariateWorkingPrior", &MCMCUpdater::setMultivariateWorkingPrior)
		.def("setReferenceDistribution", &MCMCUpdater::setReferenceDistribution)
		.def("getName", &MCMCUpdater::getName, return_value_policy<copy_const_reference>())
		.def("getPriorDescr", &MCMCUpdater::getPriorDescr)
		.def("getWeight", &MCMCUpdater::getWeight)
		.def("setName", &MCMCUpdater::setName)
		.def("setWeight", &MCMCUpdater::setWeight)
		.def("getTuningParameter", &MCMCUpdater::getTuningParameter)
		.def("setTuningParameter", &MCMCUpdater::setTuningParameter)
		.def("autotune", &MCMCUpdater::autotune)
		.def("setStartingValue", &MCMCUpdater::setStartingValue)
		.def("setTree", &MCMCUpdater::setTree)
		.def("getTree", &MCMCUpdater::getTree)
		.def("setLot", &MCMCUpdater::setLot)
		.def("setPrior", &MCMCUpdater::setPrior)
		.def("setMultivarPrior", &MCMCUpdater::setMultivarPrior)
		.def("setTreeLikelihood", &MCMCUpdater::setTreeLikelihood)
		.def("setModel", &MCMCUpdater::setModel)
		.def("setChainManager", &MCMCUpdater::setChainManager)
		.def("hasSliceSampler", &MCMCUpdater::hasSliceSampler)
		.def("isParameter", &MCMCUpdater::isParameter)
		.def("isMasterParameter", &MCMCUpdater::isMasterParameter)
		.def("isHyperParameter", &MCMCUpdater::isHyperParameter)
		.def("isMove", &MCMCUpdater::isMove)
		.def("getSliceSampler", &MCMCUpdater::getSliceSampler)
		.def("getNumAttempts", &MCMCUpdater::getNumAttempts)
		.def("getNumAccepts", &MCMCUpdater::getNumAccepts)
		.def("resetDiagnostics", &MCMCUpdater::resetDiagnostics)
		.def("update", &MCMCUpdater::update)
		.def("isFixed", &MCMCUpdater::isFixed)
		.def("fixParameter", &MCMCUpdater::fixParameter)
		.def("freeParameter", &MCMCUpdater::freeParameter)
		.def("getDebugInfo", &MCMCUpdater::getDebugInfo)
		.def("setSaveDebugInfo", &MCMCUpdater::setSaveDebugInfo)
		.def("setPower", &MCMCUpdater::setPower)
		.def("getPower", &MCMCUpdater::getPower)
		.def("setStandardHeating", &MCMCUpdater::setStandardHeating)
		.def("setLikelihoodHeating", &MCMCUpdater::setLikelihoodHeating)
		.def("isStandardHeating", &MCMCUpdater::isStandardHeating)
		.def("isLikelihoodHeating", &MCMCUpdater::isLikelihoodHeating)
		.def("isNoHeating", &MCMCUpdater::isNoHeating)
		;
	class_<phycas::EdgeLenParam, bases<phycas::MCMCUpdater, AdHocDensity>,
		boost::noncopyable, boost::shared_ptr<phycas::EdgeLenParam> >("EdgeLenParam")
		.def("getSplitReprAsString", &phycas::EdgeLenParam::getSplitReprAsString)
		;
	class_<phycas::KappaParam, bases<phycas::MCMCUpdater, AdHocDensity>,
		boost::noncopyable, boost::shared_ptr<phycas::KappaParam> >("KappaParam")
		;
	class_<phycas::GTRRateParam, bases<phycas::MCMCUpdater, AdHocDensity>,
		boost::noncopyable, boost::shared_ptr<phycas::GTRRateParam> >("GTRRateParam", init<unsigned>())
		;
	class_<phycas::StateFreqParam, bases<phycas::MCMCUpdater, AdHocDensity>,
		boost::noncopyable, boost::shared_ptr<phycas::StateFreqParam> >("StateFreqParam", init<unsigned>())
		;
	class_<phycas::HyperPriorParam, bases<phycas::MCMCUpdater, AdHocDensity>,
		boost::noncopyable, boost::shared_ptr<phycas::HyperPriorParam> >("HyperPriorParam")
		;
	class_<phycas::LargetSimonMove, bases<phycas::MCMCUpdater>,
		boost::noncopyable, boost::shared_ptr<phycas::LargetSimonMove> >("LargetSimonMove")
		.def("update", &phycas::LargetSimonMove::update)
		.def("setTuningParameter", &phycas::LargetSimonMove::setTuningParameter)
		.def("getTuningParameter", &phycas::LargetSimonMove::getTuningParameter)
		.def("topologyChanged", &phycas::LargetSimonMove::topologyChanged)
		;
	class_<phycas::TreeScalerMove, bases<phycas::MCMCUpdater>,
		boost::noncopyable, boost::shared_ptr<phycas::TreeScalerMove> >("TreeScalerMove")
		.def("update", &phycas::TreeScalerMove::update)
		.def("setTuningParameter", &phycas::TreeScalerMove::setTuningParameter)
		.def("getTuningParameter", &phycas::TreeScalerMove::getTuningParameter)
		;
	class_<phycas::DirichletMove, bases<phycas::MCMCUpdater>,
		boost::noncopyable, boost::shared_ptr<phycas::DirichletMove> >("DirichletMove")
		.def("setDimension", &phycas::DirichletMove::setDimension)
		.def("setTuningParameter", &phycas::DirichletMove::setTuningParameter)
		.def("getTuningParameter", &phycas::DirichletMove::getTuningParameter)
		;
	class_<phycas::SubsetRelRatesMove, bases<phycas::DirichletMove>,
		boost::noncopyable, boost::shared_ptr<phycas::SubsetRelRatesMove> >("SubsetRelRatesMove")
		.def("setPartitionModel", &phycas::SubsetRelRatesMove::setPartitionModel)
		// .def("setSubsetProportions", &phycas::SubsetRelRatesMove::setSubsetProportions)
		.def("update", &phycas::SubsetRelRatesMove::update)
		;
	class_<phycas::RelRatesMove, bases<phycas::DirichletMove>,
		boost::noncopyable, boost::shared_ptr<phycas::RelRatesMove> >("RelRatesMove")
		.def("update", &phycas::RelRatesMove::update)
		;
	class_<phycas::StateFreqMove, bases<phycas::DirichletMove>,
		boost::noncopyable, boost::shared_ptr<phycas::StateFreqMove> >("StateFreqMove")
		.def("update", &phycas::StateFreqMove::update)
		;
	class_<TopoProbCalculator, boost::noncopyable, boost::shared_ptr<phycas::TopoProbCalculator> >("TopoProbCalculatorBase")
        .def("GetLnTopoProb", &TopoProbCalculator::GetLnTopoProb)
        ;
	class_<FocalTreeTopoProbCalculator, bases<phycas::TopoProbCalculator>, boost::noncopyable,
		boost::shared_ptr<phycas::FocalTreeTopoProbCalculator> >("FocalTreeTopoProbCalculatorBase", init<TreeShPtr>())
		.def("setDefaultEdgeLenDist", &FocalTreeTopoProbCalculator::SetDefaultEdgeLenDist)
		.def("setEdgeLenDist", &FocalTreeTopoProbCalculator::SetEdgeLenDist)
        .def("sampleTree", &FocalTreeTopoProbCalculator::SampleTree)
        .def("calcTopologyLnProb", &FocalTreeTopoProbCalculator::CalcTopologyLnProbVector)
        ;
	class_<PolytomyTopoPriorCalculator, bases<phycas::TopoProbCalculator>, boost::noncopyable,
		boost::shared_ptr<phycas::PolytomyTopoPriorCalculator> >("PolytomyTopoPriorCalculatorBase")
		.def("setNTax", &PolytomyTopoPriorCalculator::SetNTax)
		.def("getNTax", &PolytomyTopoPriorCalculator::GetNTax)
		.def("chooseRooted", &PolytomyTopoPriorCalculator::ChooseRooted)
		.def("chooseUnrooted", &PolytomyTopoPriorCalculator::ChooseUnrooted)
        .def("getLnCount", &PolytomyTopoPriorCalculator::GetLnCount)
		.def("getLnSaturatedCount", &PolytomyTopoPriorCalculator::GetLnSaturatedCount)
		.def("getLnTotalCount", &PolytomyTopoPriorCalculator::GetLnTotalCount)
		.def("getNFactorsVect", &PolytomyTopoPriorCalculator::GetNFactorsVect)
		.def("getLnCounts", &PolytomyTopoPriorCalculator::GetLnCounts)
		.def("setLnScalingFactor", &PolytomyTopoPriorCalculator::SetLnScalingFactor)
		.def("getLnScalingFactor", &PolytomyTopoPriorCalculator::GetLnScalingFactor)
		.def("getCountsVect", &PolytomyTopoPriorCalculator::GetCountsVect)
		.def("chooseResolutionClassPrior", &PolytomyTopoPriorCalculator::ChooseResolutionClassPrior)
		.def("choosePolytomyPrior", &PolytomyTopoPriorCalculator::ChoosePolytomyPrior)
		.def("setC", &PolytomyTopoPriorCalculator::SetC)
		.def("getC", &PolytomyTopoPriorCalculator::GetC)
		.def("getLnTopologyPrior", &PolytomyTopoPriorCalculator::GetLnTopologyPrior)
		.def("getLnNormalizedTopologyPrior", &PolytomyTopoPriorCalculator::GetLnNormalizedTopologyPrior)
		.def("GetLnNormConstant", &PolytomyTopoPriorCalculator::GetLnNormConstant)
		.def("getTopoPriorVect", &PolytomyTopoPriorCalculator::GetTopoPriorVect)
		.def("isResolutionClassPrior", &PolytomyTopoPriorCalculator::IsResolutionClassPrior)
		.def("isPolytomyPrior", &PolytomyTopoPriorCalculator::IsPolytomyPrior)
		.def("isRooted", &PolytomyTopoPriorCalculator::IsRooted)
		.def("isUnrooted", &PolytomyTopoPriorCalculator::IsUnrooted)
		.def("getRealizedResClassPriorsVect", &PolytomyTopoPriorCalculator::GetRealizedResClassPriorsVect)
		.def("sample", &PolytomyTopoPriorCalculator::sample)
		;
	class_<phycas::BushMove, bases<phycas::MCMCUpdater>,
		boost::noncopyable, boost::shared_ptr<phycas::BushMove> >("BushMove")
		.def("update", &phycas::BushMove::update)
		.def("addEdgeMoveProposed", &phycas::BushMove::addEdgeMoveProposed)
		.def("setEdgeLenDistMean", &phycas::BushMove::setEdgeLenDistMean)
		.def("finalize", &phycas::BushMove::finalize)
		;
	class_<phycas::EdgeMove, bases<phycas::MCMCUpdater>,
		boost::noncopyable, boost::shared_ptr<phycas::EdgeMove> >("EdgeMove")
		.def("update", &phycas::EdgeMove::update)
		.def("setTuningParameter", &phycas::EdgeMove::setTuningParameter)
		.def("getTuningParameter", &phycas::EdgeMove::getTuningParameter)
		;
	}

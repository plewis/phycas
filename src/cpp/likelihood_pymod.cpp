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

//#include "phycas/force_include.h"
#include "probability_distribution.hpp"
#include "tree_length_distribution.hpp"
#include "basic_tree.hpp"
#include "model.hpp"
#include "tree_likelihood.hpp"
#include "tip_data.hpp"
#include "internal_data.hpp"
#include "mcmc_param.hpp"
#include "mcmc_chain_manager.hpp"
//#include "phycas/src/topo_prior_calculator.hpp"
//#include "phycas/src/larget_simon_move.hpp"
//#include "phycas/src/ncat_move.hpp"
//#include "phycas/src/bush_move.hpp"
//#include "phycas/src/edge_move.hpp"
#include "sim_data.hpp"
#include "q_matrix.hpp"
#include "xlikelihood.hpp"
#include "partition_model.hpp"
#include "char_super_matrix.hpp"
#include "joint_prior_manager.hpp"
#include "numsum.hpp"
void model_pymod();
void updater_pymod();

using namespace boost::python;
using namespace phycas;

void translateXLikelihood(const XLikelihood &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

// TreeLikelihoodWrapper is necessary in order to allow this (in C++ code)
//
// startTreeViewer(t);
//
// to be translated into a call of TreeLikelihoodDerived.startTreeViewer, where
// TreeLikelihoodDerived is a Python class derived from TreeLikelihoodBase.
// (TreeLikelihoodBase is what TreeLikelihood is called within Python, search for
// TreeLikelihoodBase below.) For more background, go to
// http://www.boost.org/libs/python/doc/index.html, click on the "Reference Manual"
// link, then see the example at the bottom of the page that results from clicking
// on the "call_method.hpp" link (under the category "Function Invocation and
// Creation").

class TreeLikelihoodWrapper : public TreeLikelihood
	{
	public:
		TreeLikelihoodWrapper(PyObject * self, PartitionModelShPtr m) : TreeLikelihood(m), m_self(self)
			{
			}
		virtual ~TreeLikelihoodWrapper() {} //@POL may need to delete uMat here (see TreeLikelihood::~TreeLikelihood)

		int startTreeViewer(TreeShPtr t, std::string msg, unsigned site) const
			{
			return call_method<int,TreeShPtr,std::string>(m_self, "startTreeViewer", t, msg, site);
			}

	private:
		PyObject * const m_self;
	};

BOOST_PYTHON_MODULE(_LikelihoodExt)
{
#if defined(USING_NUMARRAY)
#error USING_NUMARRAY defined
	// these lines required by num_util
	import_array();
	numeric::array::set_module_and_type("numarray", "NDArray");
#endif

	class_<AdHocDensity, boost::noncopyable, boost::shared_ptr<AdHocDensity> >("AdHocDensityBase", no_init)
		;
	class_<phycas::JointPriorManager, boost::noncopyable, boost::shared_ptr<phycas::JointPriorManager> >("JointPriorManagerBase")
		.def("recalcLogJointPrior", &JointPriorManager::recalcLogJointPrior)
		.def("getLogJointPrior", &JointPriorManager::getLogJointPrior)
		.def("getLogTopologyPrior", &JointPriorManager::getLogTopologyPrior)
		.def("debugPriorBreakdown", &JointPriorManager::debugPriorBreakdown)
		.def("addUnivariateDistribution", &JointPriorManager::addUnivariateDistribution)
		.def("addMultivariateDistribution", &JointPriorManager::addMultivariateDistribution)
		.def("getTopoProbCalculator", &JointPriorManager::getTopoProbCalculator)
		.def("addTopologyDistribution", &JointPriorManager::addTopologyDistribution)
		.def("addTreeLengthDistribution", &JointPriorManager::addTreeLengthDistribution)
		.def("addInternalEdgelenDistribution", &JointPriorManager::addInternalEdgelenDistribution)
		.def("addExternalEdgelenDistribution", &JointPriorManager::addExternalEdgelenDistribution)
		.def("addEdgelenHyperprior", &JointPriorManager::addEdgelenHyperprior)
		.def("univariateModified", &JointPriorManager::univariateModified)
		.def("univariateModifiedNoDebugCheck", &JointPriorManager::univariateModifiedNoDebugCheck)
		.def("multivariateModified", &JointPriorManager::multivariateModified)
		.def("externalEdgeLensModified", &JointPriorManager::externalEdgeLensModified)
		.def("internalEdgeLensModified", &JointPriorManager::internalEdgeLensModified)
        .def("allEdgeLensModified", &JointPriorManager::allEdgeLensModified)
        .def("treeLengthModified", &JointPriorManager::treeLengthModified)
		.def("edgeLenHyperparamModified", &JointPriorManager::edgeLenHyperparamModified)
		.def("topologyModified", &JointPriorManager::topologyModified)
		.def("isTreeLengthPrior", &JointPriorManager::isTreeLengthPrior)
        ;
	class_<phycas::MCMCChainManager, boost::noncopyable, boost::shared_ptr<phycas::MCMCChainManager> >("MCMCChainManagerBase", init<JointPriorManagerShPtr>())
		.def("finalize", &MCMCChainManager::finalize)
		.def("recalcLnWorkingPrior", &MCMCChainManager::recalcLnWorkingPrior)
		.def("addMove", &MCMCChainManager::addMove)
		.def("addModelParam", &MCMCChainManager::addModelParam)
		.def("addEdgeLenParam", &MCMCChainManager::addEdgeLenParam)
		.def("addEdgeLenHyperparam", &MCMCChainManager::addEdgeLenHyperparam)
		.def("setEdgeLenHyperparam", &MCMCChainManager::setEdgeLenHyperparam)
		.def("praxisLocatePosteriorMode", &MCMCChainManager::praxisLocatePosteriorMode)
		.def("getLastLnLike", &MCMCChainManager::getLastLnLike)
        .def("getAllUpdaters", &MCMCChainManager::getAllUpdaters, return_value_policy<copy_const_reference>())
        .def("updateAllUpdaters", &MCMCChainManager::updateAllUpdaters)
        .def("getMoves", &MCMCChainManager::getMoves, return_value_policy<copy_const_reference>())
        .def("getModelParams", &MCMCChainManager::getModelParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenParams", &MCMCChainManager::getEdgeLenParams, return_value_policy<copy_const_reference>())
        .def("getEdgeLenHyperparams", &MCMCChainManager::getEdgeLenHyperparams)
		.def("addMCMCUpdaters", &MCMCChainManager::addMCMCUpdaters)
		.def("clear", &MCMCChainManager::clear)
		.def("refreshLastLnLike", &MCMCChainManager::refreshLastLnLike)
		.def("setRefTree", &MCMCChainManager::setRefTree)
		.def("getRefTree", &MCMCChainManager::getRefTree)
		.def("calcRFDistance", &MCMCChainManager::calcRFDistance)
		.def("getJointPriorManager", &MCMCChainManager::getJointPriorManager)
		;
	class_<std::vector<MCMCUpdaterShPtr> >("paramVec", no_init)
		.def("__iter__",  iterator<std::vector<MCMCUpdaterShPtr> >())
		;
	class_<phycas::SimData, boost::noncopyable, boost::shared_ptr<phycas::SimData> >("SimDataBase")
		.def("clear", &phycas::SimData::clear)
		.def("zeroCounts", &phycas::SimData::zeroCounts)
		.def("createMapleTuples", &phycas::SimData::createMapleTuples)
		.def("appendCountsToFile", &phycas::SimData::appendCountsToFile)
		.def("getPatterns", &phycas::SimData::getPatterns)
		.def("getNUniquePatterns", &phycas::SimData::getNUniquePatterns)
		.def("resetPatternLength", &phycas::SimData::resetPatternLength)
		.def("wipePattern", &phycas::SimData::wipePattern)
		.def("setState", &phycas::SimData::setState)
		.def("insertPattern", &phycas::SimData::insertPattern)
		.def("saveToNexusFile", &phycas::SimData::saveToNexusFile)
		.def("getPatternLength", &phycas::SimData::getPatternLength)
		.def("getPatternVectRow", &phycas::SimData::getPatternVectRow)
		.def("patternTable", &phycas::SimData::patternTable)
		.def("divideBy", &phycas::SimData::divideBy)
		.def("addDataTo", &phycas::SimData::addDataTo)
		.def("resizePatternVect", &phycas::SimData::resizePatternVect)
		.def("calct", &phycas::SimData::calct)
		.def("multBy", &phycas::SimData::multBy)
		.def("debugAppendCountsToFile", &phycas::SimData::debugAppendCountsToFile)
		.def("calctBinned", &phycas::SimData::calctBinned)
		.def("buildBinVector", &phycas::SimData::buildBinVector)
		.def("getTotalCount", &phycas::SimData::getTotalCount)
        .def("getBinnedCounts", &phycas::SimData::getBinnedCounts)
		;
	class_<PartitionModel, boost::noncopyable, boost::shared_ptr<PartitionModel> >("PartitionModelBase")
		.def("addModel", &phycas::PartitionModel::addModel)
		.def("setModelsVect", &phycas::PartitionModel::setModelsVect)
		.def("getModel", &phycas::PartitionModel::getModel)
		.def("getModelsVect", &phycas::PartitionModel::getModelsVect, return_value_policy<copy_const_reference>())
		.def("getNumStatesVect", &phycas::PartitionModel::getNumStatesVect, return_value_policy<copy_const_reference>())
		.def("getNumRatesVect", &phycas::PartitionModel::getNumRatesVect, return_value_policy<copy_const_reference>())
		.def("getNumPatternsVect", &phycas::PartitionModel::getNumPatternsVect, return_value_policy<copy_const_reference>())
		.def("getNumSitesVect", &phycas::PartitionModel::getNumSitesVect, return_value_policy<copy_const_reference>())
		.def("getTotalNumPatterns", &phycas::PartitionModel::getTotalNumPatterns)
		.def("getSubsetRelRate", &phycas::PartitionModel::getSubsetRelRate)
		.def("getNumSubsets", &phycas::PartitionModel::getNumSubsets)
		.def("getSubsetRelRatePrior", &phycas::PartitionModel::getSubsetRelRatePrior)
		.def("setSubsetRelRatePrior", &phycas::PartitionModel::setSubsetRelRatePrior)
		.def("setSubsetRelRatesVect", &phycas::PartitionModel::setSubsetRelRatesVect)
		.def("setSiteAssignments", &phycas::PartitionModel::setSiteAssignments)
		.def("setNumSitesVect", &phycas::PartitionModel::setNumSitesVect)
		.def("getSiteAssignments", &phycas::PartitionModel::getSiteAssignments, return_value_policy<copy_const_reference>())
		.def("getFreeParameterNames", &phycas::PartitionModel::getFreeParameterNames)
		.def("getAllParameterNames", &phycas::PartitionModel::getAllParameterNames)
		.def("getNumFreeParameters", &phycas::PartitionModel::getNumFreeParameters)
		.def("getUntransformedParameters", &phycas::PartitionModel::getUntransformedParameters)
		.def("getTransformedParameters", &phycas::PartitionModel::getTransformedParameters)
		.def("setTransformedParameters", &phycas::PartitionModel::setTransformedParameters)
		.def("getLogDetJacobian", &phycas::PartitionModel::getLogDetJacobian)
		.def("getSubsetProportions", &phycas::PartitionModel::getSubsetProportions)
		.def("getJointPriorManager", &phycas::PartitionModel::getJointPriorManager)
		.def("createJointPriorManager", &phycas::PartitionModel::createJointPriorManager)
		;
	class_<TreeLengthDistribution, boost::noncopyable, boost::shared_ptr<TreeLengthDistribution> >("TreeLengthDistBase")
		.def(init<double, double, double, double>())
		.def(init<const TreeLengthDistribution &>())
		.def("cloneAndSetLot", &TreeLengthDistribution::cloneAndSetLot, return_value_policy<manage_new_object>())
		.def("clone", &TreeLengthDistribution::Clone, return_value_policy<manage_new_object>())
		.def("getDistName", &TreeLengthDistribution::GetDistributionName)
		.def("__str__", &TreeLengthDistribution::GetDistributionDescription)
		.def("__repr__", &TreeLengthDistribution::GetDistributionDescription)
		.def("setLot", &TreeLengthDistribution::SetLot)
		.def("setSeed", &TreeLengthDistribution::SetSeed)
		.def("resetLot", &TreeLengthDistribution::ResetLot)
		.def("sample", &TreeLengthDistribution::Sample)
		.def("getLnPDF", &TreeLengthDistribution::GetLnPDF)
		.def("getRelativeLnPDF", &TreeLengthDistribution::GetRelativeLnPDF)
		.def("getShape", &TreeLengthDistribution::getShape)
		.def("getScale", &TreeLengthDistribution::getScale)
		.def("getExtEdgelenParam", &TreeLengthDistribution::getExtEdgelenParam)
		.def("getIntExtEdgelenRatio", &TreeLengthDistribution::getIntExtEdgelenRatio)
		;
	class_<NumSum, boost::noncopyable, boost::shared_ptr<NumSum> >("NumSumBase")
		.def("clear", &NumSum::clear)
		.def("addValueToSum", &NumSum::addValueToSum)
		.def("getLogSumOfSquares", &NumSum::getLogSumOfSquares)
		.def("getLogSum", &NumSum::getLogSum)
		.def("getNSummed", &NumSum::getNSummed)
		.def("getMaxLogValue", &NumSum::getMaxLogValue)
		;
    class_<TreeLikelihood, TreeLikelihoodWrapper, boost::noncopyable>("TreeLikelihoodBase", init<boost::shared_ptr<PartitionModel> >())
        .def("calcLogLikeAtSubstitutionSaturation", &TreeLikelihood::calcLogLikeAtSubstitutionSaturation)
		.def("getPatternCounts", &TreeLikelihood::getPatternCounts, return_value_policy<copy_const_reference>())
		.def("getCharIndexToPatternIndex", &TreeLikelihood::getCharIndexToPatternIndex, return_value_policy<copy_const_reference>())
		.def("getSiteLikelihoods", &TreeLikelihood::getSiteLikelihoods, return_value_policy<copy_const_reference>())
		.def("getSiteUF", &TreeLikelihood::getSiteUF, return_value_policy<copy_const_reference>())
		.def("storingSiteLikelihoods", &TreeLikelihood::storingSiteLikelihoods)
		.def("storeSiteLikelihoods", &TreeLikelihood::storeSiteLikelihoods)
		//.def("storeAllCLAs", &TreeLikelihood::storeAllCLAs, return_value_policy<manage_new_object>())
		.def("copyDataFromDiscreteMatrix", &TreeLikelihood::copyDataFromDiscreteMatrix)
		.def("copyDataFromSimData", &TreeLikelihood::copyDataFromSimData)
		.def("prepareForSimulation", &TreeLikelihood::prepareForSimulation)
		.def("prepareForLikelihood", &TreeLikelihood::prepareForLikelihood)
		.def("addOrphanTip", &TreeLikelihood::addOrphanTip)
		.def("addDecoratedInternalNode", &TreeLikelihood::addDecoratedInternalNode)
		.def("sumPatternCounts", &TreeLikelihood::sumPatternCounts)
		.def("replaceModel", &TreeLikelihood::replacePartitionModel)
		.def("invalidateAwayFromNode", &TreeLikelihood::invalidateAwayFromNode)
		.def("calcLnLFromNode", &TreeLikelihood::calcLnLFromNode)
		.def("calcLnL", &TreeLikelihood::calcLnL)
		.def("simulateFirst", &TreeLikelihood::simulateFirst)
		.def("simulate", &TreeLikelihood::simulate)
		.def("listPatterns", &TreeLikelihood::listPatterns)
		.def("resetNumLikelihoodEvals", &TreeLikelihood::resetNumLikelihoodEvals)
		.def("getNumLikelihoodEvals", &TreeLikelihood::getNumLikelihoodEvals)
		.def("incrementNEvals", &TreeLikelihood::incrementNumLikelihoodEvals)
		.def("addDataTo", &TreeLikelihood::addDataTo)
		.def("recalcRelativeRates", &TreeLikelihood::recalcRelativeRates)
		.def("getCategoryLowerBoundaries", &TreeLikelihood::getCategoryLowerBoundaries)
		.def("getRateMeans", &TreeLikelihood::getRateMeans, return_value_policy<copy_const_reference>())
		.def("getRateProbs", &TreeLikelihood::getRateProbs, return_value_policy<copy_const_reference>())
		.def("getListOfAllMissingSites", &TreeLikelihood::getListOfAllMissingSites, return_value_policy<copy_const_reference>())
		.def("setNoData", &TreeLikelihood::setNoData)
		.def("setHaveData", &TreeLikelihood::setHaveData)
		.def("getLikelihoodRootNodeNum", &TreeLikelihood::getLikelihoodRootNodeNum)
		.def("setUFNumEdges", &TreeLikelihood::setUFNumEdges)
		.def("bytesPerCLA", &TreeLikelihood::bytesPerCLA)
		.def("numCLAsCreated", &TreeLikelihood::numCLAsCreated)
		.def("numCLAsStored", &TreeLikelihood::numCLAsStored)
		.def("debugCheckForUncachedCLAs", &TreeLikelihood::debugCheckForUncachedCLAs)
		.def("getNPatterns", &TreeLikelihood::getNumPatterns)
		.def("getNoData", &TreeLikelihood::getNoData)
		.def("getTreeLengthPrior", &TreeLikelihood::getTreeLengthPrior)
		.def("getTreeLengthRefDist", &TreeLikelihood::getTreeLengthRefDist)
		.def("setTreeLengthPrior", &TreeLikelihood::setTreeLengthPrior)
		.def("setTreeLengthRefDist", &TreeLikelihood::setTreeLengthRefDist)
		//.def("educateTreeLenWorkingPrior", &TreeLikelihood::educateTreeLenWorkingPrior)
		//.def("finalizeTreeLenWorkingPrior", &TreeLikelihood::finalizeTreeLenWorkingPrior)
        .def("debugUncompressedDataInfo", &TreeLikelihood::debugUncompressedDataInfo)
		;
	class_<TipData, boost::noncopyable>("TipData", no_init)
		.def("parentalCLAValid", &TipData::parentalCLAValid)
		.def("parentalCLACached", &TipData::parentalCLACached)
		;
	class_<InternalData, boost::noncopyable>("InternalData", no_init)
		.def("filialCLAValid", &InternalData::filialCLAValid)
		.def("filialCLACached", &InternalData::filialCLACached)
		.def("parentalCLAValid", &InternalData::parentalCLAValid)
		.def("parentalCLACached", &InternalData::parentalCLACached)
		;
#if 0
	class_<QMatrix, boost::noncopyable>("QMatrixBase")
		.def("getDimension", &QMatrix::getDimension)
		.def("setRelativeRates", &QMatrix::setRelativeRates)
		.def("setStateFreqs", &QMatrix::setStateFreqs)
		.def("getPMatrix", &QMatrix::getPMatrix)
		.def("getQMatrix", &QMatrix::getQMatrix)
		.def("getEigenVectors", &QMatrix::getEigenVectors)
		.def("getEigenValues", &QMatrix::getEigenValues)
		;
#endif

	// This function call necessary to avoid fatal error C1204: Compiler limit: internal structure overflow
	// in VC 7.1 (see http://www.boost.org/libs/python/doc/v2/faq.html#c1204)
	model_pymod();
	updater_pymod();

	register_exception_translator<XLikelihood>(&translateXLikelihood);
}

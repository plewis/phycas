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

//#include "phycas/force_include.h"
#include "model.hpp"
#include "gtr.hpp"
#include "jc.hpp"
#include "hky.hpp"
#include "codon_model.hpp"

using namespace boost::python;
using namespace phycas;

void model_pymod()
	{
#if 1
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
	class_<phycas::Model, boost::noncopyable, boost::shared_ptr<phycas::Model> >("Model", no_init)
		.def("fixEdgeLenHyperprior", &phycas::Model::fixEdgeLenHyperprior)
		.def("freeEdgeLenHyperprior", &phycas::Model::freeEdgeLenHyperprior)
		.def("fixEdgeLengths", &phycas::Model::fixEdgeLengths)
		.def("freeEdgeLengths", &phycas::Model::freeEdgeLengths)
		.def("fixStateFreqs", &phycas::Model::fixStateFreqs)
		.def("freeStateFreqs", &phycas::Model::freeStateFreqs)
		.def("fixShape", &phycas::Model::fixShape)
		.def("freeShape", &phycas::Model::freeShape)
		.def("shapeFixed", &phycas::Model::shapeFixed)
		.def("pinvarFixed", &phycas::Model::pinvarFixed)
		.def("stateFreqsFixed", &phycas::Model::stateFreqsFixed)
		.def("edgeLenHyperParamFixed", &phycas::Model::edgeLenHyperParamFixed)
		.def("edgeLengthsFixed", &phycas::Model::edgeLengthsFixed)
		.def("isPinvarModel", &phycas::Model::isPinvarModel)
		.def("setPinvarModel", &phycas::Model::setPinvarModel)
		.def("setNotPinvarModel", &phycas::Model::setNotPinvarModel)
		.def("fixPinvar", &phycas::Model::fixPinvar)
		.def("freePinvar", &phycas::Model::freePinvar)
		.def("getPinvarPrior", &phycas::Model::getPinvarPrior)
		.def("setPinvarPrior", &phycas::Model::setPinvarPrior)
		.def("getDiscreteGammaShapePrior", &phycas::Model::getDiscreteGammaShapePrior)
		.def("setDiscreteGammaShapePrior", &phycas::Model::setDiscreteGammaShapePrior)
		.def("getTreeLengthPrior", &phycas::Model::getTreeLengthPrior)
		.def("setTreeLengthPrior", &phycas::Model::setTreeLengthPrior)
		.def("getExternalEdgeLenPrior", &phycas::Model::getExternalEdgeLenPrior)
		.def("setExternalEdgeLenPrior", &phycas::Model::setExternalEdgeLenPrior)
		.def("getInternalEdgeLenPrior", &phycas::Model::getInternalEdgeLenPrior)
		.def("setInternalEdgeLenPrior", &phycas::Model::setInternalEdgeLenPrior)
		.def("isSeparateInternalExternalEdgeLenPriors", &phycas::Model::isSeparateInternalExternalEdgeLenPriors)
		.def("setEdgeSpecificParams", &phycas::Model::setEdgeSpecificParams)
		.def("separateInternalExternalEdgeLenPriors", &phycas::Model::separateInternalExternalEdgeLenPriors)
		.def("hasEdgeLenHyperPrior", &phycas::Model::hasEdgeLenHyperPrior)
		.def("getInternalEdgelenHyperparam", &phycas::Model::getInternalEdgelenHyperparam)
		.def("getExternalEdgelenHyperparam", &phycas::Model::getExternalEdgelenHyperparam)
		.def("getEdgeLenHyperPrior", &phycas::Model::getEdgeLenHyperPrior)
		.def("setEdgeLenHyperPrior", &phycas::Model::setEdgeLenHyperPrior)
		.def("getModelName", &phycas::Model::getModelName)
		.def("getPinvar", &phycas::Model::getPinvar)
		.def("setPinvar", &phycas::Model::setPinvar)
		.def("getShape", &phycas::Model::getShape)
		.def("setShape", &phycas::Model::setShape)
		.def("getNGammaRates", &phycas::Model::getNGammaRates)
		.def("setNGammaRates", &phycas::Model::setNGammaRates)
		.def("getGammaRateProbs", &phycas::Model::getGammaRateProbs, return_value_policy<copy_const_reference>())
		.def("setAllGammaRateProbsEqual", &phycas::Model::setAllGammaRateProbsEqual)
		.def("getPMatrix", &phycas::Model::getPMatrix)
        .def("setStateFreqsUnnorm", &phycas::Model::setStateFreqsUnnorm)
        .def("setStateFreqUnnorm", &phycas::Model::setStateFreqUnnorm)
        .def("getTimeStamp", &phycas::Model::getTimeStamp)
		;
    class_<phycas::Irreversible, bases<phycas::Model> >("IrreversibleModelBase")
		.def("getModelName", &phycas::Irreversible::getModelName)
		.def("getNStates", &phycas::Irreversible::getNumStates)
		.def("getStateFreqs", &phycas::Irreversible::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("fixScalingFactor", &phycas::Irreversible::fixScalingFactor)
		.def("freeScalingFactor", &phycas::Irreversible::freeScalingFactor)
		.def("getScalingFactor", &phycas::Irreversible::getScalingFactor)
		.def("setScalingFactor", &phycas::Irreversible::setScalingFactor)
		.def("getScalingFactorPrior", &phycas::Irreversible::getScalingFactorPrior)
		.def("setScalingFactorPrior", &phycas::Irreversible::setScalingFactorPrior)
		.def("setGainOnly", &phycas::Irreversible::setGainOnly)
		.def("setLossOnly", &phycas::Irreversible::setLossOnly)
		.def("paramHeader", &phycas::Irreversible::paramHeader)
		.def("paramReport", &phycas::Irreversible::paramReport)
		;
    class_<phycas::Binary, bases<phycas::Model> >("BinaryModelBase")
		.def("getModelName", &phycas::Binary::getModelName)
		.def("getNStates", &phycas::Binary::getNumStates)
		.def("getStateFreqs", &phycas::Binary::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::Binary::setAllFreqsEqual)
		.def("fixScalingFactor", &phycas::Binary::fixScalingFactor)
		.def("freeScalingFactor", &phycas::Binary::freeScalingFactor)
		.def("getScalingFactor", &phycas::Binary::getScalingFactor)
		.def("setScalingFactor", &phycas::Binary::setScalingFactor)
		.def("getScalingFactorPrior", &phycas::Binary::getScalingFactorPrior)
		.def("setScalingFactorPrior", &phycas::Binary::setScalingFactorPrior)
		.def("fixKappa", &phycas::Binary::fixKappa)
		.def("freeKappa", &phycas::Binary::freeKappa)
		.def("getKappa", &phycas::Binary::getKappa)
		.def("setKappa", &phycas::Binary::setKappa)
		.def("getKappaPrior", &phycas::Binary::getKappaPrior)
		.def("setKappaPrior", &phycas::Binary::setKappaPrior)
		.def("paramHeader", &phycas::Binary::paramHeader)
		.def("paramReport", &phycas::Binary::paramReport)
		;
    class_<phycas::JC, bases<phycas::Model> >("JCModelBase")
		.def("getModelName", &phycas::JC::getModelName)
		.def("getNStates", &phycas::JC::getNumStates)
		.def("getStateFreqs", &phycas::JC::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::JC::setAllFreqsEqual)
		.def("paramHeader", &phycas::JC::paramHeader)
		.def("paramReport", &phycas::JC::paramReport)
		;
	class_<phycas::HKY, bases<phycas::Model> >("HKYModelBase")
		.def("getModelName", &phycas::HKY::getModelName)
		.def("fixKappa", &phycas::HKY::fixKappa)
		.def("freeKappa", &phycas::HKY::freeKappa)
		.def("getKappa", &phycas::HKY::getKappa)
		.def("setKappa", &phycas::HKY::setKappa)
		.def("getKappaPrior", &phycas::HKY::getKappaPrior)
		.def("setKappaPrior", &phycas::HKY::setKappaPrior)
        .def("getStateFreqPrior", &phycas::HKY::getStateFreqPrior)
		.def("setStateFreqPrior", &phycas::HKY::setStateFreqPrior)
        .def("getStateFreqParamPrior", &phycas::HKY::getStateFreqParamPrior)
		.def("setStateFreqParamPrior", &phycas::HKY::setStateFreqParamPrior)
		.def("setKappaFromTRatio", &phycas::HKY::setKappaFromTRatio)
		.def("calcTRatio", &phycas::HKY::calcTRatio)
		.def("getNStates", &phycas::HKY::getNumStates)
		.def("getStateFreqs", &phycas::HKY::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::HKY::setAllFreqsEqual)
		.def("setNucleotideFreqs", &phycas::HKY::setNucleotideFreqs)
		.def("paramHeader", &phycas::HKY::paramHeader)
		.def("paramReport", &phycas::HKY::paramReport)
		;
	class_<phycas::GTR, bases<phycas::Model> >("GTRModelBase")
		.def("getModelName", &phycas::GTR::getModelName)
		//.def("createParameters", &phycas::GTR::createParameters)
		//.def("calcPMat", &phycas::GTR::calcPMat)
		.def("fixRelRates", &phycas::GTR::fixRelRates)
		.def("freeRelRates", &phycas::GTR::freeRelRates)
		.def("getRelRates", &phycas::GTR::getRelRates)
		.def("setRelRates", &phycas::GTR::setRelRates)
		.def("setRelRateUnnorm", &phycas::GTR::setRelRateUnnorm)
		.def("setRelRateParamPrior", &phycas::GTR::setRelRateParamPrior)
		.def("getRelRateParamPrior", &phycas::GTR::getRelRateParamPrior)
		.def("setRelRatePrior", &phycas::GTR::setRelRatePrior)
		.def("getRelRatePrior", &phycas::GTR::getRelRatePrior)
		.def("setNucleotideFreqs", &phycas::GTR::setNucleotideFreqs)
		.def("getStateFreqs", &phycas::GTR::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setAllFreqsEqual", &phycas::GTR::setAllFreqsEqual)
        .def("getStateFreqPrior", &phycas::GTR::getStateFreqPrior)
		.def("setStateFreqPrior", &phycas::GTR::setStateFreqPrior)
		.def("setStateFreqParamPrior", &phycas::GTR::setStateFreqParamPrior)
		.def("getStateFreqParamPrior", &phycas::GTR::getStateFreqParamPrior)
		.def("paramHeader", &phycas::GTR::paramHeader)
		.def("paramReport", &phycas::GTR::paramReport)
		.def("calcTRatio", &phycas::GTR::calcTRatio)
		;
	class_<phycas::Codon, bases<phycas::Model> >("CodonModelBase")
		.def("getModelName", &phycas::Codon::getModelName)
		.def("getNStates", &phycas::Codon::getNumStates)
		.def("getStateFreqs", &phycas::Codon::getStateFreqs, return_value_policy<copy_const_reference>())
		.def("setStateFreqUnnorm", &phycas::Codon::setStateFreqUnnorm)
		.def("setAllFreqsEqual", &phycas::Codon::setAllFreqsEqual)
		.def("setNucleotideFreqs", &phycas::Codon::setNucleotideFreqs)
		.def("fixKappa", &phycas::Codon::fixKappa)
		.def("freeKappa", &phycas::Codon::freeKappa)
		.def("fixOmega", &phycas::Codon::fixOmega)
		.def("freeOmega", &phycas::Codon::freeOmega)
		.def("getKappa", &phycas::Codon::getKappa)
		.def("setKappa", &phycas::Codon::setKappa)
		.def("getOmega", &phycas::Codon::getOmega)
		.def("setOmega", &phycas::Codon::setOmega)
		.def("getNGammaRates", &phycas::Codon::getNGammaRates)
		.def("setNGammaRates", &phycas::Codon::setNGammaRates)
		.def("getKappaPrior", &phycas::Codon::getKappaPrior)
		.def("setKappaPrior", &phycas::Codon::setKappaPrior)
		.def("getOmegaPrior", &phycas::Codon::getOmegaPrior)
		.def("setOmegaPrior", &phycas::Codon::setOmegaPrior)
        .def("getStateFreqPrior", &phycas::Codon::getStateFreqPrior)
		.def("setStateFreqPrior", &phycas::Codon::setStateFreqPrior)
		.def("setStateFreqParamPrior", &phycas::Codon::setStateFreqParamPrior)
		.def("getStateFreqParamPrior", &phycas::Codon::getStateFreqParamPrior)
		.def("paramHeader", &phycas::Codon::paramHeader)
		.def("paramReport", &phycas::Codon::paramReport)
		;
	}

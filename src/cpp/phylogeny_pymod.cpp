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
#include "basic_tree.hpp"
#include "tree_manip.hpp"
#include "cond_likelihood.hpp"
#include "cond_likelihood_storage.hpp"
#include "tip_data.hpp"
#include "internal_data.hpp"
#include "xphylogeny.hpp"
#include "split.hpp"

using namespace boost::python;
using namespace phycas;

void translateXPhylogeny(const XPhylogeny &e)
	{
    // Use the Python 'C' API to set up an exception object
    PyErr_SetString(PyExc_Exception, e.what());
    }

BOOST_PYTHON_MODULE(_PhylogenyExt)
{
	// These function pointers are used to distinguish const from non-const member functions with the same name
	// see http://www.boost.org/libs/python/doc/tutorial/doc/html/python/functions.html#python.overloading
	const TipData *			(TreeNode::*GetTipDataConst)() const		= &TreeNode::GetTipData;
	const InternalData *	(TreeNode::*GetInternalDataConst)()	const	= &TreeNode::GetInternalData;

	class_<TreeNode, boost::shared_ptr<TreeNode> >("TreeNodeBase")
		.def("getSplit", &TreeNode::GetSplit, return_internal_reference<>())
		.def("getSupport", &TreeNode::GetSupport)
		.def("setSupport", &TreeNode::SetSupport)
		.def("getX", &TreeNode::GetX)
		.def("getY", &TreeNode::GetY)
		.def("setX", &TreeNode::SetX)
		.def("setY", &TreeNode::SetY)
		.def("isSelected", &TreeNode::IsSelected)
		.def("selectNode", &TreeNode::SelectNode)
		.def("unselectNode", &TreeNode::UnselectNode)
		.def("isRoot", &TreeNode::IsTipRoot)
		.def("isTip", &TreeNode::IsTip)
		.def("isInternal", &TreeNode::IsInternal)
		.def("getLeftChild", &TreeNode::GetLeftChild, return_internal_reference<>())
		.def("getRightSib", &TreeNode::GetRightSib, return_internal_reference<>())
		.def("getParent", &TreeNode::GetParent, return_internal_reference<>())
		.def("getNodeNumber", &TreeNode::GetNodeNumber)
		.def("getNodeName", &TreeNode::GetNodeName, return_value_policy<copy_const_reference>())
        .def("setNodeName", &TreeNode::SetNodeName)
		.def("getEdgeLen", &TreeNode::GetEdgeLen)
		.def("setEdgeLen", &TreeNode::SetEdgeLen)
		.def("getNextPreorder", &TreeNode::GetNextPreorder, return_internal_reference<>())
		.def("getNextPostorder", &TreeNode::GetNextPostorder, return_internal_reference<>())
		.def("getTipData", GetTipDataConst, return_internal_reference<>())
		.def("getInternalData", GetInternalDataConst, return_internal_reference<>())
		;

	class_<Tree, boost::shared_ptr<Tree> >("TreeBase")
		.def("getNNodes", &Tree::GetNNodes)
		.def("getNTips", &Tree::GetNTips)
		.def("getNObservables", &Tree::GetNObservables)
		.def("getNInternals", &Tree::GetNInternals)
		.def("getFirstPreorder", &Tree::GetFirstPreorder, return_internal_reference<>())
		.def("getFirstPostorder", &Tree::GetLastPreorder, return_internal_reference<>())
		.def("isRooted", &Tree::IsRooted)
		.def("isPolytomous", &Tree::IsPolytomous)
//		.def("setRooted", &Tree::setRooted)
//		.def("setUnrooted", &Tree::setUnrooted)
		.def("setRootedness", &Tree::setRootedness)
		.def("hasEdgeLens", &Tree::HasEdgeLens)
		.def("forget", &Tree::Forget)
		.def("clear", &Tree::Clear)
		.def("buildFromString", &Tree::BuildFromString)
		.def("edgeLenSum", &Tree::EdgeLenSum)
		.def("edgeLens", &Tree::EdgeLens)
		.def("keyToEdges", &Tree::KeyToEdges)
		.def("replaceEdgeLens", &Tree::replaceEdgeLens)
		.def("robinsonFoulds", &Tree::robinsonFoulds)
		.def("buildTreeID", &Tree::buildTreeID)
		.def("setAllEdgeLens", &Tree::SetAllEdgeLens)
		.def("debugWalkTree", &Tree::DebugWalkTree)
		.def("debugCheckTree", &Tree::DebugCheckTree)
		.def("debugListTree", &Tree::debugListTree)
		.def("rerootAtTip", &Tree::RerootAtTip)
		.def("rerootAtThisTip", &Tree::RerootAtThisTip)
		.def("rerootAtThisInternal", &Tree::RerootAtThisInternal)
		.def("stripNodeNames", &Tree::stripNodeNames)
		.def("makeNewick", &Tree::MakeNewick)
		.def("makeNewickForRefDist", &Tree::MakeNewickForRefDist)
		.def("makeNumberedNewick", &Tree::MakeNumberedNewick)
		.def("tipNumbersSetUsingNames", &Tree::TipNumbersSetUsingNames)
		.def("findTipByName", &Tree::FindTipByName)
		.def("rectifyNumbers", &Tree::RectifyNumbers)
		.def("rectifyNames", &Tree::RectifyNames)
		.def("refreshPreorder", &Tree::RefreshPreorder)
		.def("selectAllNodes", &Tree::SelectAllNodes)
		.def("unselectAllNodes", &Tree::UnselectAllNodes)
		.def("calcTotalHeight", &Tree::calcTotalHeight)
		.def("here", &Tree::DebugHere)
		.def("debugMode", &phycas::Tree::debugMode)
		.def("recalcAllSplits", &phycas::Tree::RecalcAllSplits)
		.def("ladderize", &phycas::Tree::Ladderize)
		.def("GetNInternalsAllocated", &phycas::Tree::GetNInternalsAllocated)
		.def("deroot", &phycas::Tree::deroot)
		;

	class_<TreeManip>("TreeManipBase", init<boost::shared_ptr<Tree> >())
		.def("buildTreeFromSplitVector", &TreeManip::buildTreeFromSplitVector)
		.def("starTree", &TreeManip::starTree)
		.def("randomTree", &TreeManip::randomTree)
		.def("equiprobTree", &TreeManip::equiprobTree)
		.def("buildEquiprobTree", &TreeManip::buildEquiprobTree)
		.def("setRandomEdgeLensFromTreeLengthDist", &TreeManip::setRandomEdgeLensFromTreeLengthDist)
		.def("setRandomEdgeLens", &TreeManip::setRandomEdgeLens)
		.def("setRandomInternalExternalEdgeLens", &TreeManip::setRandomInternalExternalEdgeLens)
		.def("deleteRandomInternalEdge", &TreeManip::deleteRandomInternalEdge)
		.def("addRightChild", &TreeManip::addRightChild)
		;

	class_<Split, boost::shared_ptr<Split> >("SplitBase")
        .def("getNTaxa", &Split::GetNTaxa)
        .def("setNTaxa", &Split::SetNTaxa)
        .def("copy", &Split::Copy)
        .def("reset", &Split::Reset)
		.def("setBit", &Split::SetBit)
		.def("setBits", &Split::SetBits)
		.def("unsetBit", &Split::UnsetBit)
		.def("unsetBits", &Split::UnsetBits)
		.def("isBitSet", &Split::IsBitSet)
		.def("invertSplit", &Split::InvertSplit)
		.def("calcComplexity", &Split::CalcComplexity)
		.def("countOnBits", &Split::CountOnBits)
		.def("countOffBits", &Split::CountOffBits)
		.def("equals", &Split::Equals)
		.def("cmp", &Split::Cmp)
		.def("isCompatible", &Split::IsCompatible)
		.def("subsumedIn", &Split::SubsumedIn)
		.def("getOnList", &Split::GetOnList)
		.def("getOffList", &Split::GetOffList)
        .def("getExcludedList", &Split::GetExcludedList)
		.def("setExcluded", &Split::SetExcluded)
		.def("getExcludedSymbol", &Split::GetExcludedSymbol)
		.def("getOnSymbol", &Split::GetOnSymbol)
		.def("getOffSymbol", &Split::GetOffSymbol)
		.def("setExcludedSymbol", &Split::SetExcludedSymbol)
		.def("setOnSymbol", &Split::SetOnSymbol)
		.def("setOffSymbol", &Split::SetOffSymbol)
		.def("createNewickRepresentation", &Split::CreateNewickRepresentation)
		.def("createPatternRepresentation", &Split::CreatePatternRepresentation)
		.def("createFromPattern", &Split::CreateFromPattern)
		.def("read", &Split::Read)
		.def("write", &Split::Write)
		;

    register_exception_translator<XPhylogeny>(&translateXPhylogeny);
}

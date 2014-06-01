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

#if ! defined(MCMC_PARAM_HPP)
#define MCMC_PARAM_HPP

#include <vector>									// for std::vector
#include <boost/shared_ptr.hpp>						// for boost::shared_ptr
#include <boost/weak_ptr.hpp>						// for boost::weak_ptr
#include <boost/enable_shared_from_this.hpp>		// for boost::enable_shared_from_this
#include "mcmc_updater.hpp"				// for MCMCUpdater base class
#include "split.hpp"

namespace phycas
{

class HKY;
class Codon;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the transition/transversion rate ratio parameter of models such as the HKY model. This parameter is
|	commonly symbolized using the Greek letter kappa, hence the name.
*/
class KappaParam : public MCMCUpdater
	{
	public:
						KappaParam();
						virtual ~KappaParam();

        virtual void	sendCurrValueToModel(double v);
        virtual double  getCurrValueFromModel() const;
		virtual void	setModel(ModelShPtr p);
		virtual bool	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		HKY *			hky;
		Codon *			codon;
	};


/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the nonsynonymous/synonymous rate ratio parameter of a codon model. This parameter is commonly
|	symbolized using the Greek letter omega, hence the name.
*/
class OmegaParam : public MCMCUpdater
	{
	public:
						OmegaParam();
						virtual ~OmegaParam()
							{
							//std::cerr << "OmegaParam dying..." << std::endl;
							}

        virtual void	sendCurrValueToModel(double v);
        virtual double 	getCurrValueFromModel() const;
		virtual void	setModel(ModelShPtr p);
		virtual bool	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		Codon *			codon;
	};

class GTR;

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates one of the relative rate parameters in the GTR model.
*/
class GTRRateParam : public MCMCUpdater
	{
	public:
						GTRRateParam(unsigned w);
						virtual ~GTRRateParam()
							{
							//std::cerr << "GTRRateParam dying..." << std::endl;
							}

        virtual void	sendCurrValueToModel(double v);
        virtual double  getCurrValueFromModel() const;
		virtual void	setModel(ModelShPtr p);
		virtual bool	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	private:

		GTR *			gtr;		/**< */
		unsigned		which;		/**< */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the discrete gamma shape parameter.
*/
class DiscreteGammaShapeParam : public MCMCUpdater
	{
	public:
						DiscreteGammaShapeParam();
						virtual ~DiscreteGammaShapeParam();

		void			finalizeWorkingPrior();
        virtual void	sendCurrValueToModel(double v);
        virtual double	getCurrValueFromModel() const;
		virtual bool	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates the proportion of invariable sites parameter.
*/
class PinvarParam : public MCMCUpdater
	{
	public:
						PinvarParam();
						virtual ~PinvarParam()
							{
							//std::cerr << "PinvarParam dying..." << std::endl;
							}

        virtual void	sendCurrValueToModel(double v);
        virtual double 	getCurrValueFromModel() const;
		virtual bool	update();				// override virtual from MCMCUpdater base class
		virtual double	operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates a state frequency parameter. More precisely, this class encapsulates a parameter that governs a state
|	frequency. Because all state frequencies must add to 1.0, the value of this parameter and the other parameter values
|	must be normalized so that the model has a correct set of state frequencies.
*/
class StateFreqParam : public MCMCUpdater
{
	public:
							StateFreqParam(unsigned w);
							virtual ~StateFreqParam();

        virtual void		sendCurrValueToModel(double v);
        virtual double      getCurrValueFromModel() const;
		virtual bool		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);

	private:

		unsigned			which;
};

/*----------------------------------------------------------------------------------------------------------------------
|	Structure combining a vector of doubles (`fs') for storing a fitting sample and a probability distribution shared
|	pointer to store the working prior (`wp') that is parameterized using the fitting sample. This structure is used
|	as the value for a map (data member `edge_ref_dist') in EdgeLenMasterParam objects.
*/
struct EdgeWorkingPrior
	{
	double_vect_t	fs;		/**< vector of doubles representing samples upon which the `edge_ref_dist' will be based */
	ProbDistShPtr	wp;		/**< Gamma working prior distribution */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Encapsulates an edge length parameter. It is only used in the case of a fixed topology, otherwise edge lengths are
|	updated using a `LargetSimonMove' or a `BushMove' updater.
*/
class EdgeLenParam : public MCMCUpdater
	{
	public:
                                EdgeLenParam();
                                EdgeLenParam(unsigned i);
        virtual                 ~EdgeLenParam();

		bool                    isInternalEdge();
		bool                    isExternalEdge();

        std::string             getSplitReprAsString() const;

        virtual bool            computesTreeLengthPrior() const;
        virtual void            sendCurrValueToModel(double v);
        virtual double          getCurrValueFromModel() const;
		void                    setTreeNode(TreeNode & nd);
		virtual bool            update();				// override virtual from MCMCUpdater base class
		virtual double          operator()(double k);	// override pure virtual from AdHocDensity (base class of MCMCUpdater)

	//temp private:
		TreeNode *		_my_node;
        unsigned        _edge_len_index;                /**< is the index of this slave used when sending edge length values to the master */
	};

/*----------------------------------------------------------------------------------------------------------------------
|	Represents all edge lengths in the associated Tree. An edge length master parameter does not update any particular
|	edge length (that is left up to one of the defined moves, such as a `LargetSimonMove'); instead, its	job is to
|	compute the joint prior over all edge lengths. The update() member function and operator() of an EdgeLenMasterParam
|	are not overridden because they do nothing.
*/
class EdgeLenMasterParam : public MCMCUpdater
	{
	public:
        enum                EdgeLenType
                                {
                                internal,   /**> only computes prior for internal edges */
                                external,   /**> only computes prior for external edges */
                                both        /**> computes prior for all edges */
                                };

                            EdgeLenMasterParam(EdgeLenMasterParam::EdgeLenType t = EdgeLenMasterParam::both);
							virtual ~EdgeLenMasterParam();

        virtual bool        computesTreeLengthPrior() const;

		double				recalcWorkingPrior() const;
		double				lnWorkingPriorOneEdge(const TreeNode & nd, double v) const;
		std::string 		getWorkingPriorDescr() const;

		void				useEdgeSpecificWorkingPriors(bool use_it);

	protected:

		double				edgeLength(TreeNode & nd) const;
		double				lnPriorOneEdge(TreeNode & nd) const;

    private:

        EdgeLenType         				edgeLenType;   						/**> holds the edge length type, which determines for which edge lengths the prior is computed when recalcPrior is called */

		bool								use_edge_specific_ref_dists;	/**< if true, `edge_ref_dist' will be used; otherwise, a single generic working prior will be used for all edge lengths */
		std::map<Split,EdgeWorkingPrior>	edge_ref_dist;					/**< maps splits (keys) to EdgeWorkingPrior structs (values) so that the working prior distribution can be fetched given the split corresponding to any given node in the tree */
	};

typedef std::map<Split,EdgeWorkingPrior>::iterator			WorkingPriorMapIter;
typedef std::map<Split,EdgeWorkingPrior>::const_iterator	WorkingPriorMapConstIter;
typedef std::pair<Split,EdgeWorkingPrior> 					WorkingPriorMapPair;
typedef boost::shared_ptr<EdgeLenMasterParam> EdgeLenMasterParamShPtr;

/*----------------------------------------------------------------------------------------------------------------------
|	Represents the mean of the edge length prior. This parameter is used in models employing a hyperprior governing the
|	mean of the edge length priors. Whenever the HyperPriorParam value changes, the means of the priors affecting all
|	EdgeLenParam objects in the MCMCChainManager are reset to this new value.
*/
class HyperPriorParam : public MCMCUpdater
	{
	public:
							HyperPriorParam();
							enum EdgeTypeEnum {all, internal, external};
							HyperPriorParam(EdgeTypeEnum t);
							virtual ~HyperPriorParam();

		void				setEdgeLenMasterParam(EdgeLenMasterParamShPtr p) {edgelen_master_param = p;}
        virtual void		sendCurrValueToModel(double v);
        virtual double      getCurrValueFromModel() const;
		virtual bool		update();				// override virtual from MCMCUpdater base class
		virtual double		operator()(double k);

    private:

		EdgeTypeEnum				edge_type;				/**< governs the type of edges (internal, external, all) affected by this hyperparam */
        EdgeLenMasterParamShPtr		edgelen_master_param;   /**> is the edge length master parameter whose prior this parameter controls */
	};

} // namespace phycas

//#include "phycas/src/mcmc_param.inl"

#endif

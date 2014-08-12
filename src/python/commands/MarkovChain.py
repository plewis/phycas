import os,sys,math
from phycas import *
import phycas.phylogeny as phylogeny
import phycas.probdist as probdist
import phycas.likelihood as likelihood
from LikelihoodCore import LikelihoodCore

class MarkovChain(LikelihoodCore):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The MarkovChain class encapsulates the notion of a Markov chain used
    in Bayesian analyses. The MarkovChain class has the ability to
    orchestrate a Markov chain Monte Carlo (MCMC) analysis. In
    Metropolis-coupled MCMC, the cold chain and each of the heated chains
    are MarkovChain objects. The primary addition MarkovChain adds to the
    base class LikelihoodCore are prior distributions and the
    self.heating_power data member.

    """
    def __init__(self, parent, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The MarkovChain constructor makes a copy of the supplied parent
        object, clones all of the prior ProbabilityDistribution objects in the
        supplied parent object, and sets self.heating_power to the supplied
        power.

        """
        LikelihoodCore.__init__(self, parent)

        self.parent                     = parent    # Note: self.parent is the MCMCImpl object
        #POLTMP2 self.boldness                   = 0.0
        self.heating_power              = power
        self.chain_manager              = None
        self.tree_scaler_move           = None
        self.subset_relrates_move       = None
        self.edge_move                  = None
        self.larget_simon_move          = None
        self.bush_move                  = None
        self.topo_prior_calculator      = None

        self.state_freq_moves           = []
        self.rel_rate_moves             = []
        self.all_updaters_list = None

        self.setupChain()

    def debugUpdaters(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Feel free to modify this function to check some property of each
        updater. Should only be called when actively debugging.

        """
        print '~~~~~ entering MarkovChain.debugUpdaters()...'
        for updater in self.chain_manager.getAllUpdaters():
            n = updater.getName()
            t = updater.getTree()
            print n, '-->', t.makeNewick(5)
        print '~~~~~ leaving MarkovChain.debugUpdaters()...'

    def resetNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the resetNumLikelihoodEvals function of self.likelihood. This
        resets the number of likelihood evaluations performed to 0.

        """
        # Note: likelihood data member inherited from LikelihoodCore
        return self.likelihood.resetNumLikelihoodEvals()

    def getNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls the getNumLikelihoodEvals function of self.likelihood. This
        returns the number of likelihood evaluations performed since the last
        call of resetNumLikelihoodEvals.

        """
        # Note: likelihood data member inherited from LikelihoodCore
        return self.likelihood.getNumLikelihoodEvals()

    def paramFileHeader(self, paramf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the parameter file. The
        parameter paramf should be an (already open for writing) file object.
        The parameter file corresponds to the *.p file produced by MrBayes. It
        begins with a line containing the initial random number seed, followed
        by a line of column titles for the (tab-delimited) sampled parameter
        information on the subsequent lines.

        """
        # Note: self.parent equals MCMCImpl in this case
        # Note: r data member inherited from LikelihoodCore
        paramf.write('[ID: %d]\n' % self.r.getInitSeed())
        if self.parent.opts.doing_steppingstone_sampling:
            paramf.write('Gen\tbeta\tlnL\tlnPrior')
        else:
            paramf.write('Gen\tlnL\tlnPrior')

        if self.parent.opts.doing_steppingstone_sampling and not self.parent.opts.ss_heating_likelihood:
            paramf.write('\tlnRefDens')

        # If the user has defined a reference tree, add a column for the Robinson-Foulds
        # distance between the sampled tree and the reference tree
        if self.parent.ref_tree is not None:
            paramf.write('\tdRF')

        # If using a model that allows polytomies, include a column indicating the
        # resolution class of the tree
        if self.parent.opts.allow_polytomies:
            paramf.write('\tResClass')

        paramf.write('\tTL')
        if self.parent.opts.fix_topology:
            nbrlens = self.tree.getNNodes() - 1
            for i in range(nbrlens):
                paramf.write('\tbrlen%d' % (i+1))
            self.parent.output('\nKey to the edges (preorder traversal):\n%s' % self.tree.keyToEdges())

        param_names = self.partition_model.getAllParameterNames()
        for nm in param_names:
            paramf.write('\t%s' % nm)
        #nmodels = self.partition_model.getNumSubsets()
        #if nmodels > 1:
        #    for i in range(nmodels):
        #        paramf.write('\tm_%d' % (i+1,))
        #for i in range(nmodels):
        #    m = self.partition_model.getModel(i)
        #    paramf.write(m.paramHeader())

    def treeFileHeader(self, treef):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Writes the header line at the beginning of the tree file. The
        parameter treef should be an (already open for writing) file object.
        The tree file corresponds to the *.t file produced by MrBayes, and is
        in NEXUS format. This function writes the opening "#NEXUS" line, a
        comment containing the initial random number seed, and the beginning
        of a NEXUS TREES block, including the TRANSLATE command.

        """
        treef.write('#NEXUS\n')
        treef.write('[ID: %d]\n' % self.r.getInitSeed())
        treef.write('begin trees;\n')
        if self.parent.ntax > 0:
            treef.write('\ttranslate\n')
            for i in range(self.parent.ntax):
                if self.parent.taxon_labels[i].find(' ') < 0:
                    # no spaces found in name
                    treef.write('\t\t%d %s%s\n' % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))
                else:
                    # at least one space in taxon name, so enclose name in quotes
                    treef.write("\t\t%d '%s'%s\n" % (i + 1, self.parent.taxon_labels[i], i == self.parent.ntax - 1 and ';' or ','))

    def setPower(self, power):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the heating_power data member and calls the setPower method for
        every updater so that all updaters are immediately informed that the
        power for this chain has changed.

        """
        self.heating_power = power
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(power)

    #POLTMP2 def setBoldness(self, boldness):
    #POLTMP2     #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #POLTMP2     """
    #POLTMP2     Sets the boldness data member and calls the setBoldness method for
    #POLTMP2     every updater so that all updaters are immediately informed that the
    #POLTMP2     boldness for this chain has changed. The boldness is a value from 0 to
    #POLTMP2     100 that specifies the boldness of Metropolis-Hastings moves. Setting
    #POLTMP2     the boldness has no effect on slice sampling based updaters. Each
    #POLTMP2     move class defines what is meant by boldness. The boldness is changed
    #POLTMP2     during an MCMC run in some circumstances, such as during a path
    #POLTMP2     sampling analysis where the target distribution changes during the
    #POLTMP2     run.
    #POLTMP2
    #POLTMP2     """
    #POLTMP2     self.boldness = boldness
    #POLTMP2     for updater in self.chain_manager.getAllUpdaters():
    #POLTMP2         updater.setBoldness(boldness)

    def setupChain(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The setupChain method prepares a MarkovChain object before an MCMC
        analysis is started. This method is called at the end of the
        MarkovChain constructor. The duties performed include:
        1) calling the base class setupCore method;
        2) preparing the starting tree by adding data structures needed for
        computing likelihood (transition matrices and conditional likelihood
        arrays);
        3) setting up prior distributions for model parameters;
        4) creating an MCMCManager object and adding all relevant updaters to
        it;
        and 5) makes sure each updater knows the type of heating and the
        power.

        """
        LikelihoodCore.setupCore(self)  # creates joint prior manager
        self.partition_model.createJointPriorManager()
        self.joint_prior_manager = self.partition_model.getJointPriorManager()

        from phycas import partition,model
        if self.parent.opts.partition.noData():
            self.likelihood.setNoData()

        LikelihoodCore.prepareForLikelihood(self)

        self.python_only_moves = []

        # add priors to models already added (by LikelihoodCore.setupCore) to partition_model
        modelspecs = partition.getModels()  # partition here refers to the global object (associated with the Phycas partition command)
        #print 'modelspecs has length %d:' % len(modelspecs)
        nmodels = self.partition_model.getNumSubsets()
        #if nmodels == 1:
        #    print 'partition_model contains 1 model (i.e. unpartitioned)'
        #else:
        #    print 'partition_model contains %d models' % nmodels
        self.chain_manager = likelihood.MCMCChainManager(self.joint_prior_manager)
        for i in range(nmodels):
            # get the Model (as defined in likelihood_models.cpp) associated with partition subset i
            m = self.partition_model.getModel(i)

            # get the model specification (as defined in Model.py) stored in (Python) partition object
            mspec = modelspecs[i]

            #implemented = not (self.parent.opts.fix_topology and mspec.edgelen_hyperprior is not None)
            #self.parent.phycassert(implemented, 'Cannot currently specify an edge length hyperprior and fix the topology at the same time')

            # Copy priors related to edge lengths into model for first subset
            # Note: The first subset is the only one that needs edge length priors because those are the
            # only ones actually used (at this writing, 31 Jan 2010) because both
            # tree topology and edge lengths are always linked across subsets
            if i == 0:
                #separate_edge_len_dists = mspec.separate_edgelen_hyper
                separate_edge_len_dists = (mspec.edgelen_prior is None)
                m.separateInternalExternalEdgeLenPriors(separate_edge_len_dists)
                m.setEdgeLenHyperPrior(None)
                if mspec.tree_length_prior is not None:
                    m.setTreeLengthPrior(mspec.tree_length_prior.cloneAndSetLot(self.r))    # TREE_LENGTH_PRIOR
                    self.likelihood.setTreeLengthPrior(mspec.tree_length_prior.cloneAndSetLot(self.r))
                else:
                    if mspec.external_edgelen_prior is not None:
                        #print '~~~~~> mspec.external_edgelen_prior is not None'
                        m.setExternalEdgeLenPrior(mspec.external_edgelen_prior.cloneAndSetLot(self.r))
                    if mspec.internal_edgelen_prior is not None:
                        #print '~~~~~> mspec.internal_edgelen_prior is not None'
                        m.setInternalEdgeLenPrior(mspec.internal_edgelen_prior.cloneAndSetLot(self.r))
                    if mspec.edgelen_hyperprior is not None:
                        #print '~~~~~> mspec.edgelen_hyperprior is not None'
                        m.setEdgeLenHyperPrior(mspec.edgelen_hyperprior.cloneAndSetLot(self.r))
                        if mspec.fix_edgelen_hyperparam:
                            m.fixEdgeLenHyperprior()    #@POL should be named fixEdgeLenHyperparam

            # Copy priors for this model's parameters
            if mspec.type == 'codon':
                m.setKappaPrior(mspec.kappa_prior.cloneAndSetLot(self.r))
                m.setOmegaPrior(mspec.omega_prior.cloneAndSetLot(self.r))
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
            elif mspec.type == 'gtr':
                if mspec.update_relrates_separately:
                    m.setRelRateParamPrior(mspec.relrate_param_prior.cloneAndSetLot(self.r))
                else:
                    self.parent.phycassert(mspec.relrate_prior.getDistName() == 'Dirichlet', 'mspec.relrate_prior must be of type Dirichlet')
                    m.setRelRatePrior(mspec.relrate_prior)
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    self.parent.phycassert(mspec.state_freq_prior.getDistName() == 'Dirichlet', 'mspec.state_freq_prior must be of type Dirichlet')
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
            elif mspec.type == 'hky':
                m.setKappaPrior(mspec.kappa_prior.cloneAndSetLot(self.r))
                if mspec.update_freqs_separately:
                    m.setStateFreqParamPrior(mspec.state_freq_param_prior.cloneAndSetLot(self.r))
                else:
                    m.setStateFreqPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))

            # Copy priors related to among-site rate heterogeneity
            if mspec.num_rates > 1:
                m.setDiscreteGammaShapePrior(mspec.gamma_shape_prior.cloneAndSetLot(self.r))
            if mspec.pinvar_model:
                m.setPinvarPrior(mspec.pinvar_prior.cloneAndSetLot(self.r))

            # If user specifies fixed tree topology, make each edge length a separate parameter; otherwise, use edge length master parameters
            if self.parent.opts.fix_topology:
                m.setEdgeSpecificParams(True)
            else:
                m.setEdgeSpecificParams(False)

            # Add all necessary updaters to the MCMCManager
            if nmodels == 1:
                subset_pos = -1
            else:
                subset_pos = i

            self.chain_manager.addMCMCUpdaters(
                m,                                  # substitution model
                self.tree,                          # tree
                self.likelihood,                    # likelihood calculation machinery
                self.r,                             # pseudorandom number generator
                self.parent.opts.slice_max_units,   # maximum number of slice units allowed
                self.parent.opts.slice_weight,      # weight for each parameter added
                subset_pos)                         # i is the subset (needed so that edge length params will only be added for for first subset)

            # Add subset-model-specific moves
            if not mspec.type == 'jc' and not mspec.update_freqs_separately:
                # Create a StateFreqMove to update entire state frequency vector
                sfm = likelihood.StateFreqMove()
                if mspec.type == 'codon':
                    sfm.setDimension(61)
                else:
                    sfm.setDimension(4)
                #if nmodels > 1:
                #    sfm.setName("state_freqs_%d" % (i+1,))
                #else:
                #    sfm.setName("state_freqs")
                sfm.setName("%d_state_freqs" % (i+1,))
                sfm.setWeight(self.parent.opts.state_freq_weight)
                #POLTMP2  sfm.setPosteriorTuningParam(self.parent.opts.state_freq_psi)
                sfm.setTuningParameter(self.parent.opts.state_freq_psi)
                #POLTMP2  sfm.setPriorTuningParam(self.parent.opts.state_freq_psi0)
                sfm.setTree(self.tree)
                sfm.setModel(m)
                sfm.setTreeLikelihood(self.likelihood)
                sfm.setLot(self.r)
                if m.stateFreqsFixed():
                    sfm.fixParameter()
                sfm.setMultivarPrior(mspec.state_freq_prior.cloneAndSetLot(self.r))
                self.chain_manager.addMove(sfm)
                self.state_freq_moves.append(sfm)

            if mspec.type == 'gtr' and not mspec.update_relrates_separately:
                # Create a RelRateMove to update entire relative rates vector
                rrm = likelihood.RelRatesMove()
                #if nmodels > 1:
                #    rrm.setName("relrates_%d" % (i+1,))
                #else:
                #    rrm.setName("relrates")
                rrm.setName("%d_relrates" % (i+1,))
                rrm.setWeight(self.parent.opts.rel_rate_weight)
                #POLTMP2  rrm.setPosteriorTuningParam(self.parent.opts.rel_rate_psi)
                rrm.setTuningParameter(self.parent.opts.rel_rate_psi)
                #POLTMP2  rrm.setPriorTuningParam(self.parent.opts.rel_rate_psi0)
                rrm.setTree(self.tree)
                rrm.setModel(m)
                rrm.setTreeLikelihood(self.likelihood)
                rrm.setLot(self.r)
                #if self.model.relRatesFixed():
                #    rrm.fixParameter()
                rrm.setMultivarPrior(mspec.relrate_prior.cloneAndSetLot(self.r))
                self.chain_manager.addMove(rrm)
                self.rel_rate_moves.append(rrm)

        self.likelihood.replaceModel(self.partition_model)

        if self.parent.opts.data_source is None:
            self.likelihood.setNoData() # user apparently wants to run MCMC with no data

        model0 = self.partition_model.getModel(0)

        # Create a TreeScalerMove object to handle scaling the entire tree to allow faster
        # convergence in edge lengths. This move is unusual in using slice sampling rather
        # than Metropolis-Hastings updates: most "moves" in parent are Metropolis-Hastings.
        if self.parent.opts.tree_scaler_weight > 0:
            self.tree_scaler_move = likelihood.TreeScalerMove()
            self.tree_scaler_move.setName("tree_scaler")
            self.tree_scaler_move.setWeight(self.parent.opts.tree_scaler_weight)
            #POLTMP2  self.tree_scaler_move.setPosteriorTuningParam(self.parent.opts.tree_scaler_lambda)
            self.tree_scaler_move.setTuningParameter(self.parent.opts.tree_scaler_lambda)
            #POLTMP2  self.tree_scaler_move.setPriorTuningParam(self.parent.opts.tree_scaler_lambda0)
            self.tree_scaler_move.setTree(self.tree)
            self.tree_scaler_move.setModel(model0)
            self.tree_scaler_move.setTreeLikelihood(self.likelihood)
            self.tree_scaler_move.setLot(self.r)
            if model0.edgeLengthsFixed():
                self.tree_scaler_move.fixParameter()
            self.chain_manager.addMove(self.tree_scaler_move)

        # If more than one partition subset, add a SubsetRelRate move to modify the
        # vector of relative substitution rates for each subset
        if (nmodels > 1):
            self.subset_relrates_move = likelihood.SubsetRelRatesMove()
            self.subset_relrates_move.setDimension(nmodels)
            self.subset_relrates_move.setName("subset_relrates")
            if self.parent.__class__.__name__ == 'InflatedDensityRatio':
                # IDR method does not actually run a chain, just uses it to compute likelihoods and priors
                self.subset_relrates_move.setWeight(0)
                #POLTMP2  self.subset_relrates_move.setPosteriorTuningParam(1.0)
                self.subset_relrates_move.setTuningParameter(1.0)
                #POLTMP2  self.subset_relrates_move.setPriorTuningParam(1.0)
            else:
                self.subset_relrates_move.setWeight(self.parent.opts.subset_relrates_weight)
                #POLTMP2  sPOLTMP elf.subset_relrates_move.setPosteriorTuningParam(self.parent.opts.subset_relrates_psi)
                elf.subset_relrates_move.setTuningParameter(self.parent.opts.subset_relrates_psi)
                #POLTMP2  self.subset_relrates_move.setPriorTuningParam(self.parent.opts.subset_relrates_psi0)
            self.subset_relrates_move.setTree(self.tree)
            self.subset_relrates_move.setModel(None)    # the model data member is ignored in this case; instead, the partition model stores the parameters
            self.subset_relrates_move.setPartitionModel(self.partition_model)
            self.subset_relrates_move.setTreeLikelihood(self.likelihood)
            self.subset_relrates_move.setLot(self.r)
            #subset_proportions = partition.getSubsetProportions()
            #self.subset_relrates_move.setSubsetProportions(subset_proportions)
            if partition.fix_subset_relrates:
                self.subset_relrates_move.fixParameter()
            else:
                # only assign a prior distribution if subset relative rates are not fixed
                self.subset_relrates_move.freeParameter()
                param_list = tuple([1.0]*nmodels)
                if partition.subset_relrates_prior is None:
                    # User did NOT provide a subset relative rates prior
                    c = self.partition_model.getSubsetProportions()
                    d = probdist.RelativeRateDistribution(param_list, c.getSubsetProportions())
                    #d.setSubsetProportions(self.partition_model.getSubsetProportions())
                    self.subset_relrates_move.setMultivarPrior(d.cloneAndSetLot(self.r))
                    self.partition_model.setSubsetRelRatePrior(d.cloneAndSetLot(self.r))

                    jpm = self.chain_manager.getJointPriorManager()
                    jpm.addMultivariateDistribution("subset_relrates", d.cloneAndSetLot(self.r), param_list)
                else:
                    # User DID provide a subset relative rates prior
                    self.parent.phycassert(partition.subset_relrates_prior.getDistName() == 'RelativeRateDistribution', 'partition.subset_relrates_prior must be of type RelativeRateDistribution')
                    self.parent.phycassert(partition.subset_relrates_prior.getNParams() == nmodels, 'partition.subset_relrates_prior has dimension %d, but there are %d subsets in the partition. Try setting partion.subset_relrates_prior = None to get default flat Dirichlet prior of the appropriate dimension' % (partition.subset_relrates_prior.getNParams(), nmodels))
                    # partition.subset_relrates_prior.setCoefficients(subset_proportions)
                    partition.subset_relrates_prior.setSubsetProportions(self.partition_model.getSubsetProportions())
                    self.subset_relrates_move.setMultivarPrior(partition.subset_relrates_prior.cloneAndSetLot(self.r))
                    self.partition_model.setSubsetRelRatePrior(partition.subset_relrates_prior.cloneAndSetLot(self.r))

                    jpm = self.chain_manager.getJointPriorManager()
                    jpm.addMultivariateDistribution("subset_relrates", partition.subset_relrates_prior.cloneAndSetLot(self.r), param_list)
            self.chain_manager.addMove(self.subset_relrates_move)

        if not self.parent.opts.fix_topology:
            # Create a LargetSimonMove object to handle Metropolis-Hastings
            # updates to the tree topology and edge lengths
            self.larget_simon_move = likelihood.LargetSimonMove()
            self.larget_simon_move.setName("larget_simon_local")
            self.larget_simon_move.setWeight(self.parent.opts.ls_move_weight)
            #POLTMP2  self.larget_simon_move.setPosteriorTuningParam(self.parent.opts.ls_move_lambda)
            self.larget_simon_move.setTuningParameter(self.parent.opts.ls_move_lambda)
            #POLTMP2  self.larget_simon_move.setPriorTuningParam(self.parent.opts.ls_move_lambda0)
            self.larget_simon_move.setTree(self.tree)
            self.larget_simon_move.setModel(model0)
            self.larget_simon_move.setTreeLikelihood(self.likelihood)
            self.larget_simon_move.setLot(self.r)
            #POLTMP2  self.larget_simon_move.setLambda(self.parent.opts.ls_move_lambda)
            if model0.edgeLengthsFixed():
                self.larget_simon_move.fixParameter()
            self.chain_manager.addMove(self.larget_simon_move)


        # If requested, create a BushMove object to allow polytomous trees
        if self.parent.opts.allow_polytomies:
            # Create a BushMove object
            self.bush_move = likelihood.BushMove()

            # Set up BushMove object
            self.bush_move.setName("bush_move")
            self.bush_move.setWeight(self.parent.opts.bush_move_weight)
            self.bush_move.setTree(self.tree)
            self.bush_move.setModel(model0)
            self.bush_move.setTreeLikelihood(self.likelihood)
            self.bush_move.setLot(self.r)
            self.bush_move.setEdgeLenDistMean(self.parent.opts.bush_move_edgelen_mean)
            #self.bush_move.viewProposedMove(self.parent.bush_move_debug)
            if model0.edgeLengthsFixed():
                self.bush_move.fixParameter()
            self.bush_move.finalize()

            self.chain_manager.addMove(self.bush_move)

        if not self.parent.opts.fix_topology:
            if self.parent.opts.allow_polytomies:
                topo_prior_calculator = likelihood.PolytomyTopoPriorCalculatorBase()
                topo_prior_calculator.chooseUnrooted()
                topo_prior_calculator.setC(self.parent.opts.topo_prior_C)
                if self.parent.opts.polytomy_prior:
                    topo_prior_calculator.choosePolytomyPrior()
                else:
                    topo_prior_calculator.chooseResolutionClassPrior()
                self.chain_manager.getJointPriorManager().addTopologyDistribution("tree_topology", topo_prior_calculator, self.tree);
            else:
                self.chain_manager.getJointPriorManager().addTopologyDistribution("tree_topology", likelihood.TopoProbCalculatorBase(), self.tree);

        self.chain_manager.finalize()

        # Calculate relative rates if among-site rate heterogeneity is part of the model
        # Currently, the unnormalized rates are stored in the model; the function call
        # below normalizes the rate means and probabilities so that they will be accurately
        # recorded in the first line of the parameter file
        self.likelihood.recalcRelativeRates()

        # Make sure each updater knows the heating power and heating type
        # and set boldness to 0 for each updater (0% boldness is the default,
        # with more boldness used only in steppingstone sampling where bolder
        # moves are needed as the influence of the likelihood is progressively
        # diminished)
        for updater in self.chain_manager.getAllUpdaters():
            updater.setPower(self.heating_power)
            #POLTMP2 updater.setBoldness(0.0)
            if self.parent.opts.ss_heating_likelihood:
                updater.setLikelihoodHeating()
            else:
                updater.setStandardHeating()

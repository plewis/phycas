import os,sys,math
from phycas import *
import phycas.phylogeny as phylogeny
import phycas.probdist as probdist
import phycas.likelihood as likelihood

def cloneDistribution(d):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This function clones (deep copies) a ProbabilityDistribution object.
    Cloning is needed because while separate chains in an MCMC analysis
    need to begin with the same random number seed, it is not good for
    them to share a pseudorandom number generator if (for example)
    separate chains are run in different threads or on different
    processors.

    """
    if d == None:
        return None
    else:
        return d.clone()

class LikelihoodCore(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    The LikelihoodCore class comprises only those features of a
    MarkovChain that are needed to compute likelihoods. Thus, one could
    construct a stand-alone LikelihoodCore object if all you wanted to do
    was calculate the likelihood of a tree under a particular model. The
    LikelihoodCore class serves as the base class for MarkovChain.

    """
    def __init__(self, parent):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The constructor for the LikelihoodCore class stores the supplied
        parent object in the data member self.parent. It also creates a Tree
        object and stores it in self.tree and a pseudorandom number generator
        (Lot object), storing it in self.r. Finally, it clones the starting
        edge length distribution, storing it in self.starting_edgelen_dist.

        """
        self.parent                 = parent    # e.g. LikeImpl, MCMCImpl, etc.
        self.likelihood             = None
        self._tree                  = None
        self.r                      = self.parent._getLot()
        self.partition_model        = None
        self.joint_prior_manager    = None

    def getTree(self):
        if self._tree is None:
            self.setupCore()
        return self._tree

    def setTree(self, t):
        #self.likelihood.storeAllCLAs(self._tree)    # avoids a lot of deallocation and reallocation
        if t and not isinstance(t, phylogeny.Tree):
            tr = phylogeny.Tree()
            self._tree = t.buildTree(tr)
        else:
            self._tree = t

    def delTree(self, t):
        if self._tree:
            del self._tree
            self._tree = None

    tree = property(getTree, setTree, delTree)

    def createModel(self, model_spec):
        """
        This method creates an actual model (e.g. JCModel defined in likelihood/_Model.py) from a model_spec
        (e.g. Model defined in Phycas/Model.py). The model_spec is created by the user, who gives it a
        type (jc, hky, gtr, codon, etc.), but the model_spec cannot be used to compute likelihoods.

        """
        if model_spec.type == 'codon':
            self.parent.phycassert(model_spec.num_rates == 1, 'Cannot currently use gamma rate heterogeneity within codon model')
            self.parent.phycassert(not model_spec.pinvar_model, 'Cannot currently use invariable sites model within codon model')
            m = likelihood.CodonModel()
            m.setKappa(model_spec.kappa)
            if model_spec.fix_kappa:
                m.fixKappa()
            m.setOmega(model_spec.omega)
            if model_spec.fix_omega:
                m.fixOmega()
            if not model_spec.update_freqs_separately:
                ndimensions = model_spec.state_freq_prior.getNParams()
                self.parent.phycassert(ndimensions == 61, 'state_freq_prior should be 61-dimensional, but instead has %d dimensions' % ndimensions)
            self.parent.phycassert(model_spec.state_freqs, 'state_freqs is None, but should be a list containing 61 (unnormalized) relative codon frequencies')
            self.parent.phycassert(len(model_spec.state_freqs) == 61, 'state_freqs should be a list containing exactly 61 codon frequencies; instead, it contains %d values' % len(model_spec.state_freqs))
            for c in range(61):
                self.parent.phycassert(model_spec.state_freqs[c] >= 0.0, 'state_freqs[%d] cannot be negative (%f was specified)' % (c,model_spec.state_freqs[c]))
            m.setStateFreqsUnnorm(model_spec.state_freqs)
            if model_spec.fix_freqs:
                m.fixStateFreqs()
        elif model_spec.type in ['gtr','hky']:
            if model_spec.type == 'gtr':
                m = likelihood.GTRModel()
                m.setRelRates(model_spec.relrates)
                if model_spec.fix_relrates:
                    m.fixRelRates()
            else:
                m = likelihood.HKYModel()
                m.setKappa(model_spec.kappa)
                if model_spec.fix_kappa:
                    m.fixKappa()
            self.parent.phycassert(model_spec.state_freqs, 'state_freqs is None, but should be a list containing 4 (unnormalized) relative base frequencies')
            self.parent.phycassert(len(model_spec.state_freqs) == 4, 'state_freqs should be a list containing exactly 4 base frequencies; instead, it contains %d values' % len(model_spec.state_freqs))
            self.parent.phycassert(model_spec.state_freqs[0] >= 0.0, 'state_freqs[0] cannot be negative (%f was specified)' % model_spec.state_freqs[0])
            self.parent.phycassert(model_spec.state_freqs[1] >= 0.0, 'state_freqs[1] cannot be negative (%f was specified)' % model_spec.state_freqs[1])
            self.parent.phycassert(model_spec.state_freqs[2] >= 0.0, 'state_freqs[2] cannot be negative (%f was specified)' % model_spec.state_freqs[2])
            self.parent.phycassert(model_spec.state_freqs[3] >= 0.0, 'state_freqs[3] cannot be negative (%f was specified)' % model_spec.state_freqs[3])
            m.setNucleotideFreqs(model_spec.state_freqs[0], model_spec.state_freqs[1], model_spec.state_freqs[2], model_spec.state_freqs[3])  #POL should be named setStateFreqs?
            if model_spec.fix_freqs:
                m.fixStateFreqs()
        elif model_spec.type in ['loss','gain']:
            m = likelihood.IrreversibleModel()
            m.setScalingFactor(model_spec.scaling_factor)
            if model_spec.fix_scaling_factor:
                m.fixScalingFactor()
            else:
                m.freeScalingFactor()
            if model_spec.type == 'loss':
                m.setLossOnly()
            else:
                m.setGainOnly()
        elif model_spec.type in ['binary']:
            m = likelihood.BinaryModel()
            m.setScalingFactor(model_spec.scaling_factor)
            if model_spec.fix_scaling_factor:
                m.fixScalingFactor()
            else:
                m.freeScalingFactor()
            m.setKappa(model_spec.kappa)
            if model_spec.fix_kappa:
                m.fixKappa()
            else:
                m.freeKappa()
        else:
            m = likelihood.JCModel()

        # If rate heterogeneity is to be assumed, add it to the model here
        # Note must defer setting up pattern specific rates model until we know number of patterns
        if model_spec.num_rates > 1:
            m.setNGammaRates(model_spec.num_rates)
            m.setShape(model_spec.gamma_shape)
            if model_spec.fix_shape:
                m.fixShape()
        else:
            m.setNGammaRates(1)

        if model_spec.pinvar_model:
            m.setPinvarModel()
            m.setPinvar(model_spec.pinvar)
            if model_spec.fix_pinvar:
                m.fixPinvar()
        else:
            m.setNotPinvarModel()

        if model_spec.fix_edgelens:
            m.fixEdgeLengths()

        return m

    def setupCore(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The setupCore function does the following based on information stored
        in self.parent: 1) sets the random number seed of the local pseudo-
        random number generator (self.r); 2) creates a substitution model,
        storing it in self.model; 3) creates a TreeLikelihood object, storing
        it in self.likelihood; and 4) builds the starting tree, storing it in
        self.tree.

        """
        # Set seed if user has supplied one
        self.r = self.parent._getLot()
        #self.starting_edgelen_dist.setLot(self.r)

        from phycas import partition,model
        # Create the object
        self.partition_model = likelihood.PartitionModelBase()

        modelspecs = partition.getModels()
        if len(modelspecs) < 1:
            # model is an interface model object (i.e. Model defined in Phycas/Model.py),
            # not an actual model object (e.g. JCModel defined in likelihood/_Model.py).
            # Need to create a real model using the interface specification
            # and add the real model to partition_model.
            mod = self.createModel(model)
            self.partition_model.addModel(mod)
        else:
            for m in modelspecs:
                # m is an interface model object (i.e. Model defined in Phycas/Model.py),
                # not an actual model object (e.g. JCModel defined in likelihood/_Model.py).
                # Need to create a real model using the interface specification
                # and add the real model to partition_model.
                mod = self.createModel(m)
                self.partition_model.addModel(mod)

        # Copy the sitemodel vector (that contains for each site the index of the model assigned to that site)
        self.partition_model.setSiteAssignments(partition.getSiteModelVector())

        # Copy the number and proportion of sites in each subset of the partition
        self.partition_model.setNumSitesVect(partition.getSubsetSizes())

        # If user has specified subset relative rates in the partition interface, transfer
        # those now to partition_model
        if partition.subset_relrates is not None:
            num_subset_relrates = len(partition.subset_relrates)
            num_subsets = len(modelspecs)
            self.parent.phycassert(num_subset_relrates == num_subsets, 'Length of partition.subset_relrates list (%d) should equal the number of subsets defined (%d)' % (num_subset_relrates, num_subsets))
            self.partition_model.setSubsetRelRatesVect(partition.subset_relrates)

        self.likelihood = likelihood.TreeLikelihood(self.partition_model)
        # self.likelihood.setLot(self.r)
        self.likelihood.setUFNumEdges(self.parent.opts.uf_num_edges)
        if self.parent.data_matrix:
            #print '~!~!~!~!~! calling copyDataFromDiscreteMatrix !~!~!~!~!~' # temporary
            self.likelihood.copyDataFromDiscreteMatrix(self.parent.data_matrix, partition.getSiteModelVector())

        # Build the starting tree
        self.tree = self.parent.getStartingTree()

    def prepareForSimulation(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds data structures (conditional likelihood arrays and transition
        probability matrices) to nodes in tree to facilitate simulating data.
        Unlike prepareForLikelihood function, does not add data to the tips.

        """
        self.likelihood.prepareForSimulation(self.tree)

    def simulate(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a simulation, returning a new SimData object containing the
        simulated data.

        """
        sim_data = likelihood.SimData()
        self.parent.phycassert(self.parent.nchar > 0, 'nchar must be greater than zero in order to perform simulations')
        self.likelihood.simulateFirst(sim_data, self.tree, self.r, self.parent.nchar)
        return sim_data

    def prepareForLikelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds data structures (conditional likelihood arrays and transition
        probability matrices) to nodes in tree to facilitate likelihood
        calculations. Also adds the data to the tip nodes.

        """
        self.likelihood.prepareForLikelihood(self.tree)

    def calcLnLikelihood(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes and returns the log-likelihood.

        """
        return self.likelihood.calcLnL(self.tree)

    def debugCheckForUncachedCLAs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Scans tree for cached CLAs. If any are found, returns True. If no CLAs
        are currently cached, returns False. False is the expected result,
        because there should only be cached CLAs in the tree if we are in the
        middle of performing a Metropolis-Hastings move. If a move has been
        proposed, but not yet accepted or rejected, cached CLAs provide a way
        to return to the pre-move state after rejection without having to
        recalculate the likelihood.

        """
        return self.likelihood.debugCheckForUncachedCLAs(self.tree)


from _LikelihoodExt import *

class TreeLikelihood(TreeLikelihoodBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Likelihood objects contain methods for both manipulating trees (i.e.
    they are able to decorate trees with the data structures needed for
    computing likelihoods) and for computing likelihoods on trees once
    the trees have been equipped with these data structures. The following
    example illustrates how to open a nexus data file containing both
    data and trees, and compute the likelihood for the first tree. The
    reason True was supplied as the second argument to buildFromString
    is that tip nodes in tree descriptions read from a file start at 0,
    not 1. The True second argument informs buildFromString that the tip
    nodes are 0-based rather than the default 1-based.

    >>> from phycas import *
    >>> import time
    >>> reader = ReadNexus.NexusReader()
    >>> reader.readFile(getPhycasTestData('nyldna4.nex'))
    >>> data_matrix =  reader.getLastDiscreteMatrix()
    >>> model = Likelihood.HKYModel()
    >>> model.setStateFreqUnnorm(0, 1.0)
    >>> model.setStateFreqUnnorm(1, 2.0)
    >>> model.setStateFreqUnnorm(2, 3.0)
    >>> model.setStateFreqUnnorm(3, 4.0)
    >>> for f in model.getStateFreqs():
    ...     print "%.2f" % f
    0.10
    0.20
    0.30
    0.40
    >>> model.setKappaFromTRatio(2.0)
    >>> print "%.5f" % model.getKappa()
    4.36364
    >>> partition_model = Likelihood.PartitionModelBase()
    >>> partition_model.addModel(model)
    >>> likelihood = Likelihood.TreeLikelihood(partition_model)
    >>> likelihood.copyDataFromDiscreteMatrix(data_matrix, partition.getSiteModelVector())
    >>> for t in reader.getTrees():
    ...     tree = Phylogeny.Tree(t)
    ...     likelihood.prepareForLikelihood(tree)
    ...     lnL = likelihood.calcLnL(tree)
    ...     print '%.5f' % lnL
    -7812.79213
    
    """
    
    def __init__(self, model):
        """
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        Initializes the TreeLikelihoodBase object, establishing the model to
        be used.
        
        """
        TreeLikelihoodBase.__init__(self, model)

    def copyDataFromDiscreteMatrix(self, data_matrix, partition_info):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Copies data from a discrete data matrix object (such as that returned
        by the function ReadNexus.NexusReader.getLastDiscreteMatrix). The contents of
        data_matrix are compressed into data patterns and their counts.
        
        """
        TreeLikelihoodBase.copyDataFromDiscreteMatrix(self, data_matrix.raw_supermatrix, partition_info)

    def copyDataFromSimData(self, sim_data):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Copies data from a simulated data matrix (see SimData for details).
        
        """
        TreeLikelihoodBase.copyDataFromSimData(self, sim_data)

    def prepareForSimulation(self, tree):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a basic TipData object for all tip nodes and calls the
        function allocateInternalData for all internal nodes in the supplied
        tree. This function (prepareForSimulation) does not require a data
        matrix and does not add observed data to TipData structures. Silently
        returns if root node already has a TipData structure. This is taken
        to indicate that prepareForLikelihood or prepareForSimulation was
        previously called for this tree, in which case all data structures
        needed for simulation are already present.
        
        """
        TreeLikelihoodBase.prepareForSimulation(self, tree)

    def prepareForLikelihood(self, tree):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Prepares the specified tree for likelihood calculations by adding
        data structures to nodes for storing transition probability matrices
        and conditional likelihood arrays.
        
        """
        TreeLikelihoodBase.prepareForLikelihood(self, tree)

    def replaceModel(self, new_model):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Replaces the substitution model used in computing likelihoods with
        new_model.
        
        """
        TreeLikelihoodBase.replaceModel(self, new_model)

    def calcLnL(self, tree):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculates the log-likelihood using the current model. All transition
        probabilities and conditional likelihood arrays are recalculated, if
        necessary.
        
        """
        return TreeLikelihoodBase.calcLnL(self, tree)

    def calcLnLFromNode(self, nd):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculates the log-likelihood using the current model. All transition
        probabilities and conditional likelihood arrays are recalculated, if
        necessary.
        
        """
        return TreeLikelihoodBase.calcLnLFromNode(self, nd)

    def simulateFirst(self, sim_data, tree, lot, nchar):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Recalculates transition probabilities (or calculates them for the
        first time) in tree, then simulates data for nchar characters using
        the supplied pseudorandom number generator lot. The data is stored in
        the sim_data object. After simulateFirst is called to ensure that the
        transition probabilities are (re)computed, the function simulate can
        be called many more times to generate additional data sets using the
        same transition probabilities.
        
        """
        TreeLikelihoodBase.simulateFirst(self, sim_data, tree, lot, nchar)
    
    def simulate(self, sim_data, tree, lot, nchar):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Assuming that transition probabilities are present and calculated
        correctly in tree, simulates data for nchar characters using the
        supplied pseudorandom number generator lot. The data is stored in
        the sim_data object. After the function simulateFirst is called to
        ensure that the transition probabilities are (re)computed, simulate
        can be called many more times to generate additional data sets using
        the same transition probabilities. Only call this function after
        first calling simulateFirst at least once, otherwise the transition
        probabilities will contain garbage.
        
        """
        TreeLikelihoodBase.simulate(self, sim_data, tree, lot, nchar)
    
    def listPatterns(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string containing a list of all patterns in the form of a 
        table. The first column is the pattern index. The second column is the
        index of the partition subset to which the pattern belongs. The third
        column is the count (number of sites that had this pattern). Finally,
        the pattern itself is presented as a sequence of integer global state
        codes in the final column.
        
        """
        return TreeLikelihoodBase.listPatterns(self, False)

    def resetNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the number of likelihood evaluations to 0. The function
        getNumLikelihoodEvals can be called to obtain the number of times the
        likelihood has been evaluated since resetNumLikelihoodEvals was last
        called.
        
        """
        TreeLikelihoodBase.resetNumLikelihoodEvals(self)
    
    def getNumLikelihoodEvals(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reports the number of times the likelihood has been evaluated since
        the function resetNEvals was last called.
        
        """
        return TreeLikelihoodBase.getNumLikelihoodEvals(self)
    
    def addDataTo(self, other_sim_data):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds data currently stored to the patterns already stored in the
        supplied SimData object. Assumes that the patterns in this object 
        are based on the same number of taxa as patterns in other_sim_data.
        Adds patterns and counts only, i.e. does not preserve the identities
        of the original sites that represent a particular pattern.
        
        """
        TreeLikelihoodBase.addDataTo(self, other_sim_data)
    
    def recalcRelativeRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If the number of rates or the gamma shape parameter is changed in the
        model, this function should be called before using either the mean
        rates or the rate category boundaries. These quantities are calculated
        correctly when a TreeLikelihood object is constructed, but if either
        the number of rate categories or the shape parameter is modified in
        the model subsequently, this will not be reflected in the tuples
        returned by either getRateMeans or getCategoryLowerBoundaries until
        this function is called. For an example, see the getRateMeans
        function.
        
        """
        TreeLikelihoodBase.recalcRelativeRates(self)

    def getCategoryLowerBoundaries(self, subset_index):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns tuple comprising the lower bounds for all rate categories. In
        the discrete gamma model of among-site rate heterogeneity, the
        continuous gamma distribution is divided into ncat categories (where
        ncat is determined by the model). Normally, each category cuts out
        a probability mass equal to 1/ncat, but the model can specify unequal
        category probabilities. The tuple returned by this function specifies
        the x-coordinates of the lower boundary of each category. The only
        boundary not specified is thus the upper boundary of the last rate
        cateogory, but this is always infinity. For an example, see the
        getRateMeans function.
        
        """
        return TreeLikelihoodBase.getCategoryLowerBoundaries(self, subset_index)
    
    def getRateMeans(self, subset_index):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns tuple comprising the mean rate for all rate categories. In
        the discrete gamma model of among-site rate heterogeneity, the
        continuous gamma distribution is divided into ncat categories (where
        ncat is determined by the model). Normally, each category cuts out
        a probability mass equal to 1/ncat, but the model can specify unequal
        category probabilities. The tuple returned by this function specifies
        the mean rate in each category. These mean rates are used as
        representative relative rates in likelihood models designed to allow
        rates to vary across sites. The likelihood for a single nucleotide
        site under such a model with ncat = 3 could be written (crudely) as
        follows:

        L = Pr(data) = Pr(r0) Pr(data|r0) + Pr(r1) Pr(data|r1)
                          + Pr(r2) Pr(data|r2)

        Normally, Pr(r0) = Pr(r1) = Pr(r2) = 1/3, so this simplifies to

        L = (1/3)[Pr(data|r0) + Pr(data|r1) + Pr(data|r2)]

        The quantities r0, r1 and r3 are the mean rates returned by this
        function. Note that under this model, the likelihood must be computed
        three times, once for each relative rate, so discrete gamma models
        require roughly ncat times more computational effort than models
        that assume rate homogeneity. Here is an example producing output
        similar to the gammaplot command in PAUP* (without the plot). It
        creates a JC model with 4 rate categories (setting the number of rate
        categories to something greater than 1 triggers the discrete gamma
        version of the model). The shape parameter of the underlying
        continuous gamma distribution determines the amount of rate
        heterogeneity, with small values (e.g. 0.1) corresponding to high rate
        heterogeneity and large values (e.g. 2) corresponding to low rate
        heterogeneity. The tables produced show the lower bound for each rate
        category and the representative (mean) rate for that category. For
        example, in the case of gamma shape = 0.2, the fourth category extends
        from 0.894296 to infinity and the mean rate in this category is
        3.582035. Note that for gamma shape = 0.2, the highest rate is
        3.582035/0.000531 = 6745.8 times faster than the slowest rate. When
        the gamma shape parameter is set to 100, however, the highest rate is
        only 1.129806/0.875906 = 1.3 times faster than the slowest.
        
        >>> from phycas import *
        >>> 
        >>> m = Likelihood.JCModel()
        >>>
        >>> # First example shows high rate heterogeneity
        >>> 
        >>> m.setShape(0.2)
        >>> m.setNGammaRates(4)
        >>> partition_model = Likelihood.PartitionModelBase()
        >>> partition_model.addModel(m)
        >>> t = Likelihood.TreeLikelihood(partition_model)
        >>> t.setNoData()
        >>> mean  = t.getRateMeans(0)
        >>> lower = t.getCategoryLowerBoundaries(0)
        >>> print 'No. categories = %d' % (m.getNGammaRates())
        No. categories = 4
        >>> print 'Gamma shape = %.5f' % (m.getShape())
        Gamma shape = 0.20000
        >>> for i in range(4):
        ...     print '%12.6f %12.6f' % (lower[i], mean[i])
            0.000000     0.000531
            0.003188     0.033775
            0.103732     0.383658
            0.894296     3.582035
        >>>
        >>> # Second example shows low rate heterogeneity
        >>> 
        >>> m.setShape(100.0)
        >>> t.recalcRelativeRates()
        >>> mean  = t.getRateMeans(0)
        >>> lower = t.getCategoryLowerBoundaries(0)
        >>> print 'Gamma shape = %.5f' % (m.getShape())
        Gamma shape = 100.00000
        >>> for i in range(4):
        ...     print '%12.6f %12.6f' % (lower[i], mean[i])
            0.000000     0.875906
            0.930858     0.964739
            0.996669     1.029549
            1.065511     1.129806
            
        """
        return TreeLikelihoodBase.getRateMeans(self, subset_index)

    def getNPatterns(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of data patterns currently stored.
        
        """
        return TreeLikelihoodBase.getNPatterns(self)
    
    def getNoData(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if running with no data, False if data has been stored.
        
        """
        return TreeLikelihoodBase.getNoData(self)
    
    def usePatternSpecificRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Causes pattern-specific rates model to be used instead of a discrete
        gamma rate heterogeneity model. In the pattern-specific rates model,
        each pattern has its own relative rate and every site that shows that
        pattern is assumed to evolve at that pattern-specific relative rate.
        
        """
        TreeLikelihoodBase.setUsePatternSpecificRates(self, True)
    
    def doNotUsePatternSpecificRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Causes discrete gamma rate heterogeneity model to be used instead of
        a pattern-specific rates model. In the discrete gamma model, each
        site is assumed to evolve at one of a number of relative rates,
        whereas in the pattern-specific rates model, each pattern has its own
        relative rate parameter and every site that shows that pattern is
        assumed to evolve at that pattern-specific relative rate.
        
        """
        TreeLikelihoodBase.setUsePatternSpecificRates(self, False)
        
    def setUFNumEdges(self, nedges):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets number of edges that must be traversed before action is taken to
        prevent underflow. Likelihood calculations will be faster for higher
        values of nedges, but the tradeoff is that higher values of nedges
        also mean greater probability of underflow.
        
        """
        TreeLikelihoodBase.setUFNumEdges(self, nedges)
        
    def startTreeViewer(self, t, s, i):
        import phycas.TreeViewer
        tv = phycas.TreeViewer.TreeViewer(tree=t, msg=s, site=i)
        tv.setLikelihoodRoot(TreeLikelihoodBase.getLikelihoodRootNodeNum(self))
        return tv.run()

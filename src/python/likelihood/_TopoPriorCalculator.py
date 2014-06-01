from _LikelihoodExt import *

class PolytomyTopoPriorCalculator(PolytomyTopoPriorCalculatorBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This class can compute one of two possible topology priors: the
    polytomy prior and the resolution class prior. Both of these priors
    were described in the following paper:
    
    Lewis, P. O., M. T. Holder and K. E. Holsinger. 2005. Polytomies and
    Bayesian Phylogenetic Inference. Systematic Biology 54(2): 241-253.

    If a "resolution class" prior is used in an MCMC analysis with no
    data, so that the posterior distribution equals the prior
    distribution, the number of tree topologies sampled having m - 1
    internal nodes (termed the m - 1 "resolution class") will be
    approximately C times greater than the number of trees having m
    internal nodes (the m resolution class). The value of C can be set
    using the function setC.

    On the other hand, the "polytomy" prior assigns a tree topology with
    m - 1 internal nodes C times more prior probability than a tree
    topology with m internal nodes. There may be, however, more possible
    trees with m - 1 internal nodes, so the number of possible trees in
    each resolution class is an important consideration and may, in some
    cases, appear to overturn the effects of the polytomy prior. For
    example, if a polytomy prior is used in an MCMC analysis with no data,
    the number of tree topologies sampled having m - 1 internal nodes will
    be approximately T_{n, m - 1}*C/T_{n, m} times greater than the number
    of trees having m internal nodes. The values T_{n, m - 1} and T_{n, m}
    equal the number of tree topologies having m - 1 internal nodes and
    the number of tree topologies having m internal nodes, respectively.
    Thus, if T_{n, m - 1} is smaller than T_{n, m} by a factor greater
    than C, the resolution class m - 1 will have fewer tree topologies in
    it than the resolution class m, even though C was chosen to favor the
    m - 1 class. Note that this result is based on running the MCMC
    analysis without data so that the posterior equals the prior. If run
    with data, the information coming into the posterior through the
    likelihood can often overwhelm what appears to be a strong topology
    prior.

    Here is an example designed to replicate Table 2, p. 248, in the Lewis,
    Holder and Holsinger (2005) paper:

    >>> import math
    >>> from phycas import *
    >>> 
    >>> tpc = Likelihood.PolytomyTopoPriorCalculator()
    >>> tpc.setNTax(17)
    >>> tpc.chooseUnrooted()
    >>> 
    >>> # The columns in the table are stored as lists c1, c2, ..., c7
    >>> # This first column is for the resolution classes
    >>> c1 = [i+1 for i in range(15)]
    >>> 
    >>> # compute counts of trees in different resolution classes
    >>> c2 = tpc.getLnCounts()
    >>> 
    >>> # compute polytomy prior with C = 1
    >>> tpc.setC(1.0)
    >>> tpc.choosePolytomyPrior()
    >>> c3 = tpc.getRealizedResClassPriorsVect()
    >>> 
    >>> # compute polytomy prior with C = e
    >>> tpc.setC(math.exp(1.0))
    >>> c4 = tpc.getRealizedResClassPriorsVect()
    >>> 
    >>> # compute resolution class prior with C approx. 1.1
    >>> tpc.setC(1.0/0.9)
    >>> tpc.chooseResolutionClassPrior()
    >>> c5 = tpc.getRealizedResClassPriorsVect()
    >>> 
    >>> # compute resolution class prior with C = 2.0
    >>> tpc.setC(2.0)
    >>> c6 = tpc.getRealizedResClassPriorsVect()
    >>> 
    >>> # compute resolution class prior with C = 2.0
    >>> tpc.setC(10.0)
    >>> c7 = tpc.getRealizedResClassPriorsVect()
    >>> 
    >>> # create the table as a list of strings, then join them and output all at once
    >>> table2 = []
    >>> s = '%19s %25s %38s' % ('', '--- Polytomy prior --', '----- Resolution class prior -----')
    >>> table2.append(s)
    >>> s = '%6s %12s %12s %12s %12s %12s %12s' % ('m', 'log10(Tnm)', 'C=1', 'C=e', 'C=1.1', 'C=2', 'C=10')
    >>> table2.append(s)
    >>> for m in c1:
    ...     s = '%6d %12.6f %12.6f %12.6f %12.6f %12.6f %12.6f' % (
    ...        m,
    ...        c2[m]/math.log(10.0),
    ...        math.exp(c3[m] - c3[0]),
    ...        math.exp(c4[m] - c4[0]),
    ...        math.exp(c5[m] - c5[0]),
    ...        math.exp(c6[m] - c6[0]),
    ...        math.exp(c7[m] - c7[0]))
    ...     table2.append(s)
    >>> print '\n'.join(table2)
                            --- Polytomy prior --     ----- Resolution class prior -----
         m   log10(Tnm)          C=1          C=e        C=1.1          C=2         C=10
         1     0.000000     0.000000     0.000000     0.125927     0.500015     0.900000
         2     4.816361     0.000000     0.000000     0.113335     0.250008     0.090000
         3     7.801126     0.000000     0.000001     0.102001     0.125004     0.009000
         4    10.002016     0.000000     0.000047     0.091801     0.062502     0.000900
         5    11.729312     0.000002     0.000931     0.082621     0.031251     0.000090
         6    13.118317     0.000055     0.008388     0.074359     0.015625     0.000009
         7    14.240733     0.000730     0.040907     0.066923     0.007813     0.000001
         8    15.138396     0.005766     0.118896     0.060231     0.003906     0.000000
         9    15.836311     0.028761     0.218166     0.054208     0.001953     0.000000
        10    16.348281     0.093490     0.260894     0.048787     0.000977     0.000000
        11    16.679054     0.200235     0.205562     0.043908     0.000488     0.000000
        12    16.824003     0.279569     0.105584     0.039517     0.000244     0.000000
        13    16.765475     0.244321     0.033945     0.035566     0.000122     0.000000
        14    16.460717     0.121117     0.006190     0.032009     0.000061     0.000000
        15    15.791711     0.025954     0.000488     0.028808     0.000031     0.000000
        
    """

    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The constructor takes no arguments and simply initializes the private
        base class data members.
        
        """
        PolytomyTopoPriorCalculatorBase.__init__(self)
        
    def setNTax(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of taxa to n. Assumes that n is greater than 1 if
        using rooted trees and greater than 2 if using unrooted trees.
        
        """
        PolytomyTopoPriorCalculatorBase.setNTax(self, n)
        
    def getNTax(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of taxa on which priors and topology counts are
        currently based.
        
        """
        return PolytomyTopoPriorCalculatorBase.getNTax(self)

    def chooseRooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Forces recalculation of counts and priors for rooted trees if they had
        previously been based on unrooted trees. There are more rooted than
        unrooted trees for the same number of taxa, so this setting is
        important when asking questions that require knowledge of the numbers
        of possible trees.
        
        """
        PolytomyTopoPriorCalculatorBase.chooseRooted(self)
        
    def chooseUnrooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Forces recalculation of counts and priors for unrooted trees if they 
        had previously been based on rooted trees. There are more rooted than
        unrooted trees for the same number of taxa, so this setting is
        important when asking questions that require knowledge of the numbers
        of possible trees.
        
        """
        PolytomyTopoPriorCalculatorBase.chooseUnrooted(self)
        
    def getCount(self, n, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of trees having n taxa and m internal nodes. Calls
        RecalcCountsAndPriors function if n is not equal to the current
        number of taxa. Assumes m is greater than 0. For rooted trees, assumes
        m is less than the number of taxa. For unrooted trees, assumes m less
        than the number of taxa minus one. 
        
        """
        return PolytomyTopoPriorCalculatorBase.getCount(self, n, m)

    def getSaturatedCount(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of saturated (i.e. fully-resolved and thus having
        as many internal nodes as possible) trees of n taxa.
        
        """
        return PolytomyTopoPriorCalculatorBase.getSaturatedCount(self, n)

    def getTotalCount(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the total number of trees (including all resolution classes
        from star to saturated trees) for n taxa.
        
        """
        return PolytomyTopoPriorCalculatorBase.getTotalCount(self, n)

    def getCountsVect(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns copy of the counts vector, which contains in its first
        element the total number of tree topologies having any number of
        internal nodes, and which contains in its mth element (m > 0) the
        number of tree topologies having exactly m internal nodes.
        
        """
        return PolytomyTopoPriorCalculatorBase.getCountsVect(self)

    def chooseResolutionClassPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Forces recalculation of priors if polytomy priors had previously been 
        computed. If a resolution class prior is used in an MCMC analysis with
        no data, the number of tree topologies sampled having m - 1 internal
        nodes will be approximately C times greater than the number of trees
        having m internal nodes (the value C can be set using the function
        setC).
        
        """
        PolytomyTopoPriorCalculatorBase.chooseResolutionClassPrior(self)
        
    def choosePolytomyPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Forces recalculation of priors if resolution class priors had
        previously been computed. The polytomy prior assigns a tree topology
        with m - 1 internal nodes C times more prior probability than a tree
        topology with m internal nodes. There may be, however, more possible
        trees with m - 1 internal nodes, so the number of possible trees in
        each resolution class is an important consideration and may, in some
        cases, appear to overturn the effects of the polytomy prior. For
        example, if a polytomy prior is used in an MCMC analysis with no data,
        the number of tree topologies sampled having m - 1 internal nodes will
        be approximately T_{n, m - 1}*C/T_{n, m} times greater than the number
        of trees having m internal nodes (the value C can be set using the
        function setC). The values T_{n, m - 1} and T_{n, m} equal the number
        of tree topologies having m - 1 internal nodes and the number of tree
        topologies having m internal nodes, respectively. Thus, if T_{n, m - 1}
        is smaller than T_{n, m} by a factor greater than C, the resolution
        class m - 1 will have fewer tree topologies in it than the resolution
        class m, even though C was chosen to favor the m - 1 class. Note that
        this result is based on running the MCMC analysis without data so that
        the posterior equals the prior. If run with data, the information
        coming into the posterior through the likelihood can often overwhelm
        what appears to be a strong topology prior.
        
        """
        PolytomyTopoPriorCalculatorBase.choosePolytomyPrior(self)
        
    def setC(self, c):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value of C (which determines the amount by which less
        resolved tree topologies are favored over more resolved tree
        topologies) to c. Assumes c is greater than 0.0. Values of c greater
        than 1 favor less resolved tree topologies; values of c less than 1
        favor more resolved tree topologies. Specifying c = 1.0 and calling
        the choosePolytomyPrior function produces a prior that is flat across
        all possible topologies, from the star tree to fully-resolved (or
        saturated) tree topologies.
        
        """
        PolytomyTopoPriorCalculatorBase.setC(self, c)
        
    def getC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current value of C, which determines the amount by which
        less resolved tree topologies are favored over more resolved tree
        topologies. See the documentation for setC for more information about
        the effects of the value C.
        
        """
        return PolytomyTopoPriorCalculatorBase.getC(self)

    def getLnTopologyPrior(self, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the natural logarithm of the unnormalized topology prior. This
        represents the resolution class prior (or polytomy prior) if the
        chooseResolutionClassPrior (choosePolytomyPrior) function is called.
        The value of m must be greater than 0 and less than or equal to the
        number of internal nodes in a fully resolved tree (which equals the
        number of taxa minus 1 for rooted trees and the number of taxa minus
        2 for unrooted trees).
        
        """
        return PolytomyTopoPriorCalculatorBase.getLnTopologyPrior(self, m)

    def getLnNormalizedTopologyPrior(self, m):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the natural logarithm of the normalized topology prior. This
        represents the resolution class prior (polytomy prior) if the function
        chooseResolutionClassPrior (choosePolytomyPrior) is called. The value
        of m must be greater than 0 and less than or equal to the number of
        internal nodes in a fully resolved tree (which equals the number of
        taxa minus 1 for rooted trees and the number of taxa minus 2 for
        unrooted trees). The log of the normalized topology prior is obtained
        as getLnTopologyPrior(m) minus getLnNormConstant().
        
        """
        return PolytomyTopoPriorCalculatorBase.getLnNormalizedTopologyPrior(self, m)

    def getLnNormConstant(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the natural logarithm of the normalizing constant for the
        topology prior.
        
        """
        return PolytomyTopoPriorCalculatorBase.getLnNormConstant(self)

    def getTopoPriorVect(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a copy of the vector of priors, which contains in its mth
        element the unnormalized prior for tree topologies having exactly m
        internal nodes. The 0th element of this vector holds the normalizing
        constant. The functions getLnTopologyPrior and getLnNormConstant can
        be used to obtain the same information.
        
        """
        return PolytomyTopoPriorCalculatorBase.getTopoPriorVect(self)

    def isResolutionClassPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the prior currently calculated corresponds to the
        resolution class prior described in Lewis, Holder and Holsinger
        (2005).
        
        """
        return PolytomyTopoPriorCalculatorBase.isResolutionClassPrior(self)

    def isPolytomyPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the prior currently calculated corresponds to the
        polytomy prior described in Lewis, Holder and Holsinger (2005).
        
        """
        return PolytomyTopoPriorCalculatorBase.isPolytomyPrior(self)

    def isRooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the counts and prior currently calculated are correct
        for rooted trees.
        
        """
        return PolytomyTopoPriorCalculatorBase.isRooted(self)

    def isUnrooted(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the counts and prior currently calculated are correct
        for unrooted trees.
        
        """
        return PolytomyTopoPriorCalculatorBase.isUnrooted(self)

    def getRealizedResClassPriorsVect(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Constructs a vector of realized resolution class priors from the
        vector of priors already stored. The mth element of the returned
        vector is set to T_{n,m}*p_{n,m} for m > 0, where T_{n,m} is the
        number of trees for n taxa and m internal nodes, and p_{n,m} is the
        prior currently being used in MCMC analyses for tree topologies
        having n taxa and m internal nodes. The first element of the returned
        vector holds the normalization constant (sum of all other elements).
        This function is not efficient because it is intended only to be used
        for providing information to the user on request. Table 2, p. 248, in
        the "Polytomies and Bayesian Phylogenetic Inference" paper (Lewis,
        P. O., M. T. Holder and K. E. Holsinger. 2005. Systematic Biology
        54(2): 241-253) presented (normalized) values from this vector. The
        main example given for the PolytomyTopoPriorCalculator class shows how the
        the values from Table 2 can be generated.
        
        """
        return PolytomyTopoPriorCalculatorBase.getRealizedResClassPriorsVect(self)


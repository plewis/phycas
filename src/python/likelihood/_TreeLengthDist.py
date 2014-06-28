from phycas.probdist._PyDistributionBase import PyDistributionBase
from _LikelihoodExt import *

class TreeLengthDist(TreeLengthDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the Rannala, Zhu, and Yang (2012. MBE 29:325-335) joint 
    edge length/tree length distribution.

    Notes:
      - the tree length distribution has four parameters:
        1. Gamma tree length distribution shape parameter
        2. Gamma tree length distribution scale parameter
        3. Dirichlet parameter used for external edge lengths
        4. the ratio of internal/external edge lengths
      - the Gamma distribution used here has mean = shape/scale, which
        differs from the definition of the scale parameter used 
        everywhere else in Phycas (mean = shape*scale)

    """
    def __init__(self, shape, scale, external_edgelen_param, internal_external_edgelen_ratio):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the shape and scale parameters of the Gamma distribution of
        tree length (mean = shape/scale), the Dirichlet parameter governing
        external edge lengths, and the ratio of internal/external edge
        lengths.

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1.0, 0.1, 1.0, 0.8)

        """
        TreeLengthDistBase.__init__(self, shape, scale, external_edgelen_param, internal_external_edgelen_ratio)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this TreeLengthDist object.
        
        """
        #return TreeLengthDistBase.clone(self)
        return TreeLengthDist(self.getShape(), self.getScale(), self.getExtEdgelenParam(), self.getIntExtEdgelenRatio())
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'TreeLengthDistribution'
        
        """
        return TreeLengthDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another 
        TreeLengthDist object identical to this one. e.g.,

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1.0, 0.1, 1.0, 0.8)
        >>> print d.__str__()
        TreeLengthDist(1.00000, 0.10000, 1.00000, 0.80000)
        
        """
        return TreeLengthDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another 
        TreeLengthDist object identical to this one. e.g.,

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1.0, 0.1, 1.0, 0.8)
        >>> print d.__repr__()
        TreeLengthDist(1.00000, 0.10000, 1.00000, 0.80000)
        
        """
        return TreeLengthDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> from phycas.likelihood import *
        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = TreeLengthDist(1.0, 0.1, 1.0, 0.8)
        >>> d2 = TreeLengthDist(1.0, 0.05, 1.0, 1.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample(4)
        0.000177184566
        >>> print "%.12f" % d2.sample(4)
        0.930716753542
        
        In this example, both d1 and d2 use the same random number generator
        when their sample functions are called. Generating all random numbers
        from a single random number generator has the advantage of allowing
        an entire analysis to be recreated by simply specifying the original
        seed (here, 1357) to the master Lot object. If setLot is not used,
        each distribution object maintains its own random number generator,
        which is initialized using the system clock at the time the object is
        created. This behavior is convenient, but makes it very hard to
        replicate results.

        """
        return TreeLengthDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1.0, 0.05, 1.0, 1.0)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample(4)
        0.140064780761
        >>> print "%.12f" % d.sample(4)
        8.204670625640
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample(4)
        0.140064780761

        """
        return TreeLengthDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return TreeLengthDistBase.resetLot(self)
        
    def sample(self, ntax):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the tree length distribution by
        this TreeLengthDist object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1,.1,1,1)
        >>> d.setSeed(97531)
        >>> samples = [d.sample(4) for i in range(3)]
        >>> for x in samples:
        ...    print round(x,6)
        8.302218
        0.490945
        13.733702
        
        """
        return TreeLengthDistBase.sample(self, ntax)
        
    def getLnPDF(self, t):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function for the supplied tree t.
        Returns the natural logarithm of the density. e.g.,

        >>> from phycas.likelihood import *
        >>> d = TreeLengthDist(1, .1, 1, 1)
        >>> tree = Tree()
        >>> tree.buildFromString('(a:.1,b:.1,(c:.1,d:.1):.1')
        >>> ln_pdf = d.getLnPDF(tree)
        >>> print round(ln_pdf, 9)
        -2.291759469
        
        """
        return TreeLengthDistBase.getLnPDF(self, t)
        
    def getRelativeLnPDF(self, t):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Simply returns getLnPDF(t).
        
        """
        return TreeLengthDistBase.getRelativeLnPDF(self, t)
        
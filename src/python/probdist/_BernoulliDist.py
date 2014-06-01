from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Bernoulli(BernoulliDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate bernoulli probability distribution. The
    bernoulli probability distribution governs the probability of success
    on a single trial when the probability of success is p and the
    probability of failure is 1-p.

    """
    def __init__(self, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the probability p of a success. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getMean()
        0.1

        """
        BernoulliDistBase.__init__(self, p)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Bernoulli distribution.
        
        """
        return BernoulliDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns True because the bernoulli distribution is discrete.
        
        """
        return BernoulliDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Bernoulli'
        
        """
        return BernoulliDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        Bernoulli object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.__str__()
        Bernoulli(0.10000)
        
        """
        return BernoulliDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        Bernoulli object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.__repr__()
        Bernoulli(0.10000)
        
        """
        return BernoulliDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Bernoulli(0.5)
        >>> d2 = Bernoulli(0.2)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> [d1.sample(), d2.sample()]
        [1.0, 0.0]
        
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
        return BernoulliDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> d.setSeed(135)
        >>> d.sample()
        1.0
        >>> d.sample()
        0.0
        >>> d.setSeed(135)
        >>> d.sample()
        1.0

        """
        return BernoulliDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return BernoulliDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getMean()
        0.1

        """
        return BernoulliDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getVar()
        0.09

        """
        return BernoulliDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getStdDev()
        0.3

        """
        return BernoulliDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the sum of probabilities up to an including the state
        specified.

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getCDF(0)
        0.9
        >>> print d.getCDF(1)
        1.0

        """
        return BernoulliDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the bernoulli distribution
        specified by this Bernoulli object. Python list comprehensions
        can be used to store many simulated samples for use in subsequent
        calculations.

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> d.setSeed(97531)
        >>> print [d.sample()] * 3
        [0.0, 0.0, 0.0]
        
        """
        return BernoulliDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability of the supplied value x. Returns the
        natural logarithm of the probability of x. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getLnPDF(0)
        -0.105360515658
        >>> print d.getLnPDF(1)
        -2.30258509299
        
        """
        return BernoulliDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For Bernoulli, this function is identical to the getLnPDF
        function. For more complicated density functions, the
        getRelativeLnPDF returns a value proportional to the density function
        but does not calculate terms that are constant to save computation.
        
        """
        return BernoulliDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. This distribution is
        determined entirely by the mean, so the var argument is ignored. The
        reason this function requires both mean and variance is for
        compatibility with functions of the same name in other distributions.

        >>> from phycas.ProbDist import *
        >>> d = Bernoulli(0.1)
        >>> print d.getMean()
        0.1
        >>> print d.getVar()
        0.09
        >>> d.setMeanAndVariance(0.5, 0.25)
        >>> print d.getMean()
        0.5
        >>> print d.getVar()
        0.25
        
        """
        return BernoulliDistBase.setMeanAndVariance(self, mean, var)
    

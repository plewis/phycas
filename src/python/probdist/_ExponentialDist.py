from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Exponential(ExponentialDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate exponential probability distribution.

    Notes:
      - an exponential distribution with hazard parameter h has mean 1/h
      - the mean equals the standard deviation in exponential distributions
      - an exponential distribution with hazard parameter h is equivalent to
        a gamma distribution having shape 1 and scale parameter 1/h

    """
    def __init__(self, hazard):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the hazard parameter (i.e. inverse of the mean) when
        initializing an Exponential object. e.g.,

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getMean()
        0.5
        
        """
        ExponentialDistBase.__init__(self, hazard)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Exponential distribution.
        
        """
        return ExponentialDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the exponential distribution is
        continuous.
        """
        return ExponentialDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Exponential'
        """
        return ExponentialDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        Exponential object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.__str__()
        Exponential(2.00000)
        
        """
        return ExponentialDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        Exponential object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.__repr__()
        Exponential(2.00000)
        
        """
        return ExponentialDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Exponential(3.0)
        >>> d2 = Exponential(2.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print "%.12f" % d1.sample()
        0.003559060607
        >>> print "%.12f" % d2.sample()
        0.343362432309
        
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
        return ExponentialDistBase.setLot(self, lot)
        
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
        >>> b = Exponential(2)
        >>> b.setSeed(135)
        >>> print "%.12f" % b.sample()
        0.000528559199
        >>> print "%.12f" % b.sample()
        0.708585892848
        >>> b.setSeed(135)
        >>> print "%.12f" % b.sample()
        0.000528559199

        """
        return ExponentialDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return ExponentialDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getMean()
        0.5

        """
        return ExponentialDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getVar()
        0.25

        """
        return ExponentialDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getStdDev()
        0.5

        """
        return ExponentialDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 0.0
        up to the value x.

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getCDF(1.5)
        0.950212931632

        The integral of the exponential density is 1 - exp(-lambda*x), where
        lambda = 2 and x = 1.5 in this case. The CDF is thus 1 - exp(-3),
        which equals 0.950212931632.

        """
        return ExponentialDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the exponential distribution
        specified by this Exponential object. Python list comprehensions
        can be used to store many simulated samples for use in subsequent
        calculations.

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> b.setSeed(97531)
        >>> print "%.12f" % b.sample()
        0.720509647019
        >>> print "%.12f" % b.sample()
        0.006044079697
        >>> print "%.12f" % b.sample()
        1.429544116814
        
        """
        return ExponentialDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getLnPDF(1.5)
        -2.30685281944

        The density is lambda*exp(-lambda*x), where lambda=2 and x=1.5. The
        log density is thus -lambda*x + log(lambda) = -(2)(1.5) + log(2),
        which equals -2.30685281944.

        """
        return ExponentialDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability density function at the supplied
        value x. Returns the natural logarithm of the relative density at x.
        Use this function if speed is important but normalization is not, say
        in MCMC calculations. Use getLnPDF instead if you need to have a
        correctly normalized density value (i.e. from a density function that
        integrates to 1.0)

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getRelativeLnPDF(1.5)
        -3.0

        The density is lambda*exp(-lambda*x), where lambda=2 and x=1.5. The
        log density is thus -lambda*x + log(lambda) = -(2)(1.5) + log(2).
        The second term is not needed if only the relative PDF is required
        because it does not contain x, so the return value is just the first
        term, -3.0.
        
        """
        return ExponentialDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var=0.0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this distribution. This distribution is
        determined entirely by the mean, so no variance need be provided.
        The reason this function even has a variance argument is for
        compatibility with functions of the same name in other distributions.

        >>> from phycas.ProbDist import *
        >>> b = Exponential(2)
        >>> print b.getMean()
        0.5
        >>> print b.getVar()
        0.25
        >>> b.setMeanAndVariance(5, 0)
        >>> print b.getMean()
        5.0
        >>> print b.getVar()
        25.0
        
        """
        return ExponentialDistBase.setMeanAndVariance(self, mean, var)

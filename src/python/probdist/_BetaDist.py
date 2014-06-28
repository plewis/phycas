from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class Beta(BetaDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate beta probability distribution.

    Notes:
      - the Beta distribution has two parameters, a and b
      - a Beta distribution in which a = b is symmetric around 0.5
      - a Beta distribution in which a = b = 1 is identical to a
        Uniform(0,1) distribution
      - the mean of a Beta distribution is a/(a+b) and its variance
        is a*b/[(a+b)^2*(a+b+1)]

    """
    def __init__(self, a, b):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the a and b parameters of the Beta object. e.g.,

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getMean()
        0.4

        """
        BetaDistBase.__init__(self, a, b)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Beta distribution.

        >>> from phycas.probdist import *
        >>> d1 = Beta(2, 3)
        >>> print d1.getMean()
        0.4
        >>> d2 = d1.clone()
        >>> d2.setMeanAndVariance(0.5, 0.125)
        >>> print 'd1 mean = %.1f' % d1.getMean()
        d1 mean = 0.4
        >>> print 'd2 mean = %.1f' % d2.getMean()
        d2 mean = 0.5
        
        """
        return BetaDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the Beta distribution is continuous.
        
        """
        return BetaDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Beta'
        
        """
        return BetaDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Beta
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.__str__()
        Beta(2.00000, 3.00000)
        
        """
        return BetaDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another Beta
        object identical to this one. e.g.,

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.__repr__()
        Beta(2.00000, 3.00000)
        
        """
        return BetaDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = Beta(1.0, 1.0)
        >>> d2 = Beta(2.0, 3.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print round(d1.sample(), 8), round(d2.sample(), 8)
        0.01530994 0.35861783
        
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
        return BetaDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the random number generator of this distribution object
        using the supplied seed. Note that if you have called setLot before
        this point, calling setSeed is pointless because you have already
        replaced the random number generator for which you are setting the
        seed! If you have already called setLot, you probably want to call
        the setSeed function of that Lot ojbect. e.g.,

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.011620973198
        >>> print "%.12f" % d.sample()
        0.021080389846
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.011620973198

        """
        return BetaDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return BetaDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample is called to generate samples
        from this distribution).

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getMean()
        0.4

        """
        return BetaDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample is called to generate
        samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getVar()
        0.04

        """
        return BetaDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getStdDev()
        0.2

        """
        return BetaDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 0.0
        up to the value x.

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getCDF(0.35)
        0.43701875

        """
        return BetaDistBase.getCDF(self, x)
        
    def getQuantile(self, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given p, the integral of the probability density function from 0.0
        up to a value x, returns x.

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print '%.5f' % d.getQuantile(0.35)
        0.30096

        """
        return BetaDistBase.getQuantile(self, p)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the Beta distribution specified by
        this Beta object. Python list comprehensions can be used to
        store many simulated samples for use in subsequent calculations.

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> d.setSeed(97531)
        >>> print round(d.sample(), 8)
        0.8555905
        >>> print round(d.sample(), 8)
        0.59366384
        >>> print round(d.sample(), 8)
        0.37280974
        
        """
        return BetaDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getLnPDF(0.35)
        0.573518693104

        The Beta density is

        Gamma(a+b) p^(a-1) (1-p)^(b-1)
        ------------------------------
             Gamma(a) Gamma(b)
        
        where a = 2, b = 3, p = 0.35, and Gamma is the Gamma function, not
        the Gamma probability distribution. The density is thus

        Gamma(5) (0.35)^(2-1) (0.65)^(3-1)   4! (0.35) (0.65)^2
        ---------------------------------- = ------------------ = 1.7745
                Gamma(2) Gamma(3)                 1! 2!

        and the log density is thus 0.57351869310441410713636637798507
        Oops - looks like getLnPDF returns only relative PDF
        
        """
        return BetaDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the relative probability density function at the supplied
        value x. Returns the natural logarithm of the relative density at x.
        Use this function if speed is important but normalization is not, say
        in MCMC calculations. Use getLnPDF instead if you need to have a
        correctly normalized density value (i.e. from a density function that
        integrates to 1.0)

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getRelativeLnPDF(0.35)
        -1.91138795668

        The Beta density is

        Gamma(a+b) p^(a-1) (1-p)^(b-1)
        ------------------------------
             Gamma(a) Gamma(b)
        
        where a = 2, b = 3, p = 0.35, and Gamma is the Gamma function, not
        the Gamma probability distribution. The relative density requires
        only the two terms containing p in the numerator, so the natural
        logarithm of the relative density is

        (a-1) log(p) + (b-1) log(1-p)
        
        For the example given, this equals (2-1) log(0.35) + (3-1) log(0.65)
        = log(0.35 + 2 log(0.65) = -1.9113879566835862030933431018538
        """
        return BetaDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this Beta distribution. The values of
        parameters a and b can be determined from mean and var as follows:
           a = phi*mean
           b = phi*(1 - mean)
           phi = [mean*(1 - mean)/var] - 1
        So if mean = 0.25, var = 0.0275, we have:
           phi = 0.25*0.75/0.0275 - 1 = 5.81818
           a = 5.81818*0.25 = 1.45455
           b = 5.81818*0.75 = 4.36364
        Check:
           mean = a/(a + b)
                = 1.45455/(1.45455 + 4.36364)
                = 0.25
           var  = ab/{[(a + b)^2]*(a + b + 1)}
                = 6.347132562/(33.8513348761*6.81819)
                = 0.0275

        >>> from phycas.probdist import *
        >>> d = Beta(2, 3)
        >>> print d.getMean()
        0.4
        >>> print d.getVar()
        0.04
        >>> d.setMeanAndVariance(0.25, 0.0275)
        >>> print d.getMean()
        0.25
        >>> print d.getVar()
        0.0275
        
        """
        return BetaDistBase.setMeanAndVariance(self, mean, var)
    

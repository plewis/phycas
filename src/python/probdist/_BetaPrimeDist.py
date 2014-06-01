from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class BetaPrime(BetaPrimeDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate beta prime probability distribution with 
    shape parameters a and b. The mean is a/(b-1) and the variance is
    a*(a+b+1)/[(b-2)*(b-1)^2].

    """
    def __init__(self, a, b):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Specify the a and b parameters of the BetaPrime object. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.getMean()
        1.0

        """
        BetaPrimeDistBase.__init__(self, a, b)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this BetaPrime distribution.

        >>> from phycas.ProbDist import *
        >>> d1 = BetaPrime(2, 3)
        >>> print d1.getMean()
        1.0
        >>> d2 = d1.clone()
        >>> d2.setMeanAndVariance(0.5, 0.025)
        >>> print 'd1 mean = %.3f' % d1.getMean()
        d1 mean = 1.000
        >>> print 'd2 mean = %.3f' % d2.getMean()
        d2 mean = 0.500
        
        """
        return BetaPrimeDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the BetaPrime distribution is continuous.
        
        """
        return BetaPrimeDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'BetaPrime'
        
        """
        return BetaPrimeDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another BetaPrime
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.__str__()
        BetaPrime(2.00000, 3.00000)
        
        """
        return BetaPrimeDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another BetaPrime
        object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.__repr__()
        BetaPrime(2.00000, 3.00000)
        
        """
        return BetaPrimeDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Substitutes a different random number generator to use when drawing
        samples. e.g.,

        >>> g = Lot()
        >>> g.setSeed(1357)
        >>> d1 = BetaPrime(1.0, 1.0)
        >>> d2 = BetaPrime(2.0, 3.0)
        >>> d1.setLot(g)
        >>> d2.setLot(g)
        >>> print round(d1.sample(), 8), round(d2.sample(), 8)
        0.01554798 0.55913284
        
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
        return BetaPrimeDistBase.setLot(self, lot)
        
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
        >>> d = BetaPrime(2, 3)
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.011757608046
        >>> print "%.12f" % d.sample()
        0.021534342175
        >>> d.setSeed(135)
        >>> print "%.12f" % d.sample()
        0.011757608046

        """
        return BetaPrimeDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Resets the random number generator to point to the local Lot object.
        Because the local Lot object is used by default, this function need
        only be called if setLot has previously been called to specify an
        external random number generator.

        """
        return BetaPrimeDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the mean of the distribution. This is the theoretical mean
        (i.e., it will not change if sample() is called to generate samples
        from this distribution). The formula for the mean is a/(b-1). Note
        that the mean of this distribution is undefined if b <= 1.

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.getMean()
        1.0

        """
        return BetaPrimeDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the variance of the distribution. This is the theoretical
        variance (i.e., it will not change if sample() is called to generate
        samples from this distribution). The formula for the variance is 
        a*(a+b-1)/[(b-2)*(b-1)^2]. Note that the variance of this distribution 
        is undefined if b <= 2.

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.getVar()
        2.0

        """
        return BetaPrimeDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the standard deviation of the distribution. This is the
        theoretical standard deviation (i.e., it will not change if sample is
        called to generate samples from this distribution).

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print '%.3f' % d.getStdDev()
        1.414

        """
        return BetaPrimeDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the cumulative distribution function at the supplied value
        x. This is the integral of the probability density function from 0.0
        up to the value x.

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(1, 1)
        >>> print '%.3f' % d.getCDF(1.0)
        0.500

        """
        return BetaPrimeDistBase.getCDF(self, x)
        
    def getQuantile(self, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Given p, the integral of the probability density function from 0.0
        up to a value x, returns x.

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(1, 1)
        >>> print '%.3f' % d.getQuantile(0.5)
        1.000

        """
        return BetaPrimeDistBase.getQuantile(self, p)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Draws a single sampled value from the Beta Prime distribution 
        specified by this BetaPrime object. Python list comprehensions can be 
        used to store many simulated samples for use in subsequent 
        calculations.

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> d.setSeed(97531)
        >>> print round(d.sample(), 8)
        5.92475202
        >>> print round(d.sample(), 8)
        1.46101652
        >>> print round(d.sample(), 8)
        0.59441251
        
        """
        return BetaPrimeDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Evaluates the probability density function at the supplied value x.
        Returns the natural logarithm of the density at x. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print '%.6f' % d.getLnPDF(0.35)
        -0.065438

        The BetaPrime density is

             Gamma(a+b) p^(a-1) 
        ------------------------------
        Gamma(a) Gamma(b) (p+1)^(a+b)
        
        where a = 2, b = 3, p = 0.35, and Gamma is the Gamma function, not
        the Gamma probability distribution. The density is thus

           Gamma(5) (0.35)^(2-1)       4! (0.35)
        -------------------------- = -------------- = 0.9366567
        Gamma(2) Gamma(3) (1.35)^5   1! 2! (1.35)^5

        and the log density is thus -0.06543843696237
        
        """
        return BetaPrimeDistBase.getLnPDF(self, x)
        
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
        >>> d = BetaPrime(2, 3)
        >>> print '%.6f' % d.getRelativeLnPDF(0.35)
        -2.250240

        The Beta Prime density is

             Gamma(a+b) p^(a-1) 
        ------------------------------
        Gamma(a) Gamma(b) (p+1)^(a+b)
        
        where a = 2, b = 3, p = 0.35, and Gamma is the Gamma function, not
        the Gamma probability distribution. The relative density requires
        only the two terms containing p in the numerator, so the natural
        logarithm of the relative density is

        (a-1) log(p) - (a+b) log(1+p)
        
        For the example given, this equals (2-1) log(0.35) + (2+3-1) log(1.35)
        = log(0.35) + 4*log(1.35) = -2.250240494300031
        """
        return BetaPrimeDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the mean and variance of this Beta Prime distribution. The values
        of parameters a and b can be determined from mean and var as follows:
           a = (phi + 1)*mean
           b = phi + 2
           phi = mean*(1 + mean)/var
        So if mean = 1, var = 2, we have:
           phi = (1)(2)/2 = 1
           a = 2
           b = 3
        Check:
           mean = a/(b-1)
                = 2/2
                = 1
           var  = a(a+b-1)/[(b-2)*(b-1)^2]
                = 2*4/[1*4]
                = 2

        >>> from phycas.ProbDist import *
        >>> d = BetaPrime(2, 3)
        >>> print d.getMean()
        1.0
        >>> print d.getVar()
        2.0
        >>> d.setMeanAndVariance(1.0, 2.0)
        >>> print d.getMean()
        1.0
        >>> print d.getVar()
        2.0
        
        """
        return BetaPrimeDistBase.setMeanAndVariance(self, mean, var)
    

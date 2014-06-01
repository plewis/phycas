from _PyDistributionBase import PyDistributionBase
from _ProbDistExt import *

class ImproperUniform(ImproperUniformDistBase, PyDistributionBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Represents the univariate improper uniform probability distribution.
    It specifies a distribution that is uniform from 0.0 to infinity.
    This is an improper distribution because the area under the density
    function is infinite (not unity, which is the case for a proper
    probability distribution.

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates an improper uniform distribution. This distribution has no
        parameters, so use an empty set of parentheses after the name when
        creating one. 

        >>> from phycas.ProbDist import *
        >>> d = ImproperUniform()
        >>> print d.getRelativeLnPDF(1.5)
        0.0

        """
        ImproperUniformDistBase.__init__(self)
        
    def clone(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a copy of this Improper Uniform distribution.
        
        """
        return ImproperUniformDistBase.clone(self)
        
    def isDiscrete(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Always returns False because the improper uniform distribution is
        continuous.
        
        """
        return ImproperUniformDistBase.isDiscrete(self)
        
    def getDistName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the string 'Improper Uniform'
        
        """
        return ImproperUniformDistBase.getDistName(self)
        
    def __str__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        ImproperUniform object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = ImproperUniform()
        >>> print d.__str__()
        ImproperUniform()
        
        """
        return ImproperUniformDistBase.__str__(self)
        
    def __repr__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a string that could be used to initialize another
        ImproperUniform object identical to this one. e.g.,

        >>> from phycas.ProbDist import *
        >>> d = ImproperUniform()
        >>> print d.__repr__()
        ImproperUniform()
        
        """
        return ImproperUniformDistBase.__repr__(self)
        
    def setLot(self, lot):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because sampling from this
        distribution is undefined.
        
        """
        return ImproperUniformDistBase.setLot(self, lot)
        
    def setSeed(self, seed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because sampling from this
        distribution is undefined.

        """
        return ImproperUniformDistBase.setSeed(self, seed)
        
    def resetLot(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because sampling from this
        distribution is undefined.

        """
        return ImproperUniformDistBase.resetLot(self)
        
    def getMean(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined mean.

        """
        return ImproperUniformDistBase.getMean(self)
        
    def getVar(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined variance.

        """
        return ImproperUniformDistBase.getVar(self)
        
    def getStdDev(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined standard deviation.

        """
        return ImproperUniformDistBase.getStdDev(self)
        
    def getCDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined cumulative density function.

        """
        return ImproperUniformDistBase.getCDF(self, x)
        
    def sample(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because sampling from this
        distribution is undefined.
        
        """
        return ImproperUniformDistBase.sample(self)
        
    def getLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined probability density function.

        """
        return ImproperUniformDistBase.getLnPDF(self, x)
        
    def getRelativeLnPDF(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the value 0.0 for every supplied value of x. The relative
        probability density can be any constant; the value 1.0 is chosen here,
        the logarithm of which is 0.0.

        >>> from phycas.ProbDist import *
        >>> d = ImproperUniform()
        >>> print d.getRelativeLnPDF(1.5)
        0.0
        >>> print d.getRelativeLnPDF(2.5)
        0.0

        """
        return ImproperUniformDistBase.getRelativeLnPDF(self, x)
        
    def setMeanAndVariance(self, mean, var):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is defined to be parallel with the other probability
        distributions, but using it makes no sense because this distribution
        has an undefined mean and variance.
        
        """
        return UniformDistBase.setMeanAndVariance(self, mean, var)
    

from _ProbDistExt import *

class Lot(LotBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Lot is the pseudorandom number generator used in Phycas. It is based
    on code written by J. Monahan, Statistics Department, North Carolina
    State University, and published in the first (1990) edition of Bruce
    Weir's Genetic Data Analysis (Sinauer, Sunderland, Massachusetts).
    Based on Schrage. 1979. ACM Trans. Math. Software 5:132-138.
    Translated to C++ by Paul O. Lewis, Dec. 10, 1992.

    This class was called Lot because the noun lot is defined as "an
    object used in deciding something by chance" according to The New
    Merriam-Webster Dictionary.

    """
    def __init__(self, seed = 0):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Lot object initialized using the specified psuedorandom
        number seed. If seed is the default value (0), the actual seed used
        is taken from the system clock.

        >>> from phycas.probdist import *
        >>> r1 = Lot(3157)
        >>> print r1.getInitSeed()
        3157
        >>> r2 = Lot() # this would seed r2 using the system clock

        """
        LotBase.__init__(self, seed)

    def setSeed(self, rnseed):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the pseudorandom number seed. Useful if you want to regenerate
        a sequence of pseudorandom numbers identical to a sequence that you
        have previously generated.

        >>> from phycas.probdist import *
        >>> lot = Lot()
        >>> lot.setSeed(1357)
        >>> lot.getSeed()
        1357

        """
        return LotBase.setSeed(self, rnseed)

    def getSeed(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current pseudorandom number seed. Call this function
        and store the result before generating a sequence of pseudorandom
        numbers if you want to later regenerate the same sequence.

        """
        return LotBase.getSeed(self)

    def getInitSeed(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the pseudorandom number seed set upon construction of this
        Lot object using the system clock or, if setSeed has been called,
        returns the seed specified in the most recent call to setSeed.

        >>> from phycas.probdist import *
        >>> lot = Lot()
        >>> lot.setSeed(1357)
        >>> lot.getInitSeed()
        1357

        """
        return LotBase.getInitSeed(self)

    def uniform(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a pseudorandom number between 0.0 and 1.0.

        """
        #from phycas import Phycas
        return LotBase.uniform(self)

    def sampleUInt(self, upper_bound):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a pseudorandom integer i such that 0 <= i < upper_bound.

        """
        return LotBase.sampleUInt(self, upper_bound)


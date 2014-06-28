from _LikelihoodExt import *

class IrreversibleModel(IrreversibleModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a simple irreversible model, which assumes
    tree is rooted and root state is either absent (0) or present (1).
    Assumes by default that 1 is the root state, but this can be 
    changed via the setRootState function.
    
    """
    def getScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns current value of scaling parameter. The branch
        lengths are multiplied by this scaling factor to allow the 
        irreversible character to evolve at a different rate than the
        characters used to generate branch lengths.
        
        """
        return IrreversibleModelBase.getScalingFactor(self)
    
    def setScalingFactor(self, sf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets value of scaling parameter to supplied value sf. The branch
        lengths are multiplied by this scaling factor to allow the 
        irreversible character to evolve at a different rate than the
        characters used to generate branch lengths.
        
        """
        IrreversibleModelBase.setScalingFactor(self, sf)
    
    def fixScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fixes value of the scaling_factor to its current value.
        
        """
        IrreversibleModelBase.fixScalingFactor(self)
    
    def freeScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Frees scaling factor so that it can be modified during MCMC analyses.
        
        """
        IrreversibleModelBase.freeScalingFactor(self)

    def setGainOnly(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Turns the model into a model in which only gains are allowed.
            
            """
        IrreversibleModelBase.setGainOnly(self)
    
    def setLossOnly(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Turns the model into a model in which only losses are allowed.
            
            """
        IrreversibleModelBase.setLossOnly(self)
    
    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Returns the number of states (always 2 for this model).
            
            """
        return IrreversibleModelBase.getNStates(self)
    
    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple comprising the 2 state frequencies. Either 
        (1.0,0.0) or (0.0,1.0) for this model depending on whether 0 or 1
        is the root state, respectively.
        
        >>> import phycas.likelihood
        >>> model = phycas.likelihood.IrreversibleModel()
        >>> print model.getStateFreqs()
        (0.0,1.0)
        
        """
        return IrreversibleModelBase.getStateFreqs(self)
    
#    def setAllFreqsEqual(self):
#        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
#        """
#        Frequencies cannot be set equal for this model, so this will generate
#        and exception.
#        
#        """
#        return IrreversibleModelBase.setAllFreqsEqual(self)
#    
    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of relative rate categories.
        
        """
        return IrreversibleModelBase.getNGammaRates(self)
    
    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        IrreversibleModelBase.setNGammaRates(self, n)
    
    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return IrreversibleModelBase.getRateProbs(self)
    
    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        IrreversibleModelBase.setAllRateProbsEqual(self)
    
class BinaryModel(BinaryModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a simple 2-state Markov model with scaling_factor
    parameter, which whem multiplied by the edge length determines
    the forward rate (state 0 -> 1), and a forward/reverse rate ratio,
    rho, which when multiplied by both scaling_factor and the edge length 
    determines the reverse rate.
    
    """
    def getScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns current value of the scaling_factor parameter. The branch
        lengths are multiplied by this scaling factor to allow the 
        binary character to evolve at a different rate than the
        characters used to generate branch lengths.
        
        """
        return BinaryModelBase.getScalingFactor(self)
    
    def setScalingFactor(self, sf):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets value of scaling_factor parameter to the supplied value sf. The
        edge lengths are multiplied by this scaling factor to allow the 
        binary character to evolve at a different rate than the
        characters used to generate branch lengths.
        
        """
        BinaryModelBase.setScalingFactor(self, sf)
    
    def fixScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fixes value of the scaling_factor to its current value.
        
        """
        BinaryModelBase.fixScalingFactor(self)
    
    def freeScalingFactor(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Frees scaling_factor so that it can be modified during MCMC
        analyses.
        
        """
        BinaryModelBase.freeScalingFactor(self)
    
    def getKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns current value of the forward/reverse rate ratio, kappa.        
        """
        return BinaryModelBase.getKappa(self)

    def setKappa(self, k):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets value of forward/reverse rate ratio, kappa, to the supplied value 
        k.
        
        """
        BinaryModelBase.setKappa(self, k)

    def fixKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fixes value of the forward/reverse rate ratio kappa to its current
        value.
        
        """
        BinaryModelBase.fixKappa(self)

    def freeKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Frees forward/reverse rate ratio kappa so that it can be modified during
        MCMC analyses.
        
        """
        BinaryModelBase.freeKappa(self)

    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of states (always 2 for this model).
        
        """
        return BinaryModelBase.getNStates(self)
    
    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple comprising the 2 state frequencies. 
        
        >>> import phycas.likelihood
        >>> model = phycas.likelihood.BinaryModel()
        >>> print model.getStateFreqs()
        (0.5,0.5)
        
        """
        return BinaryModelBase.getStateFreqs(self)
    
    def setAllFreqsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets both frequencies to 0.5.
        
        """
        return BinaryModelBase.setAllFreqsEqual(self)
    
    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of relative rate categories.
        
        """
        return BinaryModelBase.getNGammaRates(self)
    
    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        BinaryModelBase.setNGammaRates(self, n)
    
    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return BinaryModelBase.getRateProbs(self)
    
    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        BinaryModelBase.setAllRateProbsEqual(self)
    
class JCModel(JCModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
        Encapsulates the Jukes and Cantor (1969) substitution model, which
        assumes base frequencies are equal and all types of substitutions
        occur at the same rate.
        
        Literature Cited:
        
        Jukes, T. H., and C. R. Cantor. 1969. Evolution of protein molecules.
        Pages 21-132 in Mammalian Protein Metabolism (H. N. Munro, ed.)
        Academic Press, New York.    
        
        """
    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Returns the number of states (always 4 for this model).
            
            """
        return JCModelBase.getNStates(self)
    
    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Returns a tuple comprising the 4 state frequencies. Always
            (0.25, 0.25, 0.25, 0.25) for this model.
            
            >>> import phycas.likelihood
            >>> model = phycas.likelihood.JCModel()
            >>> print model.getStateFreqs()
            (0.25, 0.25, 0.25, 0.25)
            
            """
        return JCModelBase.getStateFreqs(self)
    
    def setAllFreqsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Sets all four state frequencies to 0.25. Superfluous for this model,
            but included for conformity
            
            """
        return JCModelBase.setAllFreqsEqual(self)
    
    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Returns the current number of relative rate categories.
            
            """
        return JCModelBase.getNGammaRates(self)
    
    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Sets the number of relative rate categories to n (n should be greater
            than zero).
            
            """
        return JCModelBase.setNGammaRates(self, n)
    
    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Returns a list each element of which is the probability that any given
            site falls in its particular rate category.
            
            """
        return JCModelBase.getRateProbs(self)
    
    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
            Sets all rate probabilities to the inverse of the number of rate
            categories.
            
            """
        return JCModelBase.setAllRateProbsEqual(self)
    
class HKYModel(HKYModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates the Hasegawa-Kishino-Yano (1985) substitution model,
    which allows unequal base frequencies and transition-type
    substitutions to occur at a different rate than transversion-type
    substitutions. The transition/transversion rate ratio is termed kappa,
    whereas the probability of any transition divided by the probability
    of any transversion is known as tratio. This model can be made equal
    to the Felsenstein (1981) model by setting kappa to 1.0, and to the
    Kimura (1980) 2-parameter model by making the base frequencies all
    0.25. This model is identical to the Jukes-Cantor (1969) model if
    the base frequencies are equal and kappa is 1.0.

    Literature Cited:
    
    Felsenstein, J. 1981. Evolutionary trees from DNA sequences:  a
    maximum likelihood approach. Journal of Molecular Evolution
    17: 368-376.

    Hasegawa, M., H. Kishino, and T. Yano. 1985. Dating of the human-ape
    splitting by a molecular clock of mitochondrial DNA. Journal of
    Molecular Evolution 22: 160-174.

    Jukes, T. H., and C. R. Cantor. 1969. Evolution of protein molecules.
    Pages 21-132 in Mammalian Protein Metabolism (H. N. Munro, ed.)
    Academic Press, New York.

    Kimura, M. 1980. A simple method for estimating evolutionary rate of
    base substitutions through comparative studies of nucleotide sequences.
    Journal of Molecular Evolution 16:111-120.    
    
    """
    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of states (always 4 for this model).
        
        """
        return HKYModelBase.getNStates(self)

    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple comprising the 4 state frequencies.
        
        """
        return HKYModelBase.getStateFreqs(self)

    def setStateFreqParam(self, i, value):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets frequency parameter for state i to value. use i = 0 for A,
        i = 1 for C, i = 2 for G and i = 3 for T. Note that value can be any
        non-negative number; there is no need to ensure that it is between
        0.0 and 1.0 (although there is nothing wrong with providing normalized
        frequencies). The four frequency parameters are normalized for use in
        all calculations involving base frequencies. Thus, specifying 1, 2, 3,
        and 4 for the four frequency parameters will result in the relative
        base frequencies being set to 0.1, 0.2, 0.3 and 0.4.
        
        """
        return HKYModelBase.setStateFreqUnnorm(self, i, value)

    def setAllFreqsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all four state frequencies to 0.25.
        
        """
        return HKYModelBase.setAllFreqsEqual(self)

    def setNucleotideFreqs(self, freqA, freqC, freqG, freqT):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the four state frequencies to the values provided, which should
        all be greater than or equal to 0.0.
        
        """
        assert freqA >= 0.0 and freqC >= 0.0 and freqG >= 0.0 and freqT >= 0.0
        return HKYModelBase.setNucleotideFreqs(self, freqA, freqC, freqG, freqT)

    def getKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current value of the transition-transversion rate
        ratio kappa.
        
        """
        return HKYModelBase.getKappa(self)

    def setKappa(self, k):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value of the transition-transversion rate ratio kappa.
        
        """
        return HKYModelBase.setKappa(self, k)

    def setKappaFromTRatio(self, tratio):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value of the transition-transversion rate ratio kappa given
        a supplied value of tratio, the transition-transversion ratio. The
        rate ratio kappa is related to tratio as follows (where piA means the
        frequency of base A, piC means the frequency of base C, etc.

                tratio (piA + piG) (piC + piT)
        kappa = -------------------------------
                    (piA piG + piC piT)

        Thus, if piA = piC = piG = piT = 0.25, kappa is twice the tratio
        because there are twice as many kinds of transversion-type
        substitutions compared to transition-type substitutions.
        """
        return HKYModelBase.setKappaFromTRatio(self, tratio)

    def calcTRatio(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculates the transition/transversion ratio (tratio) given the
        transition/transversion rate ratio (kappa) and the relative base
        frequencies. Here are the details of the calculation (for brevity,
        piA symbolizes the frquency of base A, piC the frequency of base C,
        etc.
        
        Parameters: b = transversion rate, k = kappa, dt = infinitesimal time
        The probability that a transition from base A to base G occurs over
        time dt is the probability of starting in state A (piA) times the
        probability of a transition to base G (piG k b dt).
        
        Pr(any transition | dt) = (piA piG k b dt) + (piC piT k b dt)
          + (piG piA k b dt) + (piT piC k b dt)
          = 2 k b dt (piA piG + piC piT)
        
        Pr(any transversion | dt) = (piA piC b dt) + (piA piT b dt)
          + (piC piA b dt) + (piC piG b dt) + (piG piC b dt) + (piG piT b dt)
          + (piT piA b dt) + (piT piG b dt)
          = 2 b dt (piA + piG) (piC + piT)

                 2 k b dt (piA piG + piC piT)     k (piA piG + piC piT)
        tratio = ------------------------------ = -----------------------
                 2 b dt (piA + piG) (piC + piT)   (piA + piG) (piC + piT)
        
        """
        return HKYModelBase.calcTRatio(self)

    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of relative rate categories.
        
        """
        return HKYModelBase.getNGammaRates(self)

    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return HKYModelBase.setNGammaRates(self, n)

    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return HKYModelBase.getRateProbs(self)

    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return HKYModelBase.setAllRateProbsEqual(self)

class GTRModel(GTRModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates the General Time Reversible substitution model, first
    described by Rodriguez et al. (1990). This model allows unequal base
    frequencies as well as six different relative rates corresponding to
    the substitution classes A <-> C, A <-> G, A <-> T, C <-> G, C <-> T
    and G <-> T. Constrained versions of this model can be made equal to
    the HKY85, K80, F81 and JC models.

    Literature Cited:

    Rodrigues, F., J. L. Oliver, A. Marin, and J. R. Medina. 1990. The
    general stochastic model of nucleotide substitution. Journal of
    Theoretical Biology 142: 485-501.
    
    """
    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of states (always 4 for this model).
        
        """
        return GTRModelBase.getNStates(self)

    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple comprising the 4 state frequencies.
        
        """
        return GTRModelBase.getStateFreqs(self)

    def setStateFreqParam(self, i, value):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets frequency parameter for state i to value. use i = 0 for A,
        i = 1 for C, i = 2 for G and i = 3 for T. Note that value can be any
        non-negative number; there is no need to ensure that it is between
        0.0 and 1.0 (although there is nothing wrong with providing normalized
        frequencies). The four frequency parameters are normalized for use in
        all calculations involving base frequencies. Thus, specifying 1, 2, 3,
        and 4 for the four frequency parameters will result in the relative
        base frequencies being set to 0.1, 0.2, 0.3 and 0.4.
        
        """
        return GTRModelBase.setStateFreqUnnorm(self, i, value)

    def setAllFreqsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all four state frequencies to 0.25.
        
        """
        return GTRModelBase.setAllFreqsEqual(self)

    def setNucleotideFreqs(self, freqA, freqC, freqG, freqT):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the four state frequencies to the values provided, which should
        all be greater than or equal to 0.0.
        
        """
        assert freqA >= 0.0 and freqC >= 0.0 and freqG >= 0.0 and freqT >= 0.0
        return GTRModelBase.setNucleotideFreqs(self, freqA, freqC, freqG, freqT)

    def getRelRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list comprising the six relative rates.
        
        """
        return GTRModelBase.getRelRates(self)

    def setRelRates(self, rr):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all six relative rates. Supply the rates in a list or tuple of
        length 6. The rates should be ordered as follows: rAC, rAG, rAT, rCG,
        rCT, rGT. All rates should be greater than zero. They do not need to
        be normalized in any way.
        
        """
        return GTRModelBase.setRelRates(self, rr)

    def calcTRatio(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calculates the transition/transversion ratio (tratio) given the six
        relative rates and the relative base frequencies. Here are the
        details of the calculation.

        Symbols:
            rAC = rate of A <-> C   piA = frequency of base A
            rAG = rate of A <-> G   piC = frequency of base C
            rAT = rate of A <-> T   piG = frequency of base G
            rCG = rate of C <-> G   piT = frequency of base T
            rCT = rate of C <-> T
            rGT = rate of G <-> T   dt = infinitesimal time
            
        The probability that a transition from base A to base G occurs over
        time dt is the probability of starting in state A (piA) times the
        probability of a transition to base G (piG rAG dt).
        
        Pr(any transition | dt) = (piA piG rAG dt) + (piC piT rCT dt)
          + (piG piA rAG dt) + (piT piC rCT dt)
          = 2 dt (piA piG rAG + piC piT rCT)
        
        Pr(any transversion | dt) = (piA piC rAC dt) + (piA piT rAT dt)
          + (piC piA rAC dt) + (piC piG rCG dt) + (piG piC rCG dt)
          + (piG piT rGT dt) + (piT piA rAT dt) + (piT piG rGT dt)
          = 2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)

                            2 dt (piA piG rAG + piC piT rCT)
        tratio = ------------------------------------------------------------
                 2 dt (piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT)
        
                            piA piG rAG + piC piT rCT
               = -----------------------------------------------------
                 piA piC rAC + piA piT rAT + piC piG rCG + piG piT rGT

        Example:
        >>> from phycas import *
        >>> model = likelihood.GTRModel()
        >>> model.setRelRates([1.0, 4.0, 1.0, 1.0, 4.0, 1.0])
        >>> model.setNucleotideFreqs(0.25, 0.25, 0.25, 0.25)
        >>> print model.calcTRatio()
        2.0
        
        """
        return GTRModelBase.calcTRatio(self)

    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of relative rate categories.
        
        """
        return GTRModelBase.getNGammaRates(self)

    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return GTRModelBase.setNGammaRates(self, n)

    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list each element of which is the probability that any given
        site falls in its particular rate category.
        
        """
        return GTRModelBase.getRateProbs(self)

    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return GTRModelBase.setAllRateProbsEqual(self)

class CodonModel(CodonModelBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapsulates a codon-based substitution model. This model estimates
    codon frequencies as well as the transition/transversion rate ratio
    (kappa) and the nonsynonymous/synonymous rate ratio (omega).
    
    """
    def getModelName(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the name of this model.
        
        """
        return CodonModelBase.getModelName(self)

    def getNStates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of states, always 61 for this model (64 possible
        triplets minus the three stop codons, TAA, TAG and TGA).
        
        """
        return CodonModelBase.getNStates(self)

    def getStateFreqs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a tuple comprising the 61 state frequencies (64 possible
        triplets minus the three stop codons, TAA, TAG and TGA).
        
        """
        return CodonModelBase.getStateFreqs(self)

    def setStateFreqsUnnorm(self, values):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        For all i, sets the frequency parameter for state i to values[i]. The
        codon states are assumed to be in this order: AAA, AAC, AAG, AAT, 
        ACA, ACC, ACG, ACT, ..., TTT (with the exception of the three stop
        codons, TAA, TAG and TGA, which are never considered. Note that value
        can be any non-negative number; there is no need to ensure that it is
        between 0.0 and 1.0 (although there is nothing wrong with providing 
        normalized frequencies). The 61 frequency parameters are normalized
        for use in all calculations involving state frequencies. Thus, 
        specifying 10 for each of the 61 state frequency parameters will
        result in the relative state frequencies being all set to 1/61.
        
        """
        CodonModelBase.setStateFreqsUnnorm(self, values)

    def setStateFreqUnnorm(self, i, value):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets frequency parameter for state i to value. The codon states are
        in this order: AAA, AAC, AAG, AAT, ACA, ACC, ACG, ACT, ..., TTT (with
        the exception of the three stop codons, TAA, TAG and TGA, which are
        never considered. Note that value can be any non-negative number;
        there is no need to ensure that it is between 0.0 and 1.0 (although
        there is nothing wrong with providing normalized frequencies). The
        61 frequency parameters are normalized for use in all calculations
        involving state frequencies. Thus, specifying 10 for each of the 61
        state frequency parameters will result in the relative state
        frequencies being all set to 1/61.
        
        """
        return CodonModelBase.setStateFreqUnnorm(self, i, value)

    def setAllFreqsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all 61 state frequencies to 1/61 = 0.01639.
        
        """
        return CodonModelBase.setAllFreqsEqual(self)

    def setNucleotideFreqs(self, freqA, freqC, freqG, freqT):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the 61 codon state frequencies to the values expected based on
        the four base frequencies provided (all four values provided should
        greater than or equal to 0.0, but do not need to sum to 1.0). The
        frequency of a codon is, almost, the product of the three component
        base frequencies: the "almost" qualification being needed because the
        three stop codons are not included, so each codon frequency must be
        corrected by dividing by the sum of the 61 non-stop three-nucleotide
        products. For example, if the specified base freqencies were 0.1, 0.2,
        0.3 and 0.4, then the frequency of the ACT codon would be the product
        (0.1)*(0.2)*(0.4) = 0.008, divided by the sum of all 61 such products,
        which in this case is 0.972, yielding 0.008/0.971 = 0.00823. Note
        state frequencies set using this function will be obliterated unless
        state frequencies are fixed using the fixStateFreqs method.
        
        """
        assert freqA >= 0.0 and freqC >= 0.0 and freqG >= 0.0 and freqT >= 0.0, 'supplied nucleotide frequency parameters must all be > 0'
        return CodonModelBase.setNucleotideFreqs(self, freqA, freqC, freqG, freqT)

    def fixOmega(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fixes the value of omega, the nonsynonymous/synonymous rate ratio, to
        prevent it from being updated or estimated in the future.
        
        """
        return CodonModelBase.fixOmega(self)

    def freeOmega(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Frees the value of omega, the nonsynonymous/synonymous rate ratio,
        allowing it to be updated or estimated in the future.
        
        """
        return CodonModelBase.freeOmega(self)

    def getOmega(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current value of omega, the nonsynonymous/synonymous rate
        ratio.
        
        """
        return CodonModelBase.getOmega(self)

    def setOmega(self, w):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value of omega, the nonsynonymous/synonymous rate ratio.
        
        """
        return CodonModelBase.setOmega(self, w)

    def getOmegaPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the prior distribution of omega, the nonsynonymous/synonymous
        rate ratio.
        
        """
        return CodonModelBase.getOmegaPrior(self)

    def setOmegaPrior(self, w):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the prior distribution for omega, the nonsynonymous/synonymous
        rate ratio.
        
        """
        CodonModelBase.setOmegaPrior(self, w)

    def fixKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fixes the value of kappa, the transition/transversion rate ratio, to
        prevent it from being updated or estimated in the future.
        
        """
        return CodonModelBase.fixKappa(self)

    def freeKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Frees the value of kappa, the transition/transversion rate ratio,
        allowing it to be updated or estimated in the future.
        
        """
        return CodonModelBase.freeKappa(self)

    def getKappa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current value of kappa, the transition/transversion rate
        ratio.
        
        """
        return CodonModelBase.getKappa(self)

    def setKappa(self, k):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value of kappa, the transition/transversion rate ratio.
        
        """
        CodonModelBase.setKappa(self, k)

    def getKappaPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the prior distribution for kappa, the transition/transversion
        rate ratio.
        
        """
        return CodonModelBase.getKappaPrior(self)

    def setKappaPrior(self, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the prior distribution for kappa, the transition/transversion
        rate ratio.
        
        """
        CodonModelBase.setKappaPrior(self, p)

    def getStateFreqParamPrior(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current prior distribution for the parameters governing
        state frequencies.
        
        """
        return CodonModelBase.getStateFreqParamPrior(self)

    def setStateFreqParamPrior(self, p):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets prior distribution for the parameters governing state
        frequencies. These parameters represent unnormalized state
        frequencies, so a Gamma distribution is appropriate.
        
        """
        CodonModelBase.setStateFreqParamPrior(self, p)

    def getNGammaRates(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of relative rate categories.
        
        """
        return CodonModelBase.getNGammaRates(self)

    def setNGammaRates(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of relative rate categories to n (n should be greater
        than zero).
        
        """
        return CodonModelBase.setNGammaRates(self, n)

    def getRateProbs(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list in which each element is the probability that any given
        site falls in its particular rate category.
        
        """
        return CodonModelBase.getRateProbs(self)

    def setAllRateProbsEqual(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all rate probabilities to the inverse of the number of rate
        categories.
        
        """
        return CodonModelBase.setAllRateProbsEqual(self)


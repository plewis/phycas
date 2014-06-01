#TODO  Is this being used???

from _LikelihoodExt import *

class MCMCChainManager(MCMCChainManagerBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Sorry, no documentation yet.
    
    """
    def __init__(self, jpm):
        """
        Sorry, no documentation yet.
        
        """
        MCMCChainManagerBase.__init__(self, jpm)
        
    def finalize(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.finalize(self)

    def addMove(self, move):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addMove(self, move)
    
    def addModelParam(self, param):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addModelParam(self, param)
    
    def addEdgeLenParam(self, param):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addEdgeLenParam(self, param)
    
    def addEdgeLenHyperparam(self, param):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.addEdgeLenHyperparam(self, param)
    
    def getLastLnLike(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getLastLnLike(self)

    def getMoves(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getMoves(self)

    def getModelParams(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getModelParams(self)

    def getEdgeLenParams(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getEdgeLenParams(self)

    def getEdgeLenHyperparams(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getEdgeLenHyperparams(self)

    def clear(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.clear(self)

    def getAllUpdaters(self):
        """
        Sorry, no documentation yet.
        
        """
        return MCMCChainManagerBase.getAllUpdaters(self)


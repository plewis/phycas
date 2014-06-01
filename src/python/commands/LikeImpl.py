import os,sys,math,random
from phycas import *
from MCMCManager import LikelihoodCore
from phycas.Utilities.PhycasCommand import *
from phycas.ReadNexus import NexusReader
from phycas.Utilities.CommonFunctions import CommonFunctions

class LikeImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    To be written.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the LikeImpl object by assigning supplied phycas object
        to a data member variable.
        
        """
        CommonFunctions.__init__(self, opts)
        self.starting_tree         = None
        self.taxon_labels          = None
        self.data_matrix           = None
        self.ntax                  = None
        self.nchar                 = None
        self.reader                = NexusReader()
        self.npatterns          = []        # Will hold the actual number of patterns for each subset after data file has been read

    def _loadData(self, matrix):
        self.data_matrix = matrix
        if matrix is None:            
            self.taxon_labels = []
            self.ntax = 0
            self.nchar = 0 # used for Gelfand-Ghosh simulations only
        else:
            self.taxon_labels = matrix.taxa
            self.ntax = self.data_matrix.getNTax()
            self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")


    def getStartingTree(self):
        if self.starting_tree is None:
            try:
                tr_source = self.opts.tree_source
                tr_source.setActiveTaxonLabels(self.taxon_labels)
                i = iter(tr_source)
                self.starting_tree = i.next()
            except:
                self.stdout.error("A tree could not be obtained from the tree_source")
                raise
        return self.starting_tree
        
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current tree and current
        model.
        
        """
        
        ds = self.opts.data_source
        mat = ds and ds.getMatrix() or None
        self.phycassert(self.opts.data_source is not None, "specify data_source before calling like()")
        self._loadData(mat)
        
        self.starting_tree =  self.getStartingTree()
        if self.opts.preorder_edgelens is not None:
            self.starting_tree.replaceEdgeLens(self.opts.preorder_edgelens)
            print '@@@@@@@@@@ self.starting_tree.makeNewick() =',self.starting_tree.makeNewick()
        core = LikelihoodCore(self)
        core.setupCore()
        core.prepareForLikelihood()
        if self.opts.store_site_likes:
            core.likelihood.storeSiteLikelihoods(True)
            self.opts.pattern_counts = None
            self.opts.char_to_pattern = None
            self.opts.site_likes = None
            self.opts.site_uf = None
        else:
            core.likelihood.storeSiteLikelihoods(False)
        lnL = core.calcLnLikelihood()
        if self.opts.store_site_likes:
            self.opts.pattern_counts = core.likelihood.getPatternCounts()
            self.opts.char_to_pattern = core.likelihood.getCharIndexToPatternIndex()
            self.opts.site_likes = core.likelihood.getSiteLikelihoods()
            self.opts.site_uf = core.likelihood.getSiteUF()
        return lnL

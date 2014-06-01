import os,sys,math,random
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class CPOImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Needs to be written.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes MCMCImpl object by assigning supplied phycas object
        to a data member variable.
        
        """
        CommonFunctions.__init__(self, opts)
        self.opts = opts
        self.sitelikef = None

    def siteLikeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the site log-likelihood file.
        
        """
        self.phycassert(self.sitelikef is None, 'Attempt made to open CPOImpl.sitelikef, but it is already open!')
        sitelnl_file_spec = self.opts.out.sitelike
        try:
            self.sitelikef = sitelnl_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open site log-likelihood file (%s) failed.' % self.opts.out.sitelike.filename

        if self.sitelikef:
            print 'Site log-likelihood file was opened successfully'
            #mcmc.sitelikef = self.sitelikef

    def siteLikeFileClose(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the site log-likelihood file.
        
        """
        self.phycassert(self.sitelikef is not None, 'Attempt made to close CPOImpl.sitelikef, but it is not open!')
        self.sitelikef.close()
        self.sitelikef = None

        
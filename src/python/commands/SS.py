from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas import mcmc,partition,model

class SS(PhycasCommand):
    def __init__(self): 
        args = (   ("nstones",          20,     "The number of stepping stones (ratios) used to estimate the marginal likelihood; for example, if this value is 5, then beta will take on these values: 0.8, 0.6, 0.4, 0.2, 0.0", IntArgValidate(min=1)),
                   ("shape1",           1.0,    "The first shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from 0 to 1.", FloatArgValidate(greaterthan=0.0)),
                   ("shape2",           1.0,    "The second shape parameter of the distribution used to determine the beta values to be sampled. This distribution is, confusingly, a Beta distribution. Thus, if both shape1 and shape2 are set to 1, beta values will be chosen at uniform intervals from 0 to 1.", FloatArgValidate(greaterthan=0.0)),
                   ("ncycles",          500,    "The number of sampling cycles per stone (ratio).", IntArgValidate()),
                   ("sample_every",       1,    "The current tree topology and model parameter values will be sampled after this many cycles have elapsed since the last sample was taken", IntArgValidate(min=0)),
                   ("report_every",     100,    "A progress report will be displayed after this many cycles have elapsed since the last progress report", IntArgValidate(min=0)),
                   ("refdist_is_prior", False,  "If True, the prior will be used as the reference distribution; if False, the reference distribution definition will be that specified in refdistfile", BoolArgValidate),
                   ("refdistfile",      None,   "The file containing a reference distribution definition. If None, the stepping-stone analysis will use the prior as the reference distribution.", FileExistsOrNoneValidate)
                )
                   #("xcycles", 0, "The number of extra cycles (above and beyond mcmc.ncycles) that will be spent exploring the posterior (additional posterior cycles help stepping stone analyses formulate an effective reference distribution).", IntArgValidate()),
                   #("maxbeta", 1.0, "The first beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   #("minbeta", 0.0, "The last beta value that will be sampled.", FloatArgValidate(min=0.0, max=1.0)),
                   #("minsample", 10, "Minimum sample size needed to create a split-specific edge length reference distribution.", IntArgValidate(min=0)),
                   #("ti", False, "If True, the marginal likelihood will be estimated using thermodynamic integration and the stepping stone method with reference distribution equal to the prior; if False (the default), the generalized stepping stone method with reference distribution approximating the posterior will be used (this greatly improves the accuracy of the stepping stone method and is strongly recommended).", BoolArgValidate),
        # Specify output options
        #self.__dict__["hidden"] = True # hide from main phycas help list of commands until working
        #o = PhycasCommandOutputOptions()
        #o.__dict__["_help_order"] = ["sss"]
        PhycasCommand.__init__(self, args, "ss", "Performs stepping stone method for purposes of estimating the marginal likelihood of the current model.")

        # The data member below is hidden from the user because it overrides something that users should not be able to override 
        #self.__dict__["override_fixed_topology_restriction"] = False
        
        # The data members added below are hidden from the user because they are set when the mcmc command runs
        self.__dict__["sampled_likes"] = None
        self.__dict__["sampled_betas"] = None
        
        # The data members added below are hidden from the user because they are for developer use only
        self.__dict__["refdist_definition_file"] = None # specify file name of file containing reference distributions (one per line in same order that they are outut in the log file)

    def hidden():
        """ 
        Overrides the PhycasCommand.hidden method to keep SS's name from being displayed 
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False
        
    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a SS object.
        """
        cf = CommonFunctions(self)
        if self.refdist_is_prior and self.refdistfile is not None:
            cf.phycassert(False, 'cannot specify ss.refdist_is_prior = True and also specify a reference distribution using ss.refdistfile')
        if not self.refdist_is_prior and self.refdistfile is None:
            cf.phycassert(False, 'ss.refdist_is_prior is False but you have not specified a reference distribution (i.e. ss.refdistfile is None)')
        cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for the stepping-stone method')
        cf.phycassert(mcmc.allow_polytomies == False or self.ti == True, "mcmc.allow_polytomies must be False to use the generalized stepping-stone method (we're working on relaxing this requirement)")
        #if not self.override_fixed_topology_restriction:
        #    cf.phycassert(mcmc.fix_topology == True, "mcmc.fix_topology must be True to use the stepping-stone method (we're working on relaxing this requirement)")
        models = partition.getModels()
        if len(models) > 0:
            m = models[0]
        else:
            m = model 
        if m.edgelen_hyperprior:
            cf.warning("the generalized stepping-stone method has only been tested on models in which model.edgelen_hyperprior = None (proceed at your own risk)")
        
    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        prev_report_every = mcmc.report_every
        mcmc.report_every = self.report_every
        prev_sample_every = mcmc.sample_every
        mcmc.sample_every = self.sample_every
        mcmc.doing_steppingstone_sampling = True
        mcmc.ssobj = self
        mcmc.ss_heating_likelihood = False
        if self.refdist_is_prior:
            mcmc.ss_heating_likelihood = True   # performing TI and non-generalized SS
        mcmc()
        mcmc.ss_heating_likelihood = False
        mcmc.doing_steppingstone_sampling = False
        mcmc.ssobj = None
        self.sampled_betas = mcmc.ss_sampled_betas
        self.sampled_likes = mcmc.ss_sampled_likes
        mcmc.report_every = prev_report_every
        mcmc.sample_every = prev_sample_every

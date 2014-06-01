from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Commands.CPOImpl import CPOImpl
from phycas import mcmc

class CPO(PhycasCommand):
    def __init__(self):
        args = (("patterns_only", False, "If True, each row of the sitelike output file will contain sampled log-likelihoods for each pattern, with the first row of the file holding the counts for each pattern. If False, the rows of the sitelike file will contain the log-likelihoods for each site in the order in which the sites occur in the data file (generally produces a larger file).", BoolArgValidate),)
        # Specify output options
        o = PhycasCommandOutputOptions()
        o.__dict__["_help_order"] = ["sitelike"]
        p = TextOutputSpec(prefix='sitelike', suffix=".txt", help_str="The text file in which all sampled site log-likelihood values are saved.")
        o.__dict__["sitelike"] = p
        PhycasCommand.__init__(self, args, "cpo", "Performs a Conditional Predictive Ordinate (CPO) analysis to determine the relative fit of the model to individual sites/characters.", o)

        # The data members added below are hidden from the user because they are set when the mcmc command runs
        #self.__dict__["sampled_likes"] = None
        #self.__dict__["sampled_betas"] = None

        self.__dict__["sitelikef"] = None

    def hidden():
        """
        Overrides the PhycasCommand.hidden method to keep CPO's name from being displayed
        in the list of classes displayed when users type help. Delete this function, or
        change its return value to False, when it is ready to be advertised.
        """
        return False

    hidden = staticmethod(hidden)

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of a CPO object.
        """
        #cf = CommonFunctions(self)
        #cf.phycassert(mcmc.ncycles > 0, 'mcmc.ncycles cannot be less than 1 for path sampling')

    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        c = copy.deepcopy(self)
        cpo_impl = CPOImpl(c)
        cpo_impl.siteLikeFileOpen()
        if cpo_impl.sitelikef is not None:
            mcmc.sitelikef = cpo_impl.sitelikef
            mcmc.saving_sitelikes = True
            mcmc()
            mcmc.saving_sitelikes = False
            cpo_impl.siteLikeFileClose()
        else:
            print 'Could not run the cpo command because the sitelike file could not be opened'

from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions
from phycas import model, partition, randomtree, P
from phycas.commands.MCMCImpl import MCMCImpl
#from phycas.probdist import Beta, Exponential, InverseGamma
import copy

class MCMC(PhycasCommand):
    def __init__(self):
        args = tuple(
                PhycasCommand._getRNGOptions() +
                [("burnin",                   1000,    "The number of update cycles to ignore before sampling begins. If burnin=1000 and ncycles=10000, the total number of cycles will be 11000.", IntArgValidate(min=0)),
                ("ncycles",                  10000,    "The number of update cycles (a cycle is analogous to, but different than, a 'generation' in MrBayes; Phycas does in one cycle what MrBayes does in about 100 generations for a simple model such as JC)", IntArgValidate(min=0)),
                ("sample_every",               100,    "The current tree topology and model parameter values will be sampled after this many cycles have elapsed since the last sample was taken", IntArgValidate(min=0)),
                ("report_every",               100,    "A progress report will be displayed after this many cycles have elapsed since the last progress report", IntArgValidate(min=0)),
                ("verbose",                   True,    "You will get more output if True, less output if False", BoolArgValidate),
                ("quiet",                    False,    "If True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well", BoolArgValidate),
                ("model",                    model,    "Specifies the model to use. By default, uses the predefined model object. Type model.help to set the settings for this model."),
                ("partition",            partition,    "Specifies the partition to use. By default, uses the predefined partition object."),
                ("data_source",       P.characters,    "The DataSource that provides the data, if any, to be used in the MCMC analysis. Should be a DataSource object", DataSourceValidate),
                ("starting_tree_source",      None,    "A TreeCollection that will serve as the source of the starting tree topology. If a string is passed in, it is interpreted as a the path to a file with trees.", TreeSourceValidate),
                ("tree_topology",             None,    "Unused unless starting_tree_source is 'usertree', in which case this should be a standard newick string representation of the tree topology; e.g. '(A:0.01,B:0.071,(C:0.013,D:0.021):0.037)'"),
                ("fix_topology",             False,    "If True, an EdgeMove move will be substituted for the LargetSimonMove, so edge lengths will be updated by slice sampling but the topology will remain unchanged during an MCMC analysis", BoolArgValidate),
                ("ls_move_lambda",             0.2,    "Sets the minimum value of the tuning parameter for the LargetSimonMove Metropolis-Hastings move. This value corresponds to a boldness vlaue of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=0.01)),
                ("ls_move_weight",             100,    "Larget-Simon moves will be performed this many times per cycle", IntArgValidate(min=0)),
                ("ls_move_debug",            False,    "If set to true, TreeViewer will popup on each Larget-Simon move update showing edges affected by the proposed move", BoolArgValidate),
                ("edge_move_lambda",           0.2,    "Sets the minimum value of the tuning parameter for the EdgeMove Metropolis-Hastings move. This value corresponds to a boldness vlaue of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=0.01)),
                ("edge_move_weight",             0,    "Only used if fix_topology is True. Makes sense to set this to some multiple of the number of edges since each EdgeMove affects a single randomly-chosen edge ", IntArgValidate(min=0)),
                ("nchains",                      1,    "The number of Markov chains to run simultaneously. One chain serves as the cold chain from which samples are drawn, the other chains are heated to varying degrees and serve to enhance mixing in the cold chain.", IntArgValidate(min=1)),
                ("rel_rate_weight",              1,    "Updates of GTR relative rates will occur this many times per cycle if relative rates are being updated jointly", IntArgValidate(min=0)),
                ("rel_rate_psi",            2000.0,    "Sets the maximum value of the tuning parameter for the RelRatesMove Metropolis-Hastings move. This value corresponds to a boldness value of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=1.0)),
                ("subset_relrates_weight",       1,    "Updates of th vector of partition subset relative rates will occur this many times per cycle", IntArgValidate(min=0)),
                ("subset_relrates_psi",     2000.0,    "Sets the maximum value of the tuning parameter for the SubsetRelRatesMove Metropolis-Hastings move. This value corresponds to a boldness value of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=1.0)),
                ("state_freq_weight",            1,    "Updates of state frequencies will occur this many times per cycle if state frequencies are being updated jointly", IntArgValidate(min=0)),
                ("state_freq_psi",          2000.0,    "Sets the maximum value of the tuning parameter for the RelRatesMove Metropolis-Hastings move. This value corresponds to a boldness value of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=1.0)),
                ("tree_scaler_lambda",         0.5,    "Sets the minimum value of the tuning parameter for the TreeScalerMove Metropolis-Hastings move. This value corresponds to a boldness vlaue of 0.0 and is the value used for normal analyses.", FloatArgValidate(min=0.01)),
                ("tree_scaler_weight",           1,    "Whole-tree scaling will be performed this many times per cycle", IntArgValidate(min=0)),
                ("allow_polytomies",         False,    "If True, do Bush moves in addition to Larget-Simon moves; if False, do Larget-Simon moves only", BoolArgValidate),
                ("polytomy_prior",            True,    "If True, use polytomy prior; if False, use resolution class prior", BoolArgValidate),
                ("topo_prior_C",               2.0,    "Specifies the strength of the prior (C = 1 is flat prior; C > 1 favors less resolved topologies)", FloatArgValidate(min=0.01)),
                ("bush_move_edgelen_mean",     1.0,    "Specifies mean of exponential edge length generation distribution used by BushMove when new edges are created", FloatArgValidate(min=0.01)),
                ("bush_move_weight",           100,    "Bush moves will be performed this many times per cycle if", IntArgValidate(min=0)),
                ("bush_move_debug",          False,    "If set to true, TreeViewer will pop up on each Bush move update showing edges affected by the proposed move", BoolArgValidate),
                ("slice_weight",                 1,    "Slice sampled parameters will be updated this many times per cycle", IntArgValidate(min=0)),
                ("slice_max_units",           1000,    "Max. number of units used in slice sampling", IntArgValidate(min=0)),
                ("target_accept_rate",         0.3,    "An attempt will be made to achieve this acceptance rate (proportion) during the burn-in phase. Not applicable if burnin is set to 0.", FloatArgValidate(min=0,max=1)),
                ("report_efficiency_every",    100,    "Number of cycles between reports of updater efficiency.", IntArgValidate(min=0)),
                ("adapt_simple_param",         0.5,    "Slice sampler adaptation parameter", FloatArgValidate(min=0.01)),
                ("min_heat_power",             0.5,    "Power of the hottest chain when nchains > 1", FloatArgValidate(min=0.01)),
                ("heat_vector",               None,    "List of heating powers, one of which should be 1.0 (default value None causes this vector to be generated using min_heat_pwer)"),
                ("uf_num_edges",                50,    "Number of edges to traverse before taking action to prevent underflow", IntArgValidate(min=1)),
                ("ntax",                         0,    "To explore the prior, set to some positive value. Also set data_source to None", IntArgValidate(min=0)),
                ("ndecimals",                    8,    "Number of decimal places used for sampled parameter values", IntArgValidate(min=1)),
                ("save_sitelikes",           False,    "Saves file of site log-likelihoods (name determined by mcmc.out.sitelikes) that sump command can use in computing conditional predictive ordinates", BoolArgValidate),
                ])

        # Specify output options
        o = PhycasCommandOutputOptions()

        o.__dict__["_help_order"] = ["log", "trees", "params", "sitelikes", "pwkparams"]

        logf_spec = TextOutputSpec(prefix='mcmcoutput', help_str="The file specified by this setting saves the console output generated by mcmc() (set prefix to None to avoid creating this file).")
        o.__dict__["log"] = logf_spec

        t = TextOutputSpec(prefix='trees', suffix=".t", help_str="The nexus tree file in which all sampled tree topologies are saved. This file is equivalent to the MrBayes *.t file.")
        o.__dict__["trees"] = t

        p = TextOutputSpec(prefix='params', suffix=".p", help_str="The text file in which all sampled parameter values are saved. This file is equivalent to the MrBayes *.p file.")
        o.__dict__["params"] = p

        p = TextOutputSpec(prefix='sitelikes', suffix=".txt", help_str="The text file in which all sampled site log-likelihood values are saved.")
        o.__dict__["sitelikes"] = p

        p = TextOutputSpec(prefix='pwkparams', suffix=".txt", help_str="The text file in which all parameter values, topology index, and edge lengths are saved. Originally created so that the marginal likelihood could be estimated using the PWK method outside of Phycas.")
        o.__dict__["pwkparams"] = p

        PhycasCommand.__init__(self, args, "mcmc", "The mcmc command is used to conduct a Bayesian Markov chain Monte Carlo analysis.", o)

        # The data members added below should be hidden from the user because they are for use by phycas developers.
        # The roundabout way of introducing these data members is necessary because PhycasCommand.__setattr__ tries
        # to prevent users from adding new data members (to prevent accidental misspellings from causing problems)
        self.__dict__["debugging"] = False
        self.__dict__["draw_directly_from_prior"] = True       # If True, MCMCImpl.explorePrior function is used to draw samples directly from the prior during path sampling, which dramatically improves mixing compared to using MCMC proposals to explore the prior
        self.__dict__["reference_tree_source"] = None          # A TreeCollection that will serve as the source of the reference tree topology. The reference tree should represent the best tree topology known for the current model and data set. Specifying a reference tree makes it possible to determine if and when SAMC finds that tree topology. If a string is passed in, it is interpreted as a the path to a file with trees.
        # self.__dict__["cpo_patterns_only"] = False           # Not currently used. If True, each row of the sitelike output file will contain sampled log-likelihoods for each pattern, with the first row of the file holding the counts for each pattern. If False, the rows of the sitelike file will contain the log-likelihoods for each site in the order in which the sites occur in the data file (generally produces a larger file).

        self.__dict__["debug_allow_arbitrary_heating_powers"] = True    # if True, allow for example, heat_vector = [1.5]. Heating powers must still be positive, but are not required to be in descending order starting with 1.0 unless this is set to False

        # The data members added below are hidden from the user because they are set by the ss command
        self.__dict__["doing_steppingstone_sampling"] = False
        self.__dict__["ss_heating_likelihood"] = False
        #self.__dict__["ss_single_edgelen_ref_dist"] = False
        self.__dict__["ssobj"] = None

        # The data members added below are hidden from the user because they are set by the pwk command
        self.__dict__["doing_pwk"] = False

        # The data members added below are hidden from the user because they are experimental
        self.__dict__["pwkf"] = None            # PWKTMP
        self.__dict__["save_pwk_file"] = False  # PWKTMP # saves topology (index), edge lengths and parameters to pwkparams.txt file if True

        # The data members added below are hidden from the user because they are used internally when the users specifies a filename for cpofile
        self.__dict__["saving_sitelikes"] = False
        self.__dict__["sitelikef"] = None

        # The data members added below are hidden from the user because they are set when the mcmc command runs
        self.__dict__["ss_sampled_likes"] = None
        self.__dict__["ss_sampled_betas"] = None

    def checkSanity(self):
        """
        Place asserts in this function that should be checked before anything substantive
        is done during a call of self.run().
        """
        #cf = CommonFunctions(self)
        pass

    def __call__(self, **kwargs):
        self.set(**kwargs)
        self.checkSanity()
        #if self.fix_topology:
        #    self.ss_single_edgelen_ref_dist = False
        #else:
        #    self.ss_single_edgelen_ref_dist = True
        c = copy.deepcopy(self)
        mcmc_impl = MCMCImpl(c)

        if self.doing_pwk:
            return mcmc_impl
        else:
            if self.save_sitelikes:
                mcmc_impl.siteLikeFileOpen()
                if mcmc_impl.sitelikef:
                    self.saving_sitelikes = True
                else:
                    self.saving_sitelikes = False
                    print 'Could not open the sitelike file'

            mcmc_impl.setSiteLikeFile(self.sitelikef)

            mcmc_impl.run()

            self.ss_sampled_betas = mcmc_impl.ss_sampled_betas
            self.ss_sampled_likes = mcmc_impl.ss_sampled_likes

            if self.saving_sitelikes:
                mcmc_impl.siteLikeFileClose()

            mcmc_impl.unsetSiteLikeFile()

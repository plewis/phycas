###
### THIS FILE IS OBSOLETE AND SHOULD BE REMOVED!
###
import os, sys, math, threading, types, copy
import LikelihoodCore
import MarkovChain
import MCMCManager
from phycas.Conversions import *
from phycas.Likelihood import *
#from phycas.PDFGen import *
from phycas.Phylogeny import *
from phycas.ProbDist import *
from phycas.ReadNexus import *

# see http://mail.python.org/pipermail/python-list/2002-January/121376.html
import inspect

class Phycas(object):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Performs Bayesian phylogenetic MCMC analyses. The tree topology and
    edge lengths are updated via the Metropolis-Hastings algorithm (using
    the Larget-Simon LOCAL move without a molecular clock). Slice
    sampling is used to update all model parameters except edge lengths.

    For examples of how to use Phycas, see the Examples folder:
    phycas/Phycas/Phycas.py   <-- you are here
    phycas/Examples           <-- here is the Examples directory

    See the __init__ function below for variables that can be modified
    before Phycas is run.

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes the object with default values for all settings. If these
        defaults do not suit, they can be changed before mcmc() is called.

        """
        self.quiet                  = False     # If True, output will only be sent to the log file if open (see below); if False, output will be sent to the console as well

        # Variables associated with the brownian command
        self.brownian_input_tree_file    = None           # Set to the name of the input tree file. This setting should not be None at the time the brownian method is called.

        # Variables associated with Gelfand-Ghosh calculation
        self.gg_outfile             = 'gg.txt'  # File in which to save gg results (use None to not save results)
        self.gg_nreps               = 1         # The number of replicate simulations to do every MCMC sample
        self.gg_kvect               = [1.0]     # Vector of k values to use when computing Gm and Dm
        self.gg_save_postpreds      = False     # If True, all posterior predictive data sets will be saved
        self.gg_postpred_prefix     = 'pp'      # Prefix to use for posterior predictive dataset filenames (only used if gg_save_postpreds is True)
        self.gg_burnin              = 1         # Number of starting samples to skip when computing Gelfand-Ghosh measures
        self.gg_pfile               = None      # Name of parameter file to use for Gelfand-Ghosh calculations
        self.gg_tfile               = None      # Name of tree file to use for Gelfand-Ghosh calculations
        self.gg_bin_patterns        = False     # If True, patterns will be classified into 7 bins, corresponding to 'A only', 'C only', 'G only', 'T only', 'any 2 states', 'any 3 states' and 'any 4 states'. Gelfand-Ghosh statistics will be computed on this vector of counts instead of the complete vector of pattern counts. Can only be used for DNA/RNA data.
        self.gg_bincount_filename   = None      # If not None, and if gg_bin_patterns is True, the binned counts for the original dataset and all posterior predictive data sets will be saved to a file by this name

        # ***** IT IS BEST NOT TO CHANGE ANYTHING BELOW HERE *****
        self.debugging              = False      # If set to True expect lots of debug output (e.g. data pattern table)
        self.data_matrix            = None
        self.file_name_data_stored  = None
        self.file_name_trees_stored = None
        self.do_marginal_like       = False
        #(commented out to see if really being used anywhere) self.mcmc_manager           = MCMCManager.MCMCManager(self)
        self.heat_vector            = None      # Leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using self.nchains and self.heating_lambda
        self.stopwatch              = StopWatch()
        self.sim_model_tree         = None      # Will hold the model tree used by simulateDNA
        self.starting_tree          = None      # Will contain description of actual starting tree used
        self.warn_tip_numbers       = False     # True only if tip numbers were not able to be created using the tip names in the tree description (always False if starting_tree_source == 'random' because BuildTreeFromString is not called in this case)
        self.ntax                   = 0         # Will hold the actual number of taxa after data file read
        self.nchar                  = 0         # Will hold the actual number of characters after data file has been read
        self.npatterns              = 0         # Will hold the actual number of patterns after data file has been read
        self.taxon_labels           = []        # Will hold taxon labels from data file or default names if self.data_source equals None
        self.paramf                 = None
        self.treef                  = None
        self.tree_file_name         = ''        # Will hold tree file name (see openParameterAndTreeFiles)
        self.param_file_name        = ''        # Will hold parameter file name (see openParameterAndTreeFiles)
        self.tmp_simdata            = SimData()
        self.gg_Pm                  = 0.0       # Penalty component (same for all k)
        self.gg_Gm                  = []        # Vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []        # Vector of overall measures (one for each k in gg_kvect)
        self.reader                 = NexusReader()
        self.logf                   = None
        self._logFileName           = None
        self.addition_sequence      = []        # List of taxon numbers for addition sequence

        self.stored_tree_defs       = None
        self.ss_delta_beta          = 0.0
        self.doing_steppingstone_sampling = False
        self.path_sample            = None
        self.psf                    = None
        #self.pdf_splits_to_plot     = None
        self.param_file_name        = None
        self.tree_file_name         = None
        self.nsamples               = None
        self.ss_beta                = 1.0
        self.wangang_sampled_betas  = None
        self.wangang_sampled_likes  = None
        self.nsamples               = 0

        # make a copy of the vector of keys from __dict__ so that we can detect (in check_settings)
        # whether the user has accidentally introduced a new variable (by misspelling, for example)
        self.dict_keys = copy.copy(self.__dict__.keys())

    # see http://mail.python.org/pipermail/python-list/2002-January/121376.html
    def source_line():
        return inspect.getouterframes(inspect.currentframe())[1][2]

    def check_settings(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        At the end of the __init__ function, a variable dict_keys is
        created consisting of a copy of all keys in __dict__. Now,
        we expect len(__dict__) to be one greater than len(dict_keys)
        because __dict__ now includes an entry for dict_keys. If
        len(__dict__) is longer than this, then probably the user
        misspelled one of the variable names, thus accidentally
        adding another entry to __dict__. This gives us an
        opportunity to catch this kind of mistake.

        """
        if len(self.__dict__) > len(self.dict_keys) + 1:
            for k in self.__dict__.keys():
                if k not in self.dict_keys and not k == 'dict_keys':
                    print 'Error:',k,'is not a valid Phycas setting'

                    f = open('dist_keys.txt', 'w')
                    for kk in self.dict_keys:
                        f.write('%s\n' % kk)
                    f.close()

                    f = open('__dict__.txt', 'w')
                    for kk in self.__dict__.keys():
                        f.write('%s\n' % kk)
                    f.close()

                    sys.exit(0)

    def setEdgelenPrior(self, dist):
        self.internal_edgelen_prior = self.external_edgelen_prior = dist

    def getEdgelenPrior(self):
        self.phycassert(self.internal_edgelen_prior is self.external_edgelen_prior, "There are separate distributions for internal and external edge lengths")
        return self.internal_edgelen_prior

    edgelen_prior = property(getEdgelenPrior, setEdgelenPrior)

    def readTreesFromFile(self):
        if not self.file_name_trees_stored or (self.tree_file_name != self.file_name_trees_stored):
            self.reader.readFile(self.tree_file_name)
            self.taxon_labels = self.reader.getTaxLabels()  # shouldn't overwrite taxon_labels stored previously
            self.stored_tree_defs = self.reader.getTrees()
            self.phycassert(len(self.stored_tree_defs) > 0, 'expecting a trees block defining at least one tree in the nexus data file %s' % self.tree_file_name)
            self.file_name_trees_stored = self.tree_file_name    # prevents rereading same tree file later

    def likelihoods(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes the log-likelihood based on the current model of all trees
        in the file whose name is stored in tree_file_name.

        """
        self.check_settings()
        self.phycassert(len(self.tree_file_name) > 0, 'specify tree_file_name before calling the likelihoods function')
        self.readTreesFromFile()
        self.phycassert(self.data_source == 'file', "set data_source to 'file' and specify data_file_name before calling the likelihoods function")
        self.readDataFromFile()
        for t, topology in enumerate(self.stored_tree_defs):
            self.starting_tree = topology
            core = MCMCManager.LikelihoodCore(self)
            core.setupCore(True)    # specify zero_based_tips = True because topology came from file
            core.prepareForLikelihood()
            print 'Setting all edge lengths to 0.1'
            core.tree.setAllEdgeLens(0.1)
            print 'length of tree %d = %.5f' % (t,core.tree.edgeLenSum())
            print 'log-likelihood of tree %d = %.5f' % (t, core.calcLnLikelihood())

    def gg(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes Gelfand-Ghosh on a pre-existing MCMC sample defined in the
        files self.gg_pfile and self.gg_tfile.

        """
        self.check_settings()
        self.phycassert(self.gg_pfile, 'gg_pfile cannot be None if gg function called')
        self.phycassert(self.gg_tfile, 'gg_pfile cannot be None if gg function called')
        import GGImpl
        gelfand_ghosh = GGImpl.GelfandGhosh(self)
        self.gg_Pm, self.gg_Gm, self.gg_Dm = gelfand_ghosh.run()
        return (self.gg_Pm, self.gg_Gm, self.gg_Dm)

    def brownian(self):
        self.check_settings()
        import BrownianImpl
        brownian = BrownianImpl.Brownian(self)
        brownian.run()

    def calcDistances(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Computes a matrix of pairwise distances between sequences.

        """
        self.check_settings()
        self.phycassert(self.data_source == 'file', "data_source variable must equal 'file'")

        # matrix holds the data matrix
        matrix = []
        for i in range(self.ntax):
            matrix.append(self.data_matrix.getRow(i))

        # distance holds the JC69 distance matrix
        # Note: currently assumes no missing or ambiguous data (i.e. A vs. ? treated as difference)
        distance = []
        for i in range(self.ntax):
            distance_i = []
            for j in range(self.ntax):
                diff = 0.0
                if i < j:
                    diff = 0.0
                    for k in range(self.nchar):
                        xik = matrix[i][k]
                        xjk = matrix[j][k]
                        if xik != xjk:
                            diff += 1.0
                    p = diff/float(self.nchar)
                    if p > 0.74999:
                        p = 0.74999
                    jc = -0.75*math.log(1.0 - 4.0*p/3.0)
                    distance_i.append(jc)
                    #print 'saving distance[%d][%d] = %f' % (i,j,jc)
                elif i > j:
                    #print 'accessing [%d][%d] = %f' % (j,i,distance[j][i])
                    distance_i.append(distance[j][i])
                else:
                    distance_i.append(0.0)
            distance.append(distance_i)

#         for i in range(self.ntax):
#             for j in range(self.ntax):
#                 print "%10f" % distance[i][j],
#             print

        return distance

    # by default phycassert sys.exit.
    # When debugging, it is nice to set this to True so that you can see the stack trace
    #CPPCompiledInDebug = False

if __name__ == '__main__':
    print "The Phycas.py file should be imported, not run directly. To import it,"
    print "create a python script (i.e. a file with a name ending in .py) with the"
    print "following text:"
    print
    print "from phycas import *"
    print "myphycas = Phycas()"
    print "..."
    print
    print "See examples in the Examples and Tests folders for more information, or"
    print "consult the PDF manual."
    print

import os,sys,math,random,re 
#from phycas.TreeViewer import *
from phycas import *
from phycas.Utilities.PhycasCommand import *
from phycas.Utilities.CommonFunctions import CommonFunctions

class InvalidNumberOfColumnsError(Exception):
    def __init__(self, nparts, nexpected, line_num):
        self.msg = 'Number of values (%d) on line %d inconsistent with number of column headers (%d)' % (nparts, line_num, nexpected)
    def __str__(self):
        return self.msg

class NoTreesWarning(Exception):
    def __init__(self):
        self.msg = 'trees not specified in refdist command: assuming fixed tree topology for marginal likelihood analyses'
    def __str__(self):
        return self.msg

class ScriptGenImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Creates template scripts that user can modify before running.
    
    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes ScriptGenImpl object.
        
        """
        CommonFunctions.__init__(self, opts)
        self.script_filename = None
        self.scriptf = None
        self.sampledata_filename = None
        self.samplef = None

    def _openScriptFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the output script file.
        
        """
        self.phycassert(self.scriptf is None, 'Attempt made to open script file, but it is already open!')
        self.script_filename = self.opts.out.script._getFilename()
        try:
            self.scriptf = self.opts.out.script.open(self.stdout)
        except:
            print '*** Attempt to open script file (%s) failed.' % self.script_filename

        if self.scriptf:
            print 'Script file was opened successfully'
    
    def _closeScriptFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the output script file.
        
        """
        self.phycassert(self.scriptf is not None, 'Attempt made to close script file, but it is not open!')
        self.scriptf.close()
        if self.scriptf.closed:
            print 'Script file was closed successfully'
        self.scriptf = None
        
    def _openSampleDataFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the sample data file.
        
        """
        self.phycassert(self.samplef is None, 'Attempt made to open the sample data file, but it is already open!')
        self.sampledata_filename = self.opts.out.sampledata._getFilename()
        try:
            self.samplef = self.opts.out.sampledata.open(self.stdout)
        except:
            print '*** Attempt to open the sample data file (%s) failed.' % self.sampledata_filename

        if self.samplef:
            print 'The sample data file was opened successfully'
    
    def _closeSampleDataFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the sample data file.
        
        """
        self.phycassert(self.samplef is not None, 'Attempt made to close the sample data file, but it is not open!')
        self.samplef.close()
        if self.samplef.closed:
            print 'The sample data file was closed successfully'
        self.samplef = None
        
    def model_edgelen(self):
        self.scriptf.write("# Use independent exponential priors (mean 0.1) for each edge length parameter\n")
        self.scriptf.write("model.edgelen_prior = Exponential(10.0)\n")
        if model.edgelen_hyperprior is not None:
            self.scriptf.write("model.edgelen_hyperprior = %s\n\n" % model.edgelen_hyperprior.__repr__())
        else:
            self.scriptf.write("model.edgelen_hyperprior = None\n\n")

    def add_rate_heterogenetity(self, ratehet):
        if ratehet == 1 or ratehet==3:
            self.scriptf.write("# Add proportion of invariable sites submodel\n")
            self.scriptf.write("model.pinvar_model = True\n")
            self.scriptf.write("model.pinvar_prior = Beta(1.0, 1.0)\n\n")
        else:
            self.scriptf.write("# Assume no invariable sites\n")
            self.scriptf.write("model.pinvar_model = False\n\n")
        if ratehet == 2 or ratehet==3:
            self.scriptf.write("# Add discrete gamma rate heterogeneity submodel\n")
            self.scriptf.write("model.num_rates = 4\n")
            self.scriptf.write("model.gamma_shape_prior = Exponential(1.0)\n\n")
        else:
            self.scriptf.write("# Assume rate homogeneity across sites\n")
            self.scriptf.write("model.num_rates = 1\n\n")
    
    def model_jc(self, ratehet, title = "Set up JC model"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("model.type = 'jc'\n")
        self.add_rate_heterogenetity(ratehet)
        
    def model_hky(self, ratehet, title = "Set up HKY model"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("model.type = 'hky'\n")
        self.scriptf.write("model.kappa_prior = Exponential(1.0)\n")
        self.scriptf.write("model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))\n\n")
        self.add_rate_heterogenetity(ratehet)
        
    def model_gtr(self, ratehet, title = "Set up GTR model"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("model.type = 'gtr'\n")
        self.scriptf.write("model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))\n")
        self.scriptf.write("model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))\n\n")
        self.add_rate_heterogenetity(ratehet)
        
    def model_codon(self, title = "Set up GTR model"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("model.type = 'codon'\n")
        self.scriptf.write("model.omega_prior = Exponential(1.0)\n")
        self.scriptf.write("model.state_freq_prior = Dirichlet([1.0]*61)\n\n")
        self.add_rate_heterogenetity(0)
        
    def sump_analysis(self, cpo = False, title = "Summarize the posterior distribution of model parameters"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("sump.file = 'params.p'\n")
        self.scriptf.write("sump.burnin = %d\n" % sump.burnin)
        if cpo:
            self.scriptf.write("sump.cpofile = 'sitelikes.txt'\n")
            self.scriptf.write("sump.cpo_cutoff = 0.1\n")
        self.scriptf.write("sump.out.log.prefix = 'sump-log'\n")
        self.scriptf.write("sump.out.log.mode = REPLACE\n")
        self.scriptf.write("sump()\n\n")
        
    def sumt_analysis(self, title = "Summarize the posterior distribution of tree topologies and clades"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("sumt.trees = 'trees.t'\n")
        self.scriptf.write("sumt.burnin = %d\n" % sumt.burnin)
        self.scriptf.write("sumt.tree_credible_prob = %g  # set to 1.0 if you want KL information to be estimated\n" % sumt.tree_credible_prob)
        self.scriptf.write("sumt.save_splits_pdf = %s  # if True, may get really large PDF file if tree_credible_prob = 1.0 and number of sampled splits is large\n" % sumt.save_splits_pdf)
        self.scriptf.write("sumt.save_trees_pdf = %s   # if True, may get really large PDF file if tree_credible_prob = 1.0 and number of sampled trees is large\n" % sumt.save_trees_pdf)
        self.scriptf.write("sumt.out.log.prefix = 'sumt-log'\n")
        self.scriptf.write("sumt.out.log.mode = REPLACE\n")
        self.scriptf.write("sumt.out.trees.prefix = 'sumt-trees'\n")
        self.scriptf.write("sumt.out.trees.mode = REPLACE\n")
        self.scriptf.write("sumt.out.splits.prefix = 'sumt-splits'\n")
        self.scriptf.write("sumt.out.splits.mode = REPLACE\n")
        self.scriptf.write("sumt()\n\n")
        
    def mcmc_analysis(self, execute = True, title = "Conduct a Markov chain Monte Carlo (MCMC) analysis that samples from the posterior distribution"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("mcmc.ncycles = %g\n" % mcmc.ncycles)
        self.scriptf.write("mcmc.sample_every = %g\n" % mcmc.sample_every)
        self.scriptf.write("mcmc.report_every = %g\n" % mcmc.report_every)
        self.scriptf.write("#mcmc.starting_tree_source = TreeCollection(newick='(1:.01,2:0.01,(3:0.01,4:0.01):0.01)')\n")
        self.scriptf.write("#mcmc.starting_tree_source = TreeCollection(filename='nexustreefile.tre')\n")
        self.scriptf.write("mcmc.fix_topology = False\n")
        self.scriptf.write("mcmc.allow_polytomies = False\n")
        self.scriptf.write("mcmc.bush_move_weight = 0\n")
        self.scriptf.write("mcmc.ls_move_weight = 100\n")
        #self.scriptf.write("mcmc.polytomy_prior = True\n")
        #self.scriptf.write("mcmc.topo_prior_C = 2.0\n")
        self.scriptf.write("mcmc.out.log = 'mcmcoutput.txt'\n")
        self.scriptf.write("mcmc.out.log.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.trees = 'trees.t'\n")
        self.scriptf.write("mcmc.out.trees.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.params = 'params.p'\n")
        self.scriptf.write("mcmc.out.params.mode = REPLACE\n")
        if execute:
            self.scriptf.write("mcmc()\n\n")
        else:
            self.scriptf.write("#mcmc()\n\n")

    def polytomy_analysis(self, execute = True, title = "Conduct a Markov chain Monte Carlo (MCMC) analysis that samples polytomous (as well as fully-resolved) tree topologies"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("mcmc.ncycles = %g\n" % mcmc.ncycles)
        self.scriptf.write("mcmc.sample_every = %g\n" % mcmc.sample_every)
        self.scriptf.write("mcmc.report_every = %g\n" % mcmc.report_every)
        self.scriptf.write("#mcmc.starting_tree_source = TreeCollection(newick='(1:.01,2:0.01,(3:0.01,4:0.01):0.01)')\n")
        self.scriptf.write("#mcmc.starting_tree_source = TreeCollection(filename='nexustreefile.tre')\n")
        self.scriptf.write("mcmc.fix_topology = False\n")
        self.scriptf.write("mcmc.allow_polytomies = True\n")
        self.scriptf.write("mcmc.bush_move_weight = 50\n")
        self.scriptf.write("mcmc.ls_move_weight = 50\n")
        self.scriptf.write("mcmc.polytomy_prior = True\n")
        self.scriptf.write("mcmc.topo_prior_C = 2.0\n")
        self.scriptf.write("mcmc.out.log = 'mcmcoutput.txt'\n")
        self.scriptf.write("mcmc.out.log.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.trees = 'trees.t'\n")
        self.scriptf.write("mcmc.out.trees.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.params = 'params.p'\n")
        self.scriptf.write("mcmc.out.params.mode = REPLACE\n")
        if execute:
            self.scriptf.write("mcmc()\n\n")
        else:
            self.scriptf.write("#mcmc()\n\n")
        self.sumt_analysis()

    def cpo_analysis(self, title = "Using MCMC setup above, compute conditional predictive ordinates (CPO) for each site and for the entire data set"):
        self.mcmc_analysis(False, "Set up (but do not execute) a Markov chain Monte Carlo (MCMC) analysis that samples from the posterior distribution")
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("cpo.out.sitelike.prefix = 'sitelikes'\n")
        self.scriptf.write("cpo()\n\n")
        self.sump_analysis(True, "Run sump command to summarize CPO values saved during MCMC")

    def refdist_analysis(self, title = "Estimate the reference distribution to use with Generalized SS"):
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("refdist.burnin = %g\n" % refdist.burnin)
        self.scriptf.write("refdist.params = 'params.p'\n")
        self.scriptf.write("refdist.trees = 'trees.t'\n")
        self.scriptf.write("refdist.out.refdistfile = 'refdist.txt'\n")
        self.scriptf.write("refdist.out.refdistfile.mode = REPLACE\n")
        self.scriptf.write("refdist()\n\n")
        
    def idr_analysis(self, title = "Estimate the marginal likelihood using the Inflated Density Ratio (IDR) method"):
        self.mcmc_analysis()
        self.scriptf.write("# %s\n" % title)
        self.scriptf.write("idr.burnin = 1\n")
        self.scriptf.write("idr.data_source = '%s'\n" % self.opts.datafile)
        self.scriptf.write("idr.params = 'params.p'\n")
        self.scriptf.write("idr.trees = 'trees.t'\n")
        self.scriptf.write("idr.autork = True\n")
        self.scriptf.write("idr.out.log.prefix = 'idr-log'\n")
        self.scriptf.write("idr.out.log.mode = REPLACE\n")
        self.scriptf.write("idr()\n\n")
        
    def steppingstone_analysis(self, title = "Estimate the marginal likelihood using the Generalized Stepping-stone (GSS) method"):
        self.mcmc_analysis("Conduct MCMC analysis for purpose of generating a reference distribution")

        self.scriptf.write("# %s\n" % title)
        self.refdist_analysis()

        self.scriptf.write("# Do not fix the tree topology for stepping-stone analysis\n")
        self.scriptf.write("mcmc.fix_topology = False\n\n")

        self.scriptf.write("# Ucomment and replace mcmc.starting_tree_source below with correct newick \n")
        self.scriptf.write("# tree description complete with edge lengths if mcmc.fix_topology is set to True\n")
        self.scriptf.write("#mcmc.starting_tree_source = TreeCollection(newick='(1:0.01,2:0.01,(3:0.01,4:0.01):0.01)')\n\n")
        
        self.scriptf.write("# Choose different output file names to avoid overwriting the ones used for the reference distribution\n")
        self.scriptf.write("mcmc.out.log = 'ss-output.txt'\n")
        self.scriptf.write("mcmc.out.log.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.trees = 'ss-trees.t'\n")
        self.scriptf.write("mcmc.out.trees.mode = REPLACE\n")
        self.scriptf.write("mcmc.out.params = 'ss-params.p'\n")
        self.scriptf.write("mcmc.out.params.mode = REPLACE\n")

        self.scriptf.write("# Set up and run the ss command (this will make use of many mcmc settings)\n")
        self.scriptf.write("ss.nstones = %d\n" % ss.nstones)
        self.scriptf.write("ss.ncycles = %d\n" % ss.ncycles)
        self.scriptf.write("ss.sample_every = %d\n" % ss.sample_every)
        self.scriptf.write("ss.report_every = %d\n" % ss.report_every)
        self.scriptf.write("ss.refdist_is_prior = %s\n" % (ss.refdist_is_prior and "True" or "False"))
        self.scriptf.write("ss.refdistfile = 'refdist.txt'\n")
        self.scriptf.write("ss.shape1 = %f\n" % ss.shape1)
        self.scriptf.write("ss.shape2 = %f\n" % ss.shape2)
        self.scriptf.write("ss()\n\n")
        
        self.scriptf.write("# Running sump on the param file output by the ss command will calculate the marginal likelihood estimate\n")
        self.scriptf.write("sump.file = 'ss-params.p'\n")
        self.scriptf.write("sump.burnin = 1\n")
        self.scriptf.write("sump.out.log.prefix = 'ss-sump-log'\n")
        self.scriptf.write("sump.out.log.mode = REPLACE\n")
        self.scriptf.write("sump()\n\n")
        
    def saveSampleData(self):
        self._openSampleDataFile()
        self.samplef.write("#nexus\n\n")
        self.samplef.write("begin data;\n")
        self.samplef.write("  dimensions ntax=5 nchar=100;\n")
        self.samplef.write("  format datatype=dna missing=? gap=-;\n")
        self.samplef.write("  matrix\n")
        self.samplef.write("    Picea_pungens_AF456382       aaagattacagattAacttattatactcctgaAtatcaGaccaaagatacGgatatTttggcagcattccgagtaactcctcaaccaggGgtgccGccCg\n")
        self.samplef.write("    Asplenium_nidus_AF525270     AAAGATTATCGACTGACTTATTACACCCCCGAATACAAGACCAAAGATACCGATATCTTGGCAGCTTTCCGGATGACCCCACAACCCGGAGTACCAGCTG\n")
        self.samplef.write("    Nicotiana_tabacum_J01450     aaagagtacaaattgacttattatactcctgagtaccaaaccaaggatactgatatattggcagcattccgagtaactcctcaacctggagttccacctg\n")
        self.samplef.write("    Avena_sativa_L15300          AAAGATTATAAATTGACTTACTACACCCCGGAGTATGAAACCAAGGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCTGGGGTTCCGCCGG\n")
        self.samplef.write("    Iris_unguicularis_AJ309693   AAAGATTACAGATTGACTTATTATACTCCTGATTACGAAACCAAAGATACTGATATCTTGGCAGCATTCCGAGTAACTCCTCAACCCGGAGTTCCTGCTG\n")
        self.samplef.write("    ;\n")
        self.samplef.write("end;\n")
        self._closeSampleDataFile()
    
    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates example Python script.
        
        """
        self._openScriptFile()
        self.scriptf.write("from phycas import *\n\n")
        
        if self.opts.seed > 0:
            self.scriptf.write('setMasterSeed(%d)\n\n' % self.opts.seed)
        
        if self.opts.model == 'jc':
            self.model_jc(0)
        elif self.opts.model == 'jc+i':
            self.model_jc(1)
        elif self.opts.model == 'jc+g':
            self.model_jc(2)
        elif self.opts.model == 'jc+i+g':
            self.model_jc(3)
        elif self.opts.model == 'hky':
            self.model_hky(0)
        elif self.opts.model == 'hky+i':
            self.model_hky(1)
        elif self.opts.model == 'hky+g':
            self.model_hky(2)
        elif self.opts.model == 'hky+i+g':
            self.model_hky(3)
        elif self.opts.model == 'gtr':
            self.model_gtr(0)
        elif self.opts.model == 'gtr+i':
            self.model_gtr(1)
        elif self.opts.model == 'gtr+g':
            self.model_gtr(2)
        elif self.opts.model == 'gtr+i+g':
            self.model_gtr(3)
        elif self.opts.model == 'codon':
            self.model_codon()
            
        self.model_edgelen()

        if self.opts.datafile is None:
            self.scriptf.write("mcmc.data_source = None\n")
        elif self.opts.datafile == 'sample.nex':
            self.saveSampleData()
            self.scriptf.write("mcmc.data_source = '%s'\n\n" % self.opts.datafile)
        else:
            self.scriptf.write("mcmc.data_source = '%s'\n\n" % self.opts.datafile)

        if self.opts.analysis == 'mcmc':
            self.mcmc_analysis()
            self.sump_analysis()
            self.sumt_analysis()
        if self.opts.analysis == 'poly':
            self.polytomy_analysis()
        elif self.opts.analysis == 'cpo':
            self.cpo_analysis()
        elif self.opts.analysis == 'idr':
            self.idr_analysis()
        elif self.opts.analysis == 'ss':
            self.steppingstone_analysis()

        self._closeScriptFile()

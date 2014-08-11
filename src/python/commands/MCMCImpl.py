import os,sys,math,random
from phycas import *
from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions
from MCMCManager import MCMCManager
from phycas.probdist import StopWatch
from phycas.readnexus import NexusReader

class MCMCImpl(CommonFunctions):
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
        self.models                 = None      # This variable is set in setup

        # These copied over from Phycas.py - many are not used and should be weeded out
        self.data_matrix            = None
        self.file_name_trees_stored = None
        self.do_marginal_like       = False
        self.mcmc_manager           = MCMCManager(self)
        self.heat_vector            = None      # Leave set to None unless you are implementing some ad hoc heating scheme. This vector ordinarily computed using nchains and self.heating_lambda
        self.stopwatch              = StopWatch()
        self.sim_model_tree         = None      # Will hold the model tree used by simulateDNA
        self.starting_tree          = []        # starting_tree[i] will contain reference to phylogeny.Tree object for chain i
        self.warn_tip_numbers       = False     # True only if tip numbers were not able to be created using the tip names in the tree description (always False if starting_tree_source == 'random' because BuildTreeFromString is not called in this case)
        self.ntax                   = 0         # Will hold the actual number of taxa after data file read
        self.nchar                  = 0         # Will hold the actual number of characters after data file has been read
        self.npatterns          = []        # Will hold the actual number of patterns for each subset after data file has been read
        self.taxon_labels           = []        # Will hold taxon labels from data file or default names if self.data_source equals None
        #self.sssf                  = None
        self.paramf                 = None
        self.treef                  = None
        self.sitelikef              = None
        #self.tree_file_name         = ''        # Will hold tree file name (see openParameterAndTreeFiles)
        #self.param_file_name        = ''        # Will hold parameter file name (see openParameterAndTreeFiles)
        #self.tmp_simdata            = SimData()
        self.gg_Pm                  = 0.0       # Penalty component (same for all k)
        self.gg_Gm                  = []        # Vector of goodness-of-fit components (one for each k in gg_kvect)
        self.gg_Dm                  = []        # Vector of overall measures (one for each k in gg_kvect)
        self.reader                 = NexusReader()
        #self._logFileName           = None
        self.addition_sequence      = []        # List of taxon numbers for addition sequence
        self.ref_tree               = None      # a reference tree against which the current MCMC tree can be compared (user specifies with, e.g., "mcmc.reference_tree_source = TreeCollection(file='hkyml.tre')")
        self.stored_tree_defs       = None
        self.psf                    = None
        self.pdf_splits_to_plot     = None
        #self.param_file_name        = None
        #self.tree_file_name         = None
        self.burnin                 = 0     # same as self.opts.burnin except for path sampling, when it drops to 0 after first beta value
        self.ncycles                = 0
        self.cycle_start            = None     # used in path sampling to avoid starting over the cycle count for each beta value
        self.cycle_stop             = None     # total number of cycles (used for computing time remaining)
        self.last_adaptation        = 0
        #POLTMP self.next_adaptation        = 0
        self.ss_beta                = 1.0
        self.ss_beta_index          = 0
        self.ss_sampled_betas       = None
        self.ss_sampled_likes       = None
        self.siteIndicesForPatternIndex = None

    def setSiteLikeFile(self, sitelikef):
        if sitelikef is not None:
            self.sitelikef = sitelikef

    def siteLikeFileSetup(self, coldchain):
        if self.sitelikef is not None:
            # Set up the siteIndicesForPatternIndex, which holds a list of sites for each pattern index
            # This allows us to spit out site likelihoods for each site, even though site likelihoods
            # are stored for patterns, many of which represent numerous sites
            v = coldchain.likelihood.getCharIndexToPatternIndex()
            nexcluded = 0
            self.phycassert(len(v) == self.nchar,'Number of sites returned by coldchain.likelihood.getCharIndexToPatternIndex (%d) differs from MCMCImpl.nchar (%d) in MCMCImpl.siteLikeFileSetup()' % (len(v), self.nchar))
            npatterns = coldchain.likelihood.getNPatterns()
            self.siteIndicesForPatternIndex = []
            for i in range(npatterns):
                self.siteIndicesForPatternIndex.append([])
            for i,p in enumerate(v):
                if p < npatterns:
                    self.siteIndicesForPatternIndex[p].append(i)
                else:
                    nexcluded += 1;

    def unsetSiteLikeFile(self):
        self.sitelikef = None
        self.siteIndicesForPatternIndex = None

    def adaptSliceSamplers(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Cycle through all slice samplers and adapt each one. Adaptation 
        adjusts the slice unit width of each slice sampler in an attempt to 
        bring it closer to the optimum width using experience from past 
        sampling attempts. Non-slice-sampling updaters are also tuned so that
        eventually they hover around the target acceptance rate.

        """
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getAllUpdaters():
            if p.hasSliceSampler():
                s = p.getSliceSampler()
                if s.getNumSamples() > 0:
                    s.adaptSimple(self.opts.adapt_simple_param)

    def resetUpdaterDiagnostics(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Calls resetDiagnostics() for all slice-samplers and Metropolis-Hastings
        updaters.

        """
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getAllUpdaters():
            if p.hasSliceSampler():
                s = p.getSliceSampler()
                s.resetDiagnostics()
            else:
                p.resetDiagnostics()

    def reportUpdaterEfficiency(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If not in the burnin period (cycle negative), cycle through all
        slice samplers and adapt each one. Adaptation adjusts the slice unit 
        width of each slice sampler in an attempt to bring it closer to the 
        optimum width using experience from past sampling attempts. 
        Non-slice-sampling updaters are also tuned so that eventually they
        hover around the target acceptance rate.

        """
        summary = ''
        # need to adapt all chains, not just the cold one!
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getAllUpdaters():
            p_summary = ''
            nm = p.getName()
            #print '~POL~ nm =',nm
            if p.hasSliceSampler():
                s = p.getSliceSampler()
                if s.getNumSamples() > 0:
                    # s.adaptSimple(self.opts.adapt_simple_param)
                    mode = s.getMode()
                    accept_pct = 100.0*float(s.getNumSamples())/float(s.getNumFuncEvals())
                    p_summary += ' * efficiency = %.1f%%, mode=%.5f (%s)\n' % (accept_pct, mode, nm)
            else:
                naccepts = p.getNumAccepts()
                nattempts = p.getNumAttempts()
                if nattempts > 0:
                    accept_pct = 100.0*float(naccepts)/float(nattempts)
                    if nattempts == 1:
                        p_summary += '   accepted %.1f%% of 1 attempt (%s)\n' % (accept_pct, nm)
                    else:
                        p_summary += '   accepted %.1f%% of %d attempts (%s)\n' % (accept_pct, nattempts, nm)
                else:
                    accept_pct = 0.0
            summary += p_summary

        if self.opts.verbose and summary != '':
            self.output('\nUpdater diagnostics (* = slice sampler):')
            self.output(summary)

        if self.opts.verbose and self.opts.nchains > 1:
            self.output('Chain swaps (upper: attempted swaps, lower: accepted swaps):')

            row = '%20s' % 'heating power -->'
            row += ' '.join(['%12.2f' % self.mcmc_manager.chains[k].heating_power for k in range(self.opts.nchains)])
            self.output(row)

            for ii in range(self.opts.nchains):
                row = '%20.2f' % (self.mcmc_manager.chains[ii].heating_power,)
                row += ' '.join(['%12d' % self.mcmc_manager.swap_table[ii][jj] for jj in range(self.opts.nchains)])
                self.output(row)
            self.output()

    def obsoleteUpdateAllUpdaters(self, chain, chain_index, cycle):
        # This function abandoned; functionality moved to C++ side
        # for speed reasons: see MCMCChainManager::updateAllUpdaters
        #import gc
        #gc.enable()
        assert False, 'NOT SUPPOSED TO BE HERE'
        if self.opts.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))

        if chain.all_updaters_list is None:
            chain.all_updaters_list = [(m, m.getWeight()) for m in  chain.chain_manager.getAllUpdaters()]
            chain.all_updaters_list.extend(chain.python_only_moves)
        for t in chain.all_updaters_list:
            p, w = t
            #print "param = %s (weight = %f), chain = %d" % (p.getName(), w, chain_index)
            for x in range(w):
                if self.opts.debugging:
                    p.setSaveDebugInfo(True)
                p.update()
                #if p.getName() == 'Bush move':
                #    print '  counts after bush move =',gc.get_count()
                #    print '  thresholds =',gc.get_threshold()
                if self.opts.debugging:
                    tmpf.write('%s | %s\n' % (p.getName(), p.getDebugInfo()))

        if self.opts.debugging:
            tmpf.close()

    def showTopoPriorInfo(self):
        m = self.mcmc_manager.getColdChain()
        self.output('\nTopology prior:')
        if not self.opts.allow_polytomies:
            self.output('  flat across all fully-resolved tree topologies (polytomies not allowed)')
        else:
            if self.opts.polytomy_prior:
                self.output('  Prior type: polytomy prior')
            else:
                self.output('  Prior type: resolution class prior')

            self.output('  Prior strength (C): %s' % self.opts.topo_prior_C)
            self.output('  Prior probability for each resolution class:')
            self.output('  Note: 0.00000000 does *not* mean that the prior is zero! It simply')
            self.output('        indicates that the prior is less than 0.000000005\n')
            self.output('%20s %20s' % ('internal nodes', 'prior probability'))
            self.output('%20s %20s' % ('--------------', '-----------------'))
            jpm = self.mcmc_manager.getColdChainManager().getJointPriorManager()
            topo_prior_calculator = jpm.getTopoProbCalculator()
            topo_priors = topo_prior_calculator.getRealizedResClassPriorsVect()
            for i,v in enumerate(topo_priors):
                if i == 0:
                    denom = v   # first element of vector is log of normalization constant (sum of all other elements)
                else:
                    topo_prior = math.exp(v - denom)
                    self.output('%20d %20.8f' % (i,topo_prior))
            self.output()

    def showParamInfo(self, p):
        if p.computesUnivariatePrior() or p.computesMultivariatePrior() or p.computesTreeLengthPrior():
            self.output('  Parameter name:     %s' % p.getName())
            self.output('  Prior distribution: %s' % p.getPriorDescr())
            if p.isMasterParameter():
                self.output('  Master parameter (no current value)')
            else:
                if p.computesUnivariatePrior():
                    v = p.getCurrValueFromModel()
                    self.output('  Current value:      %s' % v)
                else:
                    v = p.listCurrValuesFromModel()
                    self.output('  Current value:      %s' % ','.join(['%.5f' % x for x in v]))
            self.output('  Prior log-density:  %s' % p.getLnPrior())
            self.output()

    def siteLikeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the site log-likelihood file.

        """
        self.phycassert(self.sitelikef is None, 'Attempt made to open MCMCImpl.sitelikef, but it is already open!')
        sitelnl_file_spec = self.opts.out.sitelikes
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
        self.phycassert(self.sitelikef is not None, 'Attempt made to close MCMCImpl.sitelikef, but it is not open!')
        self.sitelikef.close()
        self.sitelikef = None

    def treeFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the tree file and writes a translate table.

        """
        tree_file_spec = self.opts.out.trees
        self.treef = None
        try:
            self.treef = tree_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open tree file (%s) failed.' % self.opts.out.trees.filename

        if self.treef:
            self.mcmc_manager.treeFileHeader(self.treef)

    def treeFileClose(self):
        self.treef.write('end;\n')
        self.treef.close()

    def paramFileOpen(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the parameter file and writes a header line.

        """
        param_file_spec = self.opts.out.params
        self.paramf = None
        try:
            self.paramf = param_file_spec.open(self.stdout)
        except:
            print '*** Attempt to open parameter file (%s) failed.' % self.opts.out.params.filename

        if self.paramf:
            self.mcmc_manager.paramFileHeader(self.paramf)
            self.paramf.write('\n')

    def paramFileClose(self):
        if self.paramf is not None:
            self.paramf.close()

    def openParameterAndTreeFiles(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates parameter and tree file names based on the data file name or the
        user-supplied prefix and opens the files

        """
        #prefix = self.getPrefix()
        #self.param_file_name = prefix + '.p'
        #self.tree_file_name = prefix + '.t'

        self.paramFileOpen()
        self.treeFileOpen()

    def _loadData(self, matrix):
        self.phycassert(matrix is not None, 'Tried to load data from a non-existant matrix')
        self.data_matrix = matrix
        self.taxon_labels = matrix.getTaxLabels()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar() # used for Gelfand-Ghosh simulations only
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")

    def cumLagrange(self, which, x, y):
        xx = x[which]
        term1 = (xx - x[0])*(y[0] + (xx - x[0])*(y[1] - y[0])/(2.0*(x[1] - x[0])))
        term2a = 2.0*(xx**2.0) - xx*x[0] - (x[0]**2.0) + 2.0*x[0]*x[1] - 3.0*xx*x[1]
        term2b = (y[2] - y[1])/(x[2] - x[1]) - (y[1] - y[0])/(x[1] - x[0])
        term2c = (xx - x[0])/(x[2] - x[0])
        term2 = term2a*term2b*term2c/6.0
        cum = term1 + term2
        return cum

    def getStartingTree(self):
        #         if self.starting_tree is None:
        #             if False:
        #                 if self.opts.starting_tree_source == 'file':
        #                     self.phycassert(self.data_source, "Specified starting_tree_source to be 'file' when data_source was None (file was not read)")
        #                     tree_defs = self.reader.getTrees()
        #                     self.phycassert(len(tree_defs) > 0, 'a trees block defining at least one tree must be stored in the nexus data file')
        #                     # Grab first tree description in the data file
        #                     # TODO allow some other tree than the first
        #                     self.starting_tree = tree_defs[0]
        #                 elif self.opts.starting_tree_source == 'usertree':
        #                     self.starting_tree = Newick(self.opts.tree_topology)
        #                 elif self.opts.starting_tree_source == 'random':
        #                     self.phycassert(self.ntax > 0, 'expecting ntax to be greater than 0')
        #                     self.starting_tree = None
        #                 else:
        #                     self.phycassert(False, "starting_tree_source should equal 'random', 'file', or 'usertree', but instead it was this: %s" % self.starting_tree_source)
        #             else:
        # If user failed to specify starting_tree_source, get starting tree from randomtree object
        # as it is currently configured
        tr_source = self.opts.starting_tree_source
        if tr_source is None:
            tr_source = randomtree()
        try:
            tr_source.setActiveTaxonLabels(self.taxon_labels)
            i = iter(tr_source)
            # self.starting_tree = i.next()
            self.starting_tree.append(i.next())
        except:
            self.stdout.error("A starting tree could not be obtained from the starting_tree_source")
            raise
        t = self.starting_tree[-1]
        num_degree_two_nodes = t.deroot()
        if num_degree_two_nodes > 0:
            self.stdout.warning("A total of %d degree-2 nodes were removed from tree defined in starting_tree_source" % num_degree_two_nodes)
        return t

    def storeRefTreeIfSupplied(self):
        cold_chain = self.mcmc_manager.getColdChain()

        # If a reference tree was specified, create that tree now
        tr_source = self.opts.reference_tree_source
        if tr_source is not None:
            try:
                tr_source.setActiveTaxonLabels(self.taxon_labels)
                i = iter(tr_source)
                self.ref_tree = i.next()
            except:
                self.stdout.error("A reference tree could not be obtained from the specified reference_tree_source")
        if self.ref_tree is not None:
            cold_chain.chain_manager.setRefTree(self.ref_tree)

    def setup(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function is for parts of the setup that should occur right before
        run() is called. Setup is deferred until this point to give the
        user a chance to change the default settings before the call to
        run(). This function does these things:
        1) reads the data and ensures that the taxon_labels list is filled
        with the correct number of taxon labels;
        2) creates a starting tree description;
        3) creates an appropriate heat_vector
        4) calls MCMCManager's createChains function to handle setup for each
        individual chain;
        5) opens the parameter and tree files; and
        6) establishes an output log file name if requested

        """
        ds = self.opts.data_source
        if ds is None:
            # User apparently wants to run without data
            self.data_matrix = None
            self.ntax = self.opts.ntax
            self.nchar = 0 # used for Gelfand-Ghosh simulations only
            self.phycassert(self.ntax > 0, 'Number of taxa (mcmc.ntax) should be > 0 if mcmc.data_source is None')
            self.taxon_labels = ['taxon%d' % (i+1,) for i in range(self.ntax)]
        else:
            mat = ds.getMatrix()
            self.phycassert(mat is not None, 'Data matrix could not be input')
            self._loadData(mat)
            self.phycassert(self.ntax > 0, 'Number of taxa in data matrix was 0')

        #print 'In MCMCImpl.py, function setup():'
        #print '  mat =', mat
        #print '  self.nchar = %d' % self.nchar

        # Next line creates a default partition if a partition was not defined by the user
        self.opts.partition.validate(self.nchar)

        # Ask for a partition report, passing self as the reporter (object that has an output function)
        #self.opts.partition.partitionReport(self)  # delayed now until after all_missing has been determined

        self.models         = [m for (n,s,m) in self.opts.partition.subset]
        self.model_names    = [n for (n,s,m) in self.opts.partition.subset]

        # Perform sanity checks on models
        for m,n in zip(self.models, self.model_names):
            #print '==> checking model %s' % n
            bad_priors = m.checkPriorSupport()
            self.phycassert(len(bad_priors) == 0, 'In model %s, prior support is incorrect for these parameters:\n%s' % (n,'  \n'.join([p for p in bad_priors])))

            if m.edgelen_prior is not None:
                # set both internal and external edge length priors to edgelen_prior
                m.internal_edgelen_prior = m.edgelen_prior
                m.external_edgelen_prior = m.edgelen_prior
            else:
                # Ensure that user has specified both internal and external edge length priors, or has specified a tree length prior
                if m.tree_length_prior is None:
                    self.phycassert(m.internal_edgelen_prior is not None, 'In model %s, internal_edgelen_prior cannot be None if edgelen_prior and tree_length_prior are both None' % n)
                    self.phycassert(m.external_edgelen_prior is not None, 'In model %s, external_edgelen_prior cannot be None if edgelen_prior and tree_length_prior are both None' % n)

            if m.edgelen_hyperprior is not None:
                # Ensure that both internal and external edgelen priors are Exponential
                if m.internal_edgelen_prior.getDistName() != 'Exponential':
                    m.internal_edgelen_prior = Exponential(1.0)
                    self.warning('In model %s, internal_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)
                if m.external_edgelen_prior.getDistName() != 'Exponential':
                    m.external_edgelen_prior = Exponential(1.0)
                    self.warning('In model %s, external_edgelen_prior reset to Exponential because edgelen_hyperprior was specified' % n)

        # Determine heating levels if multiple chains
        if self.opts.heat_vector == None:
            if self.opts.nchains == 1:
                self.heat_vector = [1.0]
            else:
                # Determine vector of powers for each chain
                self.heat_vector = []
                for i in range(self.opts.nchains):
                    # For n=5 chains (1 cold, 4 heated), min_heat_power = 0.5, we have:
                    # lambda = (1 - min_heat_power)/(min_heat_power*(n-1))
                    #        = (1 - 0.5)/(0.5*4)
                    #        = 0.25
                    # 0 1.000 = 1/1.0 cold chain explores posterior
                    # 1 0.800 = 1/1.25
                    # 2 0.667 = 1/1.50
                    # 3 0.571 = 1/1.75
                    # 4 0.500 = 1/2.00
                    z = self.opts.min_heat_power
                    n = self.opts.nchains
                    lamda = (1.0 - z)/(z*(n - 1.0))
                    temp = 1.0/(1.0 + float(i)*lamda)
                    self.heat_vector.append(temp)
        else:
            # User supplied his/her own heat_vector; perform sanity checks
            self.heat_vector = self.opts.heat_vector
            if not self.opts.debug_allow_arbitrary_heating_powers:

                self.phycassert(1.0 in self.heat_vector, 'user-supplied heat_vector does not allow for a cold chain (one power must be 1.0)')
                h = list(self.heat_vector)
                h.sort(reverse=True)
                if h != self.heat_vector:
                    self.phycassert(False, 'chain heating powers must be in decreasing order')
                self.phycassert(h[-1] > 0.0, 'all chain heating powers must be positive')

        self.mcmc_manager.createChains()
        self.storeRefTreeIfSupplied()
        self.openParameterAndTreeFiles()

        cold_chain = self.mcmc_manager.getColdChain()
        all_missing = cold_chain.likelihood.getListOfAllMissingSites()
        partition.handleAllMissingSites(all_missing)

        # Ask for a partition report, passing self as the reporter (object that has an output function)
        self.opts.partition.partitionReport(self)

        # If all-missing sites found, need to revise the subset proportions supplied to the
        # cold chain's subset_relrates_move when the cold chain was created
        # BETTER IDEA is to give subset_relrates_move a shared pointer to the partition model object
        #             so we don't have multiple copies of the subset proportions in the first place
        #if partition.getNumSubsets() > 1:
        #    subset_proportions = partition.getSubsetProportions()
        #    cold_chain.subset_relrates_move.setSubsetProportions(subset_proportions)

        if self.opts.doing_steppingstone_sampling:
            # start with posterior (ss_beta = 1) and work toward the prior (ss_beta = 0)
            self.ss_beta = 1.0
            cc = self.mcmc_manager.getColdChain()
            cc.setPower(self.ss_beta)
        self.siteLikeFileSetup(self.mcmc_manager.getColdChain())

    def beyondBurnin(self, cycle):
        c = cycle + 1
        #POLTMP return (c > self.burnin)
        return (c > 0)

    def doThisCycle(self, cycle, burnin, mod):
        if cycle < 0:
            c = cycle - burnin + 1
        else:
            c = cycle + 1
        return ((c % mod) == 0)

    def getModelIndex(self, name):
        """
        If an updater has name like 'k_rAC', this function will return the model index k-1.
        If an updater has name like 'rAC_k', this function will return the model index k-1. <-- deprecated
        If an updater has name 'rAC', returns model index 0. <-- deprecated
        """
        import re
        mo = re.search('^([0-9]+)_',name)
        if mo is not None:
            return int(mo.group(1)) - 1
        else:
            mo = re.search('_([0-9]+)$',name)
            if mo is not None:
                return int(mo.group(1)) - 1
            else:
                return 0

    ############################ exploreWorkingPrior ############################
    def exploreWorkingPrior(self, cycle):   # GENERALIZED-STEPPING-STONE
        chain_index = 0
        cold_chain = self.mcmc_manager.getColdChain()
        tm = phylogeny.TreeManip(cold_chain.tree)
        jpm = cold_chain.chain_manager.getJointPriorManager()

        nmodels = cold_chain.partition_model.getNumSubsets()
        unpartitioned = (nmodels == 1)

        tree_length_prior_specified = (cold_chain.likelihood.getTreeLengthPrior() is not None)

        #new_edge_lens = []
        #new_internal_edge_lens = None
        #new_external_edge_lens = None

        all_updaters = cold_chain.chain_manager.getAllUpdaters()
        edge_lens_need_updating = True
        for u in all_updaters:      # good candidate for moving into C++
            if u.isFixed() or not u.isPriorSteward():
                continue
            name = u.getName()
            if name.find('edgelen_hyper') > -1:
                # draw an edge length hyperparameter value that applies to all edges
                edgelen_hyperparam = u.sampleWorkingPrior()
                m = cold_chain.partition_model.getModel(0)
                cold_chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                jpm.edgeLenHyperparamModified(name, cold_chain.tree, edgelen_hyperparam)
            elif name.find('external_hyper') > -1:
                # draw an edge length hyperparameter value for external edges
                edgelen_hyperparam = u.sampleWorkingPrior()
                cold_chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                m = cold_chain.partition_model.getModel(0)
                self.phycassert(m.isSeparateInternalExternalEdgeLenPriors(), "found updater named 'external_hyper' but isSeparateInternalExternalEdgeLenPriors returns False")
                m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                jpm.edgeLenHyperparamModified(name, cold_chain.tree, edgelen_hyperparam)
            elif name.find('internal_hyper') > -1:
                # draw an edge length hyperparameter value for external edges
                edgelen_hyperparam = u.sampleWorkingPrior()
                cold_chain.chain_manager.setEdgeLenHyperparam(1, edgelen_hyperparam)
                m = cold_chain.partition_model.getModel(0)
                self.phycassert(m.isSeparateInternalExternalEdgeLenPriors(), "found updater named 'internal_hyper' but isSeparateInternalExternalEdgeLenPriors returns False")
                m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                jpm.edgeLenHyperparamModified(name, cold_chain.tree, edgelen_hyperparam)
            elif name.find('subset_relrates') > -1:                 # C++ class SubsetRelRatesMove
                rates_vector = u.sampleMultivariateWorkingPrior()
                cold_chain.partition_model.setSubsetRelRatesVect(rates_vector)
                jpm.multivariateModified(name, rates_vector)
            elif name.find('relrates') > -1:                        # C++ class RelRatesMove
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rate_vector = u.sampleMultivariateWorkingPrior()
                m.setRelRates(rate_vector)
                jpm.multivariateModified(name, rate_vector)
            elif name.find('state_freqs') > -1:                 # C++ class StateFreqMove
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                freq_vector = u.sampleMultivariateWorkingPrior()
                m.setStateFreqsUnnorm(freq_vector)
                jpm.multivariateModified(name, freq_vector)
            elif name.find('gamma_shape') > -1:                     # C++ class DiscreteGammaShapeParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_gamma_shape = u.sampleWorkingPrior()
                m.setShape(new_gamma_shape)
                jpm.univariateModified(name, new_gamma_shape)
            elif name.find('pinvar') > -1:                          # C++ class PinvarParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_pinvar = u.sampleWorkingPrior()
                m.setPinvar(new_pinvar)
                jpm.univariateModified(name, new_pinvar)
            elif name.find('kappa') > -1:                           # C++ class KappaParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_kappa = u.sampleWorkingPrior()
                m.setKappa(new_kappa)
                jpm.univariateModified(name, new_kappa)
            elif name.find('rAC') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rAC = u.sampleWorkingPrior()
                rr[0] = new_rAC
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rAC)
            elif name.find('rAG') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rAG = u.sampleWorkingPrior()
                rr[1] = new_rAG
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rAG)
            elif name.find('rAT') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rAT = u.sampleWorkingPrior()
                rr[2] = new_rAT
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rAT)
            elif name.find('rCG') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rCG = u.sampleWorkingPrior()
                rr[3] = new_rCG
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rCG)
            elif name.find('rCT') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rCT = u.sampleWorkingPrior()
                rr[4] = new_rCT
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rCT)
            elif name.find('rGT') > -1:                         # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                rr = list(m.getRelRates())
                new_rGT = u.sampleWorkingPrior()
                rr[5] = new_rGT
                m.setRelRates(rr)
                jpm.univariateModified(name, new_rGT)
            elif name.find('freqA') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_freqA = u.sampleWorkingPrior()
                m.setStateFreqUnnorm(0, new_freqA)
                jpm.univariateModified(name, new_freqA)
            elif name.find('freqC') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_freqC = u.sampleWorkingPrior()
                m.setStateFreqUnnorm(1, new_freqC)
                jpm.univariateModified(name, new_freqC)
            elif name.find('freqG') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_freqG = u.sampleWorkingPrior()
                m.setStateFreqUnnorm(2, new_freqG)
                jpm.univariateModified(name, new_freqG)
            elif name.find('freqT') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = cold_chain.partition_model.getModel(i)
                new_freqT = u.sampleWorkingPrior()
                m.setStateFreqUnnorm(3, new_freqT)
                jpm.univariateModified(name, new_freqT)
            elif (name.find('extedge') > -1) or (name.find('intedge') > -1):                         # C++ class EdgeLenParam
                self.phycassert(self.opts.fix_topology, 'There is an EdgeLenParam updater but mcmc.fix_topology is False. This is a bug: please report it to paul.lewis@uconn.edu')
                edgelen = u.sampleWorkingPrior()
                u.sendCurrValueToModel(edgelen)
                jpm.univariateModified(name, edgelen)
                edge_lens_need_updating = False
            #elif name.find('master_edgelen') > -1:                  # C++ class EdgeLenMasterParam
            #    num_edge_lens = cold_chain.tree.getNNodes() - 1
            #    new_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
            #elif name.find('external_edgelen') > -1:                # C++ class EdgeLenMasterParam
            #    num_edge_lens = cold_chain.tree.getNTips() - 1
            #    new_external_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
            #elif name.find('internal_edgelen') > -1:                # C++ class EdgeLenMasterParam
            #    num_edge_lens = cold_chain.tree.getNInternals()
            #    new_internal_edge_lens = [u.sampleWorkingPrior() for j in range(num_edge_lens)]
            elif name == 'larget_simon_local':
                # Sample a new tree and edge lengths from the topology reference distribution
                # Note: this really has nothing to do with the Larget-Simon move, but if there
                # is a larget_simon_local updater present, it means that the topology is being
                # updated and is not fixed
                self.phycassert(not self.opts.fix_topology, 'There is a larget_simon_local move updater but mcmc.fix_topology is True. This is a bug: please report it to paul.lewis@uconn.edu')
                u.sampleWorkingPrior()
                edge_lens_need_updating = False
                jpm.allEdgeLensModified(cold_chain.tree)
                jpm.topologyModified("tree_topology", cold_chain.tree)
            else:
                self.phycassert(0, 'model uses an updater (%s) that has not yet been added to MCMCImpl.exploreWorkingPrior (workaround: specify mcmc.draw_directly_from_prior = False)' % name)

        if edge_lens_need_updating and tree_length_prior_specified and cold_chain.likelihood.getTreeLengthRefDist():
            # all edge lengths must be drawn from the tree length working prior
            tm.setRandomEdgeLensFromTreeLengthDist(cold_chain.likelihood.getTreeLengthRefDist())
            jpm.treeLengthModified('tree_length', cold_chain.tree)

        #if new_internal_edge_lens is not None:
        #    # Case of a separate master edge length parameter for internals and tips
        #    self.phycassert(new_external_edge_lens is not None, 'not expecting new_external_edge_lens to be None in MCMCImpl.exploreWorkingPrior')
        #    i = 0   # indexes internals
        #    j = 0   # indexes tips
        #    for nd in cold_chain.tree:
        #        if nd.isRoot():
        #            continue
        #        elif nd.isInternal():
        #            nd.setEdgeLen(new_internal_edge_lens[i])
        #            i += 1
        #        elif nd.isTip():
        #            nd.setEdgeLen(new_external_edge_lens[j])
        #            j += 1
        #        else:
        #            self.phycassert(0, 'nd is neither a tip nor an internal node in MCMCImpl.exploreWorkingPrior')
        #if cold_chain.likelihood.getTreeLengthRefDist():
        #    # all edge lengths must be drawn from the tree length working prior
        #    tm.setRandomEdgeLensFromTreeLengthDist(cold_chain.likelihood.getTreeLengthRefDist())
        #    jpm.treeLengthModified('tree_length', cold_chain.tree)
        #
        #    #tl_ref_dist = cold_chain.likelihood.getTreeLengthRefDist()
        #    #num_taxa = cold_chain.tree.getNTips()
        #    #new_edge_lens = tl_ref_dist.sample(num_taxa)
        #    #new_external_edge_lens = new_edge_lens[:num_taxa]
        #    #new_internal_edge_lens = new_edge_lens[num_taxa:]
        #    #iInternal = 0
        #    #iExternal = 0
        #    #for nd in cold_chain.tree:
        #    #    if nd.isRoot():
        #    #        continue
        #    #    if nd.isTip() or nd.getParent().isTip():
        #    #        nd.setEdgeLen(new_external_edge_lens[iExternal])
        #    #        iExternal += 1
        #    #    else:
        #    #        nd.setEdgeLen(new_internal_edge_lens[iInternal])
        #    #        iInternal += 1
        #elif edge_lens_need_updating:
        #    # fixed topology (and hence each edge has its own length parameter)
        #    self.phycassert(len(new_edge_lens) == cold_chain.tree.getNNodes() - 1, 'new_edge_lens has %d elements but expecting %d in MCMCImpl.exploreWorkingPrior' % (len(new_edge_lens), cold_chain.tree.getNNodes() - 1))
        #    i = 0
        #    for nd in cold_chain.tree:
        #        if nd.isRoot():
        #            continue
        #        nd.setEdgeLen(new_edge_lens[i])
        #        i += 1

        # replace the model
        cold_chain.prepareForLikelihood()
        cold_chain.likelihood.replaceModel(cold_chain.partition_model)

        # recalculate the likelihood
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        cold_chain_manager.refreshLastLnLike()

        #TODO what about polytomies?

        ############################ end exploreWorkingPrior ############################

    def explorePrior(self, cycle):
        chain_index = 0
        chain = self.mcmc_manager.getColdChain()
        tm = phylogeny.TreeManip(chain.tree)
        if self.opts.debugging:
            tmpf = file('debug_info.txt', 'a')
            tmpf.write('************** cycle=%d, chain=%d\n' % (cycle,chain_index))
        edgelens_generated = False

        nmodels = chain.partition_model.getNumSubsets()
        unpartitioned = (nmodels == 1)

        tree_length_prior_specified = (chain.likelihood.getTreeLengthPrior() is not None)

        jpm = chain.chain_manager.getJointPriorManager()
        for p in chain.chain_manager.getAllUpdaters():
            if p.isFixed():
                continue
            w = p.getWeight()
            name = p.getName()
            if (name.find('edgelen_hyper') > -1) or (name.find('external_hyper') > -1) or (name.find('internal_hyper') > -1):   # C++ class HyperPriorParam
                self.phycassert(not tree_length_prior_specified, 'Cannot specify edge length hyperpriors and tree length prior simultaneously')
                if not edgelens_generated:
                    # Choose hyperparam, then use it to choose new edge lengths for a newly-created tree
                    m = chain.partition_model.getModel(0)
                    if m.isSeparateInternalExternalEdgeLenPriors():
                        # draw an edge length hyperparameter value for external edges
                        edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
                        chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                        m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                        jpm.edgeLenHyperparamModified('external_hyper', chain.tree, edgelen_hyperparam)

                        # draw an edge length hyperparameter value for internal edges
                        edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
                        chain.chain_manager.setEdgeLenHyperparam(1, edgelen_hyperparam)
                        m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                        jpm.edgeLenHyperparamModified('internal_hyper', chain.tree, edgelen_hyperparam)
                    else:
                        # draw an edge length hyperparameter value that applies to all edges
                        edgelen_hyperparam = m.getEdgeLenHyperPrior().sample()
                        chain.chain_manager.setEdgeLenHyperparam(0, edgelen_hyperparam)
                        m.getInternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                        m.getExternalEdgeLenPrior().setMeanAndVariance(1.0/edgelen_hyperparam, 0.0) # 2nd arg. (variance) ignored for exponential distributions
                        jpm.edgeLenHyperparamModified('edgelen_hyper', chain.tree, edgelen_hyperparam)
                    if self.opts.fix_topology:
                        tm.setRandomInternalExternalEdgeLengths(m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())
                    else:
                        tm.equiprobTree(chain.tree.getNTips(), chain.r, m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())
                    edgelens_generated = True
            elif name.find('master_edgelen') > -1:
                pass
            elif name.find('external_edgelen') > -1:
                pass
            elif name.find('internal_edgelen') > -1:
                pass
            elif name.find('intedge_') > -1:
                pass
            elif name.find('extedge_') > -1:
                pass
            elif name.find('subset_relrates') > -1:
                new_subset_relrate_vector = chain.partition_model.getSubsetRelRatePrior().sample()
                chain.partition_model.setSubsetRelRatesVect(new_subset_relrate_vector)
                jpm.multivariateModified(name, new_subset_relrate_vector)
            elif name.find('kappa') > -1:                           # C++ class KappaParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_kappa = m.getKappaPrior().sample()
                m.setKappa(new_kappa)
                jpm.univariateModified(name, new_kappa)
            elif name.find('omega') > -1:                           # C++ class OmegaParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_omega = m.getOmegaPrior().sample()
                m.setOmega(new_omega)
                jpm.univariateModified(name, new_omega)
            elif name.find('rAC') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rAC = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(0, new_rAC)
                jpm.univariateModified(name, new_rAC)
            elif name.find('rAG') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rAG = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(1, new_rAG)
                jpm.univariateModified(name, new_rAG)
            elif name.find('rAT') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rAT = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(2, new_rAT)
                jpm.univariateModified(name, new_rAT)
            elif name.find('rCG') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rCG = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(3, new_rCG)
                jpm.univariateModified(name, new_rCG)
            elif name.find('rCT') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rCT = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(4, new_rCT)
                jpm.univariateModified(name, new_rCT)
            elif name.find('rGT') > -1:                             # C++ class GTRRateParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_rGT = m.getRelRateParamPrior().sample()
                m.setRelRateUnnorm(5, new_rGT)
                jpm.univariateModified(name, new_rGT)
            elif name.find('freqA') > -1:                           # C++ class StateFreqParam
                new_freq_param_A = m.getStateFreqParamPrior().sample()
                m.setStateFreqParam(0, new_freq_param_A)
                jpm.univariateModified(name, new_freq_param_A)
            elif name.find('freqC') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_freq_param_C = m.getStateFreqParamPrior().sample()
                m.setStateFreqParam(1, new_freq_param_C)
                jpm.univariateModified(name, new_freq_param_C)
            elif name.find('freqG') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_freq_param_G = m.getStateFreqParamPrior().sample()
                m.setStateFreqParam(2, new_freq_param_G)
                jpm.univariateModified(name, new_freq_param_G)
            elif name.find('freqT') > -1:                           # C++ class StateFreqParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_freq_param_T = m.getStateFreqParamPrior().sample()
                m.setStateFreqParam(3, new_freq_param_T)
                jpm.univariateModified(name, new_freq_param_T)
            elif name.find('gamma_shape') > -1:                     # C++ class DiscreteGammaShapeParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_gamma_shape = m.getDiscreteGammaShapePrior().sample()
                m.setShape(new_gamma_shape)
                jpm.univariateModified(name, new_gamma_shape)
            elif name.find('pinvar') > -1:                          # C++ class PinvarParam
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                new_pinvar = m.getPinvarPrior().sample()
                m.setPinvar(new_pinvar)
                jpm.univariateModified(name, new_pinvar)
            elif name.find('relrates') > -1:                        # C++ class RelRatesMove
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                rate_vector = m.getRelRatePrior().sample()
                # Drawing values from a Dirichlet prior, but relative rates should have mean 1, not sum 1,
                # so multiply each by the value nmodels to correct this.
                m.setRelRates([6.0*x for x in rate_vector])
                jpm.multivariateModified(name, rate_vector)
            elif name.find('state_freqs') > -1:                     # C++ class StateFreqMove
                i = unpartitioned and 0 or self.getModelIndex(name)
                m = chain.partition_model.getModel(i)
                freq_vector = m.getStateFreqPrior().sample()
                m.setStateFreqsUnnorm(freq_vector)
                jpm.multivariateModified(name, freq_vector)
            elif name.find('edge_move') > -1:                       # C++ class EdgeMove
                pass
            elif name.find('larget_simon_local') > -1:              # C++ class LargetSimonMove
                pass
            elif name.find('tree_scaler') > -1:                     # C++ class TreeScalerMove
                pass
            elif name.find('bush_move') > -1:                       # C++ class BushMove
                pass    # polytomies handled further down (by randomly pruning fully-resolved equiprobable tree)
            else:
                self.phycassert(0, 'model uses an updater (%s) that has not yet been added to MCMCImpl.explorePrior (workaround: specify mcmc.draw_directly_from_prior = False)' % name)

        # If no edge length hyperprior was specified, then build the tree with edge lengths now

        if tree_length_prior_specified:
            if not edgelens_generated:
                m = chain.partition_model.getModel(0)
                if self.opts.fix_topology:
                    tm.setRandomEdgeLensFromTreeLengthDist(chain.likelihood.getTreeLengthPrior())
                else:
                    tm.buildEquiprobTree(chain.tree.getNTips(), chain.r)

            if self.opts.allow_polytomies and not self.opts.fix_topology:
                # Choose number of internal nodes
                jpm = self.mcmc_manager.getColdChainManager().getJointPriorManager()
                topo_prior_calculator = jpm.getTopoProbCalculator()
                num_internal_nodes = topo_prior_calculator.sample(chain.r)

                # Delete edges at random from tree to achieve chosen number of internal nodes
                orig_num_internal_nodes = chain.tree.getNInternals()
                num_internals_to_delete = orig_num_internal_nodes - num_internal_nodes
                for i in range(num_internals_to_delete):
                    tm.deleteRandomInternalEdge(chain.r)

            if not self.opts.fix_topology:
                tm.setRandomEdgeLensFromTreeLengthDist(chain.likelihood.getTreeLengthPrior())

            jpm.treeLengthModified('tree_length', chain.tree)
        else:
            if not edgelens_generated:
                m = chain.partition_model.getModel(0)
                if self.opts.fix_topology:
                    tm.setRandomInternalExternalEdgeLengths(m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())
                else:
                    tm.equiprobTree(chain.tree.getNTips(), chain.r, m.getInternalEdgeLenPrior(), m.getExternalEdgeLenPrior())

            if self.opts.allow_polytomies:
                # Choose number of internal nodes
                jpm = self.mcmc_manager.getColdChainManager().getJointPriorManager()
                topo_prior_calculator = jpm.getTopoProbCalculator()
                num_internal_nodes = topo_prior_calculator.sample(chain.r)

                # Delete nodes at random from tree to achieve chosen number of internal nodes
                orig_num_internal_nodes = chain.tree.getNInternals()
                num_internals_to_delete = orig_num_internal_nodes - num_internal_nodes
                for i in range(num_internals_to_delete):
                    tm.deleteRandomInternalEdge(chain.r)

            jpm.allEdgeLensModified(chain.tree)
            #jpm.externalEdgeLensModified('external_edgelen', chain.tree)
            #jpm.internalEdgeLensModified('internal_edgelen', chain.tree)

        if not self.opts.fix_topology:
            jpm.topologyModified('tree_topology', chain.tree)

        chain.prepareForLikelihood()
        chain.likelihood.replaceModel(chain.partition_model)

        if False:
            # debugging code
            chain.likelihood.storeSiteLikelihoods(True)
            from phycas.Utilities.kappa2tratio import convert
            f = chain.model.getStateFreqs()
            k = convert(chain.model.getKappa(), f[0], f[1], f[2], f[3])
            print 'cycle = %d, model = %s' % (cycle + 1, chain.model.getModelName())
            print '  lset tratio=%.5f basefreq=(%.5f %.5f %.5f) rates=gamma ncat=4 shape=%.5f;' % (k, f[0], f[1], f[2], chain.model.getShape())
            print 'taxon names:', self.opts.data_source.taxon_labels
            chain.tree.rectifyNames(self.opts.data_source.taxon_labels)
            print '  utree one = %s;' % chain.tree.makeNewick()
            print '  sum of edge lengths = %.5f' % chain.tree.edgeLenSum()

        # recalculate the likelihood
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        cold_chain_manager.refreshLastLnLike()

        if False:
            # debugging code
            counts = chain.likelihood.getPatternCounts()
            sitelikes = chain.likelihood.getSiteLikelihoods()
            print '  lnL = %.6f' % cold_chain_manager.getLastLnLike()
            sumlikes = 0.0
            for sitelike,count in zip(sitelikes, counts):
                if count > 100:
                    print '%6.0f -> %12.5f' % (count, sitelike)
                sumlikes += count*sitelike

        if self.opts.debugging:
            tmpf.close()

    #POLTMP def computeTimeRemaining(self, secs, ndone, ntotal):
    def computeTimeRemaining(self, secs, cycle_start, cycle_stop, cycle):
        if self.opts.doing_steppingstone_sampling:
            num_betas_completed = self.ss_sampled_betas.index(self.ss_beta)
            num_betas_total = len(self.ss_sampled_betas)
        else:
            num_betas_completed = 0
            num_betas_total = 1

        ndone  = (cycle_stop - cycle_start)*num_betas_completed + (cycle - cycle_start)
        ntotal = (cycle_stop - cycle_start)*num_betas_total

        if ndone < 1:
            return ''

        days_remaining = 0
        hours_remaining = 0
        secs_remaining = float(secs)*(float(ntotal)/float(ndone) - 1.0)

        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        #print '~~~ secs                =',secs
        #print '~~~ cycle_start         =',cycle_start
        #print '~~~ cycle_stop          =',cycle_stop
        #print '~~~ cycle               =',cycle
        #print '~~~ num_betas_completed =',num_betas_completed
        #print '~~~ num_betas_total     =',num_betas_total
        #print '~~~ ndone               =',ndone
        #print '~~~ ntotal              =',ntotal
        #print '~~~ ntotal/ndone        =',(float(ntotal)/float(ndone))
        #print '~~~ secs_remaining      =',secs_remaining
        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

        time_left = []
        if secs_remaining > 86400.0:
            days_remaining = math.floor(secs_remaining/86400.0)
            secs_remaining -= 86400.0*days_remaining
            if days_remaining > 0:
                if days_remaining == 1:
                    time_left.append('1 day')
                else:
                    time_left.append('%d days' % days_remaining)
        if secs_remaining > 3600.0:
            hours_remaining = math.floor(secs_remaining/3600.0)
            secs_remaining -= 3600.0*hours_remaining
            if hours_remaining > 0:
                if hours_remaining == 1:
                    time_left.append('1 hour')
                else:
                    time_left.append('%d hours' % hours_remaining)
        if secs_remaining > 60.0:
            minutes_remaining = math.floor(secs_remaining/60.0)
            secs_remaining -= 60.0*minutes_remaining
            if minutes_remaining > 0:
                if minutes_remaining == 1 and (days_remaining + hours_remaining) == 0:
                    time_left.append('less than 2 minutes')
                else:
                    time_left.append('%d minutes' % minutes_remaining)
        if len(time_left) > 0:
            return ', '.join(time_left) + ' remaining'
        elif math.floor(secs_remaining) == 1:
            return '1 second remaining'
        else:
            return '%d seconds remaining' % math.floor(secs_remaining)

    def mainMCMCLoop(self, explore_prior = False):
        levels_file_created = False #temp!
        nchains = len(self.mcmc_manager.chains)
        # print '******** nchains =',nchains
        self.last_adaptation = 0
        #POLTMP self.next_adaptation = self.opts.adapt_first

        CPP_UPDATER = True # if False, uses python obsoleteUpdateAllUpdaters

        #POLTMP for cycle in xrange(self.burnin + self.ncycles):
        cycle = self.cycle_start
        done = False
        burning_in = True
        while not done:

            if cycle == 0:
                # switching from burnin to sampling, tell chains to stop adapting themselves
                burning_in = False
                for c in self.mcmc_manager.chains:
                    c.chain_manager.setAdapting(False) #POLTMP
            assert cycle < 0 or burning_in == False, 'should not be adapting chains during sampling phases'

            # Update all updaters
            if explore_prior and self.opts.draw_directly_from_prior:
                if self.opts.ss_heating_likelihood or not self.opts.doing_steppingstone_sampling:
                    # MCMC without data, TI, or specific SS
                    self.explorePrior(cycle)
                else:
                    # generalized SS
                    self.exploreWorkingPrior(cycle)
            else:
                for i,c in enumerate(self.mcmc_manager.chains):
                    if CPP_UPDATER:
                        c.chain_manager.updateAllUpdaters()
                    else:
                        self.obsoleteUpdateAllUpdaters(c, i, cycle)

            # Attempt to swap two random chains
            if nchains > 1:
                self.mcmc_manager.attemptChainSwap(cycle)

            # Provide progress report to user if it is time
            if self.opts.verbose and self.doThisCycle(cycle, self.burnin, self.opts.report_every):
                # Refresh log-likelihood of cold chain if necessary
                if self.ss_beta == 0.0:
                    self.mcmc_manager.getColdChainManager().refreshLastLnLike()

                self.stopwatch.normalize()
                secs = self.stopwatch.elapsedSeconds()
                #POLTMP time_remaining = self.computeTimeRemaining(secs, self.cycle_start + cycle + 1, self.cycle_stop)
                time_remaining = self.computeTimeRemaining(secs, self.cycle_start, self.cycle_stop, cycle)
                if time_remaining != '':
                    time_remaining = '(' + time_remaining + ')'
                if self.opts.doing_steppingstone_sampling:
                    cold_chain_manager = self.mcmc_manager.getColdChainManager()
                    msg = 'beta = %.5f, cycle = %d, lnL = %.5f %s' % (self.ss_beta, cycle + 1, cold_chain_manager.getLastLnLike(), time_remaining)
                else:
                    if nchains == 1:
                        cold_chain_manager = self.mcmc_manager.getColdChainManager()
                        msg = 'cycle = %d, lnL = %.5f %s' % (cycle + 1, cold_chain_manager.getLastLnLike(), time_remaining)
                    else:
                        msg = 'cycle = %d, ' % (cycle + 1)
                        for k in range(nchains):
                            c = self.mcmc_manager.chains[k]
                            msg += 'lnL(%.3f) = %.5f, ' % (c.heating_power, c.chain_manager.getLastLnLike())
                        msg += '%s' % time_remaining
                self.output(msg)

            # Sample chain if it is time
            #POLTMP if self.beyondBurnin(cycle) and self.doThisCycle(cycle - self.burnin, self.opts.sample_every):
            if self.beyondBurnin(cycle) and self.doThisCycle(cycle, self.burnin, self.opts.sample_every):
                # Refresh log-likelihood(s) if necessary
                if self.ss_beta == 0.0:
                    for i,c in enumerate(self.mcmc_manager.chains):
                        # is this necessary?
                        c.chain_manager.refreshLastLnLike()

                #POLTMP self.mcmc_manager.recordSample(self.cycle_start + cycle)
                self.mcmc_manager.recordSample(cycle)

                cold_chain_manager = self.mcmc_manager.getColdChainManager()
                sampled_lnL = cold_chain_manager.getLastLnLike()
                self.ss_sampled_likes[self.ss_beta_index].append(sampled_lnL)
                self.stopwatch.normalize()

            # Adapt slice samplers if it is time
            #POLTMP if self.doThisCycle(cycle, self.burnin, self.next_adaptation):
            #POLTMP     self.adaptSliceSamplers()
            #POLTMP     self.next_adaptation += 2*(self.next_adaptation - self.last_adaptation)
            #POLTMP     self.last_adaptation = cycle + 1
            if cycle < 0:
                self.adaptSliceSamplers()
            if self.doThisCycle(cycle, self.burnin, self.opts.report_efficiency_every):
                self.reportUpdaterEfficiency()

            # Recalculate joint prior to avoid creeping round-off error
            jpm = self.mcmc_manager.getColdChainManager().getJointPriorManager()
            jpm.recalcLogJointPrior()

            cycle += 1  #POLTMP
            if cycle == self.cycle_stop:
                done = True

        #POLTMP self.cycle_start += self.burnin + self.ncycles
        #POLTMP self.cycle_start = self.cycle_stop
        #POLTMP self.cycle_stop = self.cycle_start + self.ncycles


    #def _isTreeLengthPriorBeingUsed(self):
    #    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    #    """
    #    Returns True if tree length prior is defined, False otherwise.
    #    """
    #    cold_chain = self.mcmc_manager.getColdChain()
    #    m0 = cold_chain.partition_model.getModel(0)
    #    return bool(m0.getTreeLengthPrior())

    def _createRefDistMap(self, fn):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Takes a file name and from the contents of the
        file creates a focal tree and an associative array that matches
        reference distributions (values) with parameter names (keys). The
        first line is a tree definition with embedded clade posteriors. Lines
        following the first are reference distribution definitions, one per
        line. The file fn looks like this:

        (1:1.00000,3:1.00000,((2:1.00000,(8:1.00000,5:1.00000):1.00000):1.00000,((9:1.00000,6:1.00000):1.00000,(4:1.00000,(10:1.00000,7:1.00000):1.00000):1.00000):1.00000):1.00000)
        split_-----*--*- = Gamma(10.6102,0.00119215)
        split_-********* = Gamma(154.718,0.000845338)
        split_-----*---- = Gamma(112.402,0.000775574)
        split_------*--- = Gamma(74.0985,0.000784768)
        split_-------*-- = Gamma(138.733,0.000884312)
        split_------*--* = Gamma(29.9702,0.00115743)
        split_--*------- = Gamma(87.8519,0.000947885)
        split_--------*- = Gamma(101.101,0.000771077)
        split_-*--*--*-- = Gamma(2.27518,0.00139715)
        split_---------* = Gamma(87.2539,0.00083153)
        split_---*------ = Gamma(125.99,0.000793145)
        split_----*----- = Gamma(76.5084,0.000925474)
        split_-*-******* = Gamma(14.437,0.00150977)
        split_---*-**-** = Gamma(24.95,0.00127126)
        split_---*--*--* = Gamma(0.977012,0.00117438)
        split_----*--*-- = Gamma(32.8486,0.00113388)
        split_-*-------- = Gamma(106.928,0.000846824)
        split_NA = Gamma(2.37343,0.0256837)
        tree_length = TreeLengthDist(1226.6,1183.64,12.6972,0.227628)
        1_state_freqs = Dirichlet((559.71489, 531.95324, 507.49159, 685.90371))
        1_kappa = Gamma(249.30762, 0.01510)
        """
        cold_chain = self.mcmc_manager.getColdChain()
        m0 = cold_chain.partition_model.getModel(0)
        is_tree_length_prior = bool(m0.getTreeLengthPrior())

        ref_dist_map = {}
        lines = open(fn, 'r').readlines()

        for line in lines[1:]:
            stripped = line.strip()
            if len(stripped) == 0 or stripped[0] == '#':
                continue    # skip blank lines and comment lines
            k,v = line.split('=')
            #print '_createRefDistMap: k = |%s|' % k.strip()
            if k.strip() == 'tree_length' and not is_tree_length_prior:
                # Don't create a reference distribution for tree length if not using the tree length prior
                pass
            else:
                # Create reference distribution object and add to map
                name = k.strip()
                dist = eval(v.strip())
                ref_dist_map[name] = dist

        ref_dist_tree_newick = lines[0]    # this is just the tree definition string, not a tree object
        return (ref_dist_tree_newick,ref_dist_map)

    def _addRefDistsToJointPriorManager(self, ref_dist_map):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Adds all reference distributions stored in supplied ref_dist_map to
        the joint prior manager.
        """
        cold_chain = self.mcmc_manager.getColdChain()
        jpm = cold_chain.chain_manager.getJointPriorManager()
        for k in ref_dist_map.keys():
            key_name = k + '_ref'
            distr = ref_dist_map[k]
            distr_name = distr.__class__.__name__
            if distr_name in ['Normal', 'InverseGamma', 'Gamma', 'Beta', 'BetaPrime']:
                jpm.addUnivariateDistribution(key_name, distr, distr.sample())
            elif distr_name in ['Dirichlet', 'RelativeRateDistribution']:
                jpm.addMultivariateDistribution(key_name, distr, distr.sample())
            elif distr_name == 'TreeLengthDist':
                jpm.addTreeLengthDistribution(key_name, distr, cold_chain.tree)
            else:
                self.phycassert(False, 'MCMCImpl._addRefDistsToJointPriorManager: Need to add "%s" case' % distr_name)

    def _loadReferenceDistribution(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        User has specified a file containing the reference distribution
        definitions. Read this file and call _createRefDistMap to instantiate
        and store the reference distributions. The reference tree newick
        string is parsed and a reference tree created. Edge lengths in the
        newick tree description are assumed to be posterior clade
        probabilities.
        """
        ref_dist_tree,ref_dist_map = self._createRefDistMap(self.opts.ssobj.refdistfile)
        #self._addRefDistsToJointPriorManager(ref_dist_map)

        # Build the reference tree from the description in ref_dist_tree
        focal_tree = TreeCollection(newick=ref_dist_tree).trees[0]
        ntips = focal_tree.getNObservables()
        focal_tree.recalcAllSplits(ntips)
        nd = focal_tree.getFirstPreorder()
        assert nd.isRoot(), 'the first preorder node should be the root'
        focal_root_number = nd.getNodeNumber()
        self.phycassert(focal_tree.hasEdgeLens(), 'focal tree from reference distribution must have edge lengths (which will be interpreted as split posteriors)')
        nd = nd.getNextPreorder()
        while nd:
            # Determine whether this split represents an internal or tip node
            if not (nd.isTip() or nd.getParent().isRoot()):
                split_prob = nd.getEdgeLen()
                if split_prob == 0.0:
                    split_prob = 0.001
                    nd.setEdgeLen(split_prob)
                if split_prob == 1.0:
                    split_prob = 0.999
                    nd.setEdgeLen(split_prob)
                self.phycassert(split_prob > 0.0 and split_prob < 1.0, 'Split probabilities must be in the range (0, 1)')
                s = nd.getSplit()
                if s.isBitSet(0):
                    s.invertSplit()
            nd = nd.getNextPreorder()

        # Be sure trees for all chains are rooted with the same taxon as the reference tree
        # but only if fix_topology is False (if fix_topology is True, the reference tree is
        # irrelevant and rerooting the tree will cause problems because edge-specific updaters
        # have already been created based on the current rooting)
        if not self.opts.fix_topology:
            for i,chain in enumerate(self.mcmc_manager.chains):
                t = chain.getTree()
                t.rerootAtTip(focal_root_number)

        # Instantiate tree topology probability calculator, which is used to generate trees
        # from the topology reference distribution
        topo_ref_dist_calculator = likelihood.FocalTreeTopoProbCalculatorBase(focal_tree)

        # Inform tree topology probability calculator of edge length reference distributions
        prefix = 'split_'
        default_edge_len = None
        for k, v in ref_dist_map.iteritems():
            if k.lower().startswith(prefix):
                split_rep = k[len(prefix):] # get back end of name (part after prefix)
                v.setLot(self._getLot())
                if split_rep.lower() == 'na':
                    assert(default_edge_len is None)
                    topo_ref_dist_calculator.setDefaultEdgeLenDist(v)
                else:
                    s = phylogeny.SplitBase()
                    s.createFromPattern(split_rep)
                    #print split_rep, v, s.createPatternRepresentation()
                    topo_ref_dist_calculator.setEdgeLenDist(s, v)
        return topo_ref_dist_calculator, ref_dist_map

    def _provideRefDistToUpdaters(self, cold_chain, topo_ref_dist_calculator, ref_dist_map):
        #print '#@#@#@#@#@#@# ref_dist_map.keys() -->',ref_dist_map.keys()
        if self.opts.fix_topology and 'tree_length' in ref_dist_map.keys():
            cold_chain.likelihood.setTreeLengthRefDist(ref_dist_map['tree_length'].cloneAndSetLot(self._getLot()))
        all_updaters = cold_chain.chain_manager.getAllUpdaters()
        for u in all_updaters:
            if u.isFixed():
                continue
            else:
                u.setUseWorkingPrior(True)
                updater_name = u.getName()
                #print '~o~o~o~o~o~o~o~o~ updater_name =',updater_name
                if updater_name in ['master_edgelen','external_edgelen','internal_edgelen'] and not self.opts.fix_topology:
                    # larget_simon_local takes over this role for variable topology steppingstone analyses
                    u.setUseWorkingPrior(False)
                    u.setPrior(None);
                    u.setMultivarPrior(None);
                if u.useWorkingPrior():
                    if u.computesUnivariatePrior():
                        try:
                            d = ref_dist_map[updater_name].cloneAndSetLot(self._getLot())
                        except KeyError:
                            # Probably an edge length updater with name like 'intedge_1001' (whereas keys in ref_dist_map have names like 'split_-**-')
                            self.phycassert(u.__class__.__name__ == "EdgeLenParam", 'Expecting updater to be an EdgeLenParam. This is a bug: please report this to Paul Lewis (paul.lewis@uconn.edu)')
                            split_key = 'split_%s' % u.getSplitReprAsString()
                            d = ref_dist_map[split_key].cloneAndSetLot(self._getLot())
                        u.setWorkingPrior(d)
                        #self.output('  %s = %s' % (updater_name, u.getWorkingPriorDescr()))
                    elif u.computesMultivariatePrior():
                        u.setMultivariateWorkingPrior(ref_dist_map[updater_name].cloneAndSetLot(self._getLot()))
                        #self.output('  %s = %s' % (updater_name, u.getWorkingPriorDescr()))
                    elif u.computesTopologyPrior():
                        u.setReferenceDistribution(topo_ref_dist_calculator)

    def steppingstoneMCMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a Steppingstone MCMC analysis.

        """
        nchains = len(self.mcmc_manager.chains)
        self.phycassert(self.data_matrix is not None, 'steppingstone sampling requires data')
        self.phycassert(nchains == 1, 'steppingstone sampling currently requires nchains = 1')
        chain = self.mcmc_manager.getColdChain()

        if not self.opts.ssobj.refdist_is_prior:
            topo_ref_dist_calculator, ref_dist_map = self._loadReferenceDistribution()
            self._provideRefDistToUpdaters(chain, topo_ref_dist_calculator, ref_dist_map)

        # Compute the current log-likelihood and log-prior in case first updater
        # is a move and will thus depend on these quantities being accurate
        for c in self.mcmc_manager.chains:
            c.chain_manager.setTargetAcceptanceRate(self.opts.target_accept_rate) #POLTMP
            c.chain_manager.setAdapting(True) #POLTMP
            c.chain_manager.refreshLastLnLike()
            if c.heating_power == 1.0:
                self.output('Starting log-likelihood = %s' % c.chain_manager.getLastLnLike())
                self.output('Starting log-prior = %s' % c.chain_manager.getJointPriorManager().getLogJointPrior())

        self.mcmc_manager.recordSample()

        if self.opts.ssobj.nstones > 1:
            # Set up the list ss_sampled_betas
            # Beta distribution will be divided into nstones intervals, each of which has an equal area
            #
            # Example: nstones = 5, Beta(a,1), a = 1/4
            #
            #    quantile    ------- beta -------
            #       0/5      (0/5)^(1/a) = 0.0     <-- 0
            #       1/5      (1/5)^(1/a) = 0.0016  <-- 1
            #       2/5      (2/5)^(1/a) = 0.0256  <-- 2
            #       3/5      (3/5)^(1/a) = 0.1296  <-- 3
            #       4/5      (4/5)^(1/a) = 0.4096  <-- 4
            #       5/5      (5/5)^(1/a) = 1.0     <-- not used
            #
            segment_area = 1.0/float(self.opts.ssobj.nstones)
            cum_area = 0.0
            lower_boundary = 0.0
            self.ss_sampled_betas = [0.0]
            betadist = probdist.Beta(self.opts.ssobj.shape1, self.opts.ssobj.shape2)
            for i in range(self.opts.ssobj.nstones - 1):
                cum_area += segment_area
                upper_boundary = betadist.getQuantile(cum_area)
                self.ss_sampled_betas.append(upper_boundary)
                lower_boundary = upper_boundary

            # Reverse ss_sampled_betas so that sampled beta values start high and decrease toward 0.0
            self.ss_sampled_betas.reverse()

            # Output the beta values that will be used
            self.output('%d %s chosen from a discrete\nBeta(%.5f, %.5f) distribution:' % (self.opts.ssobj.nstones, (self.opts.ssobj.nstones == 1 and 'value was' or 'values were'), self.opts.ssobj.shape1, self.opts.ssobj.shape2))
            for i,x in enumerate(self.ss_sampled_betas):
                self.output('%6d %12.5f' % (i+1,x))
            self.output('MCMC will be used to sample from each of the')
            self.output('power posteriors defined by these values.')
            self.output()
        else:
            self.ss_sampled_betas = [0.0]

        # Run the main MCMC loop for each beta value in ss_sampled_betas
        self.ss_sampled_likes = []
        ref_dist_calculated = False
        for self.ss_beta_index, self.ss_beta in enumerate(self.ss_sampled_betas):
            self.ss_sampled_likes.append([])
            chain.setPower(self.ss_beta)
            #POLTMP boldness = 100.0*(1.0 - self.ss_beta)
            #POLTMP chain.setBoldness(boldness)
            #POLTMP self.output('Setting chain boldness to %g based on beta = %g' % (boldness,self.ss_beta))
            #POLTMP self.cycle_stop = self.opts.burnin + len(self.ss_sampled_betas)*self.opts.ssobj.ncycles
            self.ncycles = self.opts.ssobj.ncycles
            self.burnin = self.opts.ssobj.burnin

            chain.chain_manager.setTargetAcceptanceRate(self.opts.target_accept_rate) #POLTMP
            chain.chain_manager.setAdapting(True) #POLTMP

            #POLTMP if self.ss_beta_index > 0:
            #POLTMP     self.burnin = 0
            #POLTMP else:
            #POLTMP     self.burnin = self.opts.burnin
            #POLTMP     self.cycle_start = 0
            #self.burnin = self.opts.burnin
            self.cycle_start = -self.burnin
            #POLTMP self.cycle_stop = len(self.ss_sampled_betas)*self.opts.ssobj.ncycles
            self.cycle_stop = self.opts.ssobj.ncycles

            if self.ss_beta == 0.0:
                self.mainMCMCLoop(explore_prior=True)
            else:
                self.mainMCMCLoop()

    def standardMCMC(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs a standard MCMC analysis.

        """
        # Compute the current log-likelihood and log-prior in case first updater
        # is a move and will thus depend on these quantities being accurate
        for c in self.mcmc_manager.chains:
            c.chain_manager.setTargetAcceptanceRate(self.opts.target_accept_rate) #POLTMP
            c.chain_manager.setAdapting(True) #POLTMP
            c.chain_manager.refreshLastLnLike()
            if c.heating_power == 1.0:
                self.output('Starting log-likelihood = %s' % c.chain_manager.getLastLnLike())
                self.output('Starting log-prior = %s' % c.chain_manager.getJointPriorManager().getLogJointPrior())

        self.mcmc_manager.recordSample()
        self.ss_sampled_likes = []
        self.ss_sampled_likes.append([])
        self.ss_beta_index = 0

        #POLTMP self.cycle_start = 0
        #POLTMP self.cycle_stop = self.opts.burnin + self.opts.ncycles
        self.cycle_start = -self.opts.burnin
        self.cycle_stop = self.opts.ncycles

        self.burnin = self.opts.burnin
        self.ncycles = self.opts.ncycles

        if self.data_matrix is None:
            self.mainMCMCLoop(explore_prior=True)
        else:
            self.mainMCMCLoop()

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs the MCMC analysis.

        """
        self.setup()

        # If user has set quiet to True, then phycas.output calls will have no effect
        self.quiet = self.opts.quiet

        nchains = len(self.mcmc_manager.chains)

        cold_chain = self.mcmc_manager.getColdChain()

        if self.opts.verbose:
            if self.data_matrix == None:
                self.output('Data source:    None (running MCMC with no data to explore prior)')
            else:
                self.output('Data source:    %s' % str_value_for_user(self.opts.data_source))
                self.output('  No. taxa:                %d' % self.ntax)
                self.output('  No. included characters: %d' % cold_chain.likelihood.sumPatternCounts())
                all_missing = cold_chain.likelihood.getListOfAllMissingSites()
                num_excluded = len(all_missing)
                if num_excluded > 0:
                    self.output('  *** Note: the following %d sites were automatically excluded because' % num_excluded)
                    self.output('  *** they exhibited completely missing data for all taxa:')
                    while len(all_missing) > 0:
                        tmp = all_missing[:10]
                        all_missing = all_missing[10:]
                        self.output('  ***   '+','.join([str(i+1) for i in tmp]))

            if nchains > 1:
                for c in range(nchains):
                    self.output('Starting tree for chain %d:  %s' % (c, self.starting_tree[c]))
            else:
                self.output('Starting tree:  %s' % str(self.starting_tree[0]))

            if self.opts.fix_topology:
                self.output('\nTree topology fixed.\n')
            else:
                self.output('\nTree topology free to vary.\n')

            if self.opts.doing_steppingstone_sampling:
                nsamples_per_step = self.opts.ssobj.ncycles//self.opts.sample_every
                if self.opts.ss_heating_likelihood:
                    self.output('\nPerforming steppingstone sampling and thermodynamic integration to estimate marginal likelihood.')
                else:
                    self.output('\nPerforming generalized steppingstone sampling to estimate marginal likelihood.')
                self.output('Likelihood will be raised to the power beta, and beta will be')
                self.output('decremented from 1.0 to 0.0 in a series of steps.')
                self.output('  No. stones:              %s' % self.opts.ssobj.nstones)
                self.output('  No. cycles per step:     %s' % self.opts.ssobj.ncycles)
                self.output('  Sample every:            %s' % self.opts.sample_every)
                self.output('  No. samples per step:    %s' % nsamples_per_step)
                self.output('\n')
            else:
                nsamples_per_step = self.opts.ncycles//self.opts.sample_every
                self.output('No. cycles:     %s' % self.opts.ncycles)
                self.output('Sample every:   %s' % self.opts.sample_every)
                self.output('No. samples:    %s' % nsamples_per_step)
            self.output('Sampled trees will be saved in %s' % str_value_for_user(self.opts.out.trees))
            self.output('Sampled parameters will be saved in %s' % str_value_for_user(self.opts.out.params))
            self.output('Using standard MCMC (i.e. no uniformized mapping)')

            if not self.warn_tip_numbers:
                self.output('Tip node numbers were set using the names in the tree description')
            else:
                self.output('Warning: tip node numbers were NOT set using the names in the tree description')

        if nchains == 1:
            self.output('Creating one chain (i.e. not using heated chains to improve mixing)')
        else:
            self.output('Creating %d chains with these temperatures:' % (nchains))
            for t in self.heat_vector:
                self.output('  %.5f %s' % (t, t == 1.0 and '(cold chain)' or ''))

        # Show starting parameter info
        self.output('\nParameter starting values and prior densities:')
        cold_chain_manager = self.mcmc_manager.getColdChainManager()
        for p in cold_chain_manager.getAllUpdaters():
            self.showParamInfo(p)

        # Show updater names
        self.output('\nHere is a list of all updaters that will be used for this analysis:')
        for p in cold_chain_manager.getAllUpdaters():
            if p.getWeight() > 0:
                p.setUseWorkingPrior(False)
                if p.isMove():
                    self.output('  %s (Metropolis-Hastings)' % p.getName())
                elif p.hasSliceSampler():
                    self.output('  %s (slice sampler)' % p.getName())
                else:
                    self.output('  %s (computes prior but does not update)' % p.getName())

        # Debugging: show data patterns
        if self.opts.debugging and not cold_chain.likelihood.getNoData():
            #cold_chain = self.mcmc_manager.getColdChain()
            s = cold_chain.likelihood.listPatterns()
            print '\nDebug Info: List of data patterns and their frequencies (see TreeLikelihood::listPatterns):'
            print s

            cold_chain.likelihood.debugUncompressedDataInfo("all-site-patterns.txt");

        # Show information about topology prior to be used
        self.showTopoPriorInfo()

        self.stopwatch.start()
        self.mcmc_manager.resetNumLikelihoodEvals()

        if self.opts.doing_steppingstone_sampling:
            self.output('\nSampling (%d cycles for each of the %d values of beta)...' % (self.opts.ssobj.ncycles, self.opts.ssobj.nstones))
        else:
            self.output('\nSampling (%d cycles)...' % self.opts.ncycles)
        if self.opts.verbose:
            print

        # Lay down first line in params file (recorded as cycle 0) containing starting values of parameters
        if self.opts.doing_steppingstone_sampling:
            self.steppingstoneMCMC()
        else:
            self.standardMCMC()

        #POLTMP self.adaptSliceSamplers()
        self.resetUpdaterDiagnostics()
        total_evals = self.mcmc_manager.getTotalEvals() #self.likelihood.getNumLikelihoodEvals()
        total_secs = self.stopwatch.elapsedSeconds()
        self.output('%d likelihood evaluations in %.5f seconds' % (total_evals, total_secs))
        if (total_secs > 0.0):
            self.output('  = %.5f likelihood evaluations/sec' % (total_evals/total_secs))

        if self.treef:
            self.treeFileClose()
        if self.paramf:
            self.paramFileClose()


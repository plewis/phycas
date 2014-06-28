import os,sys,math,random
from phycas import *
#from MCMCManager import MCMCManager
from MCMCManager import LikelihoodCore
from phycas.utilities.PhycasCommand import *
from phycas.readnexus import NexusReader
from phycas.utilities.CommonFunctions import CommonFunctions

class GGImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Computes Gelfand-Ghosh on a pre-existing MCMC sample.

    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Gathers information relevant to Gelfand-Ghosh calculations from the
        parent GG object.

        """
        CommonFunctions.__init__(self, opts)

        self.nsamples                   = 0
        self.last_error                 = None

        self.bincountsf                 = None  # filename used for bin counts file
        self.postpredf                  = None  # filename used for post pred. files

        self.starting_tree_source       = None
        self.gg_bin_patterns            = opts.gg_bin_patterns
        self.gg_num_bins                = opts.gg_num_bins
        self.pfile                      = None
        self.tfile                      = None
        self.model                      = None
        self.model_type                 = None
        self.data_matrix                = None
        self.likelihood                 = None
        self.ntax                       = 0
        self.nchar                      = 0
        self.npatterns                  = 0
        self.params                     = None
        self.trees                      = None
        self.is_invariable_sites_model  = False
        self.is_discrete_gamma_model    = False

        self.num_edgelen_hyperparams    = 0

        self.stored_trees               = None  # list of trees built from tree definitions in the trees file
        self.taxon_labels               = None  # list of taxon labels from the data matrix
        self.param_file_lines           = None  # list of lines from the params file (header lines excluded)
        self.param_headers              = None  # list of parameter headers from the param file (self.param_names is a subset of this list)
        self.tree_objects               = None  # dictionary in which keys are unique tree identifiers, and values are tree objects
        self.parameters                 = None  # dictionary in which keys are unique tree identifiers, and values are lists
                                                #   e.g. self.parameters[tree_id][j] = {'rAG': -4.2432, 'freqC': -2.3243, ...}
                                                #   where tree_id is a list of split representations uniquely identifying a particular tree
                                                #   and j is the index of the (j+1)th sample pertaining to that tree. The keys in the
                                                #   dictionary (e.g. 'rAG') are parameter names and the values are transformed parameter values
        self.edge_lengths               = None  # dictionary in which keys are unique tree identifiers, and values are lists
                                                #   e.g. self.edge_lengths[tree_id][j] = {'-**-**--': -1.2345, '-*****--':-2.14225, ...}
                                                #   where tree_id is a list of split representations uniquely identifying a particular tree
                                                #   and j is the index of the (j+1)th sample from that tree. The keys in the dictionary (e.g. '-**-**--')
                                                #   are string representations of splits and the values are log-transformed edge lengths
        self.sample_size                = None  # self.sample_size[tid] is the sample size (after burnin is removed) for tree tid
        self.log_like                   = None  # list of log likelihood values gleaned from the params file
        self.log_prior                  = None  # list of log prior values gleaned from the params file
        self.log_posterior              = None  # list of log posterior values gleaned from the params file


        # initialize quantities used in Gelfand-Ghosh calculations
        self.gg_simdata             = likelihood.SimData()     # temporary container used to hold nascent posterior predictive simulation results until they have been analyzed
        self.gg_y                   = likelihood.SimData()     # observed dataset
        self.gg_mu                  = likelihood.SimData()     # mean of all posterior predictive datasets

        self.minbins                = None
        self.gg_binned_simdata      = None   # if self.gg_bin_patterns is True, this vector of floats is used to summarize the counts in self.gg_simdata
        self.gg_binned_y            = None   # if self.gg_bin_patterns is True, this vector of floats is used instead of self.gg_y
        if self.gg_num_bins == 7:
            self.minbins                = True
            self.gg_binned_mu           = [0.0]*7
        elif self.gg_num_bins == 15:
            self.minbins                = False
            self.gg_binned_mu           = [0.0]*15
        else:
            assert False, 'gg_num_bins must be set to either 7 or 15 (%d was specified in this case)' % self.gg_num_bins

        self.gg_a                   = []            # vector of compromise actions (one for each k in kvalues)
        self.gg_npatterns           = []            # vector containing the number of patterns in each posterior predictive dataset
        self.gg_t                   = []            # vector of t values computed from posterior predictive datasets
        self.gg_t_y                 = 0.0           # t for original dataset
        self.gg_t_mean              = 0.0           # mean of t over all posterior predictive datasets
        self.gg_t_mu                = 0.0           # t of mean over all posterior predictive datasets
        self.gg_t_a                 = []            # vector of t values computed from compromise action (one for each k in kvalues)
        self.gg_Pm                  = 0.0           # penalty component (same for all k)
        self.gg_Gm                  = []            # vector of goodness-of-fit components (one for each k in kvalues)
        self.gg_Dm                  = []            # vector of overall measures (one for each k in kvalues)
        self.gg_num_post_pred_reps  = 0.0           # counts total number of posterior predictive simulations performed
        self.gg_total               = 0

        # The data members below are needed only because an MCMCManager object is created
        #self.heat_vector = [1.0]        # specifies that just one chain will be created with power 1.0

        #self.gg_bincount_filename   = opts.gg_bincount_filename
        #self.gg_postpred_prefix     = opts.gg_postpred_prefix
        #self.datafname              = phycas.data_file_name
        #self.paramfname             = phycas.gg_pfile
        #self.treefname              = phycas.gg_tfile
        #self.rnseed                 = opts.random_seed
        #self.num_rates              = phycas.num_rates

        #self.gg_nreps               = opts.nreps
        #self.gg_kvect               = opts.kvalues
        #self.gg_burnin              = opts.burnin
        #self.heat_vector            = [1.0]

        #self.ok = True
        #if self.ok:
        #    if os.path.exists(self.paramfname):
        #        self.pfile = file(self.paramfname, 'r')
        #    else:
        #        self.ok = False
        #        print 'Parameter file specified (%s) does not exist' % self.paramfname
        #        sys.exit()
        #if self.ok:
        #    if os.path.exists(self.paramfname):
        #        self.tfile = file(self.treefname, 'r')
        #    else:
        #        self.ok = False
        #        print 'Tree file specified (%s) does not exist' % self.treefname
        #        sys.exit()

    def _openBinCountsFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the file used to save bin counts for each posterior-predivtive
        data set simulated.

        """
        self.phycassert(self.bincountsf is None, 'Attempt made to open bincounts file, but it is already open!')
        try:
            self.bincountsf = self.opts.out.bincounts.open(self.stdout)
        except:
            print '*** Attempt to open bin counts file (%s) failed.' % self.opts.out.bincounts._getFilename()

        if self.bincountsf:
            print 'Bin counts file was opened successfully'

    def _closeBinCountsFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the bin counts file.

        """
        self.phycassert(self.bincountsf is not None, 'Attempt made to close bin counts file, but it is not open!')
        self.opts.out.bincounts.close()
        self.bincountsf = None

    def _openPostPredFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the file used to save simulated posterior-predivtive data sets.

        """
        self.phycassert(self.postpredf is None, 'Attempt made to open post pred. file, but it is already open!')
        try:
            self.postpredf = self.opts.out.postpred.open(self.stdout)
        except:
            print '*** Attempt to open post pred. file (%s) failed.' % self.opts.out.postpred._getFilename()

        #if self.postpredf:
        #    print 'Post pred. file was opened successfully'

    def _closePostPredFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the post pred. file.

        """
        self.phycassert(self.postpredf is not None, 'Attempt made to close post pred. file, but it is not open!')
        self.opts.out.postpred.close()
        self.postpredf = None

    def storeTrees(self, input_trees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Open specified tree file and read trees therein, storing the list of
        tree objects in self.stored_trees.

        """
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_trees = list(input_trees)   # this triggers reading the tree file
        self.taxon_labels = input_trees.taxon_labels # this line must follow the coercion of the trees to a list
        num_stored_trees = len(self.stored_trees)
        self.stdout.phycassert(num_stored_trees > 0, 'Specified tree source (%s) contained no stored trees' %  str(input_trees))

    def storeParams(self, input_params):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Open specified params file and read parameters therein, storing the
        list of lines that represent sampled parameter vectors in
        self.param_file_lines and the header (containing the names of the
        parameters) in the list self.param_headers.

        """
        burnin = self.opts.burnin
        self.param_file_lines = open(input_params, 'r').readlines()
        self.stdout.phycassert(len(self.param_file_lines) >= 3 + burnin, "File '%s' does not look like a parameter file (too few lines)")
        self.param_headers = self.param_file_lines[1].split()
        self.stdout.phycassert(self.param_headers[1] != 'beta', "File '%s' appears to be the result of a stepping-stone analysis. This method requires a sample from the posterior (not power posterior) distribution." % input_params)

    def fillParamDict(self, param_vect):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build and return a dictionary containing the parameter names (keys)
        and transformed parameter values (values) for the sample supplied in
        param_vect.

        """
        param_dict = {}

        param_list = param_vect.split()
        assert len(param_list) == len(self.param_headers), 'Number of parameters in sample (%d) of parameter file differs from the number of parameters listed in the header line (%d)' % (len(param_list).len(self.param_headers))

        log_freqA = None
        log_rAC = None
        log_like = None
        freqs_stored = {}       # e.g. freqs_stored[1] = [0.2,0.3,0.4] <-- freqs for first subset
        relrates_stored = {}
        subset_rates_stored = []
        for h,p in zip(self.param_headers,param_list):
            if 'brlen' in h:
                # branch length parameters are gleaned from trees, so ignore the (redundant) edge lengths in the params file
                pass
            elif h in 'lnL':
                self.phycassert(log_like is None, 'expecting log_like to be None in fillParamDict')
                param_dict[h] = p
                log_like = float(p)
                self.log_like.append(log_like)
            elif h in 'lnPrior':
                self.phycassert(log_like is not None, 'expecting log_like to be not None in fillParamDict')
                param_dict[h] = p
                self.log_posterior.append(float(log_like) + float(p))
                self.log_prior.append(float(p))
                log_like = None
            elif h in ['Gen', 'TL']:
                # quantities that should be saved but not transformed
                param_dict[h] = p
            elif ('freqC' in h) or ('freqG' in h) or ('freqT' in h):
                # nucleotide frequencies for C, G and T are stored (and later transformed)
                self.phycassert(float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p)))
                param_dict[h] = float(p)
                parts = h.split('_')
                self.phycassert(len(parts) == 2, 'expecting name of state frequency parameter (%s) to have an underscore separating subset index from parameter name' % h)
                subset_index = int(parts[0])
                if subset_index in freqs_stored.keys():
                    freqs_stored[subset_index].append(h)
                else:
                    freqs_stored[subset_index] = [h]
            elif ('rAG' in h) or ('rAT' in h) or ('rCG' in h) or ('rCT' in h) or ('rGT' in h):
                # GTR exchangeabilities other than rAC are stored (and later transformed)
                self.phycassert(float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p)))
                param_dict[h] = float(p)
                parts = h.split('_')
                self.phycassert(len(parts) == 2, 'expecting name of GTR exchangeability parameter (%s) to to have an underscore separating subset index from parameter name' % h)
                subset_index = int(parts[0])
                if subset_index in relrates_stored.keys():
                    relrates_stored[subset_index].append(h)
                else:
                    relrates_stored[subset_index] = [h]
            elif 'subset_rate' in h:
                # subset relative rates othr than 1_subset_rate are stored (and later transformed)
                self.phycassert(float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p)))
                param_dict[h] = float(p)
                if h[:2] != '1_':
                    subset_rates_stored.append(h)
            else:
                # everything else is log-transformed and stored
                self.phycassert(float(p) > 0.0, '%s value (%g) zero or negative' % (h, float(p)))
                param_dict[h] = math.log(float(p))

        # If subset relative rates have been stored, compute first one and then transform all
        if len(subset_rates_stored) > 0:
            #print 'subset_rates_stored =',subset_rates_stored
            proportions = partition.getSubsetProportions()
            num_subsets = len(proportions)
            self.phycassert(len(proportions) == len(subset_rates_stored) + 1, 'Expecting to find %d subset rates in parameter file, but instead found %d' % (len(proportions)-1, len(subset_rates_stored)))
            modified_rates = [0.0]*(num_subsets-1)
            for rrname in subset_rates_stored:
                parts = rrname.split('_')
                self.phycassert(len(parts) == 3, 'Expecting subset relative rate parameter name (%s) to contain an underscore character separating subset index from the "subset_rate" part of the name' % rrname)
                subset_index = int(parts[0])
                modified_rates[subset_index - 2] = proportions[subset_index - 1]*param_dict[rrname]
            #print 'modified_rates =',modified_rates
            sum_rates = sum(modified_rates)
            log_first = math.log(1.0 - sum_rates)
            for z,modrate in enumerate(modified_rates):
                logrr = math.log(modrate)
                rrname = '%d_subset_rate' % (z+2,)
                param_dict[rrname] = logrr - log_first
                #print rrname,'-->',z,'=',logrr - log_first
            #pause('debugging subset relative rates')

        # If state frequencies have been stored, compute first one and then transform all
        if len(freqs_stored.keys()) > 0:
            #print 'freqs_stored.keys() =',freqs_stored.keys()
            for k in freqs_stored.keys():
                self.stdout.phycassert(len(freqs_stored[k]) == 3 or len(freqs_stored[k]) == 60, 'expecting either 3 or 60 state frequencies in parameter file for subset %d, but instead found %d' % (k,len(freqs_stored[k])))
                sum_freqs = sum([param_dict[frqname] for frqname in freqs_stored[k]])
                log_first = math.log(1.0 - sum_freqs)
                for frqname in freqs_stored[k]:
                    logfrq = math.log(param_dict[frqname])
                    param_dict[frqname] = logfrq - log_first

        # If GTR exchangeabilities have been stored, compute first one and then transform all
        if len(relrates_stored.keys()) > 0:
            #print 'relrates_stored.keys() =',relrates_stored.keys()
            for k in relrates_stored.keys():
                self.stdout.phycassert(len(relrates_stored[k]) == 5, 'expecting 5 GTR exchangeability parameters in parameter file for subset %d, but instead found %d' % (k,len(relrates_stored[k])))
                sum_relrates = sum([param_dict[rrname] for rrname in relrates_stored[k]])
                log_first = math.log(1.0 - sum_relrates)
                for rrname in relrates_stored[k]:
                    logrr = math.log(param_dict[rrname])
                    param_dict[rrname] = logrr - log_first

        # temporary!
        #for h,v in zip(self.param_headers,param_list):
        #    print v,'-->',h
        #print 'param_dict:',param_dict

        return param_dict

    def fillTreeNodeDict(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build a dictionary (self.curr_tree_node) in which keys are string
        representations of splits, and values are the nodes of self.curr_tree.
        Assumes that self.curr_tree is set to a tree.

        """
        self.curr_tree_node = {}

        # Recalc splits for tree so that we can replace edge lengths using splits as keys
        ntips = self.curr_tree.getNObservables()
        self.curr_tree.recalcAllSplits(ntips)

        # Traverse the tree
        nd = self.curr_tree.getFirstPreorder()
        assert nd.isRoot(), 'The first preorder node should be the root'
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                # Grab the split and invert it if necessary to attain a standard polarity
                s = nd.getSplit()
                if s.isBitSet(0):
                    s.invertSplit()

                if (s.countOnBits() == 0):
                    # root node of rooted tree, ignore because this node has no branch
                    continue

                # Create a string representation of the split
                ss = s.createPatternRepresentation()

                self.curr_tree_node[ss] = nd

    def fillEdgeLenDict(self, tree):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build and return a tuple containing: (1) a dictionary containing the
        string representations (keys) and log-transformed edge lengths
        (values) for the tree definition supplied in tree_def; and (2) the
        tree id, which is a tuple of split representations for this tree that
        can be used to uniquely identify the tree topology.

        """
        edgelen_dict = {}

        #tree = phylogeny.Tree()
        #tree_def.buildTree(tree)
        tree.rectifyNames(self.taxon_labels)
        ntips = tree.getNObservables()
        tree.recalcAllSplits(ntips)

        # Traverse the tree
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'The first preorder node should be the root'
        treelen = 0.0
        assert tree.hasEdgeLens(), 'Sampled tree has no edge lengths'
        split_list = []    # this will be a list of split representations (strings) that uniquely identifies the tree
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                # Determine whether this split represents an internal or tip node
                is_tip_node = nd.isTip() or nd.getParent().isRoot()

                # Grab the edge length
                edge_len = nd.getEdgeLen()
                treelen += edge_len

                # Grab the split and invert it if necessary to attain a standard polarity
                s = nd.getSplit()
                if s.isBitSet(0):
                    s.invertSplit()

                if (s.countOnBits() == 0):
                    # root node of rooted tree, ignore because this node has no branch
                    continue

                # Create a string representation of the split
                ss = s.createPatternRepresentation()
                split_list.append(ss)

                # Create a dictionary entry with ss as key and the log-transformed edge_len as value
                edgelen_dict[ss] = math.log(edge_len)

        split_list.sort()
        return edgelen_dict, tuple(split_list)

    def harvestTreesAndParams(self, input_trees, input_params):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Build dictionaries self.parameters and self.edge_lengths from
        supplied input_trees and input_params.

        """
        # store trees and parameter vectors from specified tree and parameter files, respectively
        self.storeTrees(input_trees)
        self.storeParams(input_params)

        # self.tree_objects will hold the constructed tree object for every tree definition read from the tree file
        self.tree_objects = {}

        # self.edge_lengths will hold all information about edge lengths
        self.edge_lengths = {}

        # Initialize the data member that will hold all information about parameter samples
        self.parameters = {}

        # eliminate burnin samples as well as 2 header lines
        post_burnin_paramvects = self.param_file_lines[self.opts.burnin+2:]

        # eliminate burnin samples
        post_burnin_trees = None
        post_burnin_trees = self.stored_trees[self.opts.burnin:]
        # post_burnin_paramvects and post_burnin_trees should now be the same length
        self.phycassert(len(post_burnin_paramvects) == len(post_burnin_trees), 'lists of trees and parameters have different lengths')

        self.nsamples = len(post_burnin_trees)

        # Loop through lists of trees and parameter vectors. Afterwards, these four
        # dictionaries will each have an entry for each distinct tree topology encountered:
        # - self.sample_size: key = treeid, value = sample size
        # - self.tree_objects: key = treeid, value = tree object
        # - self.edge_lengths: key = treeid, value = list of edgelen_dict objects
        # - self.parameters: key = treeid, value = list of param_dict objects
        self.sample_size = {}
        self.log_like = []
        self.log_prior = []
        self.log_posterior = []
        for tree, param_vect in zip(post_burnin_trees,post_burnin_paramvects):
            # This dictionary will be filled with transformed parameter values: e.g. param_dict['freqC'] = log(.251/.249); here, freqA=0.249
            param_dict = self.fillParamDict(param_vect)

            # edgelen_dict will be filled with log-transformed edge lengths: e.g. edgelen_dict[s] = log(.0023)
            # where the split object s associated with the edge is used as the key
            edgelen_dict, tree_id = self.fillEdgeLenDict(tree)

            # initialize or increase the sample size for tree
            if not tree_id in self.sample_size.keys():
                self.sample_size[tree_id] = 1
            else:
                self.sample_size[tree_id] += 1

            # add tree to the self.tree_objects dictionary if not already present
            if not tree_id in self.tree_objects.keys():
                self.tree_objects[tree_id] = tree

            # add edgelen_dict to the appropriate list in the self.edge_lengths dictionary
            if tree_id in self.edge_lengths.keys():
                # this tree has already been seen, so add to list already established for this tree
                self.edge_lengths[tree_id].append(edgelen_dict)
            else:
                # this tree has not yet been seen, so start a new list
                self.edge_lengths[tree_id] = [edgelen_dict]

            # add param_dict to the appropriate list in the self.parameters dictionary
            if tree_id in self.parameters.keys():
                # this tree has already been seen, so add to list already established for this tree
                self.parameters[tree_id].append(param_dict)
            else:
                # this tree has not yet been seen, so start a new list
                self.parameters[tree_id] = [param_dict]

    def loadData(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Loads data from data_source and sets data members taxon_labels, ntax,
        nchar, and data_matrix.

        """
        ds = self.opts.data_source
        self.phycassert(ds is not None, 'Data source is not allowed to be None for an IDR analysis')
        self.data_matrix = ds.getMatrix()
        self.phycassert(self.data_matrix is not None, 'Data matrix could not be input')
        self.taxon_labels = self.data_matrix.getTaxLabels()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar()
        self.phycassert(len(self.taxon_labels) == self.ntax, "Number of taxon labels does not match number of taxa.")
        self.phycassert(self.ntax > 0, 'Number of taxa in data matrix was 0')

    def getStartingTree(self):
        self.starting_tree = self.tree_objects[self.curr_treeid]
        return self.starting_tree

    def outputHeader(self):
        self.output('Gelfand-Ghosh Analysis:')
        self.output('  data source:        %s' % self.opts.data_source.filename)
        self.output('  parameter file:     %s' % self.opts.params)
        self.output('  tree file:          %s' % self.opts.trees.filename)
        self.output('  sims/sample:        %d' % self.opts.nreps)
        self.output('  k values:')
        for kvalue in self.opts.kvalues:
            self.output('    %f' % kvalue)
        self.output('  samples skipped:    %d' % self.opts.burnin)
        self.output()

    def outputDataInfo(self):
        self.output('  Information about the data:')
        self.output('    number of taxa:         %d' % self.ntax)
        self.output('    number of characters:   %d' % self.nchar)
        self.output('    number of patterns:     %d' % self.npatterns)
        self.output()

    def outputTreesInfo(self):
        self.output('  Information about the sampled trees:')
        self.output('    number of trees skipped:  %d' % self.gg_burnin)
        self.output('    number of trees included: %d' % self.nsamples)
        self.output('    total trees in file:      %d' % (self.gg_burnin + self.nsamples))
        self.output()

    def outputModelInfo(self):
        pinvar_str = self.is_invariable_sites_model and '+I' or ''
        dgamma_str = self.is_discrete_gamma_model and '+G' or ''

        model_str = '%s%s%s' % (self.model_type, pinvar_str, dgamma_str)
        self.output('  Information about the substitution model:')
        self.output('    model name:  %s' % model_str)
        self.output()

    def getData(self):
        self.data_matrix = phycas.readData(self.datafname)[0]
        self.taxon_labels = data_reader.getTaxLabels()
        self.ntax = self.data_matrix.getNTax()
        self.nchar = self.data_matrix.getNChar()
        assert self.likelihood, 'call setupModel before calling getData'
        self.likelihood.copyDataFromDiscreteMatrix(self.data_matrix)
        self.npatterns = self.likelihood.getNPatterns()
        self.outputDataInfo()

    def getTrees(self):
        tree_reader = readnexus.NexusReader()
        tree_reader.readFile(self.treefname)
        self.trees = tree_reader.getTrees()
        self.nsamples = len(self.trees) - self.gg_burnin
        self.outputTreesInfo()

    def debugShowBinCounts(self, msg, binned_dataset, tvalue):
        if self.minbins:
            print '\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n'
            print '%s (t = %g)' % (msg,tvalue)
            print '  A\t%12.5f' % binned_dataset[0]
            print '  C\t%12.5f' % binned_dataset[1]
            print '  G\t%12.5f' % binned_dataset[2]
            print '  T\t%12.5f' % binned_dataset[3]
            print '  2\t%12.5f' % binned_dataset[4]
            print '  3\t%12.5f' % binned_dataset[5]
            print '  4\t%12.5f' % binned_dataset[6]
        else:
            print '%s (t = %g)' % (msg,tvalue)
            print '  A   \t%12.5f' % binned_dataset[0]
            print '  C   \t%12.5f' % binned_dataset[1]
            print '  G   \t%12.5f' % binned_dataset[2]
            print '  T   \t%12.5f' % binned_dataset[3]
            print '  AC  \t%12.5f' % binned_dataset[4]
            print '  AG  \t%12.5f' % binned_dataset[5]
            print '  AT  \t%12.5f' % binned_dataset[6]
            print '  CG  \t%12.5f' % binned_dataset[7]
            print '  CT  \t%12.5f' % binned_dataset[8]
            print '  GT  \t%12.5f' % binned_dataset[9]
            print '  ACG \t%12.5f' % binned_dataset[10]
            print '  ACT \t%12.5f' % binned_dataset[11]
            print '  AGT \t%12.5f' % binned_dataset[12]
            print '  CGT \t%12.5f' % binned_dataset[13]
            print '  ACGT\t%12.5f' % binned_dataset[14]

    def ggCalculate(self):
        two_n = 2.0*float(self.nchar)

        # Compute the t function for the observed dataset
        if self.gg_bin_patterns:
            self.gg_t_y = self.gg_y.calctBinned(4, self.minbins)  #self.calcBinnedT(self.gg_binned_y)
            self.gg_binned_y = self.gg_y.getBinnedCounts()

            self.debugShowBinCounts('Observed dataset', self.gg_binned_y, self.gg_t_y)
        else:
            self.gg_t_y = self.gg_y.calct(4)

        # gg_mu (or gg_binned_mu) represents the mean over all posterior predictive datasets
        # Compute the t function for this mean dataset
        if self.gg_bin_patterns:
            self.gg_t_mu = self.gg_mu.calctBinned(4, self.minbins)  # self.calcBinnedT(self.gg_binned_mu)
            self.gg_binned_mu = self.gg_mu.getBinnedCounts()

            self.debugShowBinCounts('Mean dataset', self.gg_binned_mu, self.gg_t_mu)
        else:
            self.gg_t_mu = self.gg_mu.calct(4)

        # Compute the mean of the t values computed for individual posterior
        # predictive datasets
        self.gg_t_mean = sum(self.gg_t)/float(self.gg_total)

        # Compute the penalty term. Guaranteed to be positive by Jensen's
        # inequality and the convexity of the t function.
        self.gg_Pm = two_n*(self.gg_t_mean - self.gg_t_mu)

        # Loop over k values, computing Gm and Dm for each k value in kvalues
        for k in self.opts.kvalues:
            # Create a dataset representing the compromise "action"

            a = likelihood.SimData()
            a.resizePatternVect(self.nchar);
            self.gg_mu.addDataTo(a, 1.0)
            self.gg_y.addDataTo(a, k)
            a.divideBy(k + 1.0)
            if self.gg_bin_patterns:
                #a = []
                #for b in range(self.gg_num_bins):
                #    next_a = (self.gg_binned_mu[b] + self.gg_binned_y[b])/(k + 1.0)
                #    a.append(next_a)
                #t_a = self.calcBinnedT(a)
                t_a = a.calctBinned(4, self.minbins)
                gg_binned_a = a.getBinnedCounts()
                self.debugShowBinCounts('Compromise dataset (k = %g)' % k, gg_binned_a, t_a)
                self.gg_a.append(gg_binned_a)
            else:
                #a = likelihood.SimData()
                #self.gg_mu.addDataTo(a, 1.0)
                #self.gg_y.addDataTo(a, k)
                #a.divideBy(k + 1.0)
                t_a = a.calct(4)
                self.gg_a.append(a)
            self.gg_t_a.append(t_a)

            # Compute the goodness-of-fit term
            Gkm = (float(k) + 1.0)*two_n*((self.gg_t_mu + k*self.gg_t_y)/(k + 1.0) - t_a)
            self.gg_Gm.append(Gkm)

            # Compute the overall measure
            Dkm = self.gg_Pm + Gkm
            self.gg_Dm.append(Dkm)

        self.output('Pm = %f' % self.gg_Pm)
        for i,k in enumerate(self.opts.kvalues):
            self.output('k = %f:' % k)
            self.output('  Gm = %f' % self.gg_Gm[i])
            self.output('  Dm = %f' % self.gg_Dm[i])
        self.output()

        self.output('no. patterns in original dataset   = %d' % self.gg_y.getNUniquePatterns())
        self.output('no. patterns in mean dataset       = %d' % self.gg_mu.getNUniquePatterns())
        sum_npat = 0.0
        for npat in self.gg_npatterns:
            sum_npat += float(npat)
        self.output('mean no. patterns across datasets  = %f' % (sum_npat/float(len(self.gg_npatterns))))

        self.output('t for original dataset             = %f' % self.gg_t_y)
        self.output('t for mean dataset                 = %f' % self.gg_t_mu)
        self.output('mean of t across datasets          = %f' % self.gg_t_mean)
        for i,k in enumerate(self.opts.kvalues):
            self.output('t of compromise action for k = %.1f = %f' % (k,self.gg_t_a[i]))

        ttotal = len(self.gg_t)
        assert ttotal == self.gg_total, 'mismatch between self.gg_total and len(self.gg_t)'
        tsumsq = 0.0
        for t in self.gg_t:
            tsumsq += t*t
        tvar = tsumsq - float(ttotal)*self.gg_t_mean*self.gg_t_mean
        self.output('std. dev. of t across datasets     = %f' % math.sqrt(tvar))
        for i,k in enumerate(self.opts.kvalues):
            self.output('t of compromise action (k = %6f) = %f' % (k, self.gg_t_a[i]))

    def addPaupBlock(self, fn, tree, pheaders, pvalues):
        # this function not yet tested
        # add a paup block to allow easy estimation of parameters under maximum likelihood
        simf = file(fn, 'a')
        simf.write('\nbegin trees;\n')
        simf.write('  translate\n')
        for num,name in enumerate(self.taxon_labels):
            simf.write("  %d '%s'%s\n" % (num, name, num == self.ntax - 1 and ';' or ','))
        simf.write('  ;\n')
        simf.write('  utree one = %s;\n' % tree.makeNewick())
        simf.write('end;\n')
        simf.write('\nbegin paup;\n')
        simf.write('  log file=%s.log start replace;\n' % fn)
        simf.write('  set criterion=likelihood;\n')
        simf.write('  lset nst=6 basefreq=estimate tratio=estimate rmatrix=estimate pinvar=0.0 rates=gamma shape=estimate;\n')
        simf.write('  lscores 1;\n')
        simf.write('  describe 1 / brlens=sumonly;\n')
        simf.write('  log stop;\n')
        simf.write('  quit;\n')
        simf.write('end;\n')
        simf.write('[\n')
        simf.write('\n')
        simf.write('TL                     = %f\n' % tree.edgeLenSum())
        #if self.using_hyperprior:
        #    for p in self.chain_manager.getEdgeLenHyperparams():
        #        simf.write('edge length hyperparam = %f\n' % p.getCurrValue())
        simf.write('param headers          = %s\n' % pheaders)
        simf.write('param values           = %s\n' % pvalues)
        simf.write(']\n')
        simf.close()

    def loMismo(self, first_list, second_list, exceptions):
        """
        Compares first_list to second_list and returns True if identical with the exception of the
        values in the supplied exceptions list.
        """
        # make copies of supplied lists to allow tampering with contents
        first = list(first_list[:])
        second = list(second_list[:])

        # ensure that none of the exceptions are present in either list
        for e in exceptions:
            if e in first:
                first.remove(e)
            if e in second:
                second.remove(e)

        # sort both lists
        first.sort()
        second.sort()

        # now both lists should be identical
        lo_mismo = True
        for f,s in zip(first,second):
            if f != s:
                lo_mismo = False
                self.last_error = 'at least one parameter name mismatch: "%s" vs. "%s"' % (f,s)
                break

        return lo_mismo

    def setupModel(self, tid, core, param_map, edgelen_map):
        # create vector of transformed parameter values in correct order
        param_vect = []
        #print 'parameter names and values:'
        for k in core.partition_model.getFreeParameterNames():
            v_k = param_map[k]
            #print '  ',k,'-->',v_k,math.exp(v_k)
            param_vect.append(v_k)

        core.partition_model.setTransformedParameters(param_vect, self.curr_tree)

        #temporary
        #print 'param_vect <--> debug_param_vect:'
        #debug_param_vect = core.partition_model.getTransformedParameters()
        #sum_squared_diffs = 0.0
        #for p,q in zip(param_vect, debug_param_vect):
        #    print '%g <--> %g' % (p,q)
        #    sum_squared_diffs += (p-q)*(p-q)

        # if this is not done, new shape parameter value will be ignored
        core.likelihood.replaceModel(core.partition_model)

        # now set edge lengths
        for s in edgelen_map.keys():
            #print '>>>>',s,'<<<<'
            nd = self.curr_tree_node[s]
            edge_len = math.exp(edgelen_map[s])
            nd.setEdgeLen(edge_len)
        #print self.curr_tree.makeNewick()

    def debugPAUPNexus(self, datafn, core):
        fn = 'paup-%s' % datafn
        param_names = core.partition_model.getAllParameterNames()
        param_values = core.partition_model.getUntransformedParameters()

        pmap = {}
        for n,v in zip(param_names, param_values):
            pmap[n] = v

        rAC = pmap['1_rAC']
        rAG = pmap['1_rAG']
        rAT = pmap['1_rAT']
        rCG = pmap['1_rCG']
        rCT = pmap['1_rCT']
        rGT = pmap['1_rGT']

        rAC /= rGT
        rAG /= rGT
        rAT /= rGT
        rCG /= rGT
        rCT /= rGT
        rGT /= rGT

        print 'rmatrix = (%g %g %g %g %g %g)' % (rAC, rAG, rAT, rCG, rCT, rGT)

        frqA = pmap['1_freqA']
        frqC = pmap['1_freqC']
        frqG = pmap['1_freqG']
        frqT = pmap['1_freqT']

        print 'basefrequencies = (%g %g %g %g)' % (frqA, frqC, frqG, frqT)

        shape = pmap['1_gamma_shape']
        print 'shape = %g' % shape

        pinvar = pmap['1_pinvar']
        print 'pinvar = %g' % pinvar

        tmpf = open(fn, 'w')
        tmpf.write('#nexus\n\n')
        tmpf.write('\n')
        tmpf.write('begin paup;\n')
        tmpf.write('  set crit=like storebrlens;\n')
        tmpf.write('  exe %s;\n' % datafn)
        tmpf.write('end;\n')
        tmpf.write('\n')
        tmpf.write('begin trees;\n')
        tmpf.write('  translate\n')
        for num,name in enumerate(self.curr_tree.taxon_labels):
            tmpf.write("    %d '%s'%s\n" % (num+1, name, num + 1 == self.ntax and ';' or ','))
        tmpf.write('  utree test = %s;\n' % self.curr_tree.makeNumberedNewick(12))
        tmpf.write('end;\n')
        tmpf.write('\n')
        tmpf.write('begin paup;\n')
        tmpf.write('  lset nst=6 rmatrix=(%g %g %g %g %g) basefreq=(%g %g %g) rates=gamma shape=%g pinvar=%g;\n' % (rAC, rAG, rAT, rCG, rCT, frqA, frqC, frqG, shape, pinvar))
        tmpf.write('  lscores 1 / userbrlen;\n')
        tmpf.write('end;\n')
        tmpf.close()

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Performs simulations and computes the Gelfand-Ghosh measures Pm, Gm,
        and Dm = Pm + Gm for the parameters and trees specified when this
        object was created.

        """
        self.outputHeader()

        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        self.stdout.phycassert(input_trees is not None, 'trees cannot be None when gg is called')
        self.stdout.phycassert(input_trees.__class__.__name__ == 'TreeCollection', 'trees must be a TreeCollection object')

        # Check to make sure user specified an input params file
        input_params = self.opts.params
        self.stdout.phycassert(input_params is not None and len(input_params) > 0, 'params cannot be None or empty when gg is called')

        # Store transformed edge lengths and parameter values in dictionaries self.edge_lengths and self.parameters, respectively
        print 'Harvesting trees and parameters...'
        self.harvestTreesAndParams(input_trees, input_params)

        # Build list of tuples (t,n) where t is the tree id and n is the frequency of that tree in the sample
        tree_id_samplesize_tuples = []
        for tid in self.parameters.keys():
            tree_id_samplesize_tuples.append((tid,self.sample_size[tid]))

        # Sort list so that most commonly-encountered tree topology is first
        tree_id_samplesize_tuples.sort(cmp=lambda x,y: cmp(x[1], y[1]))
        tree_id_samplesize_tuples.reverse()

        self.loadData()
        self.gg_mu.resizePatternVect(self.nchar);

        self.curr_treeid, self.n = tree_id_samplesize_tuples[0]

        if False:
            # Create chain manager and (single) chain; no MCMC being done but the chain knows how to
            # interact with the model and compute likelihoods and priors
            mcmc_manager = MCMCManager(self)
            mcmc_manager.createChains()
            cc = self.mcmc_manager.getColdChain()
            core = cc.likelihood
        core = LikelihoodCore(self)
        core.setupCore()

        # Check to make sure model hasn't changed since creating the param file
        param_names_from_model = core.partition_model.getAllParameterNames()
        param_names_from_file  = self.parameters[tid][0].keys()
        exclude_from_comparison = ['Gen','TL','lnL','lnPrior', '1_external_hyper', '1_internal_hyper', '1_edgelen_hyper']
        #print '\n*******************\nparam_names_from_model:',param_names_from_model
        #print '\n*******************\nparam_names_from_file:',param_names_from_file
        self.phycassert(self.loMismo(param_names_from_model, param_names_from_file, exclude_from_comparison), 'Model differs from the model used to generate parameter file (%s)' % self.last_error)

        # Let gg_y contain the observed pattern counts
        core.likelihood.addDataTo(self.gg_y)

        #POLOLD
        #if self.gg_bin_patterns:
        #    self.gg_binned_y = self.gg_y.getBinnedCounts()

        if self.gg_bin_patterns and self.opts.out.bincounts is not None:
            self._openBinCountsFile()

        #######################
        ###### TREE LOOP ######
        #######################
        prev_pct_done = 0.0
        prev_secs = 0.0
        self.gg_total = 0
        stopwatch = probdist.StopWatch()
        stopwatch.start()
        num_distinct_tree_topologies = len(tree_id_samplesize_tuples)
        self.output('Performing posterior-predictive simulations (%d distinct topologies):' % num_distinct_tree_topologies)
        for tnum,(tid,sample_size) in enumerate(tree_id_samplesize_tuples):
            self.n = sample_size
            #self.output('\nEvaluating tree = %d (%s)...' % (tnum+1,sample_size == 1 and '1 sample' or '%d samples' % sample_size))

            # Replace the tree in likelihood core
            core.setTree(self.tree_objects[tid])
            self.curr_tree = core.tree

            # Build up dictionary of nodes (self.curr_tree_node) in which keys are split representations and values
            # are node objects in self.curr_tree. This will be used later to set edge lengths in the tree.
            self.fillTreeNodeDict()

            # Prepare tree for likelihood calculation by allocating tip and internal data structures
            core.likelihood.prepareForLikelihood(self.curr_tree)
            #core.likelihood.prepareForSimulation(self.curr_tree)

            #########################
            ###### SAMPLE LOOP ######
            #########################
            for param_map,edgelen_map in zip(self.parameters[tid],self.edge_lengths[tid]):

                # Report progress so user doesn't give up
                pct_done = 100.0*float(self.gg_total)/float(self.opts.nreps*self.nsamples)
                if pct_done - prev_pct_done >= 10.0:
                    prev_pct_done = pct_done
                    secs = stopwatch.elapsedSeconds()
                    proportion_finished = pct_done/100.0
                    proportion_remaining = 1.0 - proportion_finished
                    eta = secs*proportion_remaining/proportion_finished
                    self.output('  %.0f%% done (%.1fs remaining)...' % (pct_done, eta))
                    prev_secs = secs
                    cum_caltbinned_secs = 0.0

                self.setupModel(tid, core, param_map, edgelen_map)

                #core.likelihood.storeSiteLikelihoods(True)
                #core.likelihood.storeAllCLAs(self.curr_tree)

                lnL = core.likelihood.calcLnL(self.curr_tree)

                #print '-->',lnL,param_map['lnL']
                #siteLnLs = core.likelihood.getSiteLikelihoods()
                #counts = core.likelihood.getPatternCounts()
                #patterns = core.likelihood.listPatterns()
                #doof = open('doof_site_likes.txt','w')
                #doof.write('Patterns:\n%s\n\n' % patterns)
                #for s,c in zip(siteLnLs,counts):
                #    doof.write('%g\t%g\n' % (s,c))
                #doof.close()
                #core.likelihood.debugUncompressedDataInfo("all-site-patterns.txt");

                for j in range(self.opts.nreps):
                    self.gg_num_post_pred_reps += 1.0

                    self.gg_simdata = core.simulate()

                    # Save the simulated data set if desired
                    if self.opts.out.postpred.__bool__():
                        self._openPostPredFile()
                        self.gg_simdata.saveToNexusFilePython(self.postpredf, self.taxon_labels, 'dna', ('a','c','g','t'))
                        #self.gg_simdata.saveToNexusFile('doof.nex', self.taxon_labels, 'dna', ('a','c','g','t'))
                        #self.debugPAUPNexus(self.opts.out.postpred._getFilename(), core)
                        self._closePostPredFile()

                    if self.gg_bin_patterns:
                        # Compute the t function for the simulated dataset

                        curr_t = self.gg_simdata.calctBinned(4,self.minbins) # 4 = number of states

                        self.gg_binned_simdata = self.gg_simdata.getBinnedCounts()

                        if self.bincountsf is not None:
                            bstr = ['%.1f' % x for x in self.gg_binned_simdata]
                            self.bincountsf.write('%s\tposterior predictive replicate\n' % '\t'.join(bstr))
                    else:
                        # Compute the t function for the simulated dataset
                        curr_t = self.gg_simdata.calct(4)

                    # Add this value of t to the list (later the mean t will be computed)
                    self.gg_t.append(curr_t)

                    # Add the number of patterns in self.gg_simdata to the gg_npatterns list
                    self.gg_npatterns.append(self.gg_simdata.getNUniquePatterns())

                    # Update running mean vector gg_mu. A running mean is maintained because
                    # it is easy for the number of counts of constant patterns to overflow
                    # if you wait until the end of the MCMC run to divide by the total.
                    # Here is how the running mean is kept. Assume there will be four numbers
                    # (a, b, c, d) averaged. Thus, the desired quantity is (a+b+c+d)/4.
                    #
                    # gg_num_post_pred_reps   self.gg_mu
                    # ------------------------------------------------------------
                    #           1             a                      = (a)/1
                    #           2             (1/2)a + b/2           = (a+b)/2
                    #           3             (2/3)[(a+b)/2] + c/3   = (a+b+c)/3
                    #           4             (3/4)[(a+b+c)/3] + d/4 = (a+b+c+d)/4
                    # ------------------------------------------------------------
                    #
                    # Note that it is ok if gg_num_post_pred_reps = 1 (in which case
                    # gg_mu is multiplied by zero). Because gg_mu is empty, multBy is
                    # a no-op in this case
                    p = 1.0/self.gg_num_post_pred_reps
                    not_p = 1.0 - p

                    self.gg_mu.multBy(not_p)
                    self.gg_simdata.multBy(p)
                    self.gg_simdata.addDataTo(self.gg_mu, 1.0)

                    # Increment count of the total number of simulated datasets created
                    # This value is used to later compute the mean t for all simulated datasets
                    # and the mean counts for all simulated data sets
                    self.gg_total += 1

        self.ggCalculate()

        if self.gg_bin_patterns and self.bincountsf is not None:
            # write observed bin counts
            bstr = ['%.1f' % x for x in self.gg_binned_y]
            self.bincountsf.write('%s\tobserved\n' % '\t'.join(bstr))

            # write mu bin counts
            bstr = ['%.1f' % x for x in self.gg_binned_mu]
            self.bincountsf.write('%s\tmu\n' % '\t'.join(bstr))

            # write compromise action bin counts
            for i,k in enumerate(self.opts.kvalues):
                bstr = ['%.1f' % x for x in self.gg_a[i]]
                self.bincountsf.write('%s\ta for k=%.1f\n' % ('\t'.join(bstr),k))

            self._closeBinCountsFile()

        return (self.gg_Pm, self.gg_Gm, self.gg_Dm);



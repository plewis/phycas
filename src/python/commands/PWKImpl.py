import os,sys,math,random,re
from phycas import *
from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions
from phycas.commands.MCMCImpl import MCMCImpl

# class VarianceZeroError(Exception):
#     def __init__(self):
#         self.msg = 'Cannot calculate autocorrelation because variance is zero'
#     def __str__(self):
#         return self.msg
#
# class VarianceUndefinedError(Exception):
#     def __init__(self):
#         self.msg = 'Cannot calculate standard deviation because sample size is less than 2'
#     def __str__(self):
#         return self.msg
#
class InvalidNumberOfColumnsError(Exception):
    def __init__(self, nparts, nexpected, line_num):
        self.msg = 'Number of values (%d) on line %d inconsistent with number of column headers (%d)' % (nparts, line_num, nexpected)
    def __str__(self):
        return self.msg

class PartitionWeightedKernelEstimator(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Estimates the marginal likelihood using a posterior sample of
    parameters contained in params and a corresponding tree sample
    contained in trees, using the Partition-Weighted Kernel (PWK)
    method.

    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes ParamSummarizer object.

        """
        CommonFunctions.__init__(self, opts)

    def transformParamVector(self, pvalues, pnames):
        p = []
        logpiA = None
        logrAC = None
        for h,x in zip(pnames,pvalues):
            m = re.match('(\d+)_(.+)',h)
            assert m is not None
            the_subset = int(m.group(1))
            the_param_name = m.group(2)
            if the_param_name in ['edgelen_hyper', 'pinvar', 'gamma_shape', 'rAC', 'rAG', 'rAT', 'rCG', 'rCT', 'rGT', 'freqA', 'freqC', 'freqG', 'freqT']:
                try:
                    logx = math.log(float(x))
                    if the_param_name == 'gamma_shape':
                        p.append(logx)
                    elif the_param_name == 'pinvar':
                        p.append(logx - math.log(1.0 - float(x)))
                    elif the_param_name == 'edgelen_hyper':
                        p.append(logx)
                    elif the_param_name == 'freqA':
                        logpiA = logx
                    elif the_param_name in ['freqC', 'freqG', 'freqT']:
                        p.append(logx - logpiA)
                    elif the_param_name == 'rAC':
                        logrAC = logx
                    elif the_param_name in ['rAG', 'rAT', 'rCG', 'rCT', 'rGT']:
                        p.append(logx - logrAC)
                except ValueError:
                    print 'log(%s) failed for %s' % (x,h)
        return p

    def obsolete_transformParameters(self, param_headers, param_lines, tree_map, skip):
        # interpret parameters and log (or log-ratio) transform, creating param_map with keys identical to tree_map
        # param_map[k] thus holds parameters sampled for the tree topology defined in tree_map[k]

        # First, generate the param_type list from the param_headers
        #   1: edge length
        #   2: GTR exchangeability
        #   3: nucleotide frequency
        #   4. gamma shape
        #   5. pinvar
        #   6. subset relative rate
        param_type = []
        for h in param_headers:
            if h == 'alpha':
                param_type.append(4)
            elif h in ['pi(C)', 'pi(G)', 'pi(T)']:
                param_type.append(3)
            elif h in ['r(A<->G)', 'r(A<->T)', 'r(C<->G)', 'r(C<->T)', 'r(G<->T)']:
                param_type.append(2)

        # Second, store post-burnin transformed parameter vectors in a list log_params
        row_start = 2 + skip  # first line ID, second line headers
        lines_not_skipped = param_lines[row_start:]
        n = len(lines_not_skipped)
        log_params = []
        for i,line in enumerate(lines_not_skipped):
            parts = line.split()
            if len(parts) != len(param_headers):
                raise InvalidNumberOfColumnsError(len(parts), len(param_headers), i + row_start + 1)
            p = []
            logpiA = None
            logrAC = None
            for h,x in zip(param_headers,parts):
                if h in ['alpha', 'pi(A)', 'pi(C)', 'pi(G)', 'pi(T)', 'r(A<->C)', 'r(A<->G)', 'r(A<->T)', 'r(C<->G)', 'r(C<->T)', 'r(G<->T)']:
                    try:
                        logx = math.log(float(x))
                    except ValueError:
                        print 'log(%s) failed for %s' % (x,h)
                    if h == 'alpha':
                        p.append(logx)
                    elif h == 'pi(A)':
                        logpiA = logx
                    elif h in ['pi(C)', 'pi(G)', 'pi(T)']:
                        p.append(logx - logpiA)
                    elif h == 'r(A<->C)':
                        logrAC = logx
                    elif h in ['r(A<->G)', 'r(A<->T)', 'r(C<->G)', 'r(C<->T)', 'r(G<->T)']:
                        p.append(logx - logrAC)
            log_params.append(p)

        print 'log_params has length %d' % len(log_params)
        print 'log_params[0]  =', log_params[0]
        print 'log_params[-1] =', log_params[-1]

        # Now, create param_map from log_params
        param_map = {}
        for k in tree_map.keys():
            indices = tree_map[k][3]
            v = []
            for i in indices:
                param_vector_i = log_params[i-1][:]
                v.append(param_vector_i)
            param_map[k] = v

        log_params = None
        return param_type, param_map

    def segregateParamsByTopology(self, param_headers, param_lines, tree_map, skip):
        # Creates map in which keys are treeIDs and values are lists of raw parameter vectors

        # Store post-burnin parameter vectors in one big list raw_params
        row_start = 2 + skip  # first line ID, second line headers
        lines_not_skipped = param_lines[row_start:]
        n = len(lines_not_skipped)
        raw_params = []
        log_likelihoods = []
        for i,line in enumerate(lines_not_skipped):
            parts = line.split()
            if len(parts) != len(param_headers):
                raise InvalidNumberOfColumnsError(len(parts), len(param_headers), i + row_start + 1)
            p = []
            for h,x in zip(param_headers,parts):
                if h[0] in ['1','2','3','4','5','6','7','8','9']:
                    p.append(float(x))
                elif h == 'lnL':
                    log_likelihoods.append(float(x))
            raw_params.append(p)

        print 'raw_params has length %d' % len(raw_params)
        print 'log_likelihoods has length %d' % len(log_likelihoods)

        # Now, create param_map from raw_params
        param_map = {}
        lnL_map = {}
        for k in tree_map.keys():
            indices = tree_map[k][3]
            lnL = []
            v = []
            for i in indices:
                param_vector_i = raw_params[i-1][:]
                v.append(param_vector_i)
                lnL.append(log_likelihoods[i-1])
            param_map[k] = v
            lnL_map[k] = lnL

        raw_params = None
        return param_map, lnL_map

    def recordNodeInMaps(self, nd, splits_in_tree, split_weights):
        # nd is a node in the tree
        # splits_in_tree is a list, and this function will add the split corresponding to nd to this list
        # split_weights maps splits to edge lengths, and this function will add a key,value pair for nd

        # Grab the split and invert it if necessary to attain a standard polarity
        s = nd.getSplit()
        #if (not self.rooted_trees) and s.isBitSet(0):
        if s.isBitSet(0):
            s.invertSplit()

        # Create a string representation of the split
        ss = s.createPatternRepresentation()

        # Add edge length to split_weights map
        split_weights[ss] = nd.getEdgeLen()

        # Add string represention of the split to the splits_in_tree list, which
        # will be used to uniquely identify the tree topology
        is_pendant_edge = nd.isTip() or nd.getParent().isRoot()
        if not is_pendant_edge:
            splits_in_tree.append(ss)

    def recordTreeInMaps(self, tree, edge_map):
        # First build tree_key and create split_weights dict which maps splits
        # to edge lengths
        splits_in_tree = []
        split_weights = {}
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'the first preorder node should be the root'
        treelen = 0.0
        has_edge_lens = tree.hasEdgeLens()
        assert has_edge_lens, 'PWK requires trees to have edge lengths'

        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                treelen += nd.getEdgeLen()
                self.recordNodeInMaps(nd, splits_in_tree, split_weights)

        splits_in_tree.sort()
        tree_key = tuple(splits_in_tree)

        # Now that we have tree_key, add edge lengths to edge_lengths list
        if tree_key in edge_map.keys():
            edge_lengths = [[] for s in range(len(split_weights))]
            for (split,edgelen) in split_weights.items():
                i = edge_map[tree_key][split]
                edge_lengths[i] = edgelen
        else:
            edge_lengths = []
            edge_map[tree_key] = {}
            for i,(split,edgelen) in enumerate(split_weights.items()):
                edge_map[tree_key][split] = i
                edge_lengths.append(edgelen)

        return treelen, tree_key, edge_lengths

    def processTrees(self, trees, skip):
        # Read trees and create tree_map:
        #  key   = tree ID
        #  value = list containing the following elements:
        #     0: number of trees with this topology (n)
        #     1: newick tree description
        #     2: mean tree length
        #     3: list (length n) of indices of samples having this topology (first is 1, not 0)
        #     4: list (length n) of lists (length number of edges) of untransformed edge lengths
        #     5: edge_map (key = split, value = index into tree_map[4] lists), ensures that edge lengths
        #           are always in the same order in the parameter vector

        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        self.stdout.phycassert(input_trees, 'trees cannot be None or empty when pwk method is called')

        num_trees = 0
        self.num_trees_considered = 0

        # key is list of splits (a tree ID), value is tuple(count, newick, treelen, 1st time seen, 2nd time seen, ...)
        tree_map = {}

        # key is list of splits (a tree ID), value is a dictionary mapping each split in the key to an index in a list
        # of edge lengths
        edge_map = {}

        # Open sumt_tfile_name and read trees therein
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_tree_defs = list(input_trees)
        self.taxon_labels = input_trees.taxon_labels # this must be kept after the coercion of the trees to a list (in case that is what triggers the readinf of the file with the taxon labels)

        num_stored_trees = len(self.stored_tree_defs)
        self.stdout.phycassert(num_stored_trees > 0, 'Specified tree source (%s) contained no stored trees' %  str(input_trees))

        # Build each tree and add the splits and tree topolgies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        self.stdout.info('Compiling list of tree topologies...')
        t = phylogeny.Tree()

        # values used for display purposes
        split_field_width = 0
        sojourn_field_width = 2 + math.floor(math.log10(float(num_stored_trees)))

        curr_tree = 0
        polytomous_trees_found = False
        for tree_def in self.stored_tree_defs:
            curr_tree += 1
            if curr_tree % (num_stored_trees//10) == 0:
                pct_done = 100.0*curr_tree/num_stored_trees
                self.stdout.info('  %.1f%% done' % pct_done)

            if num_trees < self.opts.skip:
                num_trees += 1
                continue
            num_trees += 1
            self.num_trees_considered += 1

            # Build the tree
            tree_def.buildTree(t)
            t.rectifyNames(self.taxon_labels)
            ntips = t.getNObservables()
            if t.isPolytomous():
                polytomous_trees_found = True
            if ntips > split_field_width:
                # this is necessary only if number of taxa varies from tree to tree
                split_field_width = ntips
            t.recalcAllSplits(ntips)

            # Create a tree key, which is a list of internal node splits that can be
            # used to uniquely identify a tree topology
            treelen, tree_key, edge_lengths = self.recordTreeInMaps(t, edge_map)

            # Update tree_map, which is a map with keys equal to lists of internal node splits
            # and values equal to 2-element lists containing the frequency and newick tree
            # description
            if tree_key in tree_map.keys():
                # tree topology has been seen before
                entry = tree_map[tree_key]
                entry[0] += 1           # increment count of times this tree topology has been seen
                # entry[1] is the newick tree description for the first tree with this topology
                entry[2] += treelen     # add treelen to sum of tree lengths
                entry[3].append(self.num_trees_considered)
                entry[4].append(edge_lengths)
                # entry[5] is the edge_map storing position of each edge in the edge length vectors
            else:
                # tree topology has not yet been seen
                tree_map[tree_key] = [1, tree_def, treelen, [self.num_trees_considered], [edge_lengths], edge_map[tree_key]]

        self.stdout.info('\nSummary of sampled trees:')
        self.stdout.info('-------------------------')
        self.stdout.info('Tree source: %s' %  str(input_trees))
        self.stdout.info('Total number of trees in file = %d' % num_trees)
        self.stdout.info('Number of trees considered = %d' % self.num_trees_considered)
        self.stdout.info('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        # self.stdout.info('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNObservables()))

        # Sort tree keys by number of topologies, with largest sample sizes first
        sorted_keys = sorted([(tree_map[k][0],k) for k in tree_map.keys()])
        sorted_keys.reverse()

        self.stdout.info('Tree topology summary:')
        for i,(ignored,k) in enumerate(sorted_keys):
            self.stdout.info('\n  Topology number: %d' % (i+1,))
            self.stdout.info('  Number of trees with this topology: %d' % tree_map[k][0])
            self.stdout.info('  Newick of first tree with this topology: %s' % tree_map[k][1])
            self.stdout.info('  Mean tree length: %d' % (tree_map[k][2]/tree_map[k][0],))
            self.stdout.info('  Topology ID:')
            for s in k:
                self.stdout.info('    %s' % s)
            self.stdout.info('  Indices of trees with this topology (first tree has index 1, not 0):')
            n = tree_map[k][0]
            if n > 5:
                self.stdout.info('    %s,...,%d' % (','.join(['%d' % m for m in tree_map[k][3][:5]]),tree_map[k][3][-1]))
            else:
                self.stdout.info('    %s' % ','.join(['%d' % m for m in tree_map[k][3]]))

        if self.num_trees_considered == 0:
            self.stdout.info('\npwk was not provided with any trees.')

        return tree_map

    def processParams(self, params, tree_map, skip):
        # Read trees and create param_map:
        #
        # param_map is a dictionary in which:
        #   key   = treeID
        #   value = list of parameter vectors for one tree topology
        #
        param_lines = open(params, 'r').readlines()
        if len(param_lines) < 3 + skip:
            self.output("File '%s' does not look like a parameter file (too few lines)")
        else:
            param_headers = param_lines[1].split()
        return self.segregateParamsByTopology(param_headers, param_lines, tree_map, skip)

    def replaceEdgeLengths(self, tree, edge_lengths, edge_map):
        ntips = tree.getNObservables()
        tree.recalcAllSplits(ntips)
        nd = tree.getFirstPreorder()
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                s = nd.getSplit()
                if s.isBitSet(0):
                    s.invertSplit()
                ss = s.createPatternRepresentation()
                pos = edge_map[ss]
                nd.setEdgeLen(edge_lengths[pos])

    def processOneTopology(self, lnL_map, param_map, tree_map, treeID):
        print '*** new tree topology ***'

        # create a list v of vectors each containing all parameters
        # and calculate mean parameter vector
        untransformed_param_vectors = param_map[treeID]

        untransformed_edgelens = tree_map[treeID][4]

        # n is the number of samples for this tree topology
        n = len(untransformed_param_vectors)
        self.phycassert(len(untransformed_edgelens) == n, "Number of observations differs for params and trees in PWKImple.py standardizeParamsOneTopology function")

        # p is the number of parameters
        p = len(untransformed_param_vectors[0]) + len(untransformed_edgelens[0])

        # get cold chain and obtain parameter names from it
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        param_names = cold_chain.partition_model.getAllParameterNames()

        # set tree topology
        newick = tree_map[treeID][1]
        newick.buildTree(cold_chain.tree)
        cold_chain.likelihood.prepareForLikelihood(cold_chain.tree)

        # set edge lengths
        self.replaceEdgeLengths(cold_chain.tree, untransformed_edgelens[0], tree_map[treeID][5])

        # transform parameter vector and send to model
        param_vect = self.transformParamVector(untransformed_param_vectors[0], param_names)
        cold_chain.partition_model.setTransformedParameters(param_vect, cold_chain.tree)
        cold_chain.likelihood.recalcRelativeRates()

        # sanity check
        param_values = cold_chain.partition_model.getUntransformedParameters()
        print '\nSanity check:'
        print 'model:',cold_chain.partition_model.getModel(0).getModelName()
        print 'tree:',cold_chain.tree.makeNewick()
        for name,value in zip(param_names,param_values):
            print '>>> ',name,'-->',value

        # calculate log likelihood and compare with value in param file
        print 'calculating likelihood...'
        lnL = cold_chain.likelihood.calcLnL(cold_chain.tree)
        print 'lnL =',lnL,' (',lnL_map[treeID][0],')'
        print 'n =',n
        print 'index of this tree =',tree_map[treeID][3][0]
        #print 'cold_chain.tree.makeNewick() =',cold_chain.tree.makeNewick()
        print 'newick =',newick
        raw_input('..')

        if False:
            # Create a PWKMargLike object to manage the calculations
            V = likelihood.PWKMargLikeBase(n, p)
            V.copyNewick(tree_map[treeID][1])

            # Inform V of the parameter types
            nedges = len(log_edgelens[0])
            ptypes = [1]*nedges + param_type[:]
            V.setParamTypes(ptypes)

            # Add all parameter vectors to V and
            for q,e in zip(log_params, log_edgelens):
                eq = e+q
                V.copyParamVector(eq)

            V.estimateMargLike()

    def estimateMargLike(self, params, trees):
        marglike = None
        skip = self.opts.skip

        # Open the tree file and process lines
        tree_map = self.processTrees(trees, skip)

        # Create sorted list of tree IDs with tree having most samples first
        sorted_tree_IDs = [(tree_map[k][0],k) for k in tree_map.keys()]
        sorted_tree_IDs.sort(cmp = lambda x,y: cmp(y[0],x[0]))

        # Open the parameter file and process lines, separating parameter vectors according to tree topology
        # and creating param_map in which keys are tree IDs and values are lists of parameter vectors
        param_map, lnL_map = self.processParams(params, tree_map, skip)

        # Loop over tree topologies
        print '\nLooping over tree topologies:'
        for n,k in sorted_tree_IDs:
            print '  n = ',n,', k =',k
            if n >= self.opts.minsample:
                print 'Including topology with %d samples' % n
                self.processOneTopology(lnL_map, param_map, tree_map, k)
            else:
                print 'Excluding topology with %d samples' % n

        return marglike

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the param file and decides whether to use the
        marginal_likelihood or std_summary functions to summarize the data
        therein. The marginal_likelihood function is called if the second
        header is "beta"; otherwise, std_summary is called. If a cpofile
        name has been specified, computes conditional predictive ordinates.

        """
        global pwk_mcmc_impl

        self.phycassert(len(self.opts.params) > 0, "Must specify a parameter file ('params') when invoking the pwk command")
        self.stdout.phycassert(self.opts.trees, 'trees cannot be None or empty when pwk method is called')

        # Set mcmc.doing_pwk to True just long enough to grab the MCMCImpl object (also inhibits mcmc() from running an MCMC analysis)
        mcmc.doing_pwk = True
        pwk_mcmc_impl = mcmc()
        pwk_mcmc_impl.setup()
        mcmc.doing_pwk = False

        marglike = self.estimateMargLike(self.opts.params, self.opts.trees)
        return marglike

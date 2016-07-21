import os,sys,math,random
from phycas import *
from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions

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

    def _cpoOpenRFile(self):
        if self._cpoRFile is None:
            sp = self.optsout.cpoplot
            self._cpoRFilename = sp._getFilename()
            self._cpoRFile = sp.open(self.stdout)
        return self._cpoRFile

    def transformParameters(self, param_headers, param_lines, tree_map, skip):
        # interpret parameters and log (or log-ratio) transform, creating param_map with keys identical to tree_map
        # param_map[k] thus holds parameters sampled for the tree topology defined in tree_map[k]

        # First, store post-burnin transformed parameter vectors in a list log_params
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
            indices = tree_map[k][3:]
            v = []
            for i in indices:
                param_vector_i = log_params[i-1][:]
                v.append(param_vector_i)
            param_map[k] = v

        log_params = None # hoping this will trigger garbage collection
        return param_map

    def recordNodeInMaps(self, nd, tree_key, edge_len):
        # Grab the split and invert it if necessary to attain a standard polarity
        s = nd.getSplit()
        #if (not self.rooted_trees) and s.isBitSet(0):
        if s.isBitSet(0):
            s.invertSplit()

        # Create a string representation of the split
        ss = s.createPatternRepresentation()

        # Add string represention of the split to the tree_key list, which
        # will be used to uniquely identify the tree topology
        is_pendant_edge = nd.isTip() or nd.getParent().isRoot()
        if not is_pendant_edge:
            tree_key.append(ss)

    def recordTreeInMaps(self, tree, tree_key):
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'the first preorder node should be the root'
        treelen = 0.0
        has_edge_lens = tree.hasEdgeLens()

        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                edge_len = has_edge_lens and nd.getEdgeLen() or 1.0
                treelen += edge_len
                self.recordNodeInMaps(nd, tree_key, edge_len)
        return treelen

    def processTrees(self, trees, skip):
        # Read trees and create tree_map:
        #  key   = tree ID
        #  value = list
        #     0: number of trees with this topology (n)
        #     1: newick tree description
        #     2: mean tree length
        #     3: index (starting with 1) of first sample with this topology
        #     4: index of second sample with this topology
        #     .
        #     .
        #     .
        #   n+2: index of (n)th sample with this topology
        
        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        self.stdout.phycassert(input_trees, 'trees cannot be None or empty when pwk method is called')

        num_trees = 0
        self.num_trees_considered = 0

        # key is list of splits, value is tuple(count, newick, treelen, 1st time seen, 2nd time seen, ...)
        tree_map = {}

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

            # Create a tree key, which is a list of internal node splits that can be
            # used to uniquely identify a tree topology
            tree_key = []

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

            treelen = self.recordTreeInMaps(t, tree_key)

            # Update tree_map, which is a map with keys equal to lists of internal node splits
            # and values equal to 2-element lists containing the frequency and newick tree
            # description
            tree_key.sort()
            k = tuple(tree_key)
            if k in tree_map.keys():
                # tree topology has been seen before
                entry = tree_map[k]
                entry[0] += 1           # increment count of times this tree topology has been seen
                # entry[1] is the newick tree description for the first tree with this topology
                entry[2] += treelen     # add treelen to sum of tree lengths
                entry.append(self.num_trees_considered)
            else:
                # tree topology has not yet been seen
                tree_map[k] = [1, tree_def, treelen, self.num_trees_considered]

        self.stdout.info('\nSummary of sampled trees:')
        self.stdout.info('-------------------------')
        self.stdout.info('Tree source: %s' %  str(input_trees))
        self.stdout.info('Total number of trees in file = %d' % num_trees)
        self.stdout.info('Number of trees considered = %d' % self.num_trees_considered)
        self.stdout.info('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        # self.stdout.info('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNObservables()))

        self.stdout.info('Tree topology summary:')
        for i,k in enumerate(tree_map.keys()):
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
                self.stdout.info('    %s,...,%d' % (','.join(['%d' % m for m in tree_map[k][3:8]]),tree_map[k][-1]))
            else:
                self.stdout.info('    %s' % ','.join(tree_map[k][3:]))

        if self.num_trees_considered == 0:
            self.stdout.info('\npwk was not provided with any trees.')

        return tree_map

    def processParams(self, params, tree_map, skip):
        param_lines = open(params, 'r').readlines()
        if len(param_lines) < 3 + skip:
            self.output("File '%s' does not look like a parameter file (too few lines)")
        else:
            param_headers = param_lines[1].split()
        param_map = self.transformParameters(param_headers, param_lines, tree_map, skip)
        return param_map

    def estimateMargLike(self, params, trees):
        marglike = None
        skip = self.opts.skip

        # Open the tree file and process lines
        tree_map = self.processTrees(trees, skip)

        # Create sorted list of tree IDs with tree having most samples first
        sorted_tree_IDs = [(tree_map[k][0],k) for k in tree_map.keys()]
        sorted_tree_IDs.sort(cmp = lambda x,y: cmp(y[0],x[0]))

        # Open the parameter file and process lines, separating parameter vectors according to tree topology
        # and creating param_map in which keys are tree IDs and values are lists of transformed parameter vectors
        param_map = self.processParams(params, tree_map, skip)

        #for (count,key) in sorted_tree_IDs:
        #    print '%d = %d' % (tree_map[key][0], len(param_map[key]))


        # *** BEGIN AGAIN HERE ***

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
        self.phycassert(len(self.opts.params) > 0, "Must specify a parameter file ('params') when invoking the pwk command")
        self.stdout.phycassert(self.opts.trees, 'trees cannot be None or empty when pwk method is called')
        marglike = self.estimateMargLike(self.opts.params, self.opts.trees)
        return marglike

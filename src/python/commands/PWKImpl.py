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

    def logJacobianForEdgelens(self, log_edgelens):
        return sum(log_edgelens)

    def logJacobianForModelParams(self, log_modelparams, free_param_names):
        log_jacobian = 0.0
        xchgsum = 1.0
        freqsum = 1.0

        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~begin'
        #print 'logJacobianForModelParams:'

        for h,logx in zip(free_param_names, log_modelparams):
            m = re.match('(\d+)_(.+)',h)
            assert m is not None
            the_subset = int(m.group(1))
            the_param_name = m.group(2)
            if the_param_name in ['edgelen_hyper', 'pinvar', 'gamma_shape', 'rAC', 'rAG', 'rAT', 'rCG', 'rCT', 'rGT', 'freqA', 'freqC', 'freqG', 'freqT']:
                try:
                    if the_param_name == 'gamma_shape':
                        log_jacobian += logx
                        #print '%12.5f %s' % (logx,the_param_name)
                    elif the_param_name == 'edgelen_hyper':
                        log_jacobian += logx
                        #print '%12.5f %s' % (logx,the_param_name)
                    elif the_param_name == 'pinvar':
                        # logx = Y = log(pinvar) - log(1 - pinvar)
                        # X = e^Y / (1 + e^Y)
                        # dX/dY = X (1-X)
                        # log(dX/dY) = log(X) - log(1-X) = Y - 2 log(1 + e^Y)
                        log_jacobian += logx - 2.0*math.log(1.0 + math.exp(logx))
                        #print '%12.5f %s' % (logx - 2.0*math.log(1.0 + math.exp(logx)),the_param_name)
                    elif the_param_name in ['freqC', 'freqG', 'freqT']:
                        log_jacobian += logx
                        freqsum += math.exp(logx)
                        #print '%12.5f %s' % (logx,the_param_name)
                    elif the_param_name in ['rAG', 'rAT', 'rCG', 'rCT', 'rGT']:
                        log_jacobian += logx
                        xchgsum += math.exp(logx)
                        #print '%12.5f %s' % (logx,the_param_name)
                except ValueError:
                    print 'ValueError for h = %s, logx = %g' % (h, logx)
        log_jacobian -= 4.0*math.log(freqsum)
        log_jacobian -= 6.0*math.log(xchgsum)

        #print '%12.5f 4.0*math.log(freqsum)' % (4.0*math.log(freqsum),)
        #print '%12.5f 6.0*math.log(xchgsum)' % (6.0*math.log(xchgsum),)
        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~end'

        return log_jacobian

    def transformParamVector(self, pvalues, pnames):
        # Log-transforms parameter values in pvalues (log ratio transform used for GTR exchangeabilities and nucleotide frequencies)
        # Returns the transformed parameter vector
        assert len(pvalues) == len(pnames), 'Number of parameters (%d) does not equal number of parameter names (%d) in transformParamVector' % (len(pvalues), len(pnames))
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
                        assert logpiA is not None, 'expecting freqA to come before freqC, freqG, and freqT'
                        p.append(logx - logpiA)
                    elif the_param_name == 'rAC':
                        logrAC = logx
                    elif the_param_name in ['rAG', 'rAT', 'rCG', 'rCT', 'rGT']:
                        assert logrAC is not None, 'expecting rAC to come before rAG, rAT, rCG, rCT, and rGT'
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
        # Creates param_map, lnL_map, lnP_map in which key=treeID and values are:
        # lists of raw parameter vectors (param_map), log-likelihoods (lnL_map), or
        # log joint priors (lnP_map)

        # Store post-burnin parameter vectors in one big list raw_params
        row_start = 2 + skip  # first line ID, second line headers
        lines_not_skipped = param_lines[row_start:]
        n = len(lines_not_skipped)
        raw_params = []
        log_likelihoods = []
        log_priors = []
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
                elif h == 'lnPrior':
                    log_priors.append(float(x))
            raw_params.append(p)

        print 'raw_params has length %d' % len(raw_params)
        print 'log_likelihoods has length %d' % len(log_likelihoods)

        # Now, create param_map, lnL_map, and lnP_map from raw_params
        param_map = {}
        lnL_map = {}
        lnP_map = {}
        for k in tree_map.keys():
            indices = tree_map[k][3]
            lnL = []
            lnP = []
            v = []
            for i in indices:
                param_vector_i = raw_params[i-1][:]
                v.append(param_vector_i)
                lnL.append(log_likelihoods[i-1])
                lnP.append(log_priors[i-1])
            param_map[k] = v
            lnL_map[k] = lnL
            lnP_map[k] = lnP

        raw_params = None
        return param_map, lnL_map, lnP_map

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
        #   key   = treeID (tuple of strings, each representing a split)
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

    def calcLogKernel(self, transformed_parameters, newick, edge_map, V):
        # Recalculates log kernel given a list of parameter values in `transformed_parameters'
        # (log-transformed edge lengths followed by log-transformed substitution model parameters)
        # A tree is built from the supplied TreeCollection variable `newick', `edge_map' provides
        # the index for each edge length in `transformed_parameters' given a string representation
        # of a split in the tree, and `V' provides the log Jacobian for the standardization.
        nedges = len(edge_map)

        # Get the cold chain, which has the tree and partition model and knows how to compute the likelihood
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        free_param_names = cold_chain.partition_model.getFreeParameterNames()
        assert len(transformed_parameters) == nedges + len(free_param_names), 'length of transformed_parameters vector (%d) does not equal the number of edge lengths (%d) plus the number of free parameters (%d) in calcLogKernel' % (len(transformed_parameters), nedges, len(free_param_names))

        # Build the tree using the supplied TreeCollection variable `newick'
        newick.buildTree(cold_chain.tree)
        cold_chain.likelihood.prepareForLikelihood(cold_chain.tree)

        # Set edge lengths (log-transformed edge lengths come first in transformed_parameters list)
        try:
            edge_lengths = [math.exp(log_edgelen) for log_edgelen in transformed_parameters[:nedges]]
        except OverflowError:
            print 'Overflow Error exponentiating log-transformed edge lengths:'
            for log_edgelen in transformed_parameters[:nedges]:
                print log_edgelen
            sys.exit('aborting')
        self.replaceEdgeLengths(cold_chain.tree, edge_lengths, edge_map)

        # Send parameters to model (model parameters come last in transformed_parameters list)
        cold_chain.partition_model.setTransformedParameters(transformed_parameters[nedges:], cold_chain.tree)
        cold_chain.likelihood.replaceModel(cold_chain.partition_model) # trigger cold_chain.likelihood.recalcRelativeRates() call and invalidation of CLAs

        # Calculate the three components of the log kernel: log likelihood, log prior, and log Jacobian
        lnL = cold_chain.likelihood.calcLnL(cold_chain.tree)
        lnP = cold_chain.partition_model.getJointPriorManager().getLogJointPrior()
        lnJs =  V.logJacobianForStandardization()
        lnJe = self.logJacobianForEdgelens(transformed_parameters[:nedges])
        lnJp = self.logJacobianForModelParams(transformed_parameters[nedges:], free_param_names)
        lnJ = lnJe + lnJp + lnJs

        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
        #print 'lnJ (edge lengths)       =',lnJe
        #print 'lnJ (substitution model) =',lnJp
        #print 'lnJ (standardization)    =',lnJs
        #print 'lnJ (total)              =',lnJ
        #print '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'

        return (lnL, lnP, lnJe, lnJp, lnJs, free_param_names)

    def debugCheckKernelOneTopology(self, lnL_map, lnP_map, param_map, tree_map, treeID):
        # get list of vectors each containing untransformed parameters
        untransformed_param_vectors = param_map[treeID]

        # get list of vectors each containing untransformed edge lengths
        untransformed_edgelens = tree_map[treeID][4]
        edge_map = tree_map[treeID][5]

        # n is the number of samples for this tree topology
        n = len(untransformed_param_vectors)
        self.phycassert(len(untransformed_edgelens) == n, 'Number of observations differs for params (%d) and trees (%d) in PWKImple.py debugCheckKernelOneTopology function' % (n, len(untransformed_edgelens)))
        self.phycassert(tree_map[treeID][0] == n, 'Number of parameter vectors (%d) does not match number of samples (%d) in PWKImple.py debugCheckKernelOneTopology function' % (n, tree_map[treeID][0]))

        # get cold chain and obtain parameter names from it
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        all_param_names = cold_chain.partition_model.getAllParameterNames()
        free_param_names = cold_chain.partition_model.getFreeParameterNames()

        # set tree topology
        newick = tree_map[treeID][1]
        newick.buildTree(cold_chain.tree)
        cold_chain.likelihood.prepareForLikelihood(cold_chain.tree)

        # loop through all trees sampled with this topology, calculating the log-likelihood and log-prior for each
        print '%d trees were sampled with this topology:' % n,newick
        print '%12s %12s %12s %12s %12s %12s %12s' % ('index', 'lnL', 'lnLcf', 'lnLdiff', 'lnP', 'lnPcf', 'lnPdiff')
        for i in range(n):
            index = tree_map[treeID][3][i]
            evect = untransformed_edgelens[i]
            pvect = untransformed_param_vectors[i]
            lnLcf = lnL_map[treeID][i]
            lnPcf = lnP_map[treeID][i]

            # set edge lengths
            self.replaceEdgeLengths(cold_chain.tree, evect, edge_map)

            # transform parameter vector and send to model
            param_vect = self.transformParamVector(pvect, all_param_names)
            lnJ_log_transformation = self.logJacobianForModelParams(param_vect, free_param_names)
            cold_chain.partition_model.setTransformedParameters(param_vect, cold_chain.tree)
            cold_chain.likelihood.replaceModel(cold_chain.partition_model) # trigger cold_chain.likelihood.recalcRelativeRates() call and invalidation of CLAs

            # uncomment to examine value of each parameter
            #param_values = cold_chain.partition_model.getUntransformedParameters()
            #print '\nSanity check:'
            #print 'model:',cold_chain.partition_model.getModel(0).getModelName()
            #print 'tree:',cold_chain.tree.makeNewick()
            #for name,value in zip(all_param_names,param_values):
            #    print '>>> ',name,'-->',value
            #raw_input('..')

            # calculate log likelihood and compare with value in param file
            lnL = cold_chain.likelihood.calcLnL(cold_chain.tree)
            jpm = cold_chain.partition_model.getJointPriorManager()
            lnP = jpm.getLogJointPrior()
            print '%12d %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f' % (index, lnL, lnLcf, math.fabs(lnL-lnLcf), lnP, lnPcf, math.fabs(lnP-lnPcf))

    def debugShowUntransformedAndTransformedParams(row, untransformed_param_vector, transformed_param_vector):
        print '\nrow = %d' % row
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        tmp_free_param_names = cold_chain.partition_model.getFreeParameterNames()
        for paramname,uval,tval in zip(tmp_free_param_names, untransformed_param_vector, transformed_param_vector):
            print '%20s %12.5f %12.5f' % (paramname, uval, tval)
        #raw_input('..')

    # calcVolume returns volume of hypersphere in n dimensions with radius r
    # For n=2, returns pi r^2, the area of a circle with radius r
    # For n=3, returns (4/3) pi r^3, the volume of a sphere with radius r
    def calcVolume(self, n, r):
        log_volume = 1.0*n*math.log(r) + 0.5*n*math.log(math.pi) - math.lgamma(0.5*n + 1.0)
        return math.exp(log_volume);

    # Returns representative log kernel value from shell whose midpoint is at radius r
    # V should be the likelihood.VarCovBase object used to standardize the samples
    def getRepresentativeLogKernelPaul(self, r, delta, V, newick, edge_map):
        nsamples = V.calcRepresentativeForShell(r, delta)
        log_ratio_sum = V.getLogRatioSum()
        if nsamples > 0:
            log_kernel = V.getRepresentativeLogKernel()
            #print 'average log kernel =', log_kernel

            if False:
                log_params = V.getRepresentativeParamVect()
                check_lnL, check_lnP, check_lnJ_edgelen, check_lnJ_substmodel, check_lnJ_standardization, free_param_names = self.calcLogKernel(log_params, newick, edge_map, V)
                check_lnKernel = check_lnL + check_lnP + check_lnJ_edgelen + check_lnJ_substmodel + check_lnJ_standardization
                print check_lnKernel,'<-- log kernel of chosen point'
                print free_param_names
                raw_input('..debug stop..')

            if False:
                log_likelihood  = V.getRepresentativeLogLikelihood()
                log_prior       = V.getRepresentativeLogPrior()
                log_jacobian_edgelen = V.getRepresentativeLogJacobianEdgelen()
                log_jacobian_substmodel = V.getRepresentativeLogJacobianSubstmodel()
                log_jacobian_standardization = V.getRepresentativeLogJacobianStandardization()

                log_params = V.getRepresentativeParamVect()
                check_lnL, check_lnP, check_lnJ_edgelen, check_lnJ_substmodel, check_lnJ_standardization, free_param_names = self.calcLogKernel(log_params, newick, edge_map, V)
                print
                print 'check_lnL                        =', check_lnL
                print 'log_likelihood                   =', log_likelihood
                print '--> difference                   =', math.fabs(log_likelihood-check_lnL)

                print 'check_lnP                        =', check_lnP
                print 'log_prior                        =', log_prior
                print '--> difference                   =', math.fabs(log_prior-check_lnP)

                print 'check_lnJ_edgelen                =', check_lnJ_edgelen
                print 'log_jacobian_edgelen             =', log_jacobian_edgelen
                print '--> difference                   =', math.fabs(log_jacobian_edgelen - check_lnJ_edgelen)

                print 'check_lnJ_substmodel             =', check_lnJ_substmodel
                print 'log_jacobian_substmodel          =', log_jacobian_substmodel
                print '--> difference                   =', math.fabs(log_jacobian_substmodel - check_lnJ_substmodel)

                print 'check_lnJ_standardization        =', check_lnJ_standardization
                print 'log_jacobian_standardization     =', log_jacobian_standardization
                print '--> difference                   =', math.fabs(log_jacobian_standardization - check_lnJ_standardization)

                print 'check kernel                     =', check_lnL + check_lnP + check_lnJ_edgelen + check_lnJ_substmodel + check_lnJ_standardization
                print 'log_kernel                       =', log_kernel
                print '--> difference                   =', math.fabs(log_kernel - check_lnL - check_lnP - check_lnJ_edgelen - check_lnJ_substmodel - check_lnJ_standardization)

                for i,logp in enumerate(log_params):
                    print '%15.5f' % logp,
                    if i > 8:
                        print ' <-- %s' % free_param_names[i-9]
                    else:
                        print
                raw_input('..debug stop..')
        else:
            log_params = V.getRepresentativeParamVect()
            lnL, lnP, lnJ_edgelen, lnJ_substmodel, lnJ_standardization, free_param_names = self.calcLogKernel(log_params, newick, edge_map, V)
            log_kernel = lnL + lnP + lnJ_edgelen + lnJ_substmodel + lnJ_standardization
            #print 'manufactured log kernel =', log_kernel
            #print '  lnL                 =', lnL
            #print '  lnP                 =', lnP
            #print '  lnJ_edgelen         =', lnJ_edgelen
            #print '  lnJ_substmodel      =', lnJ_substmodel
            #print '  lnJ_standardization =', lnJ_standardization

        return nsamples, log_kernel, log_ratio_sum

    # Returns representative log kernel value from shell whose midpoint is at radius r
    # V should be the likelihood.VarCovBase object used to standardize the samples
    def getRepresentativeLogKernelMing(self, r, delta, V, newick, edge_map):
        nsamples = V.calcRepresentativeForShell(r, delta)
        log_ratio_terms = V.getLogRatioVect()
        if nsamples > 0:
            log_kernel = V.getRepresentativeLogKernel()
        else:
            log_params = V.getRepresentativeParamVect()
            lnL, lnP, lnJ_edgelen, lnJ_substmodel, lnJ_standardization, free_param_names = self.calcLogKernel(log_params, newick, edge_map, V)
            log_kernel = lnL + lnP + lnJ_edgelen + lnJ_substmodel + lnJ_standardization

        return nsamples, log_kernel, log_ratio_terms

    def processOneTopologyPaul(self, lnL_map, lnP_map, param_map, tree_map, treeID):
        debug_check_likelihoods = False

        #print '*** new tree topology ***'

        # Let n be the number of trees sampled for tree topology indexed by treeID
        n = tree_map[treeID][0]

        # Get cold chain and obtain parameter names from it
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        all_param_names = cold_chain.partition_model.getAllParameterNames()
        free_param_names = cold_chain.partition_model.getFreeParameterNames()

        if debug_check_likelihoods:
            print '~~~~ begin ~~~'
            newick = tree_map[treeID][1]
            newick.buildTree(cold_chain.tree)
            cold_chain.likelihood.prepareForLikelihood(cold_chain.tree)

        # Get list of vectors each containing untransformed parameters
        untransformed_param_vectors = param_map[treeID]

        # get list of vectors each containing untransformed edge lengths
        untransformed_edgelens = tree_map[treeID][4]
        edge_map = tree_map[treeID][5]
        nedges = len(untransformed_edgelens[0])
        assert nedges == len(edge_map.items()), 'length of untransformed_edgelens (%d) differs from number of items in edge_map (%d)' % (len(untransformed_edgelens[0]), len(edge_map.items()))

        # Let p be the length of each untransformed parameter vector
        #print 'untransformed_edgelens.__class__.__name__ =',untransformed_edgelens.__class__.__name__
        #print 'len(untransformed_edgelens) =',len(untransformed_edgelens)
        p = len(untransformed_edgelens[0]) + cold_chain.partition_model.getNumFreeParameters()

        # Create a VarCovBase object to manage the standardization
        V = likelihood.VarCovBase('topo%d' % n, n, p)

        #raw_input('..debugTestStandardization starting..')
        #V.debugTestStandardization();
        #raw_input('..debugTestStandardization finished..')

        edgelen_names = [""]*nedges
        for key in edge_map.keys():
            edgelen_names[edge_map[key]] = key
        V.setParamNames(edgelen_names + list(free_param_names))

        # Add all parameter vectors to V
        row = 0
        for q,e,lnl in zip(untransformed_param_vectors, untransformed_edgelens, lnL_map[treeID]):
            # Log (or log-ratio) transform parameter vector
            transformed_param_vector = self.transformParamVector(q, all_param_names)
            lnJ_substitution_model = self.logJacobianForModelParams(transformed_param_vector, free_param_names)

            # Log transform edge lengths
            transformed_edgelen_vector = [math.log(edgelen) for edgelen in e]

            # Add log-jacobians for edge length transformations to lnJ_log_transformation
            lnJ_edge_lengths = sum(transformed_edgelen_vector)

            # Copy transformed parameter vector, log likelihood, log prior, and log jacobian to correct row of V
            all_params = transformed_edgelen_vector + transformed_param_vector
            V.copyParamVector(row, all_params, lnL_map[treeID][row], lnP_map[treeID][row], lnJ_edge_lengths, lnJ_substitution_model)

            if debug_check_likelihoods and row < 1:
                test_params = all_params[:]
                self.replaceEdgeLengths(cold_chain.tree, e, edge_map)
                cold_chain.partition_model.setTransformedParameters(transformed_param_vector, cold_chain.tree)
                cold_chain.likelihood.replaceModel(cold_chain.partition_model) # trigger cold_chain.likelihood.recalcRelativeRates() call and invalidation of CLAs
                lnL = cold_chain.likelihood.calcLnL(cold_chain.tree)
                print 'free_param_names =',free_param_names
                print 'nedges =',nedges
                print 'len(e) =',len(e)
                print 'len(transformed_param_vector) =',len(transformed_param_vector)
                print '%6d %12.5f %12.5f %12.5f' % (row,lnl,lnL,math.fabs(lnl-lnL))
                for v in transformed_edgelen_vector:
                    print v
                for v,nm in zip(transformed_param_vector,free_param_names):
                    print v,' (%s)' % nm
                print

            row += 1

        # Standardize samples by subtracting the mean vector and multiplying by the variance-covariance matrix raised to the power -0.5
        furthest_radius = V.standardizeSamples()

        if debug_check_likelihoods:
            print '~~~~ after standardization ~~~'
            # check likelihoods again using de-standardized parameter vectors
            for row in range(1):
                transformed_parameter_vector = V.destandardizeSample(row)
                print 'free_param_names =',free_param_names
                print 'nedges =',nedges
                print 'len(transformed_parameter_vector) =',len(transformed_parameter_vector)
                e = [math.exp(v) for v in transformed_parameter_vector[:nedges]]
                self.replaceEdgeLengths(cold_chain.tree, e, edge_map)
                cold_chain.partition_model.setTransformedParameters(transformed_parameter_vector[nedges:], cold_chain.tree)
                cold_chain.likelihood.replaceModel(cold_chain.partition_model) # trigger cold_chain.likelihood.recalcRelativeRates() call and invalidation of CLAs
                lnL = cold_chain.likelihood.calcLnL(cold_chain.tree)
                lnl = lnL_map[treeID][row]
                print '%6d %12.5f %12.5f %12.5f' % (row,lnl,lnL,math.fabs(lnl-lnL))
                for v in transformed_parameter_vector:
                    print v
                print
            print '~~~~ end ~~~'
            raw_input('..debug stop..')

        # Calculate delta_r, which equals half the "thickness" of each of the K shells
        K = self.opts.shells
        furthest_radius *= float(self.opts.reach)
        delta_r = 0.5*furthest_radius/K;

        debugf = open('debug_topology_details.txt', 'a')
        debugf.write('\n*** topology with %d samples *** \n' % n)
        debugf.write('furthest_radius = %.8f\n' % furthest_radius)

        # Loop through shells to compute the contribution of this topology to the overall log marginal likelihood
        log_numerator   = []
        log_denominator = []
        volume_area_k_minus_1 = 0.0
        debugf.write('%20s %20s %20s %20s %20s %20s\n' % ('shell','radius','samples','log kernel repr', 'log ratio sum', 'log volume'))
        for k in range(K):
            # determine midpoint radius and volume of this shell
            r_k = furthest_radius*(1.0*(k+1)/K - 1.0/(2.0*K))
            volume_area_k  = self.calcVolume(p, r_k + delta_r)
            log_volume_shell_k = math.log(volume_area_k - volume_area_k_minus_1)

            # See if any samples fall within this shell
            # if so, return value will be a list with one element equalling the average log kernel of these samples
            # if not, return value will be a list representing a parameter vector of a representative point
            newick = tree_map[treeID][1]
            nsamples_this_shell, log_kernel_representative, log_ratio_sum = self.getRepresentativeLogKernelPaul(r_k, delta_r, V, newick, edge_map)

            # debugging
            if n == 77695 and k == 49:
                paramvects = V.getRepresentativesForShell(10, r_k, delta_r)
                nelements = len(paramvects)
                print '~~~~~~ debugging ~~~~~'
                print 'r =',r_k
                print 'lower =',r_k-delta_r
                print 'upper =',r_k+delta_r
                print 'p =',p
                print 'nelements =',nelements
                print '%20s' % 'log(kernel)',
                for nm in edgelen_names:
                    print '%20s' % nm,
                for nm in free_param_names:
                    print '%20s' % nm,
                print
                for which in range(0, nelements, p):
                    check_lnL, check_lnP, check_lnJ_edgelen, check_lnJ_substmodel, check_lnJ_standardization, free_param_names = self.calcLogKernel(paramvects[which:which+p], newick, edge_map, V)
                    check_lnKernel = check_lnL + check_lnP + check_lnJ_edgelen + check_lnJ_substmodel + check_lnJ_standardization
                    print '%20.5f' % check_lnKernel,
                    for pp in paramvects[which:which+p]:
                        print '%20.5f' % pp,
                    print
                raw_input('..stopped at last shell of topology 77695..')

            if nsamples_this_shell > 0: #temporary!
                log_numerator.append(log_ratio_sum)
                log_denominator.append(log_kernel_representative + log_volume_shell_k)

                debugf.write('%20d %20.5f %20d %20.5f %20.5f %20.5f\n' % (k+1,r_k,nsamples_this_shell,log_kernel_representative,log_ratio_sum,log_volume_shell_k))
        debugf.close()

        log_sum_of_ratios = self.calculateSumTermsOnLogScale(log_numerator)
        log_sum_of_volumes = self.calculateSumTermsOnLogScale(log_denominator)
        log_contribution = log_sum_of_ratios - log_sum_of_volumes

        tmpf = open('debug_breakdown.txt', 'a')
        tmpf.write('%20d %20.5f %20.5f %35.5f\n' % (n, log_sum_of_ratios, log_sum_of_volumes, log_contribution))
        tmpf.close()

        return log_sum_of_ratios, log_sum_of_volumes, log_contribution

    def processOneTopologyMing(self, lnL_map, lnP_map, param_map, tree_map, treeID):
        # Let n be the number of trees sampled for tree topology indexed by treeID
        n = tree_map[treeID][0]
        #log_tree_posterior = math.log(n) - math.log(self.num_trees_considered)

        # Get cold chain and obtain parameter names from it
        cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
        all_param_names = cold_chain.partition_model.getAllParameterNames()
        free_param_names = cold_chain.partition_model.getFreeParameterNames()

        # Get list of vectors each containing untransformed parameters
        untransformed_param_vectors = param_map[treeID]

        # get list of vectors each containing untransformed edge lengths
        untransformed_edgelens = tree_map[treeID][4]
        edge_map = tree_map[treeID][5]
        nedges = len(untransformed_edgelens[0])
        assert nedges == len(edge_map.items()), 'length of untransformed_edgelens (%d) differs from number of items in edge_map (%d)' % (len(untransformed_edgelens[0]), len(edge_map.items()))

        # Let p be the length of each untransformed parameter vector
        p = len(untransformed_edgelens[0]) + cold_chain.partition_model.getNumFreeParameters()

        # Create a VarCovBase object to manage the standardization
        V = likelihood.VarCovBase('topo%d' % n, n, p)

        # Add all parameter vectors to V
        row = 0
        for q,e,lnl in zip(untransformed_param_vectors, untransformed_edgelens, lnL_map[treeID]):
            # Log (or log-ratio) transform parameter vector
            transformed_param_vector = self.transformParamVector(q, all_param_names)
            lnJ_substitution_model = self.logJacobianForModelParams(transformed_param_vector, free_param_names)

            # Log transform edge lengths
            transformed_edgelen_vector = [math.log(edgelen) for edgelen in e]

            # Add log-jacobians for edge length transformations to lnJ_log_transformation
            lnJ_edge_lengths = sum(transformed_edgelen_vector)

            # Copy transformed parameter vector, log likelihood, log prior, and log jacobian to correct row of V
            all_params = transformed_edgelen_vector + transformed_param_vector
            V.copyParamVector(row, all_params, lnL_map[treeID][row], lnP_map[treeID][row], lnJ_edge_lengths, lnJ_substitution_model)

            row += 1

        # Standardize samples by subtracting the mean vector and multiplying by the variance-covariance matrix raised to the power -0.5
        furthest_radius = V.standardizeSamples()

        # Calculate delta_r, which equals half the "thickness" of each of the K shells
        K = self.opts.shells
        furthest_radius *= float(self.opts.reach)
        delta_r = 0.5*furthest_radius/K;
        print 'furthest_radius =',furthest_radius

        # same as processOneTopologyPaul to here

        # Loop through shells to compute the contribution of this topology to the overall log marginal likelihood
        # Note that numerator and denominator refer to the formula for 1/c, so these will ultimately need to be swapped to estimate c
        log_numerator_terms   = []
        log_denominator_terms = []
        volume_area_k_minus_1 = 0.0
        print '%12s %12s %12s %s' % ('shell','radius','samples','representative log(kernel)')
        for k in range(K):
            # determine midpoint radius and volume of this shell
            r_k = furthest_radius*(1.0*(k+1)/K - 1.0/(2.0*K))
            volume_area_k  = self.calcVolume(p, r_k + delta_r)
            log_volume_shell_k = math.log(volume_area_k - volume_area_k_minus_1)

            # See if any samples fall within this shell
            # if so, return value will be a list with one element equalling the average log kernel of these samples
            # if not, return value will be a list representing a parameter vector of a representative point
            newick = tree_map[treeID][1]
            nsamples_this_shell, log_kernel_representative, log_ratio_terms = self.getRepresentativeLogKernelMing(r_k, delta_r, V, newick, edge_map)

            log_numerator_terms.extend(log_ratio_terms)
            log_denominator_terms.append(log_kernel_representative + log_volume_shell_k)

            print '%12d %12.5f %12d %.5f' % (k+1,r_k,nsamples_this_shell,log_kernel_representative)

        return (log_numerator_terms,log_denominator_terms)

    def calculateSumTermsOnLogScale(self, log_terms):
        max_log_term = max(log_terms)
        sum_of_terms = 0.0;
        for log_term in log_terms:
            sum_of_terms += math.exp(log_term - max_log_term);
        sum_on_log_scale = max_log_term + math.log(sum_of_terms)
        return sum_on_log_scale

    def estimateMargLike(self, params, trees):
        marglike = None
        skip = self.opts.skip

        # Open the tree file and process lines, creating tree_map and setting self.num_trees_considered
        # to equal the number of trees not skipped.
        #
        # Description of tree_map:
        #
        #  key   = tuple of strings (aka treeID), each string representing a single split
        #
        #  value = list containing the following elements:
        #     0: number of trees with this topology (n)
        #     1: TreeDescription (not a string, see class TreeDescription in _NexusReader.py)
        #     2: mean tree length
        #     3: list (length n) of indices of samples having this topology (start at 1, not 0)
        #     4: list (length n) of lists (length nedges) of untransformed edge lengths
        #     5: edge_map (ensures that edge lengths are always in the same order in the parameter vector)
        #
        # Description of edge_map:
        #
        #  key   = split, i.e. a string representing a single bipartition
        #
        #  value = index into tree_map[4] lists
        #
        tree_map = self.processTrees(trees, skip)

        # Create sorted list of (n,treeID) tuples with tree having most samples first.
        # The lambda function is necessary because the standard sort would put least sampled tree first.
        sorted_tree_IDs = [(tree_map[k][0],k) for k in tree_map.keys()]
        sorted_tree_IDs.sort(cmp = lambda x,y: cmp(y[0],x[0]))

        # Open the parameter file and process lines, separating parameter vectors according to tree topology
        # and creating param_map in which keys are treeIDs and values are lists of parameter vectors
        #
        # Description of param_map:
        #
        #  key   = treeID
        #
        #  value = list of parameter vectors
        #
        #
        # Description of lnL_map:
        #
        #  key   = treeID
        #
        #  value = list of log-likelihood values for one tree topology
        #
        #
        # Description of lnP_map:
        #
        #  key   = treeID
        #
        #  value = list of log-prior values for one tree topology
        #
        param_map, lnL_map, lnP_map = self.processParams(params, tree_map, skip)

        if True:
            # This version works, but is not as elegant as it could be

            open('debug_topology_details.txt', 'w').close()
            tmpf = open('debug_breakdown.txt', 'w')
            tmpf.write('%20s %20s %20s %35s\n' % ('trees','log(sum qratio)','log(sum volume)','log(sum qratio) - log(sum volume)'))
            tmpf.write('%20s %20s %20s %35s\n' % ('-------------------', '-------------------', '-------------------', '----------------------------------'))
            tmpf.close()

            # Loop over tree topologies and, for each topology having the minimum number of samples,
            # compute the term of the denominator of the PWK estimation formula corresponding to this topology
            log_denom_terms = []
            tmp_log_sum_of_ratios_list = []
            tmp_log_sum_of_volumes_list = []
            num_trees_included = 0
            for n,treeID in sorted_tree_IDs:
                if n >= self.opts.minsample:
                    #self.debugCheckKernelOneTopology(lnL_map, lnP_map, param_map, tree_map, treeID)

                    print '\n*** Including topology with %d samples ***' % n

                    # Number of trees included equals self.num_trees_considered iff sample size is
                    # at least self.opts.minsample for all tree topologies
                    num_trees_included += n

                    # Log marginal posterior probability of this tree topology
                    log_posterior_this_topology = math.log(n) - math.log(self.num_trees_considered)

                    # Compute contribution of this topology to the marginal likelihood
                    log_sum_of_ratios, log_sum_of_volumes, log_topology_contribution = self.processOneTopologyPaul(lnL_map, lnP_map, param_map, tree_map, treeID)

                    # log_sum_of_ratios, log_sum_of_volumes, log_contribution

                    tmp_log_sum_of_ratios_list.append(log_sum_of_ratios)
                    tmp_log_sum_of_volumes_list.append(log_sum_of_volumes)
                    log_denom_terms.append(log_topology_contribution)
                else:
                    print '\n*** Excluding topology with %d samples ***' % n

            log_sum_of_ratios   = self.calculateSumTermsOnLogScale(tmp_log_sum_of_ratios_list)
            log_sum_of_volumes  = self.calculateSumTermsOnLogScale(tmp_log_sum_of_volumes_list)
            log_denominator     = self.calculateSumTermsOnLogScale(log_denom_terms)

            tmpf = open('debug_breakdown.txt', 'a')
            tmpf.write('%20s %20s %20s %35s\n' % ('-------------------', '-------------------', '-------------------', '----------------------------------'))
            tmpf.write('%20d %20.5f %20.5f %35.5f\n' % (num_trees_included, log_sum_of_ratios, log_sum_of_volumes, log_denominator))
            tmpf.write('\ntotal sample size      = %d\n' % self.num_trees_considered)
            logN = math.log(self.num_trees_considered)
            tmpf.write('log(total sample size) = %.5f\n' % logN)
            tmpf.write("sum of ratios method = %.5f = %.5f - %.5f\n" % (logN - log_denominator, logN, log_denominator))
            tmpf.write("ratio of sums method = %.5f = %.5f - (%.5f - (%.5f))\n" % (logN - (log_sum_of_ratios - log_sum_of_volumes), logN, log_sum_of_ratios, log_sum_of_volumes))
            tmpf.close()

            print 'self.num_trees_considered =',self.num_trees_considered
            print 'num_trees_included        =',num_trees_included

            # Note that self.num_trees_considered cancels in the version of the PWK estimator that uses posterior probabilities in the numerator
            if False:
                print '@@@@@ Paul (1/{\cal N}) @@@@@'
                log_numerator = math.log(num_trees_included)
                log_marginal_likelihood = log_numerator - log_denominator
            else:
                print '@@@@@ Paul (1/N) @@@@@'
                log_numerator = math.log(self.num_trees_considered)
                log_marginal_likelihood = log_numerator - log_denominator

            #

        else:
            # This version is more elegant, but may not be as accurate

            # Note that numerator and denominator below apply to the 1/c estimator
            log_numerator_terms = []
            log_denominator_terms = []
            num_trees_included = 0
            for n,treeID in sorted_tree_IDs:
                if n >= self.opts.minsample:
                    print '\n*** Including topology with %d samples ***' % n

                    # Number of trees included equals self.num_trees_considered iff sample size is
                    # at least self.opts.minsample for all tree topologies
                    num_trees_included += n

                    # Compute contribution of this topology to the marginal likelihood
                    numterms, denomterms = self.processOneTopologyMing(lnL_map, lnP_map, param_map, tree_map, treeID)
                    log_numerator_terms.extend(numterms)
                    log_denominator_terms.extend(denomterms)
                else:
                    print '\n*** Excluding topology with %d samples ***' % n

            #cold_chain = pwk_mcmc_impl.mcmc_manager.getColdChain()
            #jpm = cold_chain.partition_model.getJointPriorManager()
            #log_tree_prior_prob = jpm.getLogTopologyPrior()

            log_numerator = self.calculateSumTermsOnLogScale(log_numerator_terms) - math.log(self.num_trees_considered)
            log_denominator = self.calculateSumTermsOnLogScale(log_denominator_terms)

            # Paul's modification
            #log_denominator += math.log(num_trees_included) - math.log(self.num_trees_considered)

            # reverse numer. and denom. in order to estimate c rather than 1/c
            log_marginal_likelihood = log_denominator - log_numerator
            print '@@@@@ Ming @@@@@'

            print 'self.num_trees_considered =',self.num_trees_considered
            print 'num_trees_included        =',num_trees_included
            print 'math.log(self.num_trees_considered) =',math.log(self.num_trees_considered)
            print 'math.log(num_trees_included)        =',math.log(num_trees_included)
            print 'log difference            =',math.log(self.num_trees_considered) - math.log(num_trees_included)

        return log_marginal_likelihood

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the params and trees files specified, then 
        estimates the log(marginal likelihood) using the PWK method, using
        the number of shells specified.

        """
        global pwk_mcmc_impl

        self.phycassert(len(self.opts.params) > 0, "Must specify a parameter file ('params') when invoking the pwk command")
        self.stdout.phycassert(self.opts.trees, 'trees cannot be None or empty when pwk method is called')

        # Set mcmc.doing_pwk to True just long enough to grab the MCMCImpl object (also inhibits mcmc() from running an MCMC analysis)
        mcmc.doing_pwk = True
        pwk_mcmc_impl = mcmc()
        pwk_mcmc_impl.setup()
        mcmc.doing_pwk = False

        logmarglike = self.estimateMargLike(self.opts.params, self.opts.trees)
        print '~~> log(marginal likelihood) =',logmarglike
        return logmarglike

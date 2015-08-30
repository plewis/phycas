import os,sys,math,random,re
#from phycas.treeviewer import *
from phycas import *
from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions

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

class RefDistImpl(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Creates reference distribution for Stepping-stone anayses based on
    a previous MCMC sample from the posterior distribution.

    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes RefDistImpl object.

        """
        CommonFunctions.__init__(self, opts)
        self.skip = None
        self.epsilon = None
        self.rooted = None
        self.refdistf = None
        self.rooted_trees = False
        self.out_refdist_file_name = None
        self.min_sample_size = 1

    def _openRefDistFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Opens the reference distribution definition file.

        """
        self.phycassert(self.refdistf is None, 'Attempt made to open reference distribution definition file, but it is already open!')
        self.out_refdist_file_name = self.opts.out.refdistfile._getFilename()
        try:
            self.refdistf = self.opts.out.refdistfile.open(self.stdout)
        except:
            print '*** Attempt to open reference distribution definition file (%s) failed.' % self.out_refdist_file_name

        if self.refdistf:
            print 'Reference distribution definition file was opened successfully'

    def _closeRefDistFile(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Closes the reference distribution definition file.

        """
        self.phycassert(self.refdistf is not None, 'Attempt made to close reference distribution definition file, but it is not open!')
        self.refdistf.close()
        self.refdistf = None

    def _fitLognormalRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Lognormal reference distribution for a univariate sample x.

        """
        sum_of_log_values = 0.0
        sum_of_log_squares = 0.0
        nvalid = 0
        for v in x:
            if v > 0:
                nvalid += 1
                logv = math.log(v)
                sum_of_log_values += logv
                sum_of_log_squares += logv*logv

        self.phycassert(nvalid > self.min_sample_size, 'sample size (%d) must be at least %d to create a Lognormal reference distribution' % (nvalid, self.min_sample_size))
        logmean = sum_of_log_values/float(nvalid)
        logvar = (sum_of_log_squares - float(nvalid)*logmean*logmean)/float(nvalid-1)
        logsd = math.sqrt(logvar)
        d = Lognormal(logmean, logsd)
        return d.__str__()

    def _fitInverseGammaRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates an InverseGamma reference distribution for a univariate sample
        x.

        """
        self.phycassert(len(x) > self.min_sample_size, 'sample size (%d) must be at least %d to create an InverseGamma reference distribution' % (len(x), self.min_sample_size))
        n = len(x)
        sum_of_values = 0.0
        sum_of_squares = 0.0
        for v in x:
            sum_of_values += v
            sum_of_squares += v*v

        mean = sum_of_values/float(n)
        variance = (sum_of_squares - float(n)*mean*mean)/float(n-1)
        shape = math.pow(mean, 2.0)/variance + 2.0
        inverse_scale = math.pow(mean, 3.0)/variance + mean
        d = InverseGamma(shape, 1.0/inverse_scale)
        return d.__str__()

    def _fitGammaRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Gamma reference distribution for a univariate sample x.

        """
        self.phycassert(len(x) > self.min_sample_size, 'sample size (%d) must be at least %d to create a Gamma reference distribution' % (len(x), self.min_sample_size))
        n = len(x)
        sum_of_values = 0.0
        sum_of_squares = 0.0
        for v in x:
            sum_of_values += v
            sum_of_squares += v*v

        mean = sum_of_values/float(n)
        variance = (sum_of_squares - float(n)*mean*mean)/float(n-1)
        scale = variance/mean
        self.phycassert(scale > 0.0, 'scale (%g) must be strictly positive to compute a Gamma reference distribution' % scale)
        shape = mean/scale
        d = Gamma(shape, scale)
        return d.__str__()

    def _fitBetaRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Beta reference distribution for a univariate sample x.
        Let a, b be the parameters of a Beta(a,b) and let phi = a + b.
        Note that:
            mean = a/phi
          1-mean = b/phi
        variance = a*b/[phi^2*(phi + 1)]
        Letting z = mean*(1-mean)/variance,
        phi can be estimated as z - 1
        Now, a = mean*phi and b = (1-mean)*phi

        """
        self.phycassert(len(x) > self.min_sample_size, 'sample size (%d) must be at least %d to create a Beta reference distribution' % (len(x), self.min_sample_size))
        n = len(x)
        sum_of_values = 0.0
        sum_of_squares = 0.0
        for v in x:
            sum_of_values += v
            sum_of_squares += v*v

        mean = sum_of_values/float(n)
        variance = (sum_of_squares - float(n)*mean*mean)/float(n-1)

        phi = mean*(1.0-mean)/variance - 1.0
        a = phi*mean
        b = phi*(1.0 - mean)

        d = Beta(a, b)
        return d.__str__()

    def _estimateTreeLengthPriorDirichlet(self, n_internal, mean_internal, var_internal, n_external, mean_external, var_external):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Estimates the alpha and c parameters of the compound Dirichlet tree
        length prior from the means, variances, and numbers of external and
        internal edge length proportions. The alpha parameter is the
        Dirichlet parameter common to all external edge length proportions,
        and c is the ratio of internal to external edge length proportions.
        In an unrooted tree of 30 tips, there are 2*30 - 3 = 57 edges, of
        which 30 are external and 57 - 30 = 27 are internal. In this case,
        the function should be called with n_internal = 27 and
        n_external = 30. The return value is the 2-tuple (alpha,c).

        """
        numer  = float(n_internal)*mean_internal*mean_internal*(1.0 - mean_internal)*(1.0 - mean_internal)
        numer += float(n_external)*mean_external*mean_external*(1.0 - mean_external)*(1.0 - mean_external)
        denom  = float(n_internal)*var_internal*mean_internal*(1.0 - mean_internal)
        denom += float(n_external)*var_external*mean_external*(1.0 - mean_external)
        self.phycassert(denom > 0.0, 'denom <= zero in RefDistImpl._estimateTreeLengthPriorDirichlet')
        phi = (numer/denom) - 1.0
        alpha = mean_external*phi
        c = mean_internal/mean_external
        return (alpha,c)

    def _estimateDirichletParams(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Dirichlet parameter vector from a sample x comprising n
        m-tuples, where m is the dimension of the Dirichlet distribution and
        n is the sample size. That is, for x = [(a_1,b_1,c_1),(a_2,b_2,c_2),
        ...,(a_n,b_n,c_n)], m = 3 and this function assumes that, for each
        i = 1,2,...,n, a_i + b_i + c_i = 1.0.

        """
        # Compute sums and sums-of-squares
        dim = len(x[0])
        n = len(x)
        sums = [0.0]*dim
        ss = [0.0]*dim
        for t in x:
            for k,v in enumerate(t):
                sums[k] += v
                ss[k] += v*v

        # Estimate sample means and variances
        means = [0.0]*dim
        variances = [0.0]*dim
        for i in range(dim):
            mean = sums[i]/float(n)
            means[i] = mean
            variances[i] = (ss[i] - float(n)*mean*mean)/float(n-1)

        # Now compute the Dirichlet parameters
        # Ming-Hui Chen's least squares approach estimates phi by minimizing difference between
        # means[i]*(1 - means[i])/(phi+1) and variances[i] for all parameters. Let s = variances[i]
        # and u = means[i]. Set to zero the derivative of (s - u(1-u)/(phi+1))^2 with respect to phi.
        # This yields the estimate phi + 1 = u (1-u)/s. Summing numerators and denominators separately,
        # we have phi = [sum_i u_i^2 (1 - u_i)^2]/[sum_i s_i u_i (1 - u_i)] - 1.
        #
        # Below is a LaTeX snippet that can be executed in LaTeX-it for more detail:
        # z &= \sum_i \left( s_i - \frac{u_i (1 - u_i)}{\phi + 1} \right)^2 \\
        # \frac{\partial z}{\partial \phi}  &= \sum_i \left\{ 2 \left( s_i - \frac{u_i (1 - u_i)}{\phi + 1} \right) \left( \frac{u_i (1 - u_i)}{(\phi+1)^2} \right)\right\} = 0 \\
        # \frac{2}{(\phi+1)^2} \sum_i s_i u_i (1 - u_i) &= \frac{2}{(\phi+1)^3} \sum_i u_i^2 (1 - u_i)^2 \\
        # \phi + 1 &= \frac{\sum_i u_i^2 (1 - u_i)^2}{\sum_i s_i u_i (1 - u_i)} \\
        # \phi &= \frac{\sum_i u_i^2 (1 - u_i)^2}{\sum_i s_i u_i (1 - u_i)} - 1
        #
        # Below is a description of Paul's method-of-moments estimator, which is not as accurate as
        # Ming's least-squares estimator:
        # Let a, b, c, d be the parameters of a Dirichlet(a,b,c,d).
        # Let phi = a + b + c + d.
        # Let mu_1, mu_2, mu_3, and mu_4 be the sample means.
        # Let s_1^2, s_2^2, s_3^2, and s_4^2 be the sample variances.
        # Noting that mu_1 = a/phi, mu_2 = b/phi, mu_3 = c/phi and mu_4 = d/phi
        # and noting that s_1^2 = a*(b + c + d)/[phi^2*(phi + 1)], etc., and
        # letting z = s_1^2/mu_1 + s_2^2/mu_2 + s_3^2/mu_3 + s_4^2/mu_4,
        # phi can be estimated as (3/z) - 1, where the 3 is really k-1
        # and k is the number of Dirichlet parameters. Now,
        # a = mu_1*phi, b = mu_2*phi, c = mu_3*phi and d = mu_4*phi
        numer_sum = 0.0
        denom_sum = 0.0
        for i in range(dim):
            u = means[i]
            s = variances[i]
            numer_sum += math.pow(u*(1.0 - u),2.0)
            denom_sum += s*u*(1.0 - u)

        phi = (numer_sum/denom_sum) - 1.0
        params = []
        for i in range(dim):
            curr_param = phi*means[i]
            params.append(curr_param)
        return params

    def _fitDirichletRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a Dirichlet reference distribution for a multivariate sample
        x, which is expected to be a list of tuples representing the posterior
        sample. Note: For x = [(a_1,b_1,c_1),(a_2,b_2,c_2),...,(a_n,b_n,c_n)],
        this function assumes that a_i + b_i + c_i = 1.0 for each i in [1..n].
        """
        # First compute the sample means and variances using the data stored in mv_ref_dist
        self.phycassert(len(x) > self.min_sample_size, 'sample size (%d) must be at least %d to create a Dirichlet reference distribution' % (len(x), self.min_sample_size))
        dim = len(x[0])
        self.phycassert(dim > 1, 'Dirichlet reference distribution is only appropriate for multivariate samples')

        params = self._estimateDirichletParams(x)

        d = Dirichlet(params)
        return d.__str__()

    def _fitRelRateRefDist(self, x):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a relative rate reference distribution for a multivariate
        sample x, which is expected to be a list of tuples representing the
        posterior sample. Note: For x = [(a_1,b_1,c_1),(a_2,b_2,c_2),...,
        (a_n,b_n,c_n)], and for weights p_1, p_2, and p_3 representing e.g.
        the fraction of sites in each partition, this function assumes that
        p_1*a_i + p_2*b_i + p_3*c_i = 1.0 for each i in [1..n]. That is, the
        mean relative rate is always 1.0.
        """
        # First compute the sample means and variances using the data stored in mv_ref_dist
        self.phycassert(len(x) > self.min_sample_size, 'sample size (%d) must be at least %d to create a relative rate reference distribution in RefDistImpl._fitRelRateRefDist' % (len(x), self.min_sample_size))
        dim = len(x[0])
        self.phycassert(dim > 1, 'Relative rate reference distributions are only appropriate for multivariate samples in RefDistImpl._fitRelRateRefDist')

        # x vector comprises tuples of relative rates, but to estimate Dirichlet parameters
        # the elements of each tuple must sum to 1. Transform x to Dirichlet tuples and let
        # self._estimateDirichletParams do the work from that point
        subset_proportions = partition.getSubsetProportions()
        self.phycassert(dim == len(subset_proportions), 'Vector of subset proportions has length %d, which should be identical to the number of relative rates %d in RefDistImpl._fitRelRateRefDist' % (len(subset_proportions), dim))
        y = []
        for rel_rate_tuple in x:
            dirichlet_tuple = tuple([rel_rate_tuple[i]*subset_proportions[i] for i in range(dim)])
            self.phycassert(math.fabs(sum(dirichlet_tuple) - 1.0) < 1.e-8, 'Relative rates are expected to have mean 1.0 in RefDistImpl._fitRelRateRefDist')
            y.append(dirichlet_tuple)

        params = self._estimateDirichletParams(y)

        d = RelativeRateDistribution(params, subset_proportions)
        return d.__str__()

    def _processParamSample(self, paramfname):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the parameter sample file and returns a
        reference distribution definition string.

        """
        skip = self.opts.skip
        lines = open(paramfname, 'r').readlines()
        headers = lines[1].split()

        # Sanity checks
        self.phycassert(len(lines) >= 3 + skip, "File '%s' does not look like a parameter file (only %d lines, but expecting at least %d lines)" % (paramfname, len(lines), skip+3))
        self.phycassert(headers[1] != 'beta', "Not expecting to see a column labeled 'beta' in parameter file; params should be a parameter file from an ordinary MCMC analysis, not one produced by a Stepping-stone analysis")

        # Create params dictionary such that, for example, params['lnL'] is a list of
        # all log-likelihood values (in the order in which they were sampled)
        params = {}
        for h in headers:
            params[h] = []
        row_start = 2 + skip  # first line ID, second line headers
        for i,line in enumerate(lines[row_start:]):
            parts = line.split()
            if len(parts) != len(headers):
                raise InvalidNumberOfColumnsError(len(parts), len(headers), i + row_start + 1)
            for h,x in zip(headers,parts):
                params[h].append(float(x))

        # Create a reference distribution for each parameter and add distribution description to param_refdist_definition
        param_refdist_definition = ''

        # Example param file headers for GTR model:
        # 1_edgelen_hyper
        # 1_freqA
        # 1_freqC
        # 1_freqG
        # 1_freqT
        # 1_gamma_var
        # 1_pinvar
        # 1_rAC
        # 1_rAG
        # 1_rAT
        # 1_rCG
        # 1_rCT
        # 1_rGT
        # Gen
        # TL
        # lnL
        # lnPrior

        sorted_headers = params.keys()[:]
        sorted_headers.sort()
        num_headers = len(sorted_headers)

        subset_rates = {}
        for i in range(num_headers):
            h = sorted_headers[i]
            #print '--> h = %s (%d)' % (h, len(params[h]))
            if '_freqA' in h:
                subset_index = int(re.match('(\d+)_freqA',h).group(1))
                freqA = params[h]
                freqC = params['%d_freqC' % subset_index]
                freqG = params['%d_freqG' % subset_index]
                freqT = params['%d_freqT' % subset_index]
                x = [(a/(a+c+g+t),c/(a+c+g+t),g/(a+c+g+t),t/(a+c+g+t)) for a,c,g,t in zip(freqA,freqC,freqG,freqT)]
                param_refdist_definition += '%d_state_freqs = %s\n' % (subset_index, self._fitDirichletRefDist(x))
            elif '_freqC' in h:
                pass
            elif '_freqG' in h:
                pass
            elif '_freqT' in h:
                pass
            elif '_rAC' in h:
                subset_index = int(re.match('(\d+)_rAC',h).group(1))
                rAC = params[h]
                rAG = params['%d_rAG' % subset_index]
                rAT = params['%d_rAT' % subset_index]
                rCG = params['%d_rCG' % subset_index]
                rCT = params['%d_rCT' % subset_index]
                rGT = params['%d_rGT' % subset_index]
                x = [(ac/(ac+ag+at+cg+ct+gt),ag/(ac+ag+at+cg+ct+gt),at/(ac+ag+at+cg+ct+gt),cg/(ac+ag+at+cg+ct+gt),ct/(ac+ag+at+cg+ct+gt),gt/(ac+ag+at+cg+ct+gt)) for ac,ag,at,cg,ct,gt in zip(rAC,rAG,rAT,rCG,rCT,rGT)]
                param_refdist_definition += '%d_relrates = %s\n' % (subset_index, self._fitDirichletRefDist(x))
            elif '_rAG' in h:
                pass
            elif '_rAT' in h:
                pass
            elif '_rCG' in h:
                pass
            elif '_rCT' in h:
                pass
            elif '_rGT' in h:
                pass
            elif '_pinvar' in h:
                subset_index = int(re.match('(\d+)_pinvar',h).group(1))
                x = params[h]
                param_refdist_definition += '%d_pinvar = %s\n' % (subset_index, self._fitBetaRefDist(x))
            elif '_gamma_shape' in h:
                subset_index = int(re.match('(\d+)_gamma_shape',h).group(1))
                x = params[h]
                param_refdist_definition += '%d_gamma_shape = %s\n' % (subset_index, self._fitLognormalRefDist(x))
            elif '_kappa' in h:
                subset_index = int(re.match('(\d+)_kappa',h).group(1))
                x = params[h]
                param_refdist_definition += '%d_kappa = %s\n' % (subset_index, self._fitGammaRefDist(x))
            elif '_edgelen_hyper' in h:
                subset_index = int(re.match('(\d+)_edgelen_hyper',h).group(1))
                self.phycassert(subset_index == 1, 'Expecting edge length hyperparameter to be defined only for subset 1, but found %s in RefDistImpl._processParamSample' % h)
                param_refdist_definition += 'edgelen_hyper = %s\n' % self._fitInverseGammaRefDist(params[h])
            elif '_external_hyper' in h:
                subset_index = int(re.match('(\d+)_external_hyper',h).group(1))
                self.phycassert(subset_index == 1, 'Expecting external edge length hyperparameter to be defined only for subset 1, but found %s in RefDistImpl._processParamSample' % h)
                param_refdist_definition += 'external_hyper = %s\n' % self._fitInverseGammaRefDist(params[h])
            elif '_internal_hyper' in h:
                subset_index = int(re.match('(\d+)_internal_hyper',h).group(1))
                self.phycassert(subset_index == 1, 'Expecting internal edge length hyperparameter to be defined only for subset 1, but found %s in RefDistImpl._processParamSample' % h)
                param_refdist_definition += 'internal_hyper = %s\n' % self._fitInverseGammaRefDist(params[h])
            elif '_subset_rate' in h:
                subset_index = int(re.match('(\d+)_subset_rate',h).group(1))
                subset_rates[subset_index] = params[h]
            elif 'Gen' in h:
                pass
            elif 'TL' in h:
                pass
            elif 'lnL' in h:
                pass
            elif 'lnPrior' in h:
                pass
            elif 'brlen' in h:
                pass
            else:
                self.phycassert(False, "RefDistImpl did not account for parameter named '%s'" % h)

        num_subsets = len(subset_rates.keys())
        if num_subsets > 0:
            sorted_keys = subset_rates.keys()[:]
            sorted_keys.sort()
            n = len(subset_rates[2])
            x = []
            for i in range(n):
                t = [subset_rates[k][i] for k in sorted_keys]
                x.append(tuple(t))
            param_refdist_definition += 'subset_relrates = %s\n' % self._fitRelRateRefDist(x)

        #print 'param_refdist_definition:\n',param_refdist_definition
        return param_refdist_definition

    def _recordTreeInMaps(self, tree, split_map, tree_key):
        # Each split found in any tree is associated with a list by the dictionary split_map
        #   The list is organized as follows:
        #   - element 0 is the number of times the split was seen over all sampled trees
        #       (this provides the posterior probability of this split when divided by the
        #       number of trees sampled)
        #   - element 1 holds the sum of edge lengths for this split (this provides the posterior
        #       mean edge length corresponding to this split when divided by the value in element 0)
        #   - element 2 holds the sum of squared edge lengths for this split (this allows the variance
        #       in edge lengths corresponding to this split to be computed)
        #   - elements 3... are, for internal nodes, the indices of trees in which the split was
        #       found (this part is omitted for terminal nodes, which must appear in every tree)
        #       From this tree index list we can compute the following quantities for internal splits:
        #       1. the number of sojourns for each internal split (a sojourn is a series of consecutive
        #          samples in which a split was found that is preceded and followed by at least one
        #          sample not containing the split)
        #       2. sliding window and cumulative plots, such as those produced by AWTY
        # Each distinct tree topology is associated with a similar list by the dictionary tree_map
        #   (but this occurs after this function returns):
        #   - element 0 is again the number of times the tree topology was seen over all samples
        #   - element 1 holds the newick tree description
        #   - element 2 holds the sum of tree lengths over all sampled trees of this topology
        #   - elements 3... are the indices of trees in which the topology was found. This list
        #       allows the number and extent of each sojourn to be computed
        # The return value is a 3-tuple: (tree length, list of internal edge length proportions,
        # list of external edge length proportions).
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'the first preorder node should be the root'
        treelen = 0.0
        internal_edgelens = []
        external_edgelens = []
        has_edge_lens = tree.hasEdgeLens()
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                # Determine whether this split represents an internal or tip node
                is_tip_node = nd.isTip() or nd.getParent().isRoot()
                # Grab the edge length
                edge_len = has_edge_lens and nd.getEdgeLen() or 1.0
                treelen += edge_len
                if is_tip_node:
                    external_edgelens.append(edge_len)
                else:
                    internal_edgelens.append(edge_len)
                self._recordNodeInMaps(nd, split_map, tree_key, is_tip_node, edge_len)

        # Convert stored edge lengths to edge length proportions
        internal_proportions = [edgelen/treelen for edgelen in internal_edgelens]
        external_proportions = [edgelen/treelen for edgelen in external_edgelens]

        return (treelen, internal_proportions, external_proportions)

    def _recordNodeInMaps(self, nd, split_map, tree_key, is_tip_node, edge_len):
        # Grab the split and invert it if necessary to attain a standard polarity
        s = nd.getSplit()
        if (not self.rooted_trees) and s.isBitSet(0):
            s.invertSplit()

        # Create a string representation of the split
        ss = s.createPatternRepresentation()

        # Add string represention of the split to the tree_key list, which
        # will be used to uniquely identify the tree topology
        if not is_tip_node:
            tree_key.append(ss)

        # Update the dictionary entry corresponding to this split
        if ss in split_map.keys():
            # split has been seen before
            entry = split_map[ss]
            entry[0] += 1
            entry[1] += edge_len
            entry[2] += math.pow(edge_len,2.0)
            if not is_tip_node:
                entry.append(self.num_trees_considered)
        else:
            # split not yet seen
            if is_tip_node:
                split_map[ss] = [1, edge_len, math.pow(edge_len,2.0)]
            else:
                split_map[ss] = [1, edge_len, math.pow(edge_len,2.0), self.num_trees_considered]

    def _assignEdgeLensAndSupportValues(self, tree, split_map, total_samples, edge_lengths_are_clade_posteriors = False):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        If edge_lengths_are_clade_posteriors is True, assigns edge lengths in
        supplied tree using posterior probabilities stored in split_map. The
        variable split_map associates a list with each split. The first
        element of this list is the number of trees considered in which the
        split was found, so the edge length applied to the tree is the first
        element of the split divided by total_samples. This produces a tree
        with clade probabilities for edge lengths, which is useful in
        specifying a reference distribution for trees, used in the generalized
        stepping-stone method.

        If edge_lengths_are_clade_posteriors is False (default), assigns edge
        lengths in supplied tree using posterior mean split weights stored in
        split_map. The variable split_map associates a list with each split.
        The first element of this list is the number of trees considered in
        which the split was found, and the second element is the sum of the
        edge lengths of the split over all trees in which it was found. The
        estimated edge length applied to the tree is the second element of
        the split divided by the first element.

        """
        tree.recalcAllSplits(tree.getNObservables())
        nd = tree.getFirstPreorder()
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            s = nd.getSplit()
            if (not self.rooted_trees) and s.isBitSet(0):
                s.invertSplit()
            ss = s.createPatternRepresentation()
            try:
                v = split_map[ss]
            except KeyError:
                raise ValueError('could not find edge length information for the split %s' % ss)

            # determine fudge factor that ensures 0 < support < 1
            n1 = s.countOnBits()
            n2 = s.countOffBits()
            self.phycassert(n1 + n2 == tree.getNObservables(),"total number of taxa differs from sum of taxa on both sides of split")
            if n1 == 1 or n2 == 1:
                support = 1.0
            else:
                lognt  = probdist.lnGamma(2.0*(n1 + n2) - 4.0) - float(n1 + n2 - 3)*math.log(2.0) - probdist.lnGamma(float(n1 + n2) - 2.0)
                logns1 = probdist.lnGamma(2.0*n1 - 2.0) - float(n1 - 2)*math.log(2.0) - probdist.lnGamma(float(n1) - 1.0) - lognt
                logns2 = probdist.lnGamma(2.0*n2 - 2.0) - float(n2 - 2)*math.log(2.0) - probdist.lnGamma(float(n2) - 1.0) - lognt
                ns_over_nt = math.exp(logns1 + logns2)

                # set support value
                numer = float(v[0]) + self.epsilon*ns_over_nt
                denom = float(total_samples) + self.epsilon
                support = numer/denom
            nd.setSupport(support)

            if edge_lengths_are_clade_posteriors:
                nd.setEdgeLen(support)
            else:
                edge_len = float(v[1])/float(v[0])
                nd.setEdgeLen(edge_len)

    def _createTreeLengthRefDistStr(self, tree_lengths, n_int_edges, n_ext_edges, internal_proportions, external_proportions):
        n = len(tree_lengths)
        nint = len(internal_proportions)
        next = len(external_proportions)
        if n <= 1 or nint <= 1 or next <= 1:
            return 'tree_length = <not enough trees sampled to generate>'
        fn = float(n)
        sum_tl  = 0.0
        ss_tl   = 0.0
        for tl in tree_lengths:
            sum_tl  += tl
            ss_tl   += tl*tl
        mean_tl = sum_tl/fn
        var_tl  = (ss_tl - fn*mean_tl*mean_tl)/(fn - 1.0)
        betaT   = mean_tl/var_tl
        alphaT  = betaT*mean_tl

        # Now estimate alpha (Dirichlet parameter for external edge lengths) and
        # c (ratio of internal to external edge lengths - the Dirichlet parameter
        # for internal edge lengths will thus be alpha*c).
        fn = float(nint)
        sum_prop  = 0.0
        ss_prop   = 0.0
        for p in internal_proportions:
            sum_prop += p
            ss_prop += p*p
        mean_internal = sum_prop/fn
        var_internal  = (ss_prop - fn*mean_internal*mean_internal)/(fn - 1.0)

        fn = float(next)
        sum_prop  = 0.0
        ss_prop   = 0.0
        for p in external_proportions:
            sum_prop += p
            ss_prop += p*p
        mean_external = sum_prop/fn
        var_external  = (ss_prop - fn*mean_external*mean_external)/(fn - 1.0)

        alpha, c = self._estimateTreeLengthPriorDirichlet(n_int_edges, mean_internal, var_internal, n_ext_edges, mean_external, var_external)

        return 'tree_length = TreeLengthDist(%g,%g,%g,%g)' % (alphaT, betaT, alpha, c)

    def _processTreeSample(self, treefname):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the tree sample file and returns a
        reference distribution definition string.

        """
        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        print 'input_trees is a(n)',input_trees.__class__.__name__
        if input_trees is None:
            raise NoTreesWarning()

        num_trees = 0
        self.num_trees_considered = 0
        split_map = {}
        tree_map = {}

        # Open sumt_tfile_name and read trees therein
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_tree_defs = list(input_trees)
        self.taxon_labels = input_trees.taxon_labels # this must be kept after the coercion of the trees to a list (in case that is what triggers the readinf of the file with the taxon labels)

        n_ext_edges = len(self.taxon_labels)
        n_int_edges = 2*n_ext_edges - 3

        num_stored_trees = len(self.stored_tree_defs)
        if num_stored_trees == 0:
            raise NoTreesWarning()

        # Create list to store tree lengths; used to create a tree length reference distribution
        tree_lengths = []
        internal_proportions = []
        external_proportions = []

        # Build each tree and add the splits and tree topologies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        self.stdout.info('Compiling lists of tree topologies and splits...')
        t = phylogeny.Tree()
        if self.rooted_trees:
            t.setRooted()
            n_int_edges += 1
        for tree_def in self.stored_tree_defs:
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
            #ntips = t.getNTips()
            ntips = t.getNObservables()
            t.recalcAllSplits(ntips)

            treelen, internal_props, external_props = self._recordTreeInMaps(t, split_map, tree_key)
            tree_lengths.append(treelen)
            internal_proportions += internal_props
            external_proportions += external_props

            # Update tree_map, which is a map with keys equal to lists of internal node splits
            # and values equal to lists containing (0) the frequency, (1) tree length, and
            # (2) total number of trees considered. From (0) and (2) the
            # estimated posterior probability of the tree can be calculated.
            tree_key.sort()
            k = tuple(tree_key)
            if k in tree_map.keys():
                # tree topology has been seen before
                entry = tree_map[k]
                entry[0] += 1                           # frequency
                entry[2] += treelen                     # sum of edge lengths
                entry.append(self.num_trees_considered) # total trees sampled
            else:
                # tree topology has not yet been seen
                tree_map[k] = [1, tree_def, treelen, self.num_trees_considered]

        tree_length_refdist = self._createTreeLengthRefDistStr(tree_lengths, n_int_edges, n_ext_edges, internal_proportions, external_proportions)

        #self.stdout.info('\nSummary of sampled trees:')
        #self.stdout.info('-------------------------')
        #self.stdout.info('Tree source: %s' %  str(input_trees))
        #self.stdout.info('Total number of trees in file = %d' % num_trees)
        #self.stdout.info('Number of trees considered = %d' % self.num_trees_considered)
        #self.stdout.info('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        #self.stdout.info('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNObservables()))
        if self.num_trees_considered == 0:
            raise NoTreesWarning()
        # Sort the splits from highest posterior probabilty to lowest
        split_vect = split_map.items()
        c = lambda x,y: cmp(y[1][0], x[1][0])
        split_vect.sort(c)

        num_trivial = 0
        first_below_50 = None
        for i,(k,v) in enumerate(split_vect):
            # len(v) is 2 in trivial splits because these splits are associated with tips,
            # for which the sojourn history is omitted (v[0] is frequency and v[1] is edge length sum)
            trivial_split = len(v) == 3 and True or False
            if trivial_split:
                num_trivial += 1

            # Split frequency is simply the first element of the list
            split_freq = v[0]

            # Split posterior is split_freq divided by the number of trees considered
            split_posterior = float(split_freq)/float(self.num_trees_considered)

            # Identify the index of the first split having a posterior probability
            # lower than 0.5. The list up to this point will be used to generate the
            # majority-rule consensus tree
            if not first_below_50 and split_posterior < 0.5:
                first_below_50 = i

        ## No longer used because of the requirement that focal tree be fully resolved
        ## Build 50% majority rule tree if requested
        #majrule = phylogeny.Tree()
        #if self.rooted_trees:
        #    majrule.setRooted()
        #tm = phylogeny.TreeManip(majrule)
        #majrule_splits = []
        #for k,v in split_vect[:first_below_50]:
        #    if len(v) > 3:
        #        majrule_splits.append(k)
        #
        #if len(majrule_splits) == 0:
        #    tm.starTree(num_trivial)
        #else:
        #    tm.buildTreeFromSplitVector(majrule_splits, probdist.Exponential(10))
        ##self._assignEdgeLensAndSupportValues(majrule, split_map, self.num_trees_considered)
        ##summary_short_name_list = ['majrule']
        ##summary_full_name_list = ['Majority-rule Consensus']
        ##summary_tree_list = [majrule]

        # Find MAP tree by searching through tree_map for the stored tree with highest sample frequency
        maxfrq = 0.0
        for k in tree_map.keys():
            entry = tree_map[k]
            frq = entry[0]
            if frq > maxfrq:
                maxfrq = frq
                treekey = list(k)
        reftree = phylogeny.Tree()
        tm = phylogeny.TreeManip(reftree)
        tm.buildTreeFromSplitVector(treekey, probdist.Exponential(10))

        tree_refdist_definition = ''

        # Output reference prior information
        self._assignEdgeLensAndSupportValues(reftree, split_map, self.num_trees_considered, True) # True means assign clade posteriors as edge lengths
        tree_refdist_definition += '%s\n' % reftree.makeNewickForRefDist(5)
        for k,v in split_vect[:first_below_50]:
            self.phycassert(int(v[0]) > 1, 'Not able to estimate variance of edge length reference distribution because sample size is less than 2 for split %s' % k)
            num_edgelens = float(v[0])
            mean_edgelen = float(v[1])/num_edgelens
            var_edgelen = (float(v[2]) - num_edgelens*math.pow(mean_edgelen, 2.0))/(num_edgelens - 1.0)
            gamma_b = var_edgelen/mean_edgelen
            gamma_a = mean_edgelen/gamma_b
            tree_refdist_definition += 'split_%s = Gamma(%g,%g)\n' % (k,gamma_a,gamma_b)

        # create reference distribution for edges that did not make it into the majority rule tree
        num_NA_edgelens = 0
        sum_NA_edgelens = 0.0
        sum_squared_NA_edgelens = 0.0
        for k,v in split_vect[first_below_50:]:
            num_NA_edgelens += float(v[0])
            sum_NA_edgelens += float(v[1])
            sum_squared_NA_edgelens += float(v[2])
        if num_NA_edgelens > 1.0:
            mean_NA_edgelen = sum_NA_edgelens/float(num_NA_edgelens)
            var_NA_edgelen = (float(sum_squared_NA_edgelens) - num_NA_edgelens*math.pow(mean_NA_edgelen, 2.0))/(num_NA_edgelens - 1.0)
            gamma_b = var_NA_edgelen/mean_NA_edgelen
            gamma_a = mean_NA_edgelen/gamma_b
            tree_refdist_definition += 'split_NA = Gamma(%g,%g)' % (gamma_a,gamma_b)
        else:
            # if only one edge length was found, use Exponential distribution with mean equal to this edge length
            rate_param = 1.0/sum_NA_edgelens
            tree_refdist_definition += 'split_NA = Exponential(%g)' % rate_param

        return tree_refdist_definition + '\n' + tree_length_refdist

    def run(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Reads the contents of the param and tree file and produces a refdist
        file.

        """
        self.skip = self.opts.skip
        self.epsilon = self.opts.epsilon
        self.rooted = self.opts.rooted
        self.out_refdist_file_name = self.opts.out.refdistfile._getFilename()

        param_refdist_definition = self._processParamSample(self.opts.params)

        try:
            tree_refdist_definition = self._processTreeSample(self.opts.trees)
        except NoTreesWarning:
            tree_refdist_definition = ''

        self._openRefDistFile()
        self.refdistf.write('%s\n' % tree_refdist_definition)
        self.refdistf.write('%s\n' % param_refdist_definition)
        self._closeRefDistFile()

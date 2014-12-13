import os,sys,math,random
import os,sys,math,random
#from phycas.treeviewer import *
from phycas import *
from phycas import useWxPhycas
from phycas.utilities.PhycasCommand import *
from phycas.utilities.CommonFunctions import CommonFunctions
from phycas.utilities.PDFTree import PDFTree

class CCDTree(object):
    def __init__(self):
        self.nodes = {} # key is split, value is CCDNode containing a map associating ctuple keys with count values
        self.rootsplit = None
        self.H = None

    def lfact(self, n):
        """
        Returns log(n!)
        """
        if n == 0 or n == 1:
            return 0.0
        return math.lgamma(n+1)

    def lognrooted(self, y):
        fy = float(y)
        log_num_rooted_trees = 0.0
        if y > 2:
            log_num_rooted_trees = self.lfact(2.*fy-3) - (fy-2)*math.log(2.0) - self.lfact(fy-2)
        return log_num_rooted_trees

    def countbits(self, s):
        return sum([c == '*' and 1 or 0 for c in s])

    def showCCDTree(self, output):
        ntotal = 0
        if self.rootsplit is None:
            output.info('root not found')
        else:
            output.info('root split: %s' % self.rootsplit)
            v = self.nodes[self.rootsplit]
            ntotal = sum(v.children.values())

        assert ntotal > 0, 'root has no children in CCDTree.showCCDTree function'

        totalIw = 0.0

        for s in self.nodes.keys():
            v = self.nodes[s]
            n = sum(v.children.values())
            p = float(n)/ntotal
            if s == self.rootsplit:
                output.info('\n%s (root) prob = %g' % (s,p))
            else:
                output.info('\n%s prob = %g' % (s,p))

            h = 0.0
            nset = self.countbits(s)
            hp = self.lognrooted(nset)
            for ctuple in v.children.keys():
                nchild = v.children[ctuple]
                pchild = float(nchild)/ntotal
                pchild /= p
                h -= pchild*math.log(pchild)
                nsetleft  = self.countbits(ctuple[0])
                nsetright = self.countbits(ctuple[1])
                hp -= pchild*(self.lognrooted(nsetleft) + self.lognrooted(nsetright))
                output.info('%6d %.5f <-- %s' % (nchild, pchild, '|'.join(ctuple)))
            I = hp - h
            w = p
            Iw = I*w
            totalIw += Iw
            output.info('h   = %g' % h)
            output.info('hp  = %g' % hp)
            output.info('I   = %g' % I)
            output.info('w   = %g' % w)
            output.info('w*I = %g' % Iw)
        output.info('total I = %g' % totalIw)

    def calcI(self, output = None):
        assert self.rootsplit is not None, 'cannot compute H because basal split could not be identified'

        v = self.nodes[self.rootsplit]
        ntotal = sum(v.children.values())
        totalH  = 0.0
        totalHp = 0.0
        totalI  = 0.0
        results = []
        for s in self.nodes.keys():
            v = self.nodes[s]
            n = sum(v.children.values())
            p = float(n)/ntotal
            H = 0.0
            nset = self.countbits(s)
            Hp = self.lognrooted(nset)
            for ctuple in v.children.keys():
                nchild = v.children[ctuple]
                pchild = float(nchild)/ntotal
                pchild /= p
                H -= pchild*math.log(pchild)
                nsetleft  = self.countbits(ctuple[0])
                nsetright = self.countbits(ctuple[1])
                Hp -= pchild*(self.lognrooted(nsetleft) + self.lognrooted(nsetright))
            I = Hp - H
            w = p
            results.append((H, Hp, I*w, w, s))
            totalH  += H*w
            totalHp += Hp*w
            totalI  += I*w

        if output is not None:
            output.info('\nInformation partitioned by clade:')
            output.info('  %12s %12s %12s %12s %12s   %s' % ('H', 'Hp', 'I', '% max. I', 'prob', 'clade'))
            totalPct = 0.0
            for H, Hp, I, prob, clade in results:
                if I > 0.0:
                    pct = 100.0*I/totalI
                    totalPct += pct
                    output.info('  %12.5f %12.5f %12.5f %12.1f %12.5f   %s' % (H, Hp, I, pct, prob, clade))
            output.info('  %12.5f %12.5f %12.5f %12s' % (totalH, totalHp, totalI, totalPct))

        return (totalH, totalHp, totalI)

    #def _calcH(self, h0, p0, c0):
    #    H = h0
    #    for prob,ctuple in zip(p0, c0):
    #        for s in ctuple:
    #            if s in self.nodes.keys():
    #                h, p, c = self.nodes[s].calcH()
    #                H += prob*self._calcH(h, p, c)
    #    return H
    #
    #def calcH(self):
    #    assert self.rootsplit is not None, 'cannot compute H because basal split could not be identified'
    #    h, p, ctuples = self.nodes[self.rootsplit].calcH()
    #    self.H = self._calcH(h, p, ctuples)
    #    return self.H

    def updateCCDTree(self, s, ctuple, isroot):
        if isroot:
            self.rootsplit = s
        if s in self.nodes.keys():
            v = self.nodes[s]
            v.increment(ctuple)
        else:
            self.nodes[s] = CCDNode(ctuple)

class CCDNode(object):
    def __init__(self, ctuple):
        self.children = {ctuple:1}
        self.total = None
        self.h = None
        self.hp = None
        self.p = None
        self.I = None

    def increment(self, ctuple):
        if ctuple in self.children.keys():
            self.children[ctuple] += 1
        else:
            self.children[ctuple] = 1

    #def calcH(self):
    #    self.total = 0.0
    #    sum_terms = 0.0
    #    freqs = []
    #    for k in self.children.keys():
    #        n = float(self.children[k])
    #        self.total += n
    #        freqs.append(n)
    #        sum_terms += n*math.log(n)
    #    sum_terms -= self.total*math.log(self.total)
    #    self.h = -sum_terms/self.total
    #    self.p = [(f/self.total) for f in freqs]
    #    return (self.h, self.p, self.children.keys())

class TreeSummarizer(CommonFunctions):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Saves consensus tree and QQ plots in pdf files.

    """
    def __init__(self, opts):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes TreeSummarizer object.

        """
        CommonFunctions.__init__(self, opts)

        self.pdf_page_width = 8.5       # should be an option
        self.pdf_page_height = 11.0     # should be an option
        self.pdf_ladderize = 'right'    # should be an option - valid values are 'right', 'left' or None
        self.rooted_trees = opts.rooted
        self.outgroup = opts.outgroup_taxon
        if self.rooted_trees and opts.outgroup_taxon:
            self.outgroup = None
            self.warning('Specifying True for sumt_rooted is incompatible with specifying\nsumt_outgroup_taxon; I will pretend that you set sumt_outgroup_taxon to None')
        self.num_trees_considered = 0

    def assignEdgeLensAndSupportValues(self, tree, split_map, joint_split_map, total_samples, estimate_ccb_probs = True):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Assigns edge lengths in supplied tree using posterior mean split
        weights stored in split_vect. The variable split_vect associates a
        list with each split. The first element of this list is the number of
        trees considered in which the split was found, and the second element
        is the sum of the edge lengths of the split over all trees in which it
        was found. The estimated edge length applied to the tree is the
        second element of the split divided by the first element.

        Also estimates joint posterior probability using Larget's 2013 method
        (see Larget, B. 2013. The estimation of tree posterior probabilities
        using conditional clade probability distributions. Systematic Biology 62:501-511).

        """
        #tree.recalcAllSplits(tree.getNTips())
        tree.recalcAllSplits(tree.getNObservables())
        nd = tree.getFirstPreorder()
        eq_brlen = self.optsout.trees.equal_brlens
        larget_log_prob = 0.0
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
            support = float(v[0])/float(total_samples)
            nd.setSupport(support)

            if eq_brlen:
                edge_len = 1.0
            else:
                edge_len = float(v[1])/float(v[0])
            nd.setEdgeLen(edge_len)

            if estimate_ccb_probs:
                #is_tip_node = nd.isTip() or nd.getParent().isRoot()
                if not nd.isTip():
                    # Grab splits of child nodes and make tuple comprising
                    #NEWWAY s plus all of the child splits
                    #OLDWAY s plus all but one of the child splits
                    child = nd.getLeftChild()
                    slist = []
                    while child:
                        c = child.getSplit()
                        cc = c.createPatternRepresentation()
                        slist.append(cc)
                        child = child.getRightSib()
                    slist.sort()
                    #OLDWAY stuple = tuple([ss] + slist[:-1])
                    stuple = tuple([ss] + slist[:]) #NEWWAY

                    # Get count from joint split map
                    try:
                        joint_count = joint_split_map[stuple]
                    except KeyError:
                        raise ValueError('could not find conditional clade count for this key [%s]' % ','.join(stuple))
                    joint_prob = float(joint_count)/float(total_samples)
                    assert joint_prob > 0.0, 'joint_prob less than or equal to zero in TreeSummarizer.assignEdgeLensAndSupportValues (SumTImpl.py)'
                    assert support > 0.0, 'support less than or equal to zero in TreeSummarizer.assignEdgeLensAndSupportValues (SumTImpl.py)'
                    cond_log_prob = math.log(joint_prob) - math.log(support)
                    larget_log_prob += cond_log_prob

        return larget_log_prob

    def save_trees(self, short_names, full_names, trees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Saves the supplied tree names (tree_names) and newick descriptions
        (trees) to a nexus tree file named tree_file. If a file by that name
        exists, a new name will be invented by adding a random integer
        between 1 and 1 million to the end of the supplied tree file name.

        """
        # we use the trees file spec to make a pdf with the same filename prefix
        trees_spec = self.optsout.trees
        m = trees_spec.mode
        if isinstance(m, SkipExistingFileBehavior):
            return

        p = trees_spec.prefix
        if not p:
            p = trees_spec.filename
        print(str(trees_spec))

        if self.opts.save_trees_pdf:
            pdf_spec = PDFOutputSpec(p, "", "")
            try:
                pdf_spec.mode = trees_spec.mode
            except:
                pdf_spec.mode = ADD_NUMBER
            print(str(pdf_spec))

        pdf, treef = None, None
        try:
            # Open the output pdf file
            if self.opts.save_trees_pdf:
                pdf = pdf_spec.open(self.pdf_page_width, self.pdf_page_height, self.stdout)

            treef = trees_spec.open(self.taxon_labels, self.stdout)

            if (pdf is None) and (treef is None):
                return
            reroot = self.outgroup
            root_num = 0
            if reroot:
                try:
                    root_num = list(self.taxon_labels).index(self.outgroup)
                except ValueError:
                    self.stdout.info('Warning: could not root tree using specified outgroup: no tip having name "%s" could be found' % self.outgroup)
                    self.stdout.info('Here is the list of all taxon labels:')
                    for label in self.taxon_labels:
                        self.stdout.info('  %s' % label)
                    reroot = False
            for short,full,tree in zip(short_names, full_names, trees):
                # Reroot (if requested) tree and output newick description to tree file
                if reroot:
                    tree.rerootAtTip(root_num)
                if treef:
                    treef.writeTree(tree, short)
                if pdf is not None:
                    # Add taxon names to tip nodes and output to pdf file
                    tree.rectifyNames(self.taxon_labels)
                    if self.pdf_ladderize:
                        if self.pdf_ladderize == 'right':
                            tree.ladderizeRight()
                        else:
                            tree.ladderizeLeft()
                    pdf.newPage()
                    PDFTree().tree2pdf(pdf, tree, full, show_support = True)
        finally:
            if pdf is not None:
                pdf_spec.close()
            if treef is not None:
                trees_spec.close()

    def awtyPlot(self, pdf_factory, split_vect, ntrees, ndivisions = 10):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The supplied split_vect is a vector of key,value pairs, where the key
        is a split pattern string and v=value[2:] is a list of samples in
        which the split was found. The minimum possible value in v is 1, and
        the maximum possible value is ntrees. Consider the following
        oversimplified example involving a split with posterior probability
        0.6 (v has length 6 and ntrees = 10): v = [3,4,5,7,9,10]. Suppose we
        wish to plot the cumulative posterior probability of this split at
        five, evenly-spaced points. Dividing ntrees by 5 gives 2, and these
        are thus the values to be plotted:
        x     y
        --------   x = range(0, ntrees + 1, ntrees//ndivisions)
        0    0.0     = [0, 2, 4, 6, 8, 10]
        2    0.0
        4    0.5   To calculate y, let k track the number of values in v
        6    0.5   that are less than x
        8    0.5
        10   0.6
        --------
        The vector x can be calculated using the range function; however,
        computing y requires a loop:
        i  x[i]  k  v[k]  y = k/x[i] except first
        ----------------------------
        0    0   0   3           0.0 <- this one always 0.0
        1    2   0   3    0/2  = 0.0
        2    4   2   5    2/4  = 0.5
        3    6   3   7    3/6  = 0.5
        4    8   4   9    4/8  = 0.5
        5   10   6   -    6/10 = 0.6
        ----------------------------

        """
        # candidates for phycas variable status
        splits_to_plot = 1000
        ignore_uninteresting = True

        splits_plotted = 0
        trivial_ignored = 0
        uninteresting_ignored = 0
        data = []
        stride = ntrees // ndivisions
        if stride < 2:
            self.stdout.error("Insufficient number of trees (%d) to construct a splits through time plot with %d divisions" % (ntrees, ndivisions))
            return
        xvect = range(0, ntrees + 1, ntrees//ndivisions)
        for k,value in split_vect:
            if len(value) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and value[0] == ntrees:
                uninteresting_ignored += 1
                continue
            line_data = [(0.0,0.0)]
            v = value[2:]
            k = 0
            for i,x in enumerate(xvect[1:]):
                while k < len(v) and v[k] <= x:
                    k += 1
                y = float(k)/float(x)
                line_data.append((x,y))
            data.append(line_data)

            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        if splits_plotted == 0:
            self.stdout.info('No AWTY plot created because no splits with posteriors less than 1.0 were found')
            return False
        pdf = pdf_factory()
        self.stdout.info('\nSaving AWTY plot in file %s...' % self._splitsPdfFilename)
        pdf.scatterPlot(data, title = 'Split Probabilities Through Time', xinfo = (0,ntrees,10,0), yinfo = (0.0,1.0,10,1))

        if self.opts.useGUI and useWxPhycas():
            try:
                import wx
            except:
                pass
            else:
                from phycas import wxPhycas
                wxapp = wxPhycas.CreateAWTYApp()
                wxapp.scatterPlot(data, title = 'Split Probabilities Through Time', xinfo = (0,ntrees,10,0), yinfo = (0.0,1.0,10,1))
                wxapp.MainLoop()

        if trivial_ignored + uninteresting_ignored > 0:
            self.stdout.info('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))
        return True

    def sojournPlot(self, pdf_factory, split_vect, ntrees):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a plot comprising numerous horizontal lines, each representing
        the trajectory of a non-trivial, "interesting" split throughout the
        MCMC simulation. Trivial splits are those separating a terminal taxon
        from the other taxa, and are aways present. "Interesting" splits are
        those that were absent during some of the non-burn-in period. This
        plot provides a visual picture of mixing in the Markov chain with
        respect to splits.

        """
        # candidates for phycas variable status
        splits_to_plot = 1000
        ignore_uninteresting = True

        trivial_ignored = 0
        uninteresting_ignored = 0

        # First, save the sojourn vectors of all the splits to be plotted
        v = []
        splits_plotted = 0
        for k,value in split_vect:
            if len(value) == 2:
                trivial_ignored += 1
                continue
            if ignore_uninteresting and value[0] == ntrees:
                uninteresting_ignored += 1
                continue
            v.append(value[2:])
            splits_plotted += 1
            if splits_plotted == splits_to_plot:
                break

        if splits_plotted == 0:
            self.stdout.info('No sojourn plot created because all non-trivial splits had posterior probability 1.0')
            return False

        self.stdout.info('\nSaving Sojourn plot in file %s...' % self._splitsPdfFilename)

        # Now create the plot data
        data = []
        for i,s in enumerate(v):
            y = float(i + 1)
            prev_x = s[0]
            line_data = [(prev_x,y)]
            for x in s[1:]:
                if x - prev_x > 1:
                    # a gap has appeared, start a new line and close old one
                    data.append(line_data)
                    line_data = [(x,y)]
                else:
                    line_data.append((x,y))
                prev_x = x
            data.append(line_data)

        # Determine y-axis info (ymax and ydivs) such that ydivs is as close as possible to
        # 10 but divides evenly into ymax - ymin so that there are no fractional y-axis labels
        splits_per_tick = splits_plotted//10
        if splits_per_tick < 2:
            splits_per_tick = 2
        min_ticks = math.ceil(float(splits_plotted)/float(splits_per_tick))
        ymax = float(splits_per_tick)*min_ticks
        ydivs = int(ymax)//splits_per_tick
        if splits_plotted <= 10:
            ydivs = int(ymax)

        pdf = pdf_factory()
        pdf.scatterPlot(data, points = False, line_width = 3, line_cap_style = 'square', title = 'Split Sojourns', xinfo = (0,ntrees,10,0), yinfo = (0.0,float(ymax),ydivs,0))
        if trivial_ignored + uninteresting_ignored > 0:
            self.stdout.info('%d trivial and %d uninteresting splits were ignored.' % (trivial_ignored, uninteresting_ignored))

        return True

    def recordTreeInMaps(self, tree, split_lookup, split_map, joint_split_map, ccd_tree, tree_key):
        # Each split found in any tree is associated with a list by the dictionary split_map
        #   The list is organized as follows:
        #   - element 0 is the number of times the split was seen over all sampled trees
        #       (this provides the posterior probability of this split when divided by the
        #       number of trees sampled)
        #   - element 1 holds the sum of edge lengths for this split (this provides the posterior
        #       mean edge length corresponding to this split when divided by the value in element 0)
        #   - elements 2... are, for internal nodes, the indices of trees in which the split was
        #       found (this part is omitted for terminal nodes, which must appear in every tree)
        #   From the tree index list we can compute the following quantities for internal splits:
        #   1. the number of sojourns for each internal split (a sojourn is a series of consecutive
        #      samples in which a split was found that is preceded and followed by at least one
        #      sample not containing the split)
        #   2. sliding window and cumulative plots, such as those produced by AWTY
        #
        # Each split and its child splits provide the keys to joint_split_map, for which values equal
        # the number of times that joint combination of splits appears in any tree. This allows
        # calculation of the conditional clade probabilities used to estimate posterior probabilities
        # of tree topologies (see Larget, B. 2013. The estimation of tree posterior probabilities
        # using conditional clade probability distributions. Systematic Biology 62:501-511).
        #
        # Each distinct tree topology is associated with a similar list:
        #   - element 0 is again the number of times the tree topology was seen over all samples
        #   - element 1 holds the newick tree description
        #   - element 2 holds the sum of tree lengths over all sampled trees of this topology
        #   - elements 3... are the indices of trees in which the topology was found. This list
        #       allows the number and extent of each sojourn to be computed
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'the first preorder node should be the root'
        #print '>>> %s %s' % (nd.getSplit().createPatternRepresentation(),nd.isRoot() and "(root)" or "")
        #split_vec = []
        treelen = 0.0
        has_edge_lens = tree.hasEdgeLens()

        #print
        #print '********* new tree **********'
        #print tree.makeNumberedNewick()

        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                # Determine whether this split represents an internal or tip node
                #is_tip_node = nd.isTip() or nd.getParent().isRoot()
                # Grab the edge length
                edge_len = has_edge_lens and nd.getEdgeLen() or 1.0
                treelen += edge_len
                self.recordNodeInMaps(nd, split_lookup, split_map, joint_split_map, ccd_tree, tree_key, edge_len)
        return treelen

    def recordNodeInMaps(self, nd, split_lookup, split_map, joint_split_map, ccd_tree, tree_key, edge_len):
        # Grab the split and invert it if necessary to attain a standard polarity
        s = nd.getSplit()
        if (not self.rooted_trees) and s.isBitSet(0):
            s.invertSplit()

        # Create a string representation of the split
        ss = s.createPatternRepresentation()

        # Add string represention of the split to the tree_key list, which
        # will be used to uniquely identify the tree topology
        is_pendant_edge = nd.isTip() or nd.getParent().isRoot()
        if not is_pendant_edge:
            tree_key.append(ss)

        # Update the split_map entry corresponding to this split
        if ss in split_map.keys():
            # split has been seen before
            entry = split_map[ss]
            entry[0] += 1
            entry[1] += edge_len
            if not is_pendant_edge:
                entry.append(self.num_trees_considered)
        else:
            # split not yet seen
            split_lookup[ss] = s
            if is_pendant_edge:
                split_map[ss] = [1, edge_len]
            else:
                split_map[ss] = [1, edge_len, self.num_trees_considered]

        if not nd.isTip():
            # Store split in joint_split_map, or update count if already in map
            stuple = (ss,)
            if stuple in joint_split_map.keys():
                joint_split_map[stuple] += 1
            else:
                joint_split_map[stuple] = 1

            # Grab splits of child nodes and make tuple comprising split and splits of all children
            child = nd.getLeftChild()
            slist = []
            while child:
                c = child.getSplit()
                cc = c.createPatternRepresentation()
                slist.append(cc)
                child = child.getRightSib()
            slist.sort()
            stuple = tuple([ss] + slist[:])

            # Store split in joint_split_map, or update count if already in map
            if stuple in joint_split_map.keys():
                joint_split_map[stuple] += 1
            else:
                joint_split_map[stuple] = 1

            #print '~~~ %s %s' % (ss,nd.getParent().isRoot() and "(parent is root)" or "")
            #raw_input('...')
            ccd_tree.updateCCDTree(ss, tuple(slist), nd.getParent().isRoot())

    def findBestParentSplit(self, curr_stuple, joint_split_map):
        best_parent_key = None
        best_parent_value = 0.0
        for stuple in joint_split_map.keys():
            for s in stuple[1:]:
                if s == curr_stuple[0]:
                    if best_parent_key is None or joint_split_map[stuple] > best_parent_value:
                        best_parent_key = stuple
                        best_parent_value = joint_split_map[stuple]
        if best_parent_key is not None:
            return best_parent_key
        else:
            return None

    def findBestChildSplits(self, curr_stuple, joint_split_map):
        best_children = []
        for bs in curr_stuple[1:]:
            best_child_key = None
            best_child_value = 0.0
            for stuple in joint_split_map.keys():
                if len(stuple) > 1 and bs == stuple[0]:
                    if best_child_key is None or joint_split_map[stuple] > best_child_value:
                        best_child_key = stuple
                        best_child_value = joint_split_map[stuple]
            if best_child_key is not None:
                best_children.append(best_child_key)
        return best_children

    def findBestRelatedSplits(self, best_stuple, joint_split_map):
        # start by finding root by working down the tree from best_stuple until no more parents exist
        curr_stuple = best_stuple
        parent = None
        while curr_stuple is not None:
            parent = curr_stuple
            curr_stuple = self.findBestParentSplit(curr_stuple, joint_split_map)

        # now work up the tree finding the best children for each split
        best_descendants = []
        stuple_list = [parent]
        while len(stuple_list) > 0:
            new_stuple_list = []
            for curr_stuple in stuple_list:
                best_child_stuples = self.findBestChildSplits(curr_stuple, joint_split_map)
                new_stuple_list.extend(best_child_stuples)
            best_descendants.extend(new_stuple_list)
            stuple_list = new_stuple_list

        numer_term = float(joint_split_map[parent])/float(self.num_trees_considered)
        numer = numer_term
        denom = 1.0
        #numerator_terms = [numer_term]
        #denominator_terms = []
        for d in best_descendants:
            numer_term = float(joint_split_map[d])/float(self.num_trees_considered)
            #numerator_terms.append(numer_term)
            numer *= numer_term
            denom_term = float(joint_split_map[(d[0],)])/float(self.num_trees_considered)
            #denominator_terms.append(denom_term)
            denom *= denom_term

        #print 'numer =',numer
        #print 'numerator_terms =',numerator_terms
        #print 'denom =',denom
        #print 'denominator_terms =',denominator_terms
        #print 'numer/denom =',(numer/denom)

        return (numer/denom, [parent] + best_descendants)

    def calcKLupper(self, joint_split_map, log_num_possible_topol):
        # First, find largest joint split frequency in joint_split_map
        max_joint_split_freq = 0.0
        best_stuple = None
        for stuple in joint_split_map.keys():
            joint_split_freq = joint_split_map[stuple]
            if len(stuple) > 1 and (joint_split_freq > max_joint_split_freq):
                max_joint_split_freq = joint_split_freq
                best_stuple = stuple

        # Second, find set of compatible splits yielding highest overall posterior using Larget's method
        best_posterior, best_relateds = self.findBestRelatedSplits(best_stuple, joint_split_map)

        # Finally, to find upper bound for KL, assume 1/best_posterior tree topologies have
        # posterior = best_posterior, and all remaining topologies have posterior = 0.0
        # If p = best_posterior and n = total possible number of distinct tree topologies,
        #  KL (upper) = (1/p) p [log(p) - log(1/n)] + (n - 1/p)(0) [log(0) - log(1/n)]
        #             = log(p n)
        #             = log(p) + log(n)
        log_best_posterior = math.log(best_posterior)
        assert log_best_posterior > -log_num_possible_topol, 'log_best_posterior () <= -log_num_possible_topol (), which should not be possible' % (log_best_posterior,-log_num_possible_topol)
        return log_best_posterior + log_num_possible_topol

    def debugShowCCDTree(self, ccd_tree):
        self.stdout.info('\nCCD Tree:')
        ccd_tree.showCCDTree(self.stdout)
        raw_input('...')

    def debugShowJointSplitFreqs(self, joint_split_map, log_num_possible_topol):
        max_joint_split_freq = 0.0
        best_stuple = None
        print '\nList of joint splits and their frequencies:'
        for stuple in joint_split_map.keys():
            srep = '('
            srep += stuple[0]
            for s in stuple[1:]:
                srep += ','
                srep += s
            srep += ')'
            print '  ',joint_split_map[stuple],'<--',srep
            joint_split_freq = float(joint_split_map[stuple])/float(self.num_trees_considered)
            if len(stuple) > 1 and (joint_split_freq > max_joint_split_freq):
                max_joint_split_freq = joint_split_freq
                best_stuple = stuple
        print
        print 'max_joint_split_freq =', max_joint_split_freq
        print 'best_stuple =', best_stuple

        best_posterior, best_relateds = self.findBestRelatedSplits(best_stuple, joint_split_map)
        print 'best_posterior =',best_posterior
        print 'best related splits:'
        for k in best_relateds:
            print '  ',joint_split_map[k],'<--',k

        log_best_posterior = math.log(best_posterior)
        assert log_best_posterior > -log_num_possible_topol, 'log_best_posterior () <= -log_num_possible_topol (), which should not be possible' % (log_best_posterior,-log_num_possible_topol)
        print 'KL upper bound =', (log_best_posterior + log_num_possible_topol)

    def consensus(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This is the main member function of the class. It computes and outputs
        a summary of splits found in the sampled trees, a summary of distinct
        tree topologies found, and a majority rule consensus tree, which is
        saved to the file sumt_output_tree_file if this variable is not None.

        """
        # Check to make sure user specified an input tree file
        input_trees = self.opts.trees
        self.stdout.phycassert(input_trees, 'trees cannot be None or empty when sumt method is called')

        num_trees = 0
        self.num_trees_considered = 0
        split_lookup = {}   # keys are string representations of splits, values are split objects
        split_map = {}
        joint_split_map = {}
        ccd_tree = CCDTree()

        # key is list of splits, value is tuple(count, newick, treelen, 1st time seen, 2nd time seen, ...)
        tree_map = {}

        # Open sumt_tfile_name and read trees therein
        self.stdout.info('\nReading %s...' % str(input_trees))
        self.stored_tree_defs = list(input_trees)
        #print 'input_trees.__class__.__name__ =',input_trees.__class__.__name__
        #raw_input('...stored_tree_defs...')
        self.taxon_labels = input_trees.taxon_labels # this must be kept after the coercion of the trees to a list (in case that is what triggers the readinf of the file with the taxon labels)

        num_stored_trees = len(self.stored_tree_defs)
        self.stdout.phycassert(num_stored_trees > 0, 'Specified tree source (%s) contained no stored trees' %  str(input_trees))

        # Build each tree and add the splits and tree topolgies found there to the
        # dictionary of splits (split_map) and the dictionary of tree topologies
        # (tree_map), respectively
        self.stdout.info('Compiling lists of tree topologies and splits...')
        t = phylogeny.Tree()
        if self.rooted_trees:
            t.setRooted()

        # values used for display purposes
        split_field_width = 0
        sojourn_field_width = 2 + math.floor(math.log10(float(num_stored_trees)))

        curr_tree = 0
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
            if ntips > split_field_width:
                # this is necessary only if number of taxa varies from tree to tree
                split_field_width = ntips
            t.recalcAllSplits(ntips)

            treelen = self.recordTreeInMaps(t, split_lookup, split_map, joint_split_map, ccd_tree, tree_key)

            # Update tree_map, which is a map with keys equal to lists of internal node splits
            # and values equal to 2-element lists containing the frequency and newick tree
            # description
            tree_key.sort()
            k = tuple(tree_key)
            if k in tree_map.keys():
                # tree topology has been seen before
                entry = tree_map[k]
                entry[0] += 1           # increment count of times this tree topology has been seen
                entry[2] += treelen     # add treelen to sum of tree lengths
                entry.append(self.num_trees_considered)
            else:
                # tree topology has not yet been seen
                tree_map[k] = [1, tree_def, treelen, self.num_trees_considered]

        #raw_input('...stop now...')

        self.stdout.info('\nSummary of sampled trees:')
        self.stdout.info('-------------------------')
        self.stdout.info('Tree source: %s' %  str(input_trees))
        self.stdout.info('Total number of trees in file = %d' % num_trees)
        self.stdout.info('Number of trees considered = %d' % self.num_trees_considered)
        self.stdout.info('Number of distinct tree topologies found = %d' % len(tree_map.keys()))
        self.stdout.info('Number of distinct splits found = %d' % (len(split_map.keys()) - t.getNObservables()))
        if self.num_trees_considered == 0:
            self.stdout.info('\nSumT finished.')
            return

        # Sort the splits from highest posterior probabilty to lowest
        split_vect = split_map.items()
        c = lambda x,y: cmp(y[1][0], x[1][0])
        split_vect.sort(c)

        # Output summary of splits
        if sojourn_field_width < 5:
            # must be large enough to accommodate the column header 'freq.'
            sojourn_field_width = 5
        split_fmt_str = '%%%ds' % split_field_width
        sojourn_label_fmt_str = '%%%ds' % sojourn_field_width
        sojourn_fmt_str = '%%%dd' % sojourn_field_width
        self.stdout.info('\nSplit (split), split representation (pattern), frequency (freq.), posterior')
        self.stdout.info('  probability (prob.), mean edge length (weight), first sojourn start (s0),')
        self.stdout.info('  last sojourn end (sk), and number of sojourns (k):')
        split_str = split_fmt_str % 'pattern'
        freq_str = sojourn_label_fmt_str % 'freq.'
        s0_str = sojourn_label_fmt_str % 's0'
        sk_str = sojourn_label_fmt_str % 'sk'
        k_str = sojourn_label_fmt_str % 'k'
        self.stdout.info('%6s %s %s %10s %10s %s %s %s' % ('split', split_str, freq_str, 'prob.', 'weight', s0_str, sk_str, k_str))
        first_below_50 = None
        num_trivial = 0
        split_info = []
        for i,(k,v) in enumerate(split_vect):
            # len(v) is 2 in trivial splits because these splits are associated with tips,
            # for which the sojourn history is omitted (v[0] is frequency and v[1] is edge length sum)
            trivial_split = len(v) == 2 and True or False
            if trivial_split:
                num_trivial += 1

            # Split frequency is simply the first element of the list
            split_freq = v[0]

            # Split posterior is split_freq divided by the number of trees considered
            split_posterior = float(split_freq)/float(self.num_trees_considered)

            # split_weight is the sum of edge lengths v[1] divided by the split_freq
            split_weight = float(v[1])/float(split_freq)

            # Identify the index of the first split having a posterior probability
            # lower than 0.5. The list up to this point will be used to generate the
            # majority-rule consensus tree
            if not first_below_50 and split_posterior < 0.5:
                first_below_50 = i

            # Determine first sojourn (the third element of the list)
            first_sojourn_start = trivial_split and 1 or v[2]

            # Determine last sojourn (the final element of the list)
            last_sojourn_end = trivial_split and self.num_trees_considered or v[-1]

            # Determine the number of sojourns (this must be counted)
            num_sojourns = 1
            if not trivial_split:
                in_sojourn = True
                prev = v[2]
                for curr in v[3:]:
                    if curr - prev > 1:
                        num_sojourns += 1
                    prev = curr

            split_str = split_fmt_str % k
            freq_str = sojourn_fmt_str % split_freq
            s0_str = sojourn_fmt_str % first_sojourn_start
            sk_str = sojourn_fmt_str % last_sojourn_end
            k_str = sojourn_fmt_str % num_sojourns
            split_info.append((split_str, trivial_split, split_freq, split_posterior, split_weight, num_sojourns, first_sojourn_start, last_sojourn_end))
            self.stdout.info('%6d %s %s %10.5f %10.5f %s %s %s' % (i + 1, split_str, freq_str, split_posterior, split_weight, s0_str, sk_str, k_str))

        # Build 50% majority rule tree if requested
        self.stdout.info('\nSaving majority-rule consensus tree...')
        majrule = phylogeny.Tree()
        if self.rooted_trees:
            majrule.setRooted()
        tm = phylogeny.TreeManip(majrule)
        majrule_splits = []
        for k,v in split_vect[:first_below_50]:
            if len(v) > 2:
                majrule_splits.append(k)

        if len(majrule_splits) == 0:
            tm.starTree(num_trivial)
        else:
            tm.buildTreeFromSplitVector(majrule_splits, probdist.Exponential(10))

        self.assignEdgeLensAndSupportValues(majrule, split_map, joint_split_map, self.num_trees_considered, False)

        summary_short_name_list = ['majrule']
        summary_full_name_list = ['Majority-rule Consensus']
        summary_tree_list = [majrule]

        # treelen = self.recordTreeInMaps(t, split_lookup, split_map, joint_split_map, ccd_tree, tree_key)

        # Output summary of tree topologies
        self.stdout.info('\nTree topology (topology), frequency (freq.), mean tree length (TL),')
        self.stdout.info('  first sojourn start (s0), last sojourn end (sk), number of sojourns (k),')
        self.stdout.info('  posterior probability (prob.), cumulative probability (cum.):')
        freq_str = sojourn_label_fmt_str % 'freq.'
        s0_str = sojourn_label_fmt_str % 's0'
        sk_str = sojourn_label_fmt_str % 'sk'
        k_str = sojourn_label_fmt_str % 'k'
        self.stdout.info('%8s %s %10s %s %s %s %10s %10s %10s %10s %10s %10s' % ('topology', freq_str, 'TL', s0_str, sk_str, k_str, 'prob.', 'cum.', 'logccd', 'ccd', 'cumccd', 'maxccd'))
        cum_prob = 0.0
        larget_cum_prob = 0.0
        larget_max_prob = 0.0
        num_distinct_topologies = 0
        lindleyHnaive = 0.0
        KL_sum = 0.0
        KL_max = None
        KL_ntax = None
        done = False
        tree_vect = tree_map.items()
        c = lambda x,y: cmp(y[1][0], x[1][0])
        tree_vect.sort(c)
        if self.opts.tree_credible_prob <= 0.0:
            self.stdout.info('\nNo trees in credible set (tree_credible_prob = 0.0)')
        else:
            for i,(k,v) in enumerate(tree_vect):
                # k is list of splits, v is tuple(count, newick, treelen, num_trees_considered)

                if done:
                    break

                # Determine the posterior probability and cumulative posterior probability
                post_prob = float(v[0])/float(self.num_trees_considered)
                cum_prob += post_prob
                lindleyHnaive -= post_prob*math.log(post_prob)
                if cum_prob > self.opts.tree_credible_prob:
                    # stop after the current tree
                    done = True

                # Determine the average tree length of this tree topology
                avgTL = v[2]/float(v[0])

                # Determine the sampled tree that began the first sojourn (the third element of the list)
                first_sojourn_start = v[3]

                # Determine the sampled tree that ended the last sojourn (the final element of the list)
                last_sojourn_end = v[-1]

                # Determine the number of sojourns (this must be counted)
                num_sojourns = 1
                in_sojourn = True
                prev = v[3]
                for curr in v[4:]:
                    if curr - prev > 1:
                        num_sojourns += 1
                    prev = curr

                # Save the tree topology (decorated with posterior mean edge lengths) to the tree file
                t = phylogeny.Tree()
                if self.rooted_trees:
                    t.setRooted()
                v[1].buildTree(t)
                if KL_max is None:
                    KL_ntax = t.getNTips()
                    KL_max = 0.0
                    if t.isRooted():
                        KL_ntax += 1
                    x = [math.log(y) for y in range(3,2*KL_ntax-4,2)]
                    KL_max = sum(x) # log of (2n-5)!!
                    # raw_input('...')
                log_larget_post_prob = self.assignEdgeLensAndSupportValues(t, split_map, joint_split_map, self.num_trees_considered)
                t.stripNodeNames()
                summary_short_name_list.append('%d %d of %d' % (i+1,v[0], self.num_trees_considered))
                summary_full_name_list.append('Frequency %d of %d' % (v[0], self.num_trees_considered))
                summary_tree_list.append(t)

                # Output summary line for this tree topology
                freq_str = sojourn_fmt_str % v[0]
                s0_str = sojourn_fmt_str % first_sojourn_start
                sk_str = sojourn_fmt_str % last_sojourn_end
                k_str = sojourn_fmt_str % num_sojourns
                larget_post_prob = math.exp(log_larget_post_prob)
                if larget_post_prob > larget_max_prob:
                    larget_max_prob = larget_post_prob
                larget_cum_prob += larget_post_prob
                klterm = larget_post_prob*log_larget_post_prob
                KL_sum += klterm
                num_distinct_topologies += 1
                self.stdout.info('%8d %s %10.5f %s %s %s %10.5f %10.5f %10.5f %10.5f %10.5f %10.5f' % (i+1, freq_str, avgTL, s0_str, sk_str, k_str, post_prob, cum_prob, log_larget_post_prob, larget_post_prob, larget_cum_prob, larget_max_prob))

        if self.opts.tree_credible_prob <= 0.0:
            self.stdout.info('\nKL not available (re-run sumt command with sumt.tree_credible_prob = 1.0)')
        elif self.opts.tree_credible_prob < 1.0:
            self.stdout.info('\nKL should be based on all sampled trees (re-run sumt command with sumt.tree_credible_prob = 1.0)')
        else:
            effective_KLmax = larget_cum_prob*KL_max
            KL = KL_sum + effective_KLmax
            if effective_KLmax > 0.0:
                KL_pct = 100.0*KL/effective_KLmax
            else:
                KL_pct = 0.0

            # Estimate the lower bound for KL
            KLlower = KL_sum + KL_max
            if larget_cum_prob < 1.0:
                try:
                    first_term = math.log(1.0 - larget_cum_prob)
                    second_term = math.log(math.exp(KL_max) - float(num_distinct_topologies))
                    KLlower += (1.0 - larget_cum_prob)*(first_term - second_term)
                except ValueError:
                    # If here, it means that num_distinct_topologies equals total possible topologies, so
                    # second_term is undefined. This just means that larget_cum_prob should be equal to 1
                    # were it not for roundoff error, so we shouldn't have come here anyway.
                    pass
                except OverflowError:
                    first_term = math.log(1.0 - larget_cum_prob)
                    KLlower += (1.0 - larget_cum_prob)*first_term  # check this
            if KLlower < 0.0:
                KLlower = 0.0

            # Estimate the upper bound for KL
            KLupper = 0.0
            KLupper = KL_sum + KL_max
            if larget_cum_prob < 1.0:
                KLupper += (1.0 - larget_cum_prob)*(math.log(1.0 - larget_cum_prob))
            if KLupper > KL_max:
                KLupper = KL_max

            KL_near_upper = self.calcKLupper(joint_split_map, KL_max)

            # temporary!
            #self.debugShowCCDTree(ccd_tree)

            lindleyH, lindleyHp, lindleyI = ccd_tree.calcI(self.stdout)
            lindleyH0 = KL_max
            lindleyIpct = 100.0*lindleyI/KL_max
            self.stdout.info('\nLindley information based on conditional clade distribution:')
            self.stdout.info('  Posterior entropy:   %.5f' % lindleyH)
            self.stdout.info('  Prior entropy:       %.5f' % lindleyH0)
            self.stdout.info('  Information gain:    %.5f' % lindleyI)
            self.stdout.info('  %% max. information:  %.5f' % lindleyIpct)

            self.stdout.info('\nLarget coverage:')
            self.stdout.info('  Distinct tree topologies: %d' % num_distinct_topologies)
            self.stdout.info('  Posterior coverage:       %.5f' % larget_cum_prob)

            lindleyInaive = -(lindleyHnaive - lindleyH0)
            lindleyIpctnaive = 100.0*lindleyInaive/KL_max
            self.stdout.info('\nLindley information based on observed tree topology distribution:')
            self.stdout.info('  Posterior entropy:   %.5f' % lindleyHnaive)
            self.stdout.info('  Prior entropy:       %.5f' % lindleyH0)
            self.stdout.info('  Information gain:    %.5f' % lindleyInaive)
            self.stdout.info('  %% max. information:  %.5f' % lindleyIpctnaive)

            #self.stdout.info('\nTopological information content:')
            #self.stdout.info('  Number of taxa: %.5f' % KL_ntax)
            #self.stdout.info('  KL maximum (log of total number of distinct tree topologies): %.5f' % KL_max)
            #self.stdout.info('  Naive estimate based on %d distinct tree topologies:' % num_distinct_topologies)
            #self.stdout.info('    KL naive estimate:           %.5f' % (KL_sum + KL_max,))
            #self.stdout.info('    KL naive (%% of max.):       %.5f' % (100.0*(KL_sum + KL_max)/KL_max,))
            #self.stdout.info('  Bounds based on %d distinct tree topologies and %.5f cumulative probability:' % (num_distinct_topologies, larget_cum_prob))
            #self.stdout.info('    KL lower bound (%% of max.): %.5f' % (100.0*KLlower/KL_max,))  # not correct for polytomy analyses
            #self.stdout.info('    KL upper bound (%% of max.): %.5f' % (100.0*KLupper/KL_max,))  # not correct for polytomy analyses
            #if KL_near_upper < KLupper:
            #    self.stdout.info('  Near upper bound based on conditional clade probabilities:')
            #    self.stdout.info('    KL near upper bound:              %.5f' % KL_near_upper)
            #    self.stdout.info('    KL near upper bound (%% of max.): %.5f' % (100.0*KL_near_upper/KL_max,))  # not correct for polytomy analyses
            #self.stdout.info('  Notes:')
            #self.stdout.info('  1) the naive estimate uses sample frequencies of tree topologies (does not account for unsampled posterior mass);')
            #self.stdout.info('  2) the lower bound assumes all unsampled tree topologies have equal posterior probability;')
            #self.stdout.info('  3) the upper bound assumes one unsampled tree topology has all unsampled posterior mass;')
            #if KL_near_upper < KLupper:
            #    self.stdout.info('  4) the near upper bound assumes 1/p tree topologies have maximum posterior probability p,')
            #    self.stdout.info('     and all remaining topologies have posterior zero, where p is determined from sample')
            #    self.stdout.info('     conditional clade distribution')

        self.stdout.info('\nSaving distinct tree topologies...')
        self.save_trees(summary_short_name_list, summary_full_name_list, summary_tree_list)

        self._splitsPdfWriter = None
        self._splitsPdfFilename = None

        if self.opts.save_splits_pdf and bool(self.optsout.splits):
            try:
                pdf_source = lambda : self._getSplitsPDFWriter()
                awty_written = self.awtyPlot(pdf_source, split_vect, self.num_trees_considered)
                sojourn_written = self.sojournPlot(pdf_source, split_vect, self.num_trees_considered)
            finally:
                if self._splitsPdfWriter:
                    self.optsout.splits.close()

        self.stdout.info('\nSumT finished.')
        return split_info

    def _getSplitsPDFWriter(self):
        if self._splitsPdfWriter is None:
            sp = self.optsout.splits
            if sp.replaceMode() and sp.sameBaseName(self.optsout.trees):
                self.stdout.info("Splits pdf will overwrite the tree topology pdf")
            self._splitsPdfFilename = sp._getFilename()
            self._splitsPdfWriter = sp.open(11.0, 8.5, self.stdout)
        return self._splitsPdfWriter




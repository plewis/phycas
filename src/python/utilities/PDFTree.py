import math, types
from phycas.PDFGen import *
from phycas.Utilities.CommonFunctions import CommonFunctions
from phycas.Utilities.GlobalState import readFile
from phycas.Phylogeny import Tree

class PDFTree(CommonFunctions):
    def __init__(self):
        CommonFunctions.__init__(self, None)
        self.pdf_splits_to_plot = None
        # Variables associated with PDF tree drawing (used in pdftree() function)
        # The 14 standard fonts guaranteed to be available in all PDF consumer applications:
        #   Times-Roman      Helvetica             Courier             Symbol
        #   Times-Bold       Helvetica-Bold        Courier-Bold        ZapfDingbats
        #   Times-Italic     Helvetica-Oblique     Courier-Oblique
        #   Times-BoldItalic Helvetica-BoldOblique Courier-BoldOblique
        self.pdf_filename                 = 'trees.pdf'    # Set to desired name of pdf file to create
        self.pdf_edge_support_file        = None           # File containing PAUP* output with table of support values; if specified, the support values will be shown on trees plotted
        self.pdf_tip_label_font           = 'Times-Italic' # Font used for tip node names; should be one of the 14 standard fonts listed above
        self.pdf_tip_label_height         = 12             # Height in points of tip node name font
        self.pdf_plot_label_font          = 'Helvetica'    # Font used for plot axis labels; should be one of the 14 standard fonts listed above
        self.pdf_plot_label_height        = 12             # Height in points of plot axis label font
        self.pdf_title_font               = 'Helvetica'    # Font used for scalebar text; should be one of the 14 standard fonts listed above
        self.pdf_title_height             = 14             # Height in points of scalebar text font
        self.pdf_scalebar_position        = 'bottom'       # Valid values are 'top', 'bottom' or None
        self.pdf_scalebar_label_font      = 'Helvetica'    # Font used for scalebar text; should be one of the 14 standard fonts listed above
        self.pdf_scalebar_label_height    = 10             # Height in points of scalebar text font
        self.pdf_support_label_font       = 'Times-Roman'  # Font used for edge support values; should be one of the 14 standard fonts listed above
        self.pdf_support_label_height     = 8              # Height in points of edge support font
        self.pdf_support_as_percent       = True           # If True, support values will be shown as percentages (e.g. 93.1) rather than proportions (e.g. 0.931)
        self.pdf_support_decimals         = 1              # The number of decimal places shown in support values (e.g. to get 93.7, specify 1; to round up to 94, specify 0)
        self.pdf_ladderize                = 'right'        # Valid values are 'right', 'left' or None
        self.pdf_page_width               = 8.5            # Page width in inches
        self.pdf_page_height              = 11.0           # Page length in inches
        self.pdf_line_width               = 1.0            # Width of lines representing edges in the tree
        self.pdf_left_margin              = 1.0            # Left margin in inches (1 inch = 72 points)
        self.pdf_right_margin             = 1.0            # Right margin in inches (1 inch = 72 points)
        self.pdf_top_margin               = 1.0            # Top margin in inches (1 inch = 72 points)
        self.pdf_bottom_margin            = 1.0            # Bottom margin in inches (1 inch = 72 points)
        self.keep_xy_proportional         = True           # If True, vertical dimension of each tree in a collection will be kept proportional to its horizontal dimension
        self.keep_tip_labels_proportional = True           # If True, tip label height will be kept commensurate with size of tree for each tree in a printed collection (smaller trees will have smaller tip labels)
        self.pdf_treefile                 = None           # Set to tree file name if you want to make one pdf file with each tree from tree file on a separate page
        self.pdf_newick                   = None           # Set to the tree description to print if only want to save one tree to a pdf file
        self.pdf_outgroup_taxon           = None           # Set to taxon name of tip serving as the outgroup for display rooting purposes (note: at this time outgroup can consist of just one taxon)

    def pdftree(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a PDF file containing a single tree (if pdf_newick is 
        specified) or a collection of trees (if pdf_treefile is specified). 
        If a collection of trees is specified, scales all trees the same (i.e
        the scalebar is identical in size for all trees plotted).
        """
        #complex_outgroup = type(self.pdf_outgroup_taxon) in (types.ListType,types.TupleType)
        simple_outgroup = type(self.pdf_outgroup_taxon) == types.StringType
        self.phycassert(simple_outgroup, 'Phycas cannot yet deal with pdf_outgroup_taxon containing more than one outgroup taxon')
        self.phycassert((self.pdf_treefile and not self.pdf_newick) or (self.pdf_newick and not self.pdf_treefile), 'set either pdf_newick or pdf_treefile, but not both')
        
        # If pdf_edge_support_file has been specified, read splits table from the file
        # and store the splits in the pdf_splits_to_plot dictionary
        if self.pdf_edge_support_file and os.path.exists(self.pdf_edge_support_file):
            # Read splits file and store all splits found along with their frequencies
            contents_of_file = open(self.pdf_edge_support_file,'r').read()
            regex = re.compile('([*.]+)\s+([0-9.]+)', re.M)
            matches = regex.findall(contents_of_file)
            self.phycassert(matches, 'could not find any splits defined in the pdf_edge_support_file named %s' % self.pdf_edge_support_file)
            self.pdf_splits_to_plot = {}
            for p,f in matches:
                self.pdf_splits_to_plot[p] = float(f)
                
        # Fork depending on whether user wants to print just one tree (pdf_newick specified)
        # or an entire collection of trees (pdf_treefile specified)
        if self.pdf_newick:        
            # Build tree the newick description of which is in self.newick
            tree = self.pdf_newick.buildTree()
            
            if self.pdf_outgroup_taxon:
                num = tree.findTipByName(self.pdf_outgroup_taxon)
                self.phycassert(num is not None, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.pdf_outgroup_taxon)
                tree.rerootAtTip(num)
                
            if self.pdf_ladderize:
                if self.pdf_ladderize == 'right':
                    tree.ladderizeRight()
                else:
                    tree.ladderizeLeft()

            # Save tree in PDF  
            pdf = PDFGenerator(self.pdf_page_width, self.pdf_page_height)
            pdf.overwrite = True
            pdf.newPage()
            self.tree2pdf(pdf, tree)
            pdf.saveDocument(self.pdf_filename)
        else:
            # Open pdf_treefile and read trees therein
            self.tree_file_name = self.pdf_treefile
            contents = readFile(self.pdf_treefile)

            # Build each tree and determine its height
            tree = Tree()
            max_height = 0.0
            for tree_def in contents.trees:
                tree_def.buildTree(tree)
                tree.rectifyNames(contents.taxon_labels)
                if self.pdf_outgroup_taxon:
                    num = tree.findTipByName(self.pdf_outgroup_taxon)
                    self.phycassert(num is not None, 'could not root tree using specified outgroup: no tip having name "%s" could be found' % self.pdf_outgroup_taxon)
                    tree.rerootAtTip(num)
                h = tree.calcTotalHeight()
                if h > max_height:
                    max_height = h
                #tlen = tree.edgeLenSum()
                #print 'tlen =',tlen,', height =',h

            # Build each tree again and save in PDF file            
            pdf = PDFGenerator(self.pdf_page_width, self.pdf_page_height)
            pdf.overwrite = True
            for tree_def in contents.trees:
                tree_def.buildTree(tree)
                tree.rectifyNames(contents.taxon_labels)
                if self.pdf_outgroup_taxon:
                    num = tree.findTipByName(self.pdf_outgroup_taxon)
                    tree.rerootAtTip(num)
                if self.pdf_ladderize:
                    if self.pdf_ladderize == 'right':
                        tree.ladderizeRight()
                    else:
                        tree.ladderizeLeft()
                tree.rectifyNames(contents.taxon_labels)
                pdf.newPage()
                self.tree2pdf(pdf, tree, None, max_height)
            pdf.saveDocument(self.pdf_filename)
            
        # Prevent unintentional spillover
        self.pdf_splits_to_plot = None
        
    def tree2pdf(self, pdf, tree, title = None, xscalemax = 0.0, show_support = False):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Prints tree on a pdf object (instance of class PDFGenerator). If title
        is specified, the supplied string will be centered at the top of the
        page. The optional argument xscalemax represents the maximum height
        of a group of trees being printed on separate pages in the same pdf
        document. If xscalemax is left unspecified, each tree will be scaled
        to fit the page and the scalebar will be adjusted accordingly. If
        xscalemax is specified, it will be used to determine the scalebar, and
        the scalebar will remain the same size for all trees printed with the
        same xcalemax value.
        
        """
        # TODO: max_label_points should be calculated outside this function and passed in as an argument
        inch = 72.0
        spacer = 5.0
        max_label_points = 0.0
        rooted_tree = tree.isRooted()
        nodes = []

        # Perform a preorder traversal:
        # 1) for each node, set x-value to height above root (in units of edge length)
        # 2) for each tip, set y-value to tip index, with root tip being 0, and other
        #    tips being numbered from left to right
        # 3) find the length of the longest taxon label as it will be rendered in the
        #    PDF file so that the margin calculations can be made
        # 4) for each internal, just set y-value to 0.0 for now; these internal y-values
        #    will be calculated on the subsequent postorder traversal

        if self.pdf_splits_to_plot:
            tree.recalcAllSplits(tree.getNObservables())
            
        # Record information about the tip serving as the root
        nd = tree.getFirstPreorder()
        assert nd.isRoot(), 'first preorder node should be the root'
        if not rooted_tree:
            nodes.append(nd)
        subroot = nd.getLeftChild()
        height = subroot.getEdgeLen()
        nd.setX(height) 
        if self.pdf_ladderize and self.pdf_ladderize == 'left':
            last_tip_index = float(tree.getNObservables() - 1)
            nd.setY(last_tip_index) #--> Y is irrelevant if rooted
            ntips = 0.0
        else:
            nd.setY(0.0)
            if rooted_tree:
                ntips = 0.0
            else:
                ntips = 1.0
        max_height = height
        
        # Determine the width (in points) occupied by the longest taxon label 
        if self.pdf_tip_label_font and not rooted_tree:
            taxon_label = nd.getNodeName()
            label_width = float(self.pdf_tip_label_height)*pdf.calcStringWidth(self.pdf_tip_label_font, taxon_label)
            if label_width > max_label_points:
                max_label_points = label_width
        
        # Record information about the internal node serving as the subroot
        nd = nd.getNextPreorder()
        assert nd.getParent().isRoot(), 'second preorder node should be the subroot'
        nodes.append(nd)
        nd.setX(0.0)
        nd.setY(0.0)
        subroot = nd
        
        # Record information about the remaining nodes in the tree
        while True:
            nd = nd.getNextPreorder()
            if not nd:
                break
            else:
                ndpar = nd.getParent()
                nodes.append(nd)
                height = nd.getEdgeLen() + ndpar.getX()
                nd.setX(height)
                if height > max_height:
                    max_height = height
                if nd.isTip():
                    nd.setY(ntips)
                    ntips += 1.0
                    if self.pdf_tip_label_font:
                        taxon_label = nd.getNodeName()
                        label_width = float(self.pdf_tip_label_height)*pdf.calcStringWidth(self.pdf_tip_label_font, taxon_label)
                        if label_width > max_label_points:
                            max_label_points = label_width
                else:
                    nd.setY(0.0)

        # Compute length represented by scale bar. For example,
        #  xscalemax     = 0.00275
        #  log_xscalemax = -2.56
        #  ten_to_power  = 10^floor(-2.56)
        #                = 10^{-3}
        #                = 0.001
        #  scalebar      = 0.001*floor(0.00275/0.001)
        #                = 0.001*floor(2.75)
        #                = 0.002
        #  ndecimals     = -floor(-2.56)
        #                = 3.0
        if xscalemax == 0.0:
            xscalemax = max_height
        half_xscalemax = xscalemax/2.0
        log_xscalemax = math.log10(half_xscalemax)
        ten_to_power = 10**math.floor(log_xscalemax)
        scalebar = ten_to_power*math.floor(half_xscalemax/ten_to_power)
        ndecimals = -int(math.floor(log_xscalemax))
        if ndecimals < 0:
            ndecimals = 0
        format_str = '%%.%df' % (ndecimals)
        scalebar_str = format_str % scalebar
        scalebar_str_extent = float(self.pdf_scalebar_label_height)*pdf.calcStringWidth(self.pdf_scalebar_label_font, scalebar_str)
        scalebar_height = float(self.pdf_scalebar_label_height) + 2*spacer + self.pdf_line_width

        # Find xscaler (amount by which branch lengths must be multiplied to give x-coordinate)
        # and yscaler (amount by which the tip position must be multiplied to give y-coordinate).
        xheight = 0.0
        if self.pdf_tip_label_font:
            xheight = float(self.pdf_tip_label_height)*pdf.getXHeight(self.pdf_tip_label_font)
        half_xheight = xheight/2.0
        ntips = tree.getNObservables()
        label_width   = max_label_points + spacer
        right_margin  = self.pdf_right_margin*inch
        left_margin   = self.pdf_left_margin*inch
        top_margin    = self.pdf_top_margin*inch
        bottom_margin = self.pdf_bottom_margin*inch
        plot_right = self.pdf_page_width*inch
        plot_width = plot_right - left_margin - right_margin
        plot_top = self.pdf_page_height*inch
        plot_height = plot_top - top_margin - bottom_margin

        tree_width = plot_width - label_width
        tree_height = plot_height
        if self.pdf_scalebar_position:
            tree_height -= scalebar_height
        if title:
            tree_height -= 3.0*float(self.pdf_title_height)
        tree_x0 = left_margin
        tree_y0 = bottom_margin + scalebar_height

        xscaler = tree_width/xscalemax
        yscaler = tree_height/float(ntips - 1)

        #pdf.addRectangle(left_margin, bottom_margin, plot_width, plot_height, 1, 'dotted')

        if title and self.pdf_title_height > 0:
            # Draw title centered at top of page
            title_str_extent = float(self.pdf_title_height)*pdf.calcStringWidth(self.pdf_title_font, title)
            title_x = left_margin + (plot_width - title_str_extent)/2.0
            title_y = tree_y0 + tree_height + 2.0*float(self.pdf_title_height)
            pdf.addText(title_x, title_y, self.pdf_title_font, self.pdf_title_height, title)

        if self.pdf_scalebar_position:
            if self.pdf_scalebar_position == 'top':
                # Draw scalebar horizontally starting at top left corner
                scalebar_width = scalebar*xscaler
                scalebar_y = tree_x0 + tree_height - scalebar_height + spacer
                pdf.addLine(left_margin, scalebar_y, left_margin + scalebar_width, scalebar_y, self.pdf_line_width)

                # Draw scalebar text centered above the scalebar
                scalebar_x = left_margin + (scalebar_width - scalebar_str_extent)/2.0
                scalebar_y = tree_x0 + tree_height - float(self.pdf_scalebar_label_height)
                pdf.addText(scalebar_x, scalebar_y, self.pdf_scalebar_label_font, self.pdf_scalebar_label_height, scalebar_str)
            else:
                # Draw scalebar horizontally starting at bottom left corner
                scalebar_width = scalebar*xscaler
                pdf.addLine(left_margin, bottom_margin, left_margin + scalebar_width, bottom_margin, self.pdf_line_width)

                # Draw scalebar text centered above the scalebar
                scalebar_x = left_margin + (scalebar_width - scalebar_str_extent)/2.0
                scalebar_y = bottom_margin + spacer
                pdf.addText(scalebar_x, scalebar_y, self.pdf_scalebar_label_font, self.pdf_scalebar_label_height, scalebar_str)

        # add enough to left margin to center smaller trees horizontally
        left_margin += (xscaler*(xscalemax - max_height) + label_width*(1.0 - max_height/xscalemax))/2.0

        # add enough to the top margin to center smaller trees vertically
        top_margin += (tree_height*(1.0 - max_height/xscalemax))/2.0
        #top_margin += (plot_height*(1.0 - max_height/xscalemax))/2.0

        # adjust yscaler to keep vertical tree dimension proportional to its horizontal dimension
        if self.keep_xy_proportional:
            yscaler *= max_height/xscalemax

        # adjust tip label height (in points) to make size of tip labels commensurate with size of tree
        if self.keep_tip_labels_proportional:
            tip_font_points = self.pdf_tip_label_height*max_height/xscalemax
        else:
            tip_font_points = self.pdf_tip_label_height
        
        # Perform a postorder traversal:
        # 1) scale each x-value
        # 2) calculate y-value of each internal node as the average y-value of its children
        # 3) scale each y-value
        # 4) plot each edge
        # 5) plot names of tips
        # 6) for each internal node, draw shoulder from leftmost child to rightmost
        nodes.reverse()
        for nd in nodes:
            node_x = left_margin + nd.getX()*xscaler
            if nd.isTip():
                node_y = tree_y0 + tree_height - nd.getY()*yscaler
                if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                    node_y -= scalebar_height
                brlen = nd.isRoot() and xscaler*nd.getX() or xscaler*nd.getEdgeLen()
                # draw tip node name
                if self.pdf_tip_label_font:
                    pdf.addText(node_x + spacer, node_y - half_xheight, self.pdf_tip_label_font, tip_font_points, nd.getNodeName())
                # draw line representing edge leading to tip node
                pdf.addLine(node_x, node_y, node_x - brlen, node_y, self.pdf_line_width)
            else:
                nchildren = 1.0
                child = nd.getLeftChild()
                left_child = right_child = child
                childY = child.getY()
                while True:
                    child = child.getRightSib()
                    if child:
                        right_child = child
                        childY += child.getY()
                        nchildren += 1.0
                    else:
                        break
                if (not rooted_tree) and (nd is subroot):
                    if self.pdf_ladderize and self.pdf_ladderize == 'left':
                        right_child = nd.getParent()
                    else:
                        left_child = nd.getParent()
                else:
                    nd.setY(childY/nchildren)
                    node_y = tree_y0 + tree_height - childY*yscaler/nchildren
                    if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                        node_y -= scalebar_height
                    brlen = xscaler*nd.getEdgeLen()
                    # draw line representing edge leading to internal node
                    pdf.addLine(node_x, node_y, node_x - brlen, node_y, self.pdf_line_width)

                # draw line representing shoulders of internal node
                left_y = tree_y0 + tree_height - left_child.getY()*yscaler
                right_y = tree_y0 + tree_height - right_child.getY()*yscaler
                if self.pdf_scalebar_position and self.pdf_scalebar_position == 'top':
                    left_y -= scalebar_height
                    right_y -= scalebar_height
                pdf.addLine(node_x, left_y, node_x, right_y, self.pdf_line_width)

                # if specified, plot support value
                if show_support and self.pdf_splits_to_plot:
                    for p in self.pdf_splits_to_plot.keys():
                        s = Split()
                        s.setOnSymbol('*')
                        s.setOffSymbol('.')
                        s.createFromPattern(p)
                        if s.equals(nd.getSplit()):
                            support_x = node_x + spacer
                            support_y = (left_y + right_y)/2.0 - half_xheight
                            support_str = '%.1f' % self.pdf_splits_to_plot[p]
                            pdf.addText(support_x, support_y, self.pdf_support_label_font, self.pdf_support_label_height, support_str)
                            break
                elif show_support and nd is not subroot:
                    # Expecting each node's support data member to be set already
                    support_format = '%%.%df' % self.pdf_support_decimals
                    if self.pdf_support_as_percent:
                        support_str = support_format % (100.0*nd.getSupport(),)
                    else:
                        support_str = support_format % (nd.getSupport(),)
                    support_str_extent = float(self.pdf_support_label_height)*pdf.calcStringWidth(self.pdf_support_label_font, support_str)
                    support_x = node_x - (brlen + support_str_extent)/2.0
                    support_y = (left_y + right_y)/2.0 + half_xheight
                    pdf.addText(support_x, support_y, self.pdf_support_label_font, self.pdf_support_label_height, support_str)

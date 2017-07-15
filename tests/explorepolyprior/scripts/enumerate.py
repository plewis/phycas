import sys, re

#debugging = True
debugging = False

class Preorder(object):
    def __init__(self, t):
        self.tree = t
        self.pre = None
    def reset(self):
        self.pre = []
    def processNode(self, nd):
        self.pre.append(nd)
    def build(self):
        self.tree.traverse(self)
        return self.pre

class Node(object):
    def __init__(self):
        self.par = None
        self.lchild = None
        self.rsib = None
        self.number = None
        self.edgelen = 0.0
        self.split = None
    def clear(self):
        self.__init__()
    def getSplit(self):
        return split

class Tree(object):
    def __init__(self):
        self.root = None
        self.preorder = None
        self.id = None

    def ladderize(self):
        #print '######## begin ladderize ##########'
        #print 'before newick =',self.makeNewick()
        for nd in self.preorder:
            if nd.lchild is not None:
                nodelist = []
                child = nd.lchild
                while child is not None:
                    if child.lchild is None:
                        nodelist.append((1,[child.number],child))
                    else:
                        assert child.split is not None, 'Tree.ladderize expects splits to be calculated for each internal node'
                        nodelist.append((len(child.split),list(child.split),child))
                    child = child.rsib
                #print 'before'
                #for a,b,c in nodelist:
                #    print a,b,c.number
                nodelist.sort()
                #print 'after'
                #for a,b,c in nodelist:
                #    print a,b,c.number
                #raw_input('..')
                nchildren = len(nodelist)
                nd.lchild = nodelist[0][2]
                c = nodelist[0][2]
                for i in range(1,nchildren):
                    next_c = nodelist[i][2]
                    c.rsib = next_c
                    c = next_c
                c.rsib = None
        #print 'after newick =',self.makeNewick()
        #raw_input('######## end ladderize ##########')

    def calcSplits(self):
        # make sure every internal node has a valid split set
        # storing all tip descendant node numbers
        self.preorder = Preorder(self).build()
        postorder = self.preorder[:]
        postorder.reverse()
        self.id = []

        for nd in postorder:

            if nd.lchild is None:
                if nd.rsib is None:
                    nd.par.split = set([nd.number])
                else:
                    nd.par.split |= set([nd.number])

            else:
                if nd.par is not None:
                    if nd.rsib is None:
                        nd.par.split = set(nd.split)
                    else:
                        nd.par.split |= nd.split

                self.id.append(frozenset(nd.split))
        self.ladderize()

    def deepCopy(self, nd = None):
        global newroot, currnd, currpar

        if nd is None:
            nd = self.root
            newroot = Node()
            currnd = newroot
            currpar = None

        if nd.lchild is not None:
            currpar = currnd
            currnd.lchild = Node()
            currnd = currnd.lchild
            currnd.par = currpar
            self.deepCopy(nd.lchild)
        else:
            # nd is a tip
            currnd.number = nd.number

        if nd.rsib is not None:
            currnd.rsib = Node()
            currnd = currnd.rsib
            currnd.par = currpar
            self.deepCopy(nd.rsib)
        elif nd.par is not None:
            currnd = currnd.par
            currpar = currnd.par
        else:
            t = Tree()
            t.root = newroot
            t.calcSplits()
            return t

    def traverse(self, handler, nd = None):
        global internal_node_number

        if nd is None:
            nd = self.root
            handler.reset()
            internal_node_number = -1

        if nd.lchild is not None:
            nd.number = internal_node_number
            internal_node_number -= 1

        handler.processNode(nd)

        if nd.lchild is not None:
            self.traverse(handler, nd.lchild)

        if nd.rsib is not None:
            self.traverse(handler, nd.rsib)

    def buildFromNewick(self, newick):
        self.root = Node()
        self.preorder = None
        self.id = None

        # Some flags to keep track of what we did last
        Prev_Tok_LParen		= 1    # previous token was a left parenthesis ('(')
        Prev_Tok_RParen		= 2    # previous token was a right parenthesis (')')
        Prev_Tok_Colon		= 3    # previous token was a colon (':')
        Prev_Tok_Comma		= 4    # previous token was a comma (',')
        Prev_Tok_Name		= 5    # previous token was a node name (e.g. '2', 'P._articulata')
        Prev_Tok_EdgeLen	= 6	   # previous token was an edge length (e.g. '0.1', '1.7e-3')

        # Some useful flag combinations
        LParen_Valid = set([Prev_Tok_LParen, Prev_Tok_Comma])
        RParen_Valid = set([Prev_Tok_RParen, Prev_Tok_Name,  Prev_Tok_EdgeLen])
        Comma_Valid  = set([Prev_Tok_RParen, Prev_Tok_Name,  Prev_Tok_EdgeLen])
        Colon_Valid  = set([Prev_Tok_RParen, Prev_Tok_Name])
        Name_Valid   = set([Prev_Tok_RParen, Prev_Tok_LParen,  Prev_Tok_Comma])

        previous = Prev_Tok_LParen

        i = 0
        nd = self.root
        while i < len(newick):
            ch = newick[i]

            if ch == ' ':
                pass

            elif ch == ';':
                break

            elif ch == ')':
                # If nd is bottommost node, expecting left paren or semicolon, but not right paren
                assert nd.par is not None, 'Too many right parentheses at position %d in tree description' % i

                # Expect right paren only after an edge length, a node name, or another right paren
                assert previous in RParen_Valid, 'Unexpected right parenthesis at position %d in tree description' % i

                # Go down a level
                nd = nd.par
                assert nd.lchild.rsib is not None, 'Internal node has only one child at position %d in tree description' % i

                previous = Prev_Tok_RParen

            elif ch == ':':
                # Expect colon only after a node name or another right paren
                assert previous in Colon_Valid, 'Unexpected colon at position %d in tree description' % i
                previous = Prev_Tok_Colon

            elif ch == ',':
                # Expect comma only after an edge length, a node name, or a right paren
                assert nd.par is not None and previous in Comma_Valid, 'Unexpected comma at position %d in tree description' % i

                # Create the sibling
                nd.rsib = Node()
                nd.rsib.par = nd.par
                nd = nd.rsib
                previous = Prev_Tok_Comma

            elif ch == '(':
                # Expect left paren only after a comma or another left paren
                assert previous in LParen_Valid, 'Not expecting left parenthesis at position %d in tree description' % i

                # Create new node above and to the left of the current node
                assert nd.lchild is None
                nd.lchild = Node()
                nd.lchild.par = nd
                nd = nd.lchild
                previous = Prev_Tok_LParen

            elif ch == "'":
                # Encountered an apostrophe, which always indicates the start of a
                # node name (but note that node names do not have to be quoted)

                # Expect node name only after a left paren (child's name), a comma (sib's name)
                # or a right paren (parent's name)
                assert previous in Name_Valid, 'Not expecting node name at position %d in tree description' % i

                # Get the rest of the name
                name = ''
                i += 1
                while i < len(newick):
                    ch = newick[i]
                    if ch == "'":
                        break
                    elif ch in ' \n\t':
                        name += ' '
                    else:
                        name += ch
                    i += 1
                assert ch == "'", 'Expecting single quote to mark the end of node name at position %d in tree description' % i

                if nd.lchild is None:
                    nd.number = int(name)

                previous = Prev_Tok_Name

            else:
                # Expecting either an edge length or an unquoted node name
                if previous == Prev_Tok_Colon:
                    # Edge length expected (e.g. "235", "0.12345", "1.7e-3")
                    j = i
                    while i < len(newick):
                        ch = newick[i]
                        if ch == ',' or ch == ')' or ch in ' \t\n':
                            i -= 1
                            break
                        assert ch =='e' or ch == 'E' or ch =='.' or ch == '-' or ch == '+' or ch in '0123456789', 'Invalid branch length character (%c) at position %d in tree description' % i
                        i += 1

                    edge_length_str = newick[j:i+1]
                    nd.edgelen = float(edge_length_str)
                    if nd.edgelen < 1.e-10:
                        ndedgelen = 1.e-10

                    previous = Prev_Tok_EdgeLen
                else:
                    # Get the node name
                    name = ''
                    while i < len(newick):
                        ch = newick[i]
                        assert ch != '(', 'Unexpected left parenthesis inside node name at position %d in tree description' % i

                        if ch in ' \t\n' or ch == ':' or ch == ',' or ch == ')':
                            i -= 1;
                            break
                        name += ch;
                        i += 1

                    # Expect node name only after a left paren (child's name), a comma (sib's name) or a right paren (parent's name)
                    assert previous in Name_Valid, 'Unexpected node name (%s) at position %d in tree description' % (name, i)

                    if nd.lchild is None:
                        nd.number = int(name)

                    previous = Prev_Tok_Name
            if i == len(newick):
                break
            i += 1
        self.calcSplits()

    def rerootHelper(self, m, t):
        assert m is not None
        assert t is not None

        # Save nodes to which m attaches
        m_left_child = m.lchild
        m_right_sib  = m.rsib
        m_parent     = m.par

        # Starting with t, walk down tree to identify x, the child of m that is on the path from m to t
        x = t
        while x.par != m:
            x = x.par
            assert x is not None
        x_right_sib = x.rsib

        # Identify x_left_sib, the immediate left sibling of x (will be NULL if x is _left_child of m)
        x_left_sib = None
        if x != m_left_child:
            x_left_sib = m_left_child
            while x_left_sib.rsib != x:
                x_left_sib = x_left_sib.rsib
                assert x_left_sib is not None

        # identify m_left_sib, the immediate left sibling of m (will be NULL if m is root node or is _left_child of its parent)
        m_left_sib = None
        if (m_parent is not None) and (m != m_parent.lchild):
            m_left_sib = m_parent.lchild
            while m_left_sib.rsib != m:
                m_left_sib = m_left_sib.rsib
                assert m_left_sib is not None

        # Put x where m is now
        if m_parent is None:
            # m is the root node
            assert m_right_sib is None
            assert m_left_sib is None
            x.rsib = None
            x.par = None
            if x == m_left_child:
                m.lchild = x_right_sib
            else:
                x_left_sib.rsib = x_right_sib
        elif m == m_parent.lchild:
            # m is leftmost child of its parent
            x.rsib = m_right_sib
            x.par = m_parent
            m.rsib = None
            m.par = None
            m_parent.lchild = x
            if x == m_left_child:
                m.lchild = x_right_sib
            else:
                x_left_sib.rsib = x_right_sib
        else:
            # m is not leftmost child of its parent
            m_left_sib.rsib = x
            x.rsib = m_right_sib
            x.par = m_parent
            m.rsib = None
            m.par = None
            if x == m_left_child:
                m.lchild = x_right_sib
            else:
                x_left_sib.rsib = x_right_sib

        # Make m the new rightmost child of t
        m.par = t
        if t.lchild is None:
            t.lchild = m
        else:
            # Find rightmost child of t
            m_left_sib = t.lchild
            while m_left_sib.rsib:
                m_left_sib = m_left_sib.rsib

            # Make rightmost child of t the left sib of m
            m_left_sib.rsib = m

    def rerootAt(self, which):
        # for nd in self.preorder:
        #     if nd.lchild is None:
        #         if nd == self.preorder[-1]:
        #             print '%d' % nd.number
        #         else:
        #             print '%d ->' % nd.number,
        #     else:
        #         print '(%d) ->' % nd.number,
        # raw_input('debug stop before rerootAt')

        # Locate node having number which that will be the new root node
        new_root = self.root
        for nd in self.preorder:
            if nd.number == which:
                new_root = nd
        assert new_root != self.root, "no node found with number equal to %d" % which

        nd = new_root
        t = new_root
        m = new_root.par
        while nd.par is not None:
            # Begin by swapping the mover's edge length with nd's edge length
            tmp_edgelen = m.edgelen
            m.edgelen = nd.edgelen
            nd.edgelen = tmp_edgelen

            # Make the "mover" node m (along with all of its children except nd) the rightmost child of the "target" node t
            self.rerootHelper(m, t)

            # The next target node is always the previous mover, and the next mover node is always nd's parent
            t = m
            m = nd.par
        assert nd.lchild is not None
        assert nd.rsib is None
        assert nd.par is None
        assert nd.lchild.rsib is None
        self.root = nd.lchild
        nd.lchild = None
        self.calcSplits()

        # for nd in self.preorder:
        #     if nd.lchild is None:
        #         if nd == self.preorder[-1]:
        #             print '%d' % nd.number
        #         else:
        #             print '%d ->' % nd.number,
        #     else:
        #         print '(%d) ->' % nd.number,
        # raw_input('debug stop after rerootAt')

    def makeNewick(self, nd = None):
        global newick
        if nd is None:
            nd = self.root
            newick = ''

        if nd.lchild:
            newick += '('
            self.makeNewick(nd.lchild)
        else:
            newick += '%d' % nd.number

        if nd.rsib:
            newick += ','
            self.makeNewick(nd.rsib)
        elif nd.par is not None:
            newick += ')'

        return newick

    def deleteNode(self, nd):
        # nd is leftmost child
        #     nd has 2 or more siblings
        #         nd is a leaf
        #         nd is internal
        #     nd has only 1 sibling
        #         nd's parent is the root
        #             nd is a leaf
        #             nd is internal
        #         nd's parent is NOT the root
        #             nd.par is leftmost child
        #                 nd is a leaf
        #                 nd is internal
        #             nd.par is rightmost child
        #                 nd is a leaf
        #                 nd is internal
        #             nd.par is a middle child
        # nd is rightmost child
        #     nd has 2 or more siblings
        #         nd is a leaf
        #         nd is internal
        #     nd has 1 sibling
        #         nd's parent is the root
        #             nd is a leaf
        #             nd is internal
        #         nd's parent is NOT the root
        #             nd.par is the leftmost child
        #                 nd is a leaf
        #                 nd is internal
        #             nd.par is the rightmost child
        #                 nd is a leaf
        #                 nd is internal
        #             nd.par is a middle child
        # nd is a middle child
        #         nd is a leaf
        #         nd is internal
        #print '########## begin deleteNode(%d) ##########' % nd.number
        #print 'before newick = ',self.makeNewick()
        #raw_input('^^')
        deleting = nd.number
        assert nd.par is not None, 'called deleteNode on root node'
        if nd == nd.par.lchild:
            # nd is leftmost child
            assert nd.rsib is not None, 'expecting lchild to have an rsib'
            if nd.rsib.rsib is not None:
                # nd has 2 or more siblings
                # (nd is leftmost child)
                if nd.lchild is None:
                    # nd is a leaf
                    # (nd has 2 or more siblings)
                    # (nd is leftmost child)
                    #     nd--x--y       x---y
                    #      \\ | /   -->  \\ /
                    #        par          par
                    if debugging:
                        print 'del case 1a (%d)' % deleting
                    x = nd.rsib
                    #y = nd.rsib.rsib
                    nd.par.lchild = x
                    nd.clear()
                else:
                    # nd is internal
                    # (nd has 2 or more siblings)
                    # (nd is leftmost child)
                    #   a---b
                    #    \\/
                    #     nd--x--y       a-b-x-y
                    #      \\ | /   -->  \\| |/
                    #        par           par
                    if debugging:
                        print 'del case 1b (%d)' % deleting
                    a = nd.lchild
                    b = a
                    x = nd.rsib
                    while b.rsib is not None:
                        b.par = nd.par
                        b = b.rsib
                    b.par = nd.par
                    b.rsib = x
                    nd.par.lchild = a
                    nd.clear()
            else:
                # nd has only 1 sibling
                # (nd is leftmost child)
                if nd.par.par is None:
                    # nd's parent is the root
                    # (nd has only 1 sibling)
                    # (nd is leftmost child)
                    if nd.lchild is not None:
                        # nd is a leaf
                        # (nd's parent is the root)
                        # (nd has only 1 sibling)
                        # (nd is leftmost child)
                        #    u-v-w
                        #    \\|/
                        #     nd----y               u-v-w-y
                        #      \\  /                \\| |/
                        #        par = self.root      par
                        if debugging:
                            print 'del case 2a (%d)' % deleting
                        y = nd.rsib
                        u = nd.lchild
                        u.par = nd.par
                        nd.par.lchild = u
                        v = u.rsib
                        w = v
                        while w.rsib is not None:
                            w.par = nd.par
                            w = w.rsib
                        w.par = nd.par
                        w.rsib = y
                        nd.clear()
                    else:
                        # nd is internal
                        # (nd's parent is the root)
                        # (nd has only 1 sibling)
                        # (nd is leftmost child)
                        #     nd---y
                        #      \\ /              --> y = self.root
                        #       par = self.root
                        if debugging:
                            print 'del case 2b (%d)' % deleting
                        y = nd.rsib
                        y.par = None
                        tree.root = y   #bug: should be self.root = y?
                        nd.par.clear()
                        nd.clear()
                else:
                    # nd's parent is NOT the root
                    # (nd has only 1 sibling)
                    # (nd is leftmost child)
                    z = nd.par.par
                    if nd.par == z.lchild:
                        # nd.par is leftmost child
                        # (nd's parent is NOT the root)
                        # (nd has only 1 sibling)
                        # (nd is leftmost child)
                        if nd.lchild is not None:
                            # nd is a leaf
                            # (nd.par is leftmost child)
                            # (nd's parent is NOT the root)
                            # (nd has only 1 sibling)
                            # (nd is leftmost child)
                            #    x--w
                            #    \\ /
                            #     nd---y         x--w--y
                            #      \\ /          \\ | /
                            #       par--u--v  --> par--u--v
                            #         \\ | /         \\ | /
                            #            z              z
                            if debugging:
                                print 'del case 3a (%d)' % deleting
                            x = nd.lchild
                            y = nd.rsib
                            w = x
                            while w.rsib is not None:
                                w.par = nd.par
                                w = w.rsib
                            w.par = nd.par
                            w.rsib = y
                            nd.par.lchild = x
                            nd.clear()
                        else:
                            # nd is internal
                            # (nd.par is leftmost child)
                            # (nd's parent is NOT the root)
                            # (nd has only 1 sibling)
                            # (nd is leftmost child)
                            #     nd---y
                            #      \\ /
                            #       par--u--v  -->   y--u--v
                            #         \\ | /         \\ | /
                            #            z              z
                            if debugging:
                                print 'del case 3b (%d)' % deleting
                            y = nd.rsib
                            u = nd.par.rsib
                            y.rsib = u
                            y.par = z
                            z.lchild = y
                            nd.par.clear()
                            nd.clear()
                    elif nd.par.rsib is None:
                        # nd.par is rightmost child
                        # (nd's parent is NOT the root)
                        # (nd has only 1 sibling)
                        # (nd is leftmost child)
                        if nd.lchild is not None:
                            # nd is a leaf
                            # (nd.par is rightmost child)
                            # (nd's parent is NOT the root)
                            # (nd has only 1 sibling)
                            # (nd is leftmost child)
                            #      w--x
                            #      \\/
                            #       nd---y           w--x--y
                            #        \\ /            \\ | /
                            #   u--v--par  -->   u--v--par
                            #   \\ | /           \\ | /
                            #      z                z
                            if debugging:
                                print 'del case 4a (%d)' % deleting
                            y = nd.rsib
                            w = nd.lchild
                            x = w
                            while x.rsib is not None:
                                x.par = nd.par
                                x = x.rsib
                            x.rsib = y
                            x.par = nd.par
                            nd.par.lchild = w
                            nd.clear()
                        else:
                            # nd is internal
                            # (nd.par is rightmost child)
                            # (nd's parent is NOT the root)
                            # (nd has only 1 sibling)
                            # (nd is leftmost child)
                            #       nd---y
                            #        \\ /
                            #   u--v--par  -->   u--v--y
                            #   \\ | /           \\ | /
                            #      z                z
                            if debugging:
                                print 'del case 4b (%d)' % deleting
                            y = nd.rsib
                            v = z.lchild
                            while v.rsib != nd.par:
                                v = v.rsib
                            v.rsib = y
                            y.par = z
                            nd.par.clear()
                            nd.clear()
                    else:
                        # nd.par is a middle child
                        # (nd's parent is NOT the root)
                        # (nd has only 1 sibling)
                        # (nd is leftmost child)
                        #     nd---y
                        #      \\ /
                        #    u--par--v  -->   u--y--v
                        #     \\ | /          \\ | /
                        #        z               z
                        if debugging:
                            print 'del case 5 (%d)' % deleting
                        y = nd.rsib
                        u = z.lchild
                        v = nd.par.rsib
                        while u.rsib != nd.par:
                            u = u.rsib
                        u.rsib = y
                        y.rsib = v
                        y.par = z
                        nd.par.clear()
                        nd.clear()
        elif nd.rsib is None:
            # nd is rightmost child
            if nd != nd.par.lchild.rsib:
                # nd has 2 or more siblings
                if nd.lchild is None:
                    # nd is a leaf
                    # (nd has 2 or more siblings)
                    # (nd is rightmost child)
                    #      x--y--nd      x---y
                    #      \\ | /   -->   \ /
                    #        par          par
                    if debugging:
                        print '*del case 6a (%d)' % deleting
                    y = nd.par.lchild
                    while y.rsib != nd:
                        y = y.rsib
                    y.rsib = None
                    nd.clear()
                else:
                    # nd is internal
                    # (nd has 2 or more siblings)
                    # (nd is rightmost child)
                    #           a---b
                    #           \\ /
                    #      x--y--nd      x-y-a-b
                    #      \\ | /   -->   \| |/
                    #        par           par
                    if debugging:
                        print '*del case 6b (%d)' % deleting
                    a = nd.lchild
                    a.par = nd.par
                    b = a.rsib
                    y = nd.par.lchild
                    while y.rsib != nd:
                        y = y.rsib
                    y.rsib = a
                    while b is not None:
                        b.par = nd.par
                        b = b.rsib
                    nd.clear()
            else:
                # nd has 1 sibling
                # (nd is rightmost child)
                if nd.par.par is None:
                    # nd's parent is the root
                    # (nd has 1 sibling)
                    # (nd is rightmost child)
                    if nd.lchild is None:
                        # nd is a leaf
                        # (nd's parent is the root)
                        # (nd has 1 sibling)
                        # (nd is rightmost child)
                        #     y---nd
                        #      \\ /
                        #       par = self. root   -->  y = self.root
                        if debugging:
                            print 'del case 7a (%d)' % deleting
                        y = nd.par.lchild
                        y.par = None
                        self.root = y
                        nd.par.clear()
                        nd.clear()
                    else:
                        # nd is internal
                        # (nd's parent is the root)
                        # (nd has 1 sibling)
                        # (nd is rightmost child)
                        #       x----y
                        #        \  /
                        #     z---nd          z--x--y
                        #      \\ /           \\ | /
                        #       par      -->    par
                        if debugging:
                            print 'del case 7b (%d)' % deleting
                        x = nd.lchild
                        y = x.rsib
                        z = nd.par.lchild
                        z.rsib = x
                        x.par = nd.par
                        while y is not None:
                            y.par = nd.par
                            y = y.rsib
                        nd.clear()
                else:
                    # nd's parent is NOT the root
                    # (nd has 1 sibling)
                    # (nd is rightmost child)
                    z = nd.par.par
                    if nd.par == z.lchild:
                        # nd.par is the leftmost child
                        # (nd's parent is NOT the root)
                        # (nd has 1 sibling)
                        # (nd is rightmost child)
                        if nd.lchild is not None:
                            # nd is a leaf
                            # (nd.par is the leftmost child)
                            # (nd's parent is NOT the root)
                            # (nd has 1 sibling)
                            # (nd is rightmost child)
                            #        w---x
                            #        \\ /
                            #     y---nd         y--w--x
                            #      \\ /          \\ | /
                            #       par--u--v  --> par--u--v
                            #         \\ | /         \\ | /
                            #            z              z
                            if debugging:
                                print 'del case 8a (%d)' % deleting
                            y = nd.par.lchild
                            w = nd.lchild
                            x = w
                            while x.rsib is not None:
                                x.par = nd.par
                                x = x.rsib
                            y.rsib = w
                            nd.clear()
                        else:
                            # nd is internal
                            # (nd.par is the leftmost child)
                            # (nd's parent is NOT the root)
                            # (nd has 1 sibling)
                            # (nd is rightmost child)
                            #     y---nd
                            #      \\ /
                            #       par--u--v  -->   y--u--v
                            #         \\ | /         \\ | /
                            #            z              z
                            if debugging:
                                print 'del case 8b (%d)' % deleting
                            y = nd.par.lchild
                            u = nd.par.rsib
                            y.rsib = u
                            y.par = z
                            z.lchild = y
                            nd.par.clear()
                            nd.clear()
                    elif nd.par.rsib is None:
                        # nd.par is the rightmost child
                        # (nd's parent is NOT the root)
                        # (nd has 1 sibling)
                        # (nd is rightmost child)
                        if nd.lchild is None:
                            # nd is a leaf
                            # (nd.par is the rightmost child)
                            # (nd's parent is NOT the root)
                            # (nd has 1 sibling)
                            # (nd is rightmost child)
                            #     u---nd
                            #      \\ /
                            #  a--b-par  -->     a--b--u
                            #  \\ | /            \\ | /
                            #     z                 z
                            if debugging:
                                print 'del case 9a (%d)' % deleting
                            u = nd.par.lchild
                            z = nd.par.par
                            a = z.lchild
                            b = a
                            while b.rsib != nd.par:
                                b = b.rsib
                            b.rsib = u
                            u.rsib = None
                            u.par = z
                            nd.par.clear()
                            nd.clear()
                        else:
                            # nd is internal
                            # (nd.par is the rightmost child)
                            # (nd's parent is NOT the root)
                            # (nd has 1 sibling)
                            # (nd is rightmost child)
                            #      x---y--z
                            #       \\ | /
                            #     u---nd      u--v--x--y--z
                            #      \\ /        \\ |  |  | /
                            #       par  -->       par
                            #       /              /
                            if debugging:
                                print 'del case 9b (%d)' % deleting
                            x = nd.lchild
                            y = x.rsib
                            v = nd.par.lchild
                            while v.rsib != nd:
                                v = v.rsib
                            v.rsib = x
                            x.par = nd.par
                            while y.rsib is not None:
                                y.par = nd.par
                                y = y.rsib
                            nd.clear()
                    else:
                        # nd.par is a middle child
                        # (nd's parent is NOT the root)
                        # (nd has 1 sibling)
                        # (nd is rightmost child)
                        #     y---nd
                        #      \\ /
                        #    u--par--v  -->   u--y--v
                        #     \\ | /          \\ | /
                        #        z               z
                        if debugging:
                            print '*del case 10 (%d)' % deleting
                        y = nd.par.lchild
                        u = z.lchild
                        v = nd.par.rsib
                        while u.rsib != nd.par:
                            u = u.rsib
                        u.rsib = y
                        y.par = z
                        y.rsib = v
                        nd.par.clear()
                        nd.clear()
        else:
            # nd is a middle child
            if nd.lchild is None:
                # nd is a leaf
                # (nd is a middle child)
                #      x--nd--y       x---y
                #      \\  | /   -->   \ /
                #         par          par
                if debugging:
                    print '*del case 11a (%d)' % deleting
                y = nd.rsib
                x = nd.par.lchild
                while x.rsib != nd:
                    x = x.rsib
                x.rsib = y
                nd.clear()
            else:
                # nd is internal
                # (nd is a middle child)
                #        a-b-c
                #        \\|/
                #      x--nd--y       x-a-b-c-y
                #      \\  | /   -->   \| | |/
                #         par            par
                if debugging:
                    print '*del case 11b (%d)' % deleting
                x = nd.par.lchild
                while x.rsib != nd:
                    x = x.rsib
                y = nd.rsib
                a = nd.lchild
                a.par = nd.par
                x.rsib = a
                c = a.rsib
                while c.rsib is not None:
                    c.par = nd.par
                c.par = nd.par
                c.rsib = y
                nd.clear()

        self.calcSplits()
        #print 'after newick = ',self.makeNewick()
        #print '########## end deleteNode ##########'

    def addNodeTo(self, nd, ndnum):
        #  x----y     newnd--x--y
        #  \\  /  -->     \\ | /
        #    nd             nd
        if debugging:
            print 'add case 5'
        newnd = Node()
        newnd.number = ndnum
        newnd.rsib = nd.lchild
        newnd.lchild = None
        newnd.par = nd
        nd.lchild = newnd
        self.calcSplits()
        return newnd

    def addNodeBelow(self, nd, ndnum):
        newnd = Node()
        newnd.number = ndnum
        newpar = Node()
        if not nd.par:
            #   nd       newnd-----nd
            #               \\     /
            #        -->     \\   /
            #                newpar
            if debugging:
                print 'add case 1'
            newpar.lchild = newnd
            newnd.rsib = nd
            newnd.par = newpar
            nd.par = newpar
            tree.root = newpar  #bug: should be self.root = newpar?
        elif nd == nd.par.lchild:
            # nd---x---y       newnd---nd
            #  \\  |  /           \\   /
            #   \\ | /    -->    newpar---x-y
            #    \\|/                \\  / /
            #     par                  par
            if debugging:
                print 'add case 2'
            x = nd.rsib
            newnd.lchild = None
            newnd.rsib = nd
            newnd.par = newpar
            newpar.lchild = newnd
            newpar.rsib = x
            newpar.par = nd.par
            #nd.lchild =
            nd.rsib = None
            nd.par.lchild = newpar
            nd.par = newpar
        elif nd.rsib is None:
            #  x---y---nd                newnd--nd
            #  \\  |  /                     \\  /
            #   \\ | /    -->        x---y--newpar
            #    \\|/                 \\ |   /
            #     par                    par
            if debugging:
                print 'add case 3'
            y = nd.par.lchild
            while y.rsib != nd:
                y = y.rsib
                assert y is not None
            y.rsib = newpar
            newnd.lchild = None
            newnd.rsib = nd
            newnd.par = newpar
            newpar.lchild = newnd
            newpar.rsib = None
            newpar.par = nd.par
            #nd.lchild =
            nd.rsib = None
            nd.par = newpar
        else:
            #  x---nd--y               newnd--nd
            #  \\  |  /                  \\  /
            #   \\ | /    -->          x-newpar-y
            #    \\|/                   \\  |  /
            #     par                      par
            if debugging:
                print 'add case 4'
            x = nd.par.lchild
            while x.rsib != nd:
                x = x.rsib
                assert x is not None
            x.rsib = newpar
            y = nd.rsib
            newnd.lchild = None
            newnd.rsib = nd
            newnd.par = newpar
            newpar.lchild = newnd
            newpar.rsib = y
            newpar.par = nd.par
            #nd.lchild =
            nd.rsib = None
            nd.par = newpar
        self.calcSplits()
        return newnd

class Description(object):
    def __init__(self, t):
        self.tree = t
        self.newick = None
        self.nodeNumbers = None
        self.numTips = None
        self.numInternals = None
        self.treefile = open('trees.tre', 'w')
        self.treefile.write('#nexus\n\nbegin trees;\n')
    def close(self):
        self.treefile.write('end;\n')
        self.treefile.close()
    def reset(self):
        self.nodeNumbers = []
        self.numTips = 0
        self.numInternals = 0
    def processNode(self, nd):
        self.nodeNumbers.append(nd.number)
        if nd.lchild is None:
            self.numTips += 1
        else:
            self.numInternals += 1
    def describe(self, n = None, r = None):
        self.tree.traverse(self)
        self.newick = self.tree.makeNewick()

        if n is not None and r is not None and self.numTips == n and self.numInternals == r:
            self.treefile.write('tree %s = %s;\n' % ('N%d_R%d' % (self.numTips, self.numInternals), self.newick))
        elif n is not None and self.numTips == n:
            self.treefile.write('tree %s = %s;\n' % ('N%d_R%d' % (self.numTips, self.numInternals), self.newick))
        elif r is not None and self.numInternals == r:
            self.treefile.write('tree %s = %s;\n' % ('N%d_R%d' % (self.numTips, self.numInternals), self.newick))
        else:
            self.treefile.write('tree %s = %s;\n' % ('N%d_R%d' % (self.numTips, self.numInternals), self.newick))

        return (self.numTips, self.numInternals, self.newick)

def output(s = None):
    if s is None:
        print
        logf.write('\n')
    else:
        print s
        logf.write('%s\n' % s)

def reduce(t, n, r, binary_r, binary_id, binary_newick, treedict):
    # Take out each edge from t in turn and store the resulting tree's newick in treedict[r][id]['trees'].
    # Also store id in the reduced tree's treedict[r-1][id]['trees'] entry.
    global newick_lookup
    if r == 0:
        return

    # There are r-1 internal nodes with edges
    for i in range(1, r):
        tcopy = t.deepCopy()
        j = 0
        for nd in tcopy.preorder:
            if nd.par is not None and nd.lchild is not None:
                j += 1
            if j == i:
                tcopy.deleteNode(nd)
                tcopy.preorder = Preorder(tcopy).build()
                tcopy_newick = tcopy.makeNewick()
                tcopy_id = frozenset(tcopy.id)
                newick_lookup[tcopy_id] = tcopy_newick
                #output('                       %s' % tcopy_newick)
                treedict[binary_r][binary_id]['trees'].add(tcopy_id)
                treedict[r-1][tcopy_id]['trees'].add(binary_id)
                reduce(tcopy, n, r-1, binary_r, binary_id, binary_newick, treedict)
                break

def recurse(tree, scribe, nmax, create_treedict, n = None):
    global treedict, progress

    if n is None:
        progress = 0

        # start with trees having one more tip than the current tree
        tmp = [nd.number for nd in tree.preorder if nd.lchild is None]
        n = 1 + len(tmp)

        if create_treedict:
            # create dictionary to store trees having nmax tips, with trees in each resolution
            # class stored in a vector with key equal to number of internal nodes
            treedict = {}
            for r in range(1,nmax):
                treedict[r] = {}

    if n <= nmax:
        for nd in tree.preorder:
            newnd = tree.addNodeBelow(nd, n)
            n, r, newick = scribe.describe()
            if n == nmax and r == nmax-1:
                #output('* %6d %6d      %s' % (n, r, newick))
                if create_treedict:
                    treedict[r][frozenset(tree.id)] = {'newick':newick, 'count':0, 'trees':set()}
                else:
                    binary_r      = r
                    binary_id     = frozenset(tree.id)
                    binary_newick = newick
                    newick_lookup[binary_id] = binary_newick
                    progress += 1
                    #print 'Reducing %d: %s' % (progress, binary_newick)
                    reduce(tree, n, r, binary_r, binary_id, binary_newick, treedict)
            elif create_treedict:
                #output('  %6d %6d      %s' % (n, r, newick))
                if n == nmax:
                    treedict[r][frozenset(tree.id)] = {'newick':newick, 'count':0, 'trees':set()}
            recurse(tree, scribe, nmax, create_treedict, n+1)
            tree.deleteNode(newnd)
            if nd.lchild is not None:
                newnd = tree.addNodeTo(nd, n)
                n, r, newick = scribe.describe()
                if create_treedict:
                    #output('  %6d %6d      %s' % (n, r, newick))
                    if n == nmax:
                        treedict[r][frozenset(tree.id)] = {'newick':newick, 'count':0, 'trees':set()}
                recurse(tree, scribe, nmax, create_treedict, n+1)
                tree.deleteNode(newnd)

def listTrees(r, id):
    s = ''
    for treeid in list(treedict[r][id]['trees']):
        newick = newick_lookup[treeid]
        ninternals = len(treeid)
        s += ' | %d: %s' % (ninternals, newick)
        # ' | '.join(['%s' % newick_lookup[treeid] for treeid in list(treedict[r][id]['trees'])])
    return s

def resolutionClassFreq(r, id):
    freqs = [0]*rmax
    for treeid in list(treedict[r][id]['trees']):
        ninternals = len(treeid)
        freqs[ninternals-1] += 1
    s = ' '.join(['%3d' % n for n in freqs])
    return s

def processTreeFile(treefname, tree, descr, treedict, rootat):
    stuff = open(treefname, 'r').read()
    newicks = re.findall('^\s*tree\s+\S+\s*=\s*(.+?);\s*$', stuff, re.M | re.S)
    for newick in newicks:
        tree.buildFromNewick('%s;' % newick)
        tree.rerootAt(rootat)
        n, r, newick = descr.describe()
        treedict[r][frozenset(tree.id)]['count'] += 1

if __name__ == '__main__':
    # key = tree id, value = newick description
    newick_lookup = {}

    logf = open('output.txt', 'w')

    # start with 2 tip tree
    tree = Tree()
    tree.root = Node()
    tree.addNodeTo(tree.root, 2)
    tree.addNodeTo(tree.root, 1)

    #output('  %6s %6s      %s' % ('N', 'R', 'newick'))
    descr = Description(tree)
    n, r, newick = descr.describe()
    #output('  %6d %6d      %s' % (n, r, newick))

    # nmax is the number of taxa in rooted trees
    # nmax determines which file is read (e.g. nmax=6 causes 7taxa-unrooted-trees.t to
    # be processed, and generates 6taxa-rooted-output.txt)
    nmax = 6
    rmax = nmax-1

    print 'Enumerating all possible unrooted trees for %d taxa...' % (nmax+1,)

    # Running recurse to create treedict. Every multifurcating tree topology will be stored
    # in treedict under its resolution class: e.g. treedict[r][id]['newick'] stores the newick
    # string for the tree with the given id in resolution class r, with id being a set of sets
    # in which the outer set comprises r inner sets, each storing the node numbers of leaves
    # above one internal node in the tree. Other keys in treedict[r][id] include 'count',
    # which is a count of trees having this topology, and 'trees', which is the set of all
    # binary (i.e. fully-resolved) trees compatible with this tree. * = binary tree topology.
    create_treedict = True
    recurse(tree, descr, nmax, create_treedict)

    # Initialize nresclass and nsamples dictionaries that will provide a tally of the
    # number of sampled trees representing each topology (nsamples) and falling in each
    # possible resolution class (nresclass)
    nresclass = {}
    nsamples = {}
    for r in range(1,rmax+1):
        nresclass[r] = len(treedict[r].keys())
        nsamples[r] = 0

    print 'Determining "spread" by finding all binary trees compatible with each polytomous tree'
    print '   (this could be done much faster if spread was calculated directly)...'

    # Running recurse to reduce binary trees. Each internal node in each binary
    # tree will be recursively removed and the resulting non-binary tree stored
    # in the 'trees' element of the binary tree. Also, the binary tree will be
    # added to the 'trees' element of the non-binary tree. This allows one to
    # see all non-binary trees that are compatible with each binary tree and
    # all binary trees that can be generated by resolving each non-binary tree.
    create_treedict = False
    recurse(tree, descr, nmax, create_treedict)

    tree_file_name = '../%dtaxa-unrooted-trees.t' % (nmax+1,)
    print 'Processing "%s"...' % tree_file_name

    # Reads trees.t and stores each tree there in the 'counts' element of the
    # tree's treedict record. The term 'spread' is defined as the number of distinct
    # fully-resolved (i.e. binary) labeled trees compatible with a particular tree.
    # Thus, a binary tree has spread 1, while the star tree containing just 1 internal
    # node has spread equal to the total number of binary trees.
    processTreeFile(tree_file_name, tree, descr, treedict, nmax+1)

    descr.close()

    for r in range(1,rmax+1):
        #output()
        #output('The %d %s in resolution class %d:' % (nresclass[r], nresclass[r] != 1 and 'trees' or 'tree', r))
        tally_by_spread = {}
        total_by_spread = {}
        for i,id in enumerate(treedict[r].keys()):
            nsamples[r] += treedict[r][id]['count']
            #output('%6d  %s --> %6d %6d %s %s' % (i+1, treedict[r][id]['newick'], treedict[r][id]['count'], len(treedict[r][id]['trees']), resolutionClassFreq(r, id), listTrees(r, id)))
            count  = treedict[r][id]['count']
            spread = len(treedict[r][id]['trees'])
            if spread in tally_by_spread.keys():
                tally_by_spread[spread] += count
                total_by_spread[spread] += 1
            else:
                tally_by_spread[spread] = count
                total_by_spread[spread] = 1
        output()
        #output('Tally by spread:')
        output('Resolution class having %d internal nodes:' % r)
        spread_total = 0
        distinct_topologies_total = 0
        samples_total = 0
        each_total = 0.0
        for k in tally_by_spread.keys():
            spread_total += k
            distinct_topologies_total += total_by_spread[k]
            samples_total += tally_by_spread[k]
            each = float(tally_by_spread[k])/total_by_spread[k]
            each_total += each
        if each_total > 0.0:
            output('  %20s %20s %20s %20s' % ('spread', 'topologies', 'samples', 'samples/topology'))
            for spread in tally_by_spread.keys():
                kpct = 100.0*spread/spread_total
                each = float(tally_by_spread[spread])/total_by_spread[spread]
                eachpct = 100.0*each/each_total
                #output('  %6d (%6.1f%%) %6d %6d %12.2f (%6.1f%%)' % (spread, kpct, tally_by_spread[spread], total_by_spread[spread], each, eachpct))
                output('  %20d %20d %20d %20.2f' % (spread, total_by_spread[spread], tally_by_spread[spread], each))
            output('  %20s %20d %20s %20s' % ('totals -->', distinct_topologies_total, samples_total, ' '))
        else:
            output('  (no results)')

    output('\n%20s %12s %12s %12s' % ('resolution class', 'topologies', 'samples', 'proportion'))
    total_samples = sum([nsamples[r] for r in range(1,rmax+1)])
    for r in range(1,rmax+1):
        output('%20d %12d %12d %12.3f' % (r, nresclass[r], nsamples[r], float(nsamples[r])/total_samples))

    logf.close()

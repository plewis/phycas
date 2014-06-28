from _PhylogenyExt import *

class Split(SplitBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    Encapulates a split, or taxon bipartition.

    """
    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Initializes data members and sets up the Split for the 4-taxon case.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.setBits((0,1))
        >>> print s.createPatternRepresentation()
        **--

        """
        SplitBase.__init__(self)
        
    def copy(self, other):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Makes this split a copy of other.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---*')
        >>> r = Split()
        >>> r.copy(s)
        >>> print r.createPatternRepresentation()
        -*--*---*

        """
        SplitBase.copy(self, other)

    def reset(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets all bits to 0, but does not change anything else.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---*')
        >>> s.countOnBits()
        3
        >>> s.reset()
        >>> s.countOnBits()
        0
 
        """
        SplitBase.reset(self)

    def setBit(self, b):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets bit b, where 0 <= b < number of taxa. The following example sets
        the first two bits (out of 4, which is the default number of taxa).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.setBit(0)
        >>> s.setBit(1)
        >>> print s.createPatternRepresentation()
        **--

        """
        SplitBase.setBit(self, b)

    def setBits(self, bits_to_set):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets bits in tuple bits_to_set, where each value b in bits_to_set
        obeys 0 <= b < number of taxa. The following example sets the first
        two bits (out of 4, which is the default number of taxa).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.setBits((0,1))
        >>> print s.createPatternRepresentation()
        **--

        """
        SplitBase.setBits(self, bits_to_set)

    def unsetBit(self, taxon):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Unsets (clears) bit b, where 0 <= b < number of taxa. The following
        example sets all 4 bits (4 is is the default number of taxa), then
        clears bits 0 and 3.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.setBit(0)
        >>> s.setBit(1)
        >>> s.setBit(2)
        >>> s.setBit(3)
        >>> print s.createPatternRepresentation()
        ****
        >>> s.unsetBit(0)
        >>> s.unsetBit(3)
        >>> print s.createPatternRepresentation()
        -**-

        """
        SplitBase.unsetBit(self, taxon)

    def unsetBits(self, bits_to_unset):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Clears (unsets) bits in tuple bits_to_unset, where each value b in
        bits_to_unset obeys 0 <= b < number of taxa. The following example
        all 4 bits (4 is is the default number of taxa), then clears bits 0
        and 3.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.setBits((0,1,2,3))
        >>> print s.createPatternRepresentation()
        ****
        >>> s.unsetBits((0,3))
        >>> print s.createPatternRepresentation()
        -**-

        """
        SplitBase.unsetBits(self, bits_to_unset)

    def isBitSet(self, b):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Allows one to query whether a particular bit b is set, where
        0 <= b < number of taxa.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> print s.isBitSet(0)
        False
        >>> s.setBit(0)
        >>> print s.isBitSet(0)
        True

        """
        return SplitBase.isBitSet(self, b)

    def invertSplit(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Clears (unsets) all bits that are currently set, and sets all bits
        that are currently not set. Note: this function does not return a
        value.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---*')
        >>> s.invertSplit()
        >>> print s.createPatternRepresentation()
        *-**-***-

        """
        SplitBase.invertSplit(self)

    def calcComplexity(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the minimum of n and m, where n is the number of taxa on one
        side of the split and m is the number on the other side. Trivial
        splits have m = 1 or n = 1, and thus have compexity 1, whereas the
        most complex split has complexity split_ntax/2 (note that this maximum
        holds whether or not the number of taxa is even or odd).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---*')
        >>> print s.calcComplexity()
        3
        
        """
        return SplitBase.calcComplexity(self)

    def countOnBits(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of bits that are currently set.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---***-*')
        >>> print s.countOnBits()
        6

        """
        return SplitBase.countOnBits(self)

    def countOffBits(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of bits that are currently unset. This function
        will be slower than countOnBits if some bits have been excluded.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---***-*')
        >>> print s.countOffBits()
        7

        """
        return SplitBase.countOffBits(self)

    def createNewickRepresentation(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a newick-style (nested-parenthetical) tree description from
        the split. Such a tree, when displayed, would have two internal nodes
        connected by a single internal edge. Attached to one of the nodes are
        all taxa present on one side of the split, with the remaining taxa
        attached to the other internal node. In the tree description, bits
        that are on are listed first, and are 1-offset (i.e. the first bit
        is represented by 1, not 0).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---****')
        >>> print s.createNewickRepresentation()
        (2,5,9,10,11,12,(1,3,4,6,7,8))

        """
        return SplitBase.createNewickRepresentation(self, False)

    def createPatternRepresentation(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Creates a string representing the split as a series of characters. By
        default, the set bits are represented by '*' and unset bits by '-'
        although this can be changed through the use of the functions
        SetOnSymbol and SetOffSymbol, respectively.
        
        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---****')
        >>> print s.createPatternRepresentation()
        -*--*---****

        """
        return SplitBase.createPatternRepresentation(self)

    def equals(self, other_split):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if this split equals other_split, and False otherwise.
        To be equal, the two splits must have exactly the same pattern of
        set bits.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---****')
        >>> r = Split()
        >>> r.createFromPattern('-*--*---****')
        >>> print r.equals(s)
        True
        >>> t = Split()
        >>> t.createFromPattern('-*--*---***-')
        >>> print t.equals(s)
        False

        """
        return SplitBase.equals(self, other_split)

    def cmp(self, other_split):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns -1 if this split is less than other_split, 0 if this split 
        equals other_split, and 1 if this split is greater than other_split.
        To be equal, the two splits must have exactly the same pattern of
        set bits. One split is less than another if the sum of its component
        bit fields is smaller.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-*--*---****')
        >>> r = Split()
        >>> r.createFromPattern('--*-*---****')
        >>> print r.cmp(s)
        -1
        >>> t = Split()
        >>> t.createFromPattern('*---*---***-')
        >>> print t.cmp(s)
        1
        >>> print s.cmp(s)
        0

        """
        return SplitBase.cmp(self, other_split)

    def isCompatible(self, other_split):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns true if other_split is compatible with this split. The two
        splits a and b are compatible if: (1) the set of "on" bits in a is a
        subset of the set of "on" bits in b; (2) the set of "on" bits in b is
        a subset of the set of "on" bits in a; or (3) the intersection of the
        sets of "on" bits in a and b is the empty set. For example

        split a: -***---*--
        split b: ----***--*

        are compatible, because of reason (3) above. The two splits below are
        also compatible because b is a subset of (i.e. subsumed in) a:

        split a: -****-*---
        split b: --**--*---

        These two splits, on the other hand, are not compatible because the
        intersection is not empty and is also not equal to either a or b:

        split a: -***---*--
        split b: ---***---*

        >>> from phycas.phylogeny import *
        >>> a = Split()
        >>> a.createFromPattern('-***---*--')
        >>> b = Split()
        >>> b.createFromPattern('---***---*')
        >>> print a.isCompatible(b)
        False
        
        """
        return SplitBase.isCompatible(self, other_split)

    def subsumedIn(self, other_split):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns True if the intersection of this split and other_split is
        non-empty and equal to this split. See the documentation for the
        function isCompatible for more information.

        >>> from phycas.phylogeny import *
        >>> a = Split()
        >>> a.createFromPattern('-****-*---')
        >>> b = Split()
        >>> b.createFromPattern('--**--*---')
        >>> print a.subsumedIn(b)
        False
        >>> print b.subsumedIn(a)
        True

        """
        return SplitBase.subsumedIn(self, other_split, 0)

    def createFromPattern(self, pattern_string):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Documentation needs to be written...

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-****-*---')
        >>> print s.createPatternRepresentation()
        -****-*---

        """
        SplitBase.createFromPattern(self, pattern_string)

    def getOnList(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list of bits that are currently set.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-****-*---')
        >>> print s.getOnList()
        [1, 2, 3, 4, 6]

        """
        return list(SplitBase.getOnList(self))

    def getOffList(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list of bits that are currently not set.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-****-*---')
        >>> print s.getOffList()
        [0, 5, 7, 8, 9]

        """
        return list(SplitBase.getOffList(self))

    def getExcludedList(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns a list of bits that are currently excluded (not to be
        considered either on or off).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-x-')
        >>> print s.getExcludedList()
        [4, 5, 8]

        """
        return list(SplitBase.getExcludedList(self))

    def setExcluded(self, excluded_list):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Supplies a list of bits that should be excluded (i.e. not to be
        considered either on or off).

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***--*---')
        >>> s.setExcluded([4, 5, 7, 8])
        >>> print s.createPatternRepresentation()
        -***xx*xx-

        """
        SplitBase.setExcluded(self, excluded_list)

    def getOnSymbol(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the symbol used to represent "on" (i.e. set) bits in functions
        such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> print s.getOnSymbol()
        *

        """
        return SplitBase.getOnSymbol(self)

    def setOnSymbol(self, c):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the symbol used to represent "on" (i.e. set) bits in functions
        such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***--*---')
        >>> print s.createPatternRepresentation()
        -***--*---
        >>> s.setOnSymbol('^')
        >>> print s.createPatternRepresentation()
        -^^^--^---

        """
        SplitBase.setOnSymbol(self, c)

    def getOffSymbol(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the symbol used to represent "off" (i.e. unset) bits in
        functions such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> print s.getOffSymbol()
        -

        """
        return SplitBase.getOffSymbol(self)

    def setOffSymbol(self, c):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the symbol used to represent "off" (i.e. unset) bits in functions
        such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***--*---')
        >>> print s.createPatternRepresentation()
        -***--*---
        >>> s.setOffSymbol('~')
        >>> print s.createPatternRepresentation()
        ~***~~*~~~

        """
        SplitBase.setOffSymbol(self, c)

    def getExcludedSymbol(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the symbol used to represent excluded (i.e. neither set nor
        unset) bits in functions such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> print s.getExcludedSymbol()
        x

        """
        return SplitBase.getExcludedSymbol(self)

    def setExcludedSymbol(self, c):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the symbol used to represent excluded (i.e. neither set nor
        unset) bits in functions such as createPatternRepresentation.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-xx')
        >>> print s.createPatternRepresentation()
        -***xx*-xx
        >>> s.setExcludedSymbol('#')
        >>> print s.createPatternRepresentation()
        -***##*-##

        """
        SplitBase.setExcludedSymbol(self, c)

    def write(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Saves this split as a string. All information needed to restore a
        string is saved. Here are the fields in the order in which they are
        saved. Note that because the first three fields are single
        characters, no delimiter is needed. To extract just the pattern, use
        a slice (see example below).

        o the symbol used to indicate bits that are set
        o the symbol used to indicate bits that are unset
        o the symbol used to indicate bits that are excluded
        o a sequence of characters showing bits that are set, unset
          and excluded (using the three previously-defined symbols)
        
        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-xx')
        >>> split_string = s.write()
        >>> print split_string
        *-x-***xx*-xx
        >>> print split_string[3:]
        -***xx*-xx
        >>> r = Split()
        >>> r.read(split_string)
        >>> print r.createPatternRepresentation()
        -***xx*-xx

        """
        return SplitBase.write(self)

    def readDoof(self, s):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Builds this split from information read from s. The string s should
        have been previously created using the write function. See write
        function for details about how a split is stored as a string.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-xx')
        >>> split_string = s.write()
        >>> print split_string
        *-x-***xx*-xx
        >>> print split_string[3:]
        -***xx*-xx
        >>> r = Split()
        >>> r.read(split_string)
        >>> print r.createPatternRepresentation()
        -***xx*-xx

        """
        SplitBase.read(self, s)

    def getNTaxa(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the number of taxa currently supported by this Split object.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-xx')
        >>> print s.getNTaxa()
        10

        """
        return SplitBase.getNTaxa(self)

    def setNTaxa(self, n):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the number of taxa currently supported by this Split object to n.
        Note that this function clears the object (erases all information
        previously stored in it), even if n equals the number of taxa
        currently supported by the split.

        >>> from phycas.phylogeny import *
        >>> s = Split()
        >>> s.createFromPattern('-***xx*-xx')
        >>> print s.createPatternRepresentation()
        -***xx*-xx
        >>> s.setNTaxa(6)
        >>> print s.createPatternRepresentation()
        ------

        """
        SplitBase.setNTaxa(self, n)


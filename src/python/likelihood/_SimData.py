from _LikelihoodExt import *

class SimData(SimDataBase):
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    """
    This class is a container for simulated data. It stores simulated
    data internally (in the base class SimDataBase) as a pattern map
    (an associative array in which the first element is a pattern and the
    second element is a count of the number of times that pattern was
    generated). A pattern is a vector (list or tuple in Python
    terminology) of integers, each representing the character state of a
    single taxon. Before using a SimData object, call the
    resetPatternLength(ntax) function, where ntax is the number of taxa.
    One side effect of the resetPatternLength function is that it deletes
    all patterns from the pattern map, so only use this function before
    you begin inserting patterns. The resetPatternLength function sets up
    a workspace for building a single pattern, the elements of which can
    be assigned states using the setState() function. After all states
    have been assigned, call insertPattern() to store this pattern in the
    pattern map. Use wipePattern() to fill the temporary pattern with
    invalid values. If a pattern is submitted using insertPattern(), but
    not all states have been assigned, these invalid values will cause the
    insert to fail. To return SimData to its just-constructed state, use
    clear(). Finally, use saveToNexusFile() to save the stored data to a
    file in NEXUS format.

    Two examples are provided here in lieu of individual examples in the
    documentation for each function. The first example is simple: three
    patterns are built up by hand and stored in the SimData object, then
    the internal map representation is printed out.

    >>> from phycas import *
    >>> d = Likelihood.SimData()
    >>> d.resetPatternLength(4) # specify that there are 4 taxa
    >>> 
    >>> # Create the first pattern
    >>> d.setState(0, 0)      # store state 0 for taxon 0
    >>> d.setState(1, 1)      # store state 1 for taxon 1
    >>> d.setState(2, 0)      # store state 0 for taxon 2
    >>> d.setState(3, 1)      # store state 1 for taxon 3
    >>> d.insertPattern(1.0)  # store the pattern
    >>> d.insertPattern(1.0)  # store the pattern again
    >>> 
    >>> # Create the second pattern
    >>> d.wipePattern()  # clear the workspace to start another pattern
    >>> d.setState(0, 0) # store state 0 for taxon 0
    >>> d.setState(3, 1) # store state 1 for taxon 3 (can define states out of order)
    >>> d.setState(2, 2) # store state 2 for taxon 2
    >>> d.setState(1, 3) # store state 3 for taxon 1
    >>> d.insertPattern(1.0)
    >>> 
    >>> # Show the contents of the pattern map (symbols used for states are
    >>> # supplied as a tuple argument)
    >>> print d.patternTable(('a','c','g','t'))
         Count  Pattern
           2.0  acac
           1.0  atgc
    
    #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
    The second example is more complex, but perhaps closer to the common
    SimData use case. It reads a data file, creates a tree, sets up an
    HKY model, computes the likelihood for the tree using the data and
    model, and finishes by simulating a new data set using the transition
    probability matrices established during the calculation of the
    likelihood. Note that it is currently not possible to simulate data
    without first computing the likelihood for a data set already in
    memory. The likelihood calculation is used to prepare the tree for
    simulation. In the future this restriction will be removed, but
    presently the only use for simulation is for posterior predictive
    simulation, and in this case the requirement is automatically
    satisfied.

    Note: the >>> have been removed from the following example so that
    doctest does not see this as an example. The reason for this is that
    simulation has not yet been fully worked out for the case of multiple
    chains (svn > 425). -POL 24 Aug 2007
    
    from phycas import *
    # Read a data file
    phycas.data_matrix =  phycas.readData(getPhycasTestData("nyldna4.nex"))
    phycas.ntax = phycas.data_matrix.getNTax()
    taxon_names = phycas.reader.getTaxLabels()
    
    # Create a tree
    phycas.tree = Phylogeny.Tree()
    model_tree = '(0:0.1,1:0.15,(2:0.025,3:0.15):0.05)'
    phycas.tree.buildFromString(model_tree, True)
    
    # Create a model
    phycas.model = Likelihood.HKYModel()
    phycas.model.setKappa(4.0)
    phycas.model.setNucleotideFreqs(0.1, 0.2, 0.3, 0.4)
    
    # Create a likelihood object to orchestrate everything
    phycas.likelihood = Likelihood.TreeLikelihood(phycas.model)
    phycas.likelihood.copyDataFromDiscreteMatrix(phycas.data_matrix)
    
    # Prepare the tree (e.g. equip it with transition matrices)
    phycas.likelihood.prepareForLikelihood(phycas.tree)
    
    # Simulation setup
    phycas.mcmc_manager.setRandomSeedAllChains(13579)
    phycas.sim_nreps = 1
    phycas.sim_outfile = 'simout.nex'
    num_sites = 100
    
    # Simulate the data
    sim_data = Likelihood.SimData()
    phycas.likelihood.simulateFirst(sim_data, phycas.tree, phycas.r, num_sites)
    
    # Output a table showing the patterns that were simulated and the number of
    # sites exhibiting each pattern
    print sim_data.patternTable('A C G T'.split())
         Count  Pattern
           4.0  AAAA
           1.0  AAAC
           2.0  AAAG
           1.0  ATCA
           1.0  CAAA
          15.0  CCCC
           4.0  CCCT
           1.0  CCTC
           1.0  CGCC
           3.0  CGGG
           1.0  CTCC
           2.0  CTTT
           1.0  GAAA
           2.0  GCGG
           1.0  GGAA
           1.0  GGGC
          19.0  GGGG
           1.0  GGGT
           1.0  GTGA
           1.0  GTGG
           3.0  TCCC
           2.0  TCTC
           4.0  TCTT
           3.0  TGGG
           1.0  TGTT
           1.0  TTCT
           1.0  TTTC
          22.0  TTTT
            
    """

    def __init__(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        The constructor takes no arguments and simply initializes the private
        base class data members.
        
        """
        SimDataBase.__init__(self)
        
    def getNUniquePatterns(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the current number of stored patterns.
        
        """
        return SimDataBase.getNUniquePatterns(self)
        
    def clear(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Returns the SimData object to its just-constructed state. Calling this
        function will delete any data stored in the internal pattern map.
        
        >>> from phycas import *
        >>> d = Likelihood.SimData()
        >>> d.resetPatternLength(4) # specify that there are 4 taxa
        >>> d.setState(0, 0)      # store state 0 for taxon 0
        >>> d.setState(1, 1)      # store state 1 for taxon 1
        >>> d.setState(2, 0)      # store state 0 for taxon 2
        >>> d.setState(3, 1)      # store state 1 for taxon 3
        >>> d.insertPattern(1.0)  # store the pattern
        >>> print d.patternTable('A C G T'.split())
             Count  Pattern
               1.0  ACAC
        >>> d.clear()  # clear everything
        >>> print d.patternTable('A C G T'.split())
        Sorry, no patterns are stored
        
        """
        SimDataBase.clear(self)
        
    def resetPatternLength(self, ntaxa):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Establishes the length of a single pattern and sets up a temporary
        pattern internally that can be filled using ntaxa invocations of the
        function setState. Warning: this function has no effect if the last
        time it was called used the same value for ntaxa, and, more
        importantly, it causes all stored patterns to be quietly deleted if
        ntaxa is different than the value used in the previous invocation.
        Changing the pattern length also has the effect of wiping the
        temporary pattern (all states previously set will be lost).
        
        """
        SimDataBase.resetPatternLength(self, ntaxa)
        
    def wipePattern(self):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Fills the temporary workspace used to build up a new pattern with
        invalid states. These invalid states must all be replaced with valid
        states (positive integers) before insertPattern() will succeed. It is
        never necessary to call this function, but it is good practice in
        order to prevent accidentally inserting a pattern in which one or
        more states are still set to previous values.
        
        """
        SimDataBase.wipePattern(self)
        
    def setState(self, pos, state):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Sets the value in position pos of the temporary pattern workspace to
        state. After calling wipePattern(), you should call the setState
        function ntaxa times (once for every position in the temporary pattern
        workspace) before calling insertPattern() to store the pattern. Note
        that pos is a 0-based index (i.e. 0, 1, ..., ntaxa - 1). The quantity
        ntaxa is the value supplied to the resetPatternLength function the
        last time it was called.
        
        """
        SimDataBase.setState(self, pos, state)

    def insertPattern(self, weight):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        After building up a pattern in the temporary workspace using calls to
        setState, this function stores the pattern in the internal pattern
        map. The value of weight specifies by how much the count for this
        pattern should be incremented (the usual value is 1.0). If
        wipePattern() was used before calling setState, then failure to set a
        valid state for all positions in the temporary pattern workspace will
        result in an error when insertPattern() is called. This behavior is
        designed to catch accidental mistakes, but you are not required to
        call wipePattern() and in fact you may want to avoid this for
        efficiency reasons if new patterns only differ by one position, for
        example.
        
        """
        SimDataBase.insertPattern(self, weight)

    def saveToNexusFile(self, filename, taxon_names, data_type, state_symbols):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function saves the patterns currently stored to a file named
        filename. If the file already exists, it will be overwritten without
        warning. If the file does not yet exist, it will be created. The file
        so written will be a valid NEXUS data file suitable for executing in
        phylogenetic analysis software that reads the NEXUS file format. The
        supplied taxon_names will be used in the matrix command of the NEXUS
        file to specify labels for the taxa. Assumes that the number of
        elements in taxon_names equals the pattern length specified by a
        previous call to the resetPatternLength function. The data_type
        argument should be the correct NEXUS datatype (e.g. "dna", "standard")
        for the data simulated. The symbols used for the states are supplied
        in the state_symbols tuple or list. Each element of this vector
        should be a single-character string. Assumes that no state in any
        stored pattern has a value greater than or equal to the length of
        state_symbols (because states are used as indices into state_symbols.
        
        """
        SimDataBase.saveToNexusFile(self, filename, taxon_names, data_type, state_symbols)

    def saveToNexusFilePython(self, outf, taxon_names, data_type, state_symbols):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        Same as saveToNexusFile except this version is written in pure Python
        so that it can be supplied with a file rather than a file name (file
        objects are hard to pass over to C++). It represents an independent
        implementation of the C++ version.
        
        """
        pattern_length = SimDataBase.getPatternLength(self)
        total_count = int(SimDataBase.getTotalCount(self))
        assert len(state_symbols) > 0
        assert len(taxon_names) == pattern_length

        # Find length of longest string in taxon_names vector; this is used later for formatting purposes
        # The 2 is included in case apostrophes are needed when taxon names are output in the matrix command
        length_of_longest_name = 2 + max([len(nm) for nm in taxon_names])
	
        outf.write('#nexus\n\n')
        outf.write('begin data;\n')
        outf.write('  dimensions ntax=%d nchar=%d;\n' % (pattern_length,total_count))
        outf.write('  format datatype=%s;\n' % data_type)
        outf.write('  matrix\n')

        # Create a format string to use with boost::format that left-justifies (the "-" flag) the 
        # taxon names in a field of width length_of_longest_name 
        fmtstr = "    %%-%ds  " % length_of_longest_name

        for i,nm in enumerate(taxon_names):
            if ' ' in nm:
                s = "'"
                s += nm
                s += "'"
                outf.write(fmtstr % s)
            else:
                outf.write(fmtstr % nm)

            # Spit out characters in the order in which they were simulated. While this is a nice feature, 
            # it currently requires storing the data twice (sim_pattern_vect and sim_pattern_map)
            row = SimDataBase.getPatternVectRow(self, i)
            for c in row:
                outf.write(state_symbols[c])
            outf.write('\n')

        outf.write('  ;\n')
        outf.write('end;\n')

    def patternTable(self, state_symbols):
        #---+----|----+----|----+----|----+----|----+----|----+----|----+----|
        """
        This function returns a string containing a tabular representation of
        the internal pattern map. There are two columns in the table, one
        labeled "Count" and the other labeled "Pattern". The Count column
        shows the number of times its associated pattern was inserted using
        the insertPattern function. The Pattern column shows a representation
        of the pattern itself, using symbols for states provided in the
        state_symbols argument. The state_symbols argument should be a list
        or tuple supplying a single-character string to represent each
        state that might show up in any pattern. For example, any of the
        following could be used as valid state_symbols arguments for DNA
        patterns:

        state_symbols = ['A', 'C', 'G', 'T']
        state_symbols = ('A', 'C', 'G', 'T')
        state_symbols = 'a c g t'.split()
        state_symbols = 'R Y R Y'.split()
        
        """
        return SimDataBase.patternTable(self, state_symbols)

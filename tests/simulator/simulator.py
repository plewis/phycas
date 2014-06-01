from phycas import *

model_tree = '(1:0.1,2:0.1,(3:0.2,4:0.2):0.01)'

# Create a model
model.type =  'hky'
model.kappa = 4.0
model.state_freqs = [.1, .2, .3, .4]
m12 = model()
m3 = model()

partition.addSubset(subset(1,9999,3) + subset(2,9999,3), m12, 'first_second')
partition.addSubset(subset(3,9999,3), m3, 'third')
partition.subset_relrates = [1,10]
partition()

# Define the names of the taxa to use when the simulated data set is saved to a file
sim.taxon_labels = ['P. parksii', 'P. articulata', 'P._gracilis', 'P. macrophylla']

# Create a model tree
sim.tree_source = TreeCollection(newick=Newick(model_tree))

# Simulation settings
sim.random_seed = 13579
sim.file_name = 'simulated.nex'
sim.nchar = 100000  # note: ignored because partition model has more than 1 subset

if False:
    import sys,os
    if os.path.basename(sys.executable) == 'python_d.exe':
        raw_input('attach debugger to python_d process now')

# Simulate data
simulator = sim()

# Now compute the likelihood of the model tree
#like.data_source = 'file'
##like.data_source = 'simulated.nex'
##like.tree_source = TreeCollection(newick=Newick(model_tree))
##lnL = like()
##print 'lnL =',lnL

# Add a PAUP block to the simulated.nex file to make it easy to check the results
f = file('simulated.nex', 'a')
f.write('\n')
f.write('[!\n')
f.write('Relative rates of subsets:\n')
relrates_norm = partition.getNormalizedSubsetRelRates()
f.write('  firstsecond: %g\n' % relrates_norm[0])
f.write('  third:       %g\n' % relrates_norm[1])
f.write(']\n')
f.write('\nbegin paup;')
f.write('\n  set criterion=likelihood storebrlen;')
f.write('\nend;')
f.write('\n')
f.write('\nbegin trees;')
f.write('\n  translate')
for i,nm in enumerate(sim.taxon_labels):
    if nm.count(' ') > 0:
        f.write("\n    %d '%s'" % (i + 1, nm))
    else:
        f.write("\n    %d %s" % (i + 1, nm))
    if i < len(sim.taxon_labels) - 1:
        f.write(',')
    else:
        f.write(';')
f.write('\n  utree simtree = %s;' % simulator.starting_tree)
f.write('\nend;')
f.write('\n')
f.write('\nbegin paup;')
f.write('\n  charpartition codons = firstsecondpos:1-.\\3 2-.\\3, thirdpos:3-.\\3;')
f.write('\n  lset nst=2 variant=hky basefreq=estimate tratio=estimate rates=equal;')
f.write('\n  lscores 1 / nouserbrlen rates=sitespec siterates=partition:codons;')
f.write('\n  describe 1 / brlens;')
f.write('\nend;')
f.write('\n')
f.close()

# This example program simulates data on a 5-taxon tree containing one polytomy
# This data set is then analyzed using an MCMC analysis that allows polytomies
# to demonstrate that the polytomous true tree can be recovered

from phycas import *

setMasterSeed(98765)

print
print '~~~~~~~~~~~~~~~~~~~'
print 'Simulate a data set'
print '~~~~~~~~~~~~~~~~~~~'

# Define the names of the taxa to use when the simulated data set is saved to a file
sim.taxon_labels = ['P._fimbriata', 'P._parksii', 'P._articulata', 'P._gracilis', 'P._macrophylla']

# Create a model tree containing one polytomy
sim.tree_source = TreeCollection(newick='(1:0.1,2:0.1,(3:0.1,4:0.1,5:0.1):0.1)')

# Create a model
sim.model.type = 'hky'
sim.model.kappa = 4.0
sim.model.kappa_prior = Exponential(1.0)
sim.model.state_freq_prior = Dirichlet((1.00000, 1.00000, 1.00000, 1.00000))

sim.nchar = 5000
sim.file_name = 'simHKY.nex'
sim.edgelen_dist = Exponential(10.0)
sim()

# Add a MrBayes block to make it easier to summarize trees later
# A temporary measure: this functionality should be built into Phycas
dataf = file('simHKY.nex', 'a')
dataf.write('\n\nbegin mrbayes;')
dataf.write('\n  sumt file=analHKY.nex nrun=1;')
dataf.write('\nend;\n')
dataf.close()

print
print '~~~~~~~~~~~~~~~~~~~~~~'
print 'HKY analysis beginning'
print '~~~~~~~~~~~~~~~~~~~~~~'
# Set up HKY model
model.type = 'hky'
model.kappa_prior = Exponential(1.0)
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

# Assume no invariable sites
model.pinvar_model = False

# Assume rate homogeneity across sites
model.num_rates = 1

# Use independent exponential priors (mean 0.1) for each edge length parameter
model.edgelen_prior = Exponential(10.0)
model.edgelen_hyperprior = InverseGamma(2.10000, 0.90909)

mcmc.data_source = 'simHKY.nex'

# Conduct a Markov chain Monte Carlo (MCMC) analysis that samples from the posterior distribution
mcmc.ncycles = 500
mcmc.sample_every = 10
mcmc.report_every = 100
#mcmc.starting_tree_source = TreeCollection(newick='(1:.01,2:0.01,(3:0.01,4:0.01):0.01)')
#mcmc.starting_tree_source = TreeCollection(filename='nexustreefile.tre')
mcmc.fix_topology = False
mcmc.allow_polytomies = True
mcmc.polytomy_prior = True
mcmc.topo_prior_C = 2.0
mcmc.bush_move_weight = 50
mcmc.ls_move_weight = 50
mcmc.out.log = 'mcmcoutput.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = 'trees.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = 'params.p'
mcmc.out.params.mode = REPLACE
mcmc()

# Summarize the posterior distribution of model parameters
sump.file = 'params.p'
sump.burnin = 1
sump.out.log.prefix = 'sump-log'
sump.out.log.mode = REPLACE
sump()

# Summarize the posterior distribution of tree topologies and clades
sumt.trees = 'trees.t'
sumt.burnin = 1
sumt.tree_credible_prob = 0.95  # set to 1.0 if you want KL information to be estimated
sumt.save_splits_pdf = True  # if True, may get really large PDF file if tree_credible_prob = 1.0 and number of sampled splits is large
sumt.save_trees_pdf = True   # if True, may get really large PDF file if tree_credible_prob = 1.0 and number of sampled trees is large
sumt.out.log.prefix = 'sumt-log'
sumt.out.log.mode = REPLACE
sumt.out.trees.prefix = 'sumt-trees'
sumt.out.trees.mode = REPLACE
sumt.out.splits.prefix = 'sumt-splits'
sumt.out.splits.mode = REPLACE
sumt()

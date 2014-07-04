from phycas import *

# Set up GTR model
model.type = 'gtr'
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))

model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

model.pinvar_model = True
model.pinvar_prior = Beta(1.0, 1.0)
model.num_rates = 4
model.gamma_shape_prior = Exponential(1.0)

# Use independent exponential priors for each edge length parameter
model.edgelen_prior = Exponential(10.0)
model.edgelen_hyperprior = InverseGamma(2.1, 1.0/1.1)

mcmc.data_source = 'sphaero-reduced.nex'

firsts   = subset(1,5487,3)
seconds  = subset(2,5487,3)
thirds   = subset(3,5487,3)

m1 = model()
m2 = model()
m3 = model()

partition.addSubset(firsts,  m1, 'First codon positions')
partition.addSubset(seconds, m2, 'Second codon positions')
partition.addSubset(thirds,  m3, 'Third codon positions')
partition()

# Compute conditional predictive ordinates (CPO) for each site and for the entire data set
mcmc.ncycles = 201000
mcmc.sample_every = 100
mcmc.report_every = 1000
#mcmc.starting_tree_source = TreeCollection(newick='(1:.01,2:0.01,(3:0.01,4:0.01):0.01)')
mcmc.fix_topology = False
mcmc.out.log = 'mcmcoutput.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = 'trees.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = 'params.p'
mcmc.out.params.mode = REPLACE
cpo.out.sitelike.prefix = 'sitelikes'
cpo()

# Run sump command to summarize CPO values saved during MCMC
sump.file = 'params.p'
sump.burnin = 11
sump.cpofile = 'sitelikes.txt'
sump.out.log.prefix = 'sump-log'
sump.out.log.mode = REPLACE
sump()

# Summarize the posterior distribution of tree topologies and clades
sumt.trees = 'trees.t'
sumt.burnin = 11
sumt.out.log.prefix = 'sumt-log'
sumt.out.log.mode = REPLACE
sumt.out.trees.prefix = 'sumt-trees'
sumt.out.trees.mode = REPLACE
sumt.out.splits.prefix = 'sumt-splits'
sumt.out.splits.mode = REPLACE
sumt()


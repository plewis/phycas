from phycas import *

setMasterSeed(13579)

filename = getPhycasTestData('FRT2000rbcL.nex')
blob = readFile(filename)
#mcmc.data_source = 'sample.nex'
mcmc.data_source = blob.characters

# Set up GTR model
model.type = 'gtr'
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

# Add proportion of invariable sites submodel
model.pinvar_model = True
model.pinvar_prior = Beta(1.0, 1.0)

# Add discrete gamma rate heterogeneity submodel
model.num_rates = 4
model.gamma_shape_prior = Exponential(1.0)

# Use independent exponential priors (mean 0.1) for each edge length parameter
model.edgelen_prior = Exponential(1.0)
model.edgelen_hyperprior = InverseGamma(2.10000, 0.90909)

# Conduct a Markov chain Monte Carlo (MCMC) analysis that samples from the posterior distribution
mcmc.ncycles = 10000
mcmc.burnin = 1000
mcmc.target_accept_rate = 0.3
mcmc.sample_every = 100
mcmc.report_every = 100
#mcmc.starting_tree_source = TreeCollection(newick='(1:.01,2:0.01,(3:0.01,4:0.01):0.01)')
#mcmc.starting_tree_source = TreeCollection(filename='nexustreefile.tre')
mcmc.fix_topology = False
mcmc.allow_polytomies = False
mcmc.bush_move_weight = 0
mcmc.ls_move_weight = 100
mcmc.out.log = 'mcmcoutput.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = 'trees.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = 'params.p'
mcmc.out.params.mode = REPLACE
mcmc()

# Estimate the marginal likelihood using the Generalized Stepping-stone (GSS) method
# Estimate the reference distribution to use with Generalized SS
refdist.skip = 1
refdist.params = 'params.p'
refdist.trees = 'trees.t'
refdist.out.refdistfile = 'refdist.txt'
refdist.out.refdistfile.mode = REPLACE
refdist()

# Do not fix the tree topology for stepping-stone analysis
mcmc.fix_topology = False

# Ucomment and replace mcmc.starting_tree_source below with correct newick
# tree description complete with edge lengths if mcmc.fix_topology is set to True
#mcmc.starting_tree_source = TreeCollection(newick='(1:0.01,2:0.01,(3:0.01,4:0.01):0.01)')

# Choose different output file names to avoid overwriting the ones used for the reference distribution
mcmc.out.log = 'ss-output.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = 'ss-trees.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = 'ss-params.p'
mcmc.out.params.mode = REPLACE
# Set up and run the ss command (this will make use of many mcmc settings)
ss.nstones = 50
ss.ncycles = 1000
ss.sample_every = 10
ss.report_every = 100
ss.refdist_is_prior = False
ss.refdistfile = 'refdist.txt'
ss.shape1 = 1.000000
ss.shape2 = 1.000000
ss()

# Running sump on the param file output by the ss command will calculate the marginal likelihood estimate
sump.file = 'ss-params.p'
sump.skip = 1
sump.out.log.prefix = 'ss-sump-log'
sump.out.log.mode = REPLACE
sump()


# This example is designed to demonstrate that the MCMC code is working correctly by
# exploring the polytomy reference class prior.

from phycas import *

#setMasterSeed(6517649)
setMasterSeed(9467156)

# Set up GTR model
model.type = 'jc'
#model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
#model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

# Add proportion of invariable sites submodel
model.pinvar_model = False
#model.pinvar_prior = Beta(1.0, 1.0)

# Add discrete gamma rate heterogeneity submodel
model.num_rates = 1
#model.gamma_shape_prior = Exponential(1.0)

# Use independent exponential priors (mean 0.1) for each edge length parameter
model.edgelen_prior = Exponential(1.0)
model.edgelen_hyperprior = None #InverseGamma(2.10000, 0.90909)

mcmc.ntax = 7
mcmc.data_source = None

# Conduct a Markov chain Monte Carlo (MCMC) analysis that samples polytomous (as well as fully-resolved) tree topologies
mcmc.ncycles = 10000
mcmc.burnin = 100
mcmc.target_accept_rate = 0.3
mcmc.sample_every = 1
mcmc.report_every = 1000
mcmc.fix_topology = False
mcmc.allow_polytomies = True
mcmc.bush_move_weight = 50
mcmc.ls_move_weight = 50
mcmc.polytomy_prior = False
mcmc.topo_prior_C = 1.0
mcmc.out.log = '7taxa-unrooted-mcmcoutput.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = '7taxa-unrooted-trees.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = '7taxa-unrooted-params.p'
mcmc.out.params.mode = REPLACE
mcmc()

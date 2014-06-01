from phycas import *

setMasterSeed(13579)

model.type = 'hky'
model.num_rates = 1
model.pinvar_model = False

#model.gamma_shape_prior = Exponential(1.0)
model.kappa_prior = BetaPrime(1.0, 1.0)
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

model.tree_length_prior = TreeLengthDist(1.0, 0.1, 20.0, 0.05)
#model.tree_length_prior = TreeLengthDist(1.0, 0.1, 1.0, 1.0)

model.update_freqs_separately = True
model.update_relrates_separately = True

mcmc.out.log.prefix = 'rzy-explore-prior'
mcmc.out.log.mode = REPLACE

mcmc.out.trees.prefix = 'rzy-explore-prior'
mcmc.out.trees.mode = REPLACE

mcmc.out.params.prefix = 'rzy-explore-prior'
mcmc.out.params.mode = REPLACE

mcmc.nchains = 1
mcmc.ncycles = 5000
mcmc.sample_every = 1
mcmc.report_every = 100
mcmc.adapt_first = 2
mcmc.verbose = True
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_tree_source = randomtree(n_taxa=10)
mcmc.data_source = None
mcmc.ntax = 10
mcmc()

#sumt.burnin = 1
#sumt.trees = 'rzy-explore-prior.t'
#sumt.out.splits.prefix = 'splits-rzy-explore-prior'
#sumt.out.splits.mode = REPLACE
#sumt.out.trees.prefix = 'trees-rzy-explore-prior'
#sumt.out.trees.mode = REPLACE
#sumt()
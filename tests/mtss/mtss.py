# Tests stepping-stone when tree topology varies

from phycas import *

setMasterSeed(37597)

model.type = 'jc'
model.kappa_prior = Exponential(1.0)
model.num_rates = 1
#model.gamma_shape_prior = Exponential(1.0)
model.pinvar_model = False
#model.pinvar_prior = Beta(1.0, 1.0)
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

model.edgelen_prior = Exponential(10.0)
model.edgelen_hyperprior = None

prefix = 'ss-multi'
mcmc.starting_tree_source = TreeCollection(newick='(1:0.130765,2:0.001,(3:0.086676,4:0.041752):0.041753)')

mcmc.data_source = 'rokas-4taxa-25sites.nex'
mcmc.fix_topology = False
mcmc.ncycles = 5000
mcmc.sample_every = 10
mcmc.report_every = 1000
mcmc.out.log = prefix + '.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = prefix + '.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = prefix + '.p'
mcmc.out.params.mode = REPLACE
#mcmc()

#refdist.burnin = 1
#refdist.params = prefix + '.p'
#refdist.trees = prefix + '.t'
#refdist.out.refdistfile = prefix + 'ref.txt'
#refdist.out.refdistfile.mode = REPLACE
#refdist()

mcmc.out.log = prefix + 'ss.txt'
mcmc.out.log.mode = REPLACE
mcmc.out.trees = prefix + 'ss.t'
mcmc.out.trees.mode = REPLACE
mcmc.out.params = prefix + 'ss.p'
mcmc.out.params.mode = REPLACE

mcmc.ls_move_weight = 50

ss.nstones = 10
ss.ncycles = 500
ss.sample_every = 1
ss.report_every = 500
ss.refdist_is_prior = False
ss.refdistfile = prefix + 'ref.txt'
ss.shape1 = 1.0
ss.shape2 = 1.0
ss()

sump.burnin = 1
sump.file = prefix + 'ss.p'
sump.out.log = prefix + 'sump.txt'
sump.out.log.mode = REPLACE
sump()

sumt.burnin = 1
sumt.trees = prefix + 'ss.t'
sumt.out.log = prefix + 'sumt.txt'
sumt.out.log.mode = REPLACE
sumt.out.trees = prefix + 'sumt_trees.txt'
sumt.out.trees.mode = REPLACE
sumt.out.splits = prefix + 'sumt_splits.txt'
sumt.out.splits.mode = REPLACE
#sumt()


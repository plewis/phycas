# This file does a short MCMC analysis under the GTR+I+G model

from phycas import *
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

setMasterSeed(13579)

model.type = 'hky'
model.num_rates = 1
model.pinvar_model = False

model.kappa_prior = BetaPrime(1.0, 1.0)
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

model.edgelen_hyperprior = None
model.tree_length_prior = likelihood.TreeLengthDist(1.0, 0.1, 20.0, 0.05)

mcmc.out.log.prefix = 'rzy-ss'
mcmc.out.log.mode = REPLACE

mcmc.out.trees.prefix = 'rzy-ss'
mcmc.out.trees.mode = REPLACE

mcmc.out.params.prefix = 'rzy-ss'
mcmc.out.params.mode = REPLACE

mcmc.nchains = 1
mcmc.ncycles = 500
mcmc.sample_every = 1
mcmc.report_every = 50
mcmc.adapt_first = 2
mcmc.verbose = True
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_tree_source = TreeCollection(newick='(1:0.13390943,2:0.05352159,(3:0.08435058,(4:0.07149961,((5:0.06443949,8:0.12534620):0.03287487,(6:0.07127331,((7:0.04537793,10:0.07253604):0.01807743,9:0.05706492):0.03786300):0.02492948):0.02543656):0.02903737):0.02248610);')
#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
mcmc.data_source = blob.characters
mcmc.fix_topology = True

mcmc()

refdist.params = 'rzy-ss.p'
refdist.trees = 'rzy-ss.t'
refdist.out.refdistfile.prefix = 'rzy-ss-refdist'
refdist.out.refdistfile.mode = REPLACE
refdist()

ss.ncycles = 100
ss.nstones = 11
ss.shape1 = 1.0
ss.shape2 = 1.0
ss.refdist_is_prior = False
ss.refdistfile = 'rzy-ss-refdist.txt'
ss()

sump.file = 'rzy-ss.p'
sump.out.log.prefix = 'rzy-ss-sump'
sump.out.log.mode = REPLACE
sump()



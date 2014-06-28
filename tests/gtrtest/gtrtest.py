# This file does a short MCMC analysis under the GTR+I+G model

from phycas import *
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

rng = probdist.Lot()
rng.setSeed(13579)

model.type = 'gtr'
model.num_rates = 4
model.gamma_shape = 0.5
model.pinvar_model = True

model.edgelen_prior = Exponential(1.0)
#model.separate_edgelen_hyper = False

model.update_freqs_separately = True
model.update_relrates_separately = True

mcmc.out.log = 'output.txt'
mcmc.out.log.mode = REPLACE

mcmc.out.trees.prefix = 'gtr_test'
mcmc.out.trees.mode = REPLACE

mcmc.out.params.prefix = 'gtr_test'
mcmc.out.params.mode = REPLACE

mcmc.nchains = 1
mcmc.ncycles = 20
mcmc.sample_every = 5
mcmc.report_every = 1
mcmc.adapt_first = 2
mcmc.verbose = True
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
#mcmc.starting_tree_source = 'random'
mcmc.rng = rng
mcmc.data_source = blob.characters

mcmc()

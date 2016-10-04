from phycas import *

number_of_samples = 10000
sampling_freq = 1

#if __name__ == '__main__':
#     pwk.out.log      = 'logfile.txt'
#     pwk.out.log.mode = REPLACE
#     #pwk.params       = 'mcmc.run1-short.p'
#     #pwk.trees        = 'mcmc.run1-short.t'
#     #pwk.skip         = 1
#     pwk.params       = 'mcmc.run1.p'
#     pwk.trees        = 'mcmc.run1.t'
#     pwk.skip         = 51113
#     pwk()

# Perform a short MCMC analysis under the GTR+I+G model to generate samples

filename = getPhycasTestData('FRT2000rbcL.nex')
blob = readFile(filename)

#rng = probdist.Lot()
#rng.setSeed(13579)
setMasterSeed(91357)

model.type = 'gtr'
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

model.num_rates = 4
model.gamma_shape = 0.5
model.gamma_shape_prior = Exponential(1.0)

model.pinvar_model = True
model.pinvar_prior = Beta(1.0, 1.0)

model.state_freqs = [.1,.2,.3,.4]

model.edgelen_prior = Exponential(1.0)
model.edgelen_hyperprior = InverseGamma(2.10000, 0.90909)

mcmc.out.log = 'pwk_test.txt'
mcmc.out.log.mode = REPLACE

mcmc.out.trees.prefix = 'pwk_test'
mcmc.out.trees.mode = REPLACE

mcmc.out.params.prefix = 'pwk_test'
mcmc.out.params.mode = REPLACE

mcmc.nchains = 1
mcmc.target_accept_rate = 0.3
mcmc.ncycles = number_of_samples*sampling_freq
mcmc.burnin = number_of_samples*sampling_freq//100
mcmc.sample_every = sampling_freq
mcmc.report_every = number_of_samples*sampling_freq//20
mcmc.report_efficiency_every = number_of_samples*sampling_freq//20
mcmc.verbose = True
mcmc.ls_move_weight = 100
mcmc.tree_scaler_weight = 1
mcmc.slice_weight = 1
mcmc.slice_max_units = 0
mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels))
#mcmc.starting_tree_source = randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
#mcmc.starting_tree_source = 'random'
#mcmc.rng = rng
mcmc.data_source = blob.characters

# import wingdbstub
#
# if 'WINGDB_ACTIVE' in os.environ:
#     print 'Success starting debugger'
# else:
#     print 'Failed to start debugger'

mcmc()

pwk.params         = 'pwk_test.p'
pwk.trees          = 'pwk_test.t'
pwk.skip           = 1
#pwk.minsample      = 5000
#pwk.minsample      = 1525
#pwk.minsample      = 1500
#pwk.minsample      = 300
#pwk.minsample      = 211
#pwk.minsample      = 200
pwk.minsample      = 100
pwk.shells         = 50
pwk.reach          = 1.0
pwk.out.log.prefix = 'pwklog'
pwk.out.log.mode   = REPLACE
#raw_input('..attach now..')
pwk()


from phycas import *
filename = getPhycasTestData('nyldna4.nex')
blob = readFile(filename)

rng = ProbDist.Lot()
rng.setSeed(13579)

model.type                         = 'hky'
model.kappa                        = 4.0
model.fix_kappa                    = True
model.kappa_prior                  = ProbDist.Exponential(1.0)

model.state_freqs                  = [0.1, 0.2, 0.3, 0.4]
model.fix_freqs                    = True
model.state_freq_param_prior       = ProbDist.Exponential(1.0)

model.num_rates                    = 4
model.gamma_shape                  = 0.14
model.fix_shape                    = True
model.gamma_shape_prior            = ProbDist.Exponential(1.0)

model.pinvar_model                 = True
model.pinvar                       = 0.27
model.fix_pinvar                   = True
model.pinvar_prior                 = ProbDist.Beta(1.0, 1.0)

model.fix_edgelens                 = True
model.edgelen_prior                = ProbDist.Exponential(10.0)
model.edgelen_hyperparam           = 0.05
model.fix_edgelen_hyperparam       = True
model.edgelen_hyperprior           = ProbDist.InverseGamma(2.1, 0.9090909)
#model.separate_edgelen_hyper       = False


mcmc.out.log                      = 'output.txt'
mcmc.out.log.mode                 = REPLACE
mcmc.out.trees                    = 'fixed.t'
mcmc.out.trees.mode               = REPLACE
mcmc.out.params                   = 'fixed.p'
mcmc.out.params.mode              = REPLACE
mcmc.nchains                      = 1
mcmc.ncycles                      = 2500
mcmc.rng                          = rng
mcmc.data_source                   = blob.characters
mcmc.starting_tree_source         = randomtree(n_taxa=len(blob.taxon_labels), rng=rng)

if False:
    import sys,os
    if os.path.basename(sys.executable) == 'python_d.exe':
        raw_input('debug stop')

mcmc()

# This example is designed to demonstrate that the MCMC code is working correctly by
# exploring only the priors. The posterior resulting from such an analysis should
# approximate the joint prior distribution.

from phycas import *

setMasterSeed(463859)

model.num_rates                 = 4
model.gamma_shape_prior         = Exponential(10.0)

model.edgelen_hyperprior        = None
model.edgelen_prior             = Exponential(1.0)

m1 = model()

model.pinvar_model              = True
model.pinvar_prior              = Beta(0.5, 1.5)
model.gamma_shape_prior         = Gamma(1.0, 1.0)
model.gamma_shape_prior.setMeanAndVariance(0.5, 2.0)

m2 = model()

partition.addSubset(subset(1,300),    m1, 'small')
partition.addSubset(subset(301,1000), m2, 'large')
partition()

mcmc.data_source                = None
mcmc.ntax                       = 10

mcmc.starting_tree_source       = randomtree(n_taxa=10)

mcmc.out.log.prefix             = 'nodata.nex'
mcmc.out.log.mode               = REPLACE

mcmc.out.trees.prefix           = 'nodata.nex'
mcmc.out.trees.mode             = REPLACE

mcmc.out.params.prefix          = 'nodata.nex'
mcmc.out.params.mode            = REPLACE

mcmc.nchains                    = 1
mcmc.ncycles                    = 1000
mcmc.burnin                     = 100
mcmc.sample_every               = 20
mcmc.report_every               = mcmc.ncycles//100
mcmc.report_efficiency_every    = mcmc.ncycles//10
mcmc.verbose                    = True
mcmc.ls_move_weight             = 100
mcmc.slice_weight               = 1

mcmc()

# This file does a short MCMC analysis under the codon model

from phycas import *
filename = getPhycasTestData('green.nex')
blob = readFile(filename)

rng = probdist.Lot()
rng.setSeed(13579)

model.type                      = 'codon'
model.num_rates                 = 1
model.pinvar_model              = False
model.edgelen_hyperprior        = None
model.edgelen_prior             = Exponential(1.0)
model.state_freqs               = [1.0/61.0]*61
model.state_freq_prior          = Dirichlet([1.0]*61)
model.state_freq_param_prior    = Exponential(1.0)

model.update_freqs_separately   = False

mcmc.nchains                    = 1
mcmc.ncycles                    = 20
mcmc.sample_every               = 10
mcmc.report_every               = 1
mcmc.adapt_first                = 2
mcmc.verbose                    = True

mcmc.ls_move_weight             = 0
mcmc.edge_move_weight           = 0
mcmc.tree_scaler_weight         = 0
mcmc.state_freq_weight          = 1
mcmc.slice_weight               = 1

mcmc.fix_topology               = True

mcmc.slice_max_units            = 0
mcmc.starting_tree_source       = TreeCollection(newick='(8:0.56880388,(((3:0.40888265,(1:1.03799510,2:0.41917430):0.03417782):0.16416599,4:0.29333306):0.14865078,(6:0.28599164,((7:0.14870266,10:0.32973086):0.06151508,9:0.24129778):0.17828009):0.11396143):0.15762955,5:0.29601916);') # randomtree(n_taxa=len(blob.taxon_labels), rng=rng)
mcmc.rng                        = rng
mcmc.data_source                = blob.characters

mcmc.state_freq_psi             = 3000.0     # max_psi
mcmc.state_freq_psi0            = 2.0        # min_psi

mcmc.out.log                    = None

mcmc.out.log                    = 'codontest.txt'
mcmc.out.log.mode               = REPLACE

mcmc.out.params.prefix          = 'params'
mcmc.out.params.mode            = REPLACE

mcmc.out.trees.prefix           = 'trees'
mcmc.out.trees.mode             = REPLACE

if False:
    import sys,os
    if os.path.basename(sys.executable) == 'python_d.exe':
        raw_input('debug stop')

mcmc()

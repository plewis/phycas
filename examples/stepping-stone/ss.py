from phycas import *

setMasterSeed(13579)

fnprefix = 'marshall-codon'
partition_scheme = 'codon'  # 'unpartitioned', 'gene', 'codon', or 'genecodon'

# HKY+G maximum likelihood tree
user_tree_def = '(1:0.00156293,(((((((2:0.00325591,16:0.00402707):0.02205143,((20:0.00630496,29:0.00695889):0.01138659,28:0.01818355):0.00809410):0.01335766,(26:0.01188256,27:0.01250127):0.02960036):0.00967768,((((((3:0.00860238,15:0.01056681):0.00220042,4:0.01193905):0.00408679,5:0.01787676):0.00155750,((17:0.01386172,18:0.01057494):0.00315644,19:0.01190134):0.00143905):0.00324977,30:0.01220535):0.01044754,22:0.01935610):0.02073536):0.00355082,(((7:0.07257248,(25:0.07364494,(31:0.12962943,32:0.09693177):0.08359470):0.00544096):0.03392672,(((((8:0.01400035,(10:0.00393309,12:0.00684284):0.00452834):0.02117048,9:0.02142141):0.00187281,14:0.01577341):0.00403985,11:0.02125429):0.00374345,21:0.05267665):0.01590936):0.00295028,13:0.03355061):0.00389731):0.00448117,23:0.02567495):0.01624640,6:0.01444287):0.01588953,24:0.00385810);'

# Here are the charset definitions
COI           = subset(   1,  774)
COIfirst      = subset(   1,  774, 3)
COIsecond     = subset(   2,  774, 3)
COIthird      = subset(   3,  774, 3)
COII          = subset( 775, 1476)
COIIfirst     = subset( 775, 1476, 3)
COIIsecond    = subset( 776, 1476, 3)
COIIthird     = subset( 777, 1476, 3)
ATPase8       = subset(1477, 1627)
ATPase8first  = subset(1477, 1627, 3)
ATPase8second = subset(1478, 1627, 3)
ATPase8third  = subset(1479, 1627, 3)
ATPase6       = subset(1628, 2090)
ATPase6first  = subset(1629, 2090, 3)
ATPase6second = subset(1630, 2090, 3)
ATPase6third  = subset(1628, 2090, 3)
firsts        = COIfirst  + COIIfirst  + ATPase8first  + ATPase6first
seconds       = COIsecond + COIIsecond + ATPase8second + ATPase6second
thirds        = COIthird  + COIIthird  + ATPase8third  + ATPase6third

# For starting values, use MLEs based on unpartitioned GTR+G analysis of user_tree_def
relrate_MLEs    = [4.28396, 34.10179,  1.39027,  2.53545, 40.60048,  1.00000]     # see note 1 below
state_freq_MLEs = [0.376859, 0.096002, 0.094516, 0.432623]
shape_MLE       = 0.146955
pinvar_MLE      = 0.621227

# model type (jc, hky or gtr)
model.type = 'gtr'

# GTR relative rates
model.update_relrates_separately = False
model.relrates = [x/sum(relrate_MLEs) for x in relrate_MLEs]     # see note 1 below
model.relrate_prior = Dirichlet((1.0, 1.0, 1.0, 1.0, 1.0, 1.0))

# state frequencies
model.update_freqs_separately = False
model.state_freqs = state_freq_MLEs
model.state_freq_prior = Dirichlet((1.0, 1.0, 1.0, 1.0))

# discrete gamma rate heterogeneity
model.num_rates = 4
model.gamma_shape = shape_MLE
model.gamma_shape_prior = Exponential(1.0)

# proportion of invariable sites
model.pinvar_model = False

# edge lengths
model.edgelen_hyperprior = None
model.edgelen_prior = Exponential(40.0)
#model.internal_edgelen_prior = Exponential(40.0)
#model.external_edgelen_prior = Exponential(40.0)

if partition_scheme == 'unpartitioned':
    # unpartitioned
    partition()
elif partition_scheme == 'gene':
    # partition by gene
    m1 = model()
    m2 = model()
    m3 = model()
    m4 = model()
    partition.addSubset(COI,     m1, 'COI')
    partition.addSubset(COII,    m2, 'COII')
    partition.addSubset(ATPase8, m3, 'ATPase8')
    partition.addSubset(ATPase6, m4, 'ATPase6')
    partition()
elif partition_scheme == 'codon':
    # partition by codon
    m1 = model()
    m2 = model()
    m3 = model()
    partition.addSubset(firsts,  m1, 'First codon positions')
    partition.addSubset(seconds, m2, 'Second codon positions')
    partition.addSubset(thirds,  m3, 'Third codon positions')
    partition()
elif partition_scheme == 'genecodon':
    # partition by gene and codon position
    m1  = model()
    m2  = model()
    m3  = model()
    m4  = model()
    m5  = model()
    m6  = model()
    m7  = model()
    m8  = model()
    m9  = model()
    m10 = model()
    m11 = model()
    m12 = model()
    partition.addSubset(COIfirst,       m1, 'COIfirst')
    partition.addSubset(COIsecond,      m2, 'COIsecond')
    partition.addSubset(COIthird,       m3, 'COIthird')
    partition.addSubset(COIIfirst,      m4, 'COIIfirst')
    partition.addSubset(COIIsecond,     m5, 'COIIsecond')
    partition.addSubset(COIIthird,      m6, 'COIIthird')
    partition.addSubset(ATPase8first,   m7, 'ATPase8first')
    partition.addSubset(ATPase8second,  m8, 'ATPase8second')
    partition.addSubset(ATPase8third,   m9, 'ATPase8third')
    partition.addSubset(ATPase6first,  m10, 'ATPase6first')
    partition.addSubset(ATPase6second, m11, 'ATPase6second')
    partition.addSubset(ATPase6third,  m12, 'ATPase6third')
    partition()
else:
    print 'Sorry, unrecognized scheme (%s)' % partition_scheme
    sys.exit()

# read in the data
mcmc.data_source = 'marshall.nex'

# set output file names
mcmcprefix = fnprefix+'-mcmc'
mcmc.out.log.prefix = mcmcprefix
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = mcmcprefix
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = mcmcprefix
mcmc.out.params.mode = REPLACE

mcmc.edge_move_weight = 1
mcmc.ls_move_weight = 100

mcmc.state_freq_weight = 10
mcmc.state_freq_psi  = 1000.0

mcmc.rel_rate_weight = 10
mcmc.rel_rate_psi = 1000.0

mcmc.subset_relrates_weight = 10
mcmc.subset_relrates_psi  = 1000.0

mcmc.tree_scaler_weight = 1

mcmc.slice_weight = 1

mcmc.starting_tree_source = TreeCollection(newick=user_tree_def)
mcmc.fix_topology = True
mcmc.sample_every = 1
mcmc.report_every = 100

mcmc.burnin = 500
mcmc.ncycles = 10000
mcmc.sample_every = 10
mcmc.report_every = 100
mcmc.report_efficiency_every = 1000
mcmc()

refdist.burnin = 1
refdist.params = mcmcprefix+'.p'
refdist.trees = mcmcprefix+'.t'
refdist.out.refdistfile = 'refdist.txt'
refdist.out.refdistfile.mode = REPLACE
refdist()

ssprefix = fnprefix+'-ss'
mcmc.out.log.prefix = ssprefix
mcmc.out.log.mode = REPLACE
mcmc.out.trees.prefix = ssprefix
mcmc.out.trees.mode = REPLACE
mcmc.out.params.prefix = ssprefix
mcmc.out.params.mode = REPLACE

# steppingstone sampling method
ss.burnin = 100
ss.ncycles = 1000
ss.sample_every = 1
ss.report_every = 100
ss.nstones = 25
ss.shape1 = 1.0
ss.shape2 = 1.0
ss.refdist_is_prior = False
ss.refdistfile = 'refdist.txt'
ss()

sumpprefix = fnprefix+'-sump'
sump.file         = ssprefix + '.p'
sump.out.log      = sumpprefix + '.sump.txt'
sump.out.log.mode = REPLACE
sump()


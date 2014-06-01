import math
from phycas import *

outf = open('output.txt', 'w')

# Calculate likelihood under conditions in which underflow correction is used
# for a tree without polytomies

model.type = 'gtr'
model.pinvar_model = True
model.edgelen_hyperprior = None
model.state_freqs = [0.339271, 0.154491, 0.134649, 0.371589]
model.relrates    = [1.144048, 5.419204, 0.454958, 1.766404, 5.546350, 1.0]
model.gamma_shape = 0.906291 
model.pinvar      = 0.442154
model.num_rates   = 4

filename = getPhycasTestData('rbcL50.nex')
blob = readFile(filename)
#nchar = blob.characters.getMatrix().getNChar()
#partition.validate(nchar)

like.data_source = blob.characters
like.tree_source = TreeCollection(filename='gtrig.rbcL50.best.tre')
like.starting_edgelen_dist = None
like.uf_num_edges = 5
like.store_site_likes = True
lnL = like()

outf.write('Computing likelihood for non-polytomous tree:\n')
outf.write('  correct lnL = -18117.830737\n')
outf.write('          lnL = %.6f\n\n' % lnL)

#print 'Computing likelihood for non-polytomous tree:'
#print '  correct lnL = -18117.830737'
#print '          lnL = %.6f' % lnL
#print
    
if 0:
    pattern_log_likes = like.getSiteLikes()
    pattern_counts = like.getPatternCounts()
    char_to_pattern = like.getCharToPattern()
    npatterns = len(pattern_log_likes)
    assert npatterns == len(pattern_counts), 'pattern_counts is not the same length as pattern_log_likes'
    assert nchar == len(char_to_pattern), 'char_to_pattern should have length %d' % nchar
    
    patternf = open('patterns.txt','w')
    for site in range(nchar):
        pattern = char_to_pattern[site]
        site_lnL = pattern_log_likes[pattern]
        pattern_count = pattern_counts[pattern]
        patternf.write('%10d %20.6f\n' % (site+1,site_lnL))
    patternf.close()

# Calculate likelihood under conditions in which underflow correction is used
# for a tree that has polytomies

partition.resetPartition()

model.type      = 'hky'   # use the Hasegawa-Kishino-Yano (1985) model
model.pinvar_model = False
model.num_rates = 4        # add discrete gamma rate heterogeneity with 4 rate categories
model.state_freqs = [0.25, 0.25, 0.25, 0.25]    # replace with actual value 
model.kappa = 4.0                               # replace with actual value 
model.gamma_shape = 0.5                         # replace with actual value 

filename = getPhycasTestData('ShoupLewis.nex')
blob = readFile(filename)
#nchar = blob.characters.getMatrix().getNChar()
#partition.validate(nchar)

like.data_source = blob.characters
like.tree_source = TreeCollection(filename='polytomous.tre')
like.starting_edgelen_dist = None
like.uf_num_edges = 5
like.store_site_likes = True
lnL = like()

outf.write('Computing likelihood for tree with polytomies:\n')
outf.write('  correct lnL = -16403.967004359674\n')
outf.write('          lnL = %.6f\n\n' % lnL)

outf.close()

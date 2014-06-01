tree_lengths = []
lines = open('rzy-explore-prior.p', 'r').readlines()
for line in lines[3:]:
    parts = line.split()
    Gen,lnL,lnPrior,TL,kappa,freqC,freqG,freqT = parts
    tree_lengths.append(float(TL))
print 'No. tree lengths sampled =',len(tree_lengths)
outf = open('plot.R', 'w')
outf.write('tl = c(%s)\n' % ','.join(['%g' % tl for tl in tree_lengths]))
outf.write('x = seq(0,30,.01)\n')
outf.write('y = dgamma(x, shape=1, scale=10)\n')
outf.write('plot(x,y,type="l",lwd=2,xlim=c(0,30))\n')
outf.write('lines(density(tl), lty="dotted")\n')
outf.close()

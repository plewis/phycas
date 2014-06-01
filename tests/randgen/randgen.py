from phycas import *

rng = setMasterSeed(12345)

outf = open('first100.txt','w')
for i in range(100):
    #print rng.uniform()
    outf.write('%.10f\n' % rng.uniform())
outf.close()


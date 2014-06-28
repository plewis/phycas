from phycas.phylogeny import *

outf = file('out.txt', 'w')

z = '-**---****-----******------*'
outf.write('s = %s\n' % z)

s = Split();
s.createFromPattern(z)
outf.write('s = %s (after Split construction)\n' % s.createPatternRepresentation())

s.setExcluded([1,3,5,7,9])
outf.write('s = %s (after excluding 1,3,5,7,9)\n' % s.createPatternRepresentation())

on_list = s.getOnList()
outf.write('list of bits set:')
for onbit in on_list:
    outf.write(' %d' % onbit)
outf.write('\n')

off_list = s.getOffList()
outf.write('list of bits unset:')
for offbit in off_list:
    outf.write(' %d' % offbit)
outf.write('\n')

excl_list = s.getExcludedList()
outf.write('list of excluded bits:')
for exclbit in excl_list:
    outf.write(' %d' % exclbit)
outf.write('\n')

s.setOnSymbol('^')
s.setOffSymbol('~')
s.setExcludedSymbol('!')
outf.write('s = %s (on = %s, off = %s, excluded = %s)\n' % (s.createPatternRepresentation(), s.getOnSymbol(), s.getOffSymbol(), s.getExcludedSymbol()))

s.setOnSymbol('*')
s.setOffSymbol('-')
s.setExcludedSymbol('x')
outf.write('s = %s (on = %s, off = %s, excluded = %s)\n' % (s.createPatternRepresentation(), s.getOnSymbol(), s.getOffSymbol(), s.getExcludedSymbol()))

num_on_bits = s.countOnBits()
outf.write('s = %s (number of set bits is %d)\n' % (s.createPatternRepresentation(), num_on_bits))
num_off_bits = s.countOffBits()
outf.write('s = %s (number of unset bits is %d)\n' % (s.createPatternRepresentation(), num_off_bits))
r = Split()
r.copy(s)
outf.write('r = %s (copy of split s before inverting s)\n' % r.createPatternRepresentation())
s.invertSplit()
num_on_bits_inv = s.countOnBits()
outf.write('s = %s (after inverting, number of set bits now %d)\n' % (s.createPatternRepresentation(), num_on_bits_inv))
complexity = s.calcComplexity()
outf.write('Complexity is %d (minimum of %d and %d)\n' % (complexity, num_on_bits, num_off_bits))
outf.write('First bit of r %s set\n' % (r.isBitSet(0) and 'is' or 'is not'))
outf.write('First bit of s %s set\n' % (s.isBitSet(0) and 'is' or 'is not'))
outf.write('r %s equal to s\n' % (r.equals(s) and 'is' or 'is not'))
t = Split()
t.copy(r)
t.unsetBit(2)
outf.write('t = %s (copy of split r but with bit 2 unset)\n' % t.createPatternRepresentation())
r.setBit(4)
outf.write('r = %s (with bit 4 set)\n' % r.createPatternRepresentation())
outf.write('t %s subsumed in r\n' % (t.subsumedIn(r) and 'is' or 'is not'))
outf.write('r %s subsumed in t\n' % (r.subsumedIn(t) and 'is' or 'is not'))
outf.write('t %s compatible with r\n' % (t.isCompatible(r) and 'is' or 'is not'))
outf.write('r %s compatible with t\n' % (r.isCompatible(t) and 'is' or 'is not'))
outf.write('t in newick tree format: %s\n' % t.createNewickRepresentation())
t.unsetBits((15,16,17,18,19,20))
outf.write('t = %s (after unsetting bits 15 through 20)\n' % t.createPatternRepresentation())
t.setBits((21,22,23,24,25,26))
outf.write('t = %s (after setting bits 21 through 26)\n' % t.createPatternRepresentation())
tsaved = t.write()
outf.write('t = %s (string representation)\n' % tsaved)
outf.write('t = %s (minus the first 3 characters, which are the on, off and excluded symbols)\n' % tsaved[3:])
t.reset()
outf.write('t = %s (after calling reset)\n' % t.createPatternRepresentation())
t.read(tsaved)
outf.write('t = %s (after restoring from saved string)\n' % t.createPatternRepresentation())
outf.write('total number of taxa represented by t = %d' % t.getNTaxa())
t.setNTaxa(10)
outf.write('t = %s (after setting number of taxa to 10, which forces split object to be cleared)\n' % t.createPatternRepresentation())
outf.close()


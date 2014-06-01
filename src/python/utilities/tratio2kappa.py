# The tratio2kappa.py script converts the transition/transversion ratio (tratio)
# (used by PAUP*) into the transition/transversion rate ratio (kappa). The quantity 
# tratio is the expected number of all transitions divided by the expected number of all
# transversions in an infinitesimal time period, whereas kappa is the ratio of the
# instantaneous rate of transitions to the instantaneous rate of transversions.
#
# This script assumes the HKY substitution model, the instantaneous rate matrix
# of which looks like this:
#
#          A           C           G           T
#  A      ---        beta     beta*kappa     beta
#  C     beta         ---        beta     beta*kappa
#  G  beta*kappa     beta         ---        beta
#  T     beta     beta*kappa     beta         ---
#
# The expected number of all transitions per instant of time is:
#     2*beta*kappa*(piA*piG + piC*piT)
# The expected number of all transversions per instant of time is:
#     2*beta*((piA + piG)*(piC + piT))
# The ratio of these quantities is tratio:
#
#          kappa*(piA*piG + piC*piT)
# tratio = -------------------------
#           (piA + piG)*(piC + piT)
#
# kappa can be determined given tratio and the base frequencies:
#
#          tratio*(piA + piG)*(piC + piT)
# kappa = -------------------------------
#               (piA*piG + piC*piT)
# 
# Usage:
#   python tratio2kappa.py ratio piA piC piG piT
#
# Example:
#   python tratio2kappa.py 1.31 0.290 0.205 0.158 0.347
#   kappa = 2.76993510324
#
import sys,math

if __name__ == '__main__':
    if len(sys.argv) != 6:
        print 'Error: you must specify the value of tratio and all four'
        print 'nucleotide relative frequencies when invoking this script'
        print 'Usage:'
        print '  python tratio2kappa.py tratio piA piC piG piT'
        print 'Example:'
        print '  python tratio2kappa.py 1.31 0.290 0.205 0.158 0.347'
        print '  kappa = 2.76993510324'

    tratio = float(sys.argv[1])
    piA = float(sys.argv[2])
    piC = float(sys.argv[3])
    piG = float(sys.argv[4])
    piT = float(sys.argv[5])

    # Normalize frequencies if the ones supplied do not add to 1.0 +- 1.e-6
    sum = piA + piC + piG + piT
    if math.fabs(sum - 1.0) > 1.e-6:
        print 'Warning: base frequencies should sum to %f (should sum to 1.0)' % sum
        print 'The following normalized frequencies will be used instead:'
        piA = piA/sum
        piC = piC/sum
        piG = piG/sum
        piT = piT/sum
        print '  piA =',piA
        print '  piC =',piC
        print '  piG =',piG
        print '  piT =',piT

    kappa = tratio*(piA + piG)*(piC + piT)/(piA*piG + piC*piT)
    print 'kappa =',kappa

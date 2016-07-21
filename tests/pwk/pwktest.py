from phycas import *

if __name__ == '__main__':
    pwk.out.log      = 'logfile.txt'
    pwk.out.log.mode = REPLACE
    pwk.params       = 'mcmc.run1.p'
    pwk.trees        = 'mcmc.run1.t'
    pwk.skip         = 11113
    pwk()


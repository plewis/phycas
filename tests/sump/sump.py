from phycas import *

if __name__ == '__main__':
    sump.out.log      = 'logfile.txt'
    sump.out.log.mode = REPLACE
    sump.file         = 'ss-jc4-params.p'
    sump.skip         = 1
    sump()


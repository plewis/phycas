import sys
from phycas import OutputFilter, getDefaultOutFilter, probdist

output_stream = None
outputter = None

def getDefaultOutputContainer():
    return getDefaultOutputter()

def writeMessageToStdOut(msg):
    "Writes a message to standard out and adds a trailing newline"
    sys.stdout.write("%s" % msg)

def getDefaultOutputStream():
    return writeMessageToStdOut

def getDefaultOutputter():
    global outputter
    if outputter is None:
        outputter = OutputFilter(getDefaultOutFilter(), getDefaultOutputStream())
    return outputter

class CommonFunctions(object):
    def __init__(self, opts):
        self.log_file_spec = None
        self.quiet = False
        self.opts = opts
        try:
            o = self.opts.out
            self.optsout = o
            self.stdout = o.getStdOutputter()
            self.open()
        except:
            self.optsout = None
            self.stdout = getDefaultOutputter()
        self._rng = None
    def __del__(self):
        self.close()

    def open(self):
        if (self.log_file_spec is None):
            try:
                self.log_file_spec = self.optsout.log
            except:
                pass
            if self.log_file_spec:
                self.log_file_spec.openAsLog(self.stdout)
    def close(self):
        if self.log_file_spec:
            self.log_file_spec.close()
            self.log_file_spec = None
        
    def phycassert(self, assumption, msg):
        self.stdout.phycassert(assumption, msg)

    def output(self, msg=""):
        self.stdout.info(msg)
    info = output
            
    def abort(self, msg):
        self.stdout.abort(msg)
    def warning(self, msg):
        self.stdout.warning(msg)
    def error(self, msg):
        self.stdout.error(msg)
    def add_mirror(self, m):
        self.stdout.add_mirror(m)
    def remove_mirror(self, m):
        self.stdout.remove_mirror(m)
    def _getLot(self):
        l = self.opts.rng
        if l:
            self._rng = l
        elif self._rng is None:
            self._rng = probdist.Lot()
        try:            
            i = int(self.opts.random_seed)
            if i != 0:
                self._rng.setSeed(i)
                self.opts.random_seed = 0
        except:
            pass
        return self._rng


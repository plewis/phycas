#/usr/bin/env python
import os
from subprocess import Popen, PIPE

def runPhycasUpdateChecker(outstream, update_url, branch_string, revision_string):
    phycas_path = os.path.dirname(os.path.abspath(os.path.dirname(__file__)))
    phycas_path
    dot_svn_path = os.path.join(phycas_path, ".svn")
    
    underSVN = os.path.exists(dot_svn_path)
    if underSVN:
        try:
            svnStatOut = "" #Popen(["svn", "status", "-u", phycas_path], stdout=PIPE).communicate()[0]
        except:
            outstream.warning('Could not run "svn status" command on directory %d to see if it is up-to-date' % phycas_path)
            return
        if not "*" in svnStatOut:
            outstream.verbose_info("Your copy of phycas is up-to-date")
            return
        
    print "runPhycasUpdateChecker() not implemented"

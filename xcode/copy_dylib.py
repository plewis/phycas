#!/usr/bin/python

import shutil,os
f = open('doof-xcode.txt', 'w')
for k in os.environ.keys():
    f.write('%s --> %s\n' % (k, os.environ[k]))
f.close()

product     = os.environ['PRODUCT_NAME']
phycas_root = os.environ['PHYCAS_ROOT']
build_dir   = os.environ['TARGET_BUILD_DIR']
fromfile    = os.path.join(build_dir,product)+'.so'
tofile      = os.path.join(phycas_root, 'phycas', product[1:-3], product+'.so')
print "copying fromfile=",fromfile," to tofile=",tofile
shutil.copyfile(fromfile, tofile)

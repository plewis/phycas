#!/bin/bash

# Fixes the dylib install names and dependencies when creating a distribution for MacOS
# that will be installed in the standard Python 2.7 location /Library/Python/2.7/site-packages
#
# Note: this script expects...
# 1. existence of /usr/local/lib/ncl/libncl-2.1.18.dylib
# 2. existence of ../phycas
# 3. this script is <phycas dir>/scripts/fixit.sh
#
# If libncl-2.1.18.dylib does not yet exist, download NCL from http://sourceforge.net/projects/ncl
# and configure/make/install it to the standard location
#
# If ../phycas does not yet exist, build it using ../build.sh
#
# See https://blogs.oracle.com/dipol/entry/dynamic_libraries_rpath_and_mac

# Specify location of libncl-2.1.18.dylib
NCL_DYLIB_DIR=/Users/plewis/Documents/libraries/ncl/lib/ncl

# Specify location of libpython2.7.dylib
# PYTHON_DYLIB_DIR=/Users/plewis/Documents/libraries/pydbg/lib
PYTHON_DYLIB_DIR=/System/Library/Frameworks/Python.framework/Versions/2.7/lib

# Make sure only root can run this script
if [[ $EUID -ne 0 ]]; then
  echo "This script must be run as root" 1>&2
  exit 1
fi

# Save the starting directory so that we can exit in the same place
CWD=`pwd`

# Collect all dylibs in a newly-created lib directory
cd ../phycas
mkdir lib
cd conversions
mv *.dylib ../lib
#
cd ../lib
cp $NCL_DYLIB_DIR/libncl-2.1.18.dylib .
ln -s libncl-2.1.18.dylib libncl.dylib
cp $PYTHON_DYLIB_DIR/libpython2.7.dylib .

# Set the name used by all dylibs to be relative to the binary (so file) that calls them
# Thus, if
#   <python root>/Lib/site-packages/phycas/datamatrix/_DataMatrixExt.so
# needs to load libncl.dylib, @loader_path will be replaced with
#   <python root>/Lib/site-packages/phycas/datamatrix
# and the library will be looked for in
#   <python root>/Lib/site-packages/phycas/lib
# For this to work, the id of libncl.dylib must match, so in this section we set the id
# of each dylib to the path that will be used by _DataMatrixExt.so when it tries to load
# the dylib.
install_name_tool -id "@loader_path/../lib/libncl.dylib" libncl.dylib
install_name_tool -id "@loader_path/../lib/libpython2.7.dylib" libpython2.7.dylib
install_name_tool -id "@loader_path/../lib/libboost_thread.dylib" libboost_thread.dylib
install_name_tool -id "@loader_path/../lib/libboost_python.dylib" libboost_python.dylib
install_name_tool -id "@loader_path/../lib/libboost_chrono.dylib" libboost_chrono.dylib
install_name_tool -id "@loader_path/../lib/libboost_system.dylib" libboost_system.dylib

# libboost_chrono.dylib references libncl.dylib and libboost_system.dylib
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" libboost_chrono.dylib
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" libboost_chrono.dylib

# libboost_system.dylib references libncl.dylib
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" libboost_system.dylib

# libboost_thread.dylib references libncl.dylib and libboost_system
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" libboost_thread.dylib
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" libboost_thread.dylib

# libboost_python.dylib references libpython2.7.dylib
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" libboost_python.dylib

cd ../conversions
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _ConversionsExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _ConversionsExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _ConversionsExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _ConversionsExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _ConversionsExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _ConversionsExt.so

cd ../datamatrix
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _DataMatrixExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _DataMatrixExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _DataMatrixExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _DataMatrixExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _DataMatrixExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _DataMatrixExt.so

cd ../likelihood
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _LikelihoodExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _LikelihoodExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _LikelihoodExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _LikelihoodExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _LikelihoodExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _LikelihoodExt.so

cd ../phylogeny
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _PhylogenyExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _PhylogenyExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _PhylogenyExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _PhylogenyExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _PhylogenyExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _PhylogenyExt.so

cd ../probdist
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _ProbDistExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _ProbDistExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _ProbDistExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _ProbDistExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _ProbDistExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _ProbDistExt.so

cd ../readnexus
install_name_tool -change "$NCL_DYLIB_DIR/libncl-2.1.18.dylib" "@loader_path/../lib/libncl.dylib" _ReadNexusExt.so
install_name_tool -change "$PYTHON_DYLIB_DIR/libpython2.7.dylib" "@loader_path/../lib/libpython2.7.dylib" _ReadNexusExt.so
install_name_tool -change "libboost_thread.dylib" "@loader_path/../lib/libboost_thread.dylib" _ReadNexusExt.so
install_name_tool -change "libboost_python.dylib" "@loader_path/../lib/libboost_python.dylib" _ReadNexusExt.so
install_name_tool -change "libboost_chrono.dylib" "@loader_path/../lib/libboost_chrono.dylib" _ReadNexusExt.so
install_name_tool -change "libboost_system.dylib" "@loader_path/../lib/libboost_system.dylib" _ReadNexusExt.so

cd $CWD

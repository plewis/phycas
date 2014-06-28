#!/bin/bash

# This script expects certain environmental variables (BOOST_ROOT, PYTHON_ROOT,
# and NCL_INSTALL_DIR) to be defined before it is run. Examples of the necessary
# definitions are commented out below; modify appropriately for your installation.

# Specify "linux" or "clang" (for macos mavericks)
export OSTYPE="clang"

# The path to your boost installation (download from http://www.boost.org/)
export BOOST_ROOT="$HOME/Documents/libraries/boost_1_55_0"

# The path to your python installation (download from https://www.python.org/)
export PYTHON_ROOT="$HOME/Documents/libraries/pydbg/bin/python2.7"

# Modify PATH so that bjam executable can be found. This assumes bjam is just
# inside $BOOST_ROOT, which will be the case if bootstrap.sh was run: e.g.
#   cd $BOOST_ROOT
#   ./bootstrap.sh --with-python=$PYTHON_ROOT --with-toolset=clang --with-libraries=python
export PATH="$BOOST_ROOT:${PATH}"

# This is needed for bjam to find its way
export BOOST_BUILD_PATH="$BOOST_ROOT/tools/build/v2"

# Provide path to preinstalled Nexus Class Library
# (download from http://sourceforge.net/projects/ncl/)
export NCL_ALREADY_INSTALLED=1
export NCL_INSTALL_DIR="$HOME/Documents/libraries/ncl"

# This removes dynamic link libraries from previous builds
scripts/remove_dylibs.sh

# This command initiates the build
if [ $# -eq 1 ] ; then
    echo Running bjam $1...
    bjam -q $1
else
    echo Running bjam release...
    bjam -q release
fi

scripts/copy_src_python.sh


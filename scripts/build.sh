#!/bin/bash

# This script expects certain environmental variables (BOOST_ROOT, PYTHON_ROOT,
# and NCL_INSTALL_DIR) to be defined before it is run. Examples of the necessary
# definitions are commented out below; modify appropriately for your installation.

# Specify "linux" or "clang" (for MacOS Mavericks) or "windows" (for Windows)
# e.g. export OSTYPE="linux"
export OSTYPE=

# The path to your boost installation (download from http://www.boost.org/)
# e.g. export BOOST_ROOT="$HOME/Documents/libraries/modular-boost"
export BOOST_ROOT=

# The path to your python installation (download from https://www.python.org/)
# e.g. export PYTHON_ROOT="/usr/bin/python"
export PYTHON_ROOT=

# Provide path to preinstalled Nexus Class Library (http://sourceforge.net/projects/ncl/)
# The directory specified should have subdirectories include/ncl and lib/ncl
# e.g. export NCL_INSTALL_DIR="$HOME/Documents/libraries/ncl"
export NCL_INSTALL_DIR=

######## You should not need to change anything below here #######

if [ ! -e "scripts/fixit.sh" ]; then
    echo You seem to be running this script from the wrong directory
    echo Please copy $0 from scripts directory to its parent directory and run from there
    exit 1
fi

if [ ! -n "$OSTYPE" ] || ( [ "$OSTYPE" != "linux" ] && [ "$OSTYPE" != "clang" ] && [ "$OSTYPE" != "windows" ] ); then
    echo Please edit $0 and specify a valid string for OSTYPE
    exit 1
fi

if [ ! -n "$BOOST_ROOT" ] || [ ! -e "$BOOST_ROOT" ]; then
    echo Please edit $0 and specify a valid directory for BOOST_ROOT
    exit 1
fi

if [ ! -n "$PYTHON_ROOT" ] || [ ! -e "$PYTHON_ROOT" ]; then
    echo Please edit $0 and specify a valid directory for PYTHON_ROOT
    exit 1
fi

if [ ! -n "$BOOST_ROOT" ] || [ ! -e "$BOOST_ROOT" ]; then
    echo Please edit $0 and specify a valid directory for BOOST_ROOT
    exit 1
fi

if [ ! -n "$NCL_INSTALL_DIR" ] || [ ! -e "$NCL_INSTALL_DIR" ]; then
    echo Please edit $0 and specify a valid directory for NCL_INSTALL_DIR
    exit 1
fi

# Modify PATH so that bjam executable can be found. This assumes bjam is just
# inside $BOOST_ROOT, which will be the case if bootstrap.sh was run: e.g.
#   cd $BOOST_ROOT
#   ./bootstrap.sh
export PATH="$BOOST_ROOT:${PATH}"

# This is needed for bjam to find its way
export BOOST_BUILD_PATH="$BOOST_ROOT/tools/build/v2"

# Assuming NCL already installed (do not change)
export NCL_ALREADY_INSTALLED=1

# This removes the phycas directory created in a previous build
rm -rf phycas

# This command initiates the build
if [ $# -eq 1 ] ; then
    echo Running bjam $1...
    bjam -q $1
else
    echo Running bjam release...
    bjam -q release
fi

# bjam examples:
# ./bjam -j2 variant=release link=shared install

scripts/copy_src_python.sh


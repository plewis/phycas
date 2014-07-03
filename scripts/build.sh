#!/bin/bash

# This script expects certain environmental variables (BOOST_ROOT, PYTHON_ROOT,
# and NCL_INSTALL_DIR) to be defined before it is run. Examples of the necessary
# definitions are commented out below; modify appropriately for your installation.

# Specify "linux" or "clang" (for MacOS Mavericks) or "windows" (for Windows)
# e.g. export OSTYPE="linux"
export OSTYPE="PLEASE_SUPPLY"

# The path to the root of your boost installation (download from http://www.boost.org/)
# e.g. export BOOST_ROOT="$HOME/boost_1_55_0"
export BOOST_ROOT="PLEASE_SUPPLY"

# The path to your b2 (i.e bjam) executable
# e.g. export BJAM_DIR="$BOOST_ROOT/tools/build/v2/engine/bin.linuxx86_64"
export BJAM_DIR="PLEASE_SUPPLY"

# The path to your python interpreter (download from https://www.python.org/)
# e.g. export PYTHON_ROOT="/opt/python/bin/python"
export PYTHON_ROOT="PLEASE_SUPPLY"

# Provide path to preinstalled Nexus Class Library (http://sourceforge.net/projects/ncl/)
# The directory specified should have subdirectories include/ncl and lib/ncl
# e.g. export NCL_INSTALL_DIR="$HOME/nclib"
export NCL_INSTALL_DIR="PLEASE_SUPPLY"

######## You should not need to change anything below here #######

if [ ! -e "scripts/fixit.sh" ]; then
    echo You seem to be running this script from the wrong directory
    echo Please copy build.sh from scripts directory to its parent directory and run from there
    exit 1
fi

if [ "$OSTYPE" == "PLEASE_SUPPLY" ] || ( [ "$OSTYPE" != "linux" ] && [ "$OSTYPE" != "clang" ] && [ "$OSTYPE" != "windows" ] ); then
    echo Please edit $0 and specify a valid string for OSTYPE \(either \"linux\", \"clang\" or \"windows\"\)
    exit 1
fi

if [ "$BOOST_ROOT" == "PLEASE_SUPPLY" ] || [ ! -e "$BOOST_ROOT" ]; then
    echo Please edit $0 and specify a valid directory for BOOST_ROOT
    exit 1
fi

if [ "$BJAM_DIR" == "PLEASE_SUPPLY" ] || [ ! -e "$BJAM_DIR/bjam" ]; then
    echo Please edit $0 and specify a valid directory for BJAM_DIR
    exit 1
fi

if [ "$PYTHON_ROOT" == "PLEASE_SUPPLY" ] || ( [ ! -e "$PYTHON_ROOT" ] ||  [ ! -f "$PYTHON_ROOT" ]); then
    echo Please edit $0 and specify a valid path to the Python interpreter for PYTHON_ROOT \(you should be able to start Python if you execute the supplied value\)
    exit 1
fi

if [ "$NCL_INSTALL_DIR" == "PLEASE_SUPPLY" ] || [ ! -e "$NCL_INSTALL_DIR" ]; then
    echo Please edit $0 and specify a valid directory for NCL_INSTALL_DIR
    exit 1
fi

# Modify PATH so that bjam executable can be found
export PATH="$BJAM_DIR:${PATH}"

# This is needed for bjam to find boost-build.jam
export BOOST_BUILD_PATH="$BOOST_ROOT"

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
# bjam -d2 -q release
# bjam -j2 variant=release link=shared install

scripts/copy_src_python.sh

#!/bin/bash

# Xcode 6 places all *.so files in $TARGET_BUILD_DIR, which is a subdirectory
# of ~/Library/Developer/Xcode/DerivedData and thus completely outside the staging
# directory that will be tarred and gzipped for distribution. On the other hand, bjam
# creates the staging directory "phycas" at the same level as the Jamroot file. It places,
# e.g., _ConversionsExt.so inside of the "phycas/conversions" directory, _ProbDistExt.so
# inside the "phycas/probdist" directory, etc., so in this case the staging directory must
# be completed rather than completely constructed. These differences are the reason why
# the environmental variable RUNNING_FROM_XCODE is necessary.
#
# This script should be run from directory containing Jamroot (which should have a "src"
# subdirectory containing the "python" and "cpp" subdirectories), and expects either 0
# or 1 command line arguments. Usage:
#
# Case 1: This script is run from a "Run Script" Build Phase within Xcode 6. The single
#   command line argument in this case can be anything (it is the number of command line
#   arguments that the script uses, not their identity):
#
#   copy_src_python.sh xcode
#
# Case 2: This script is run from the build.sh script. In this case it is assumed that
#   both the built products and staging directories are the same, namely the "phycas"
#   directory at the same level as Jamroot:
#
#   copy_src_python.sh
#

if [[ -z "$1" ]]; then
    # No command line arguments provided
    # Assume script is being run by build.sh from the directory containing Jamroot
    RUNNING_FROM_XCODE=0
    PHYCAS_STAGE="phycas"
else
    # At least one command line argument provided
    # Assume script is being run by Xcode from the xcode project directory
    cd ..
    RUNNING_FROM_XCODE=1
    BUILT_PRODUCTS_DIR="$TARGET_BUILD_DIR"
    PHYCAS_STAGE="$PHYCAS_ROOT/phycas"
fi

echo "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
echo "$0 starting..."
echo "RUNNING_FROM_XCODE --> $RUNNING_FROM_XCODE"
#echo "BUILT_PRODUCTS_DIR --> $BUILT_PRODUCTS_DIR"
#echo "PHYCAS_STAGE       --> $PHYCAS_STAGE"

##########################################################
# Check that we are being run from the correct directory #
##########################################################
if [[ ! -e "src/python" ]]; then
    echo "$0 expected to find a directory named "src/python" inside `pwd`"
    exit 1
fi

######################################################################
# Create staging directory (any existing contents will be deleted)   #
# This is not necessary if using bjam because bjam creates it for us #
######################################################################
if [[ $RUNNING_FROM_XCODE -eq 1 && -e "$PHYCAS_STAGE" ]]; then
    echo "Deleting existing staging directory, $PHYCAS_STAGE"
    rm -rf $PHYCAS_STAGE
    mkdir $PHYCAS_STAGE
fi

################################
# Copy master __init__.py file #
################################
cp src/python/__init__.py $PHYCAS_STAGE

########################################
# Copy files for conversions extension #
########################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/conversions
    cp $BUILT_PRODUCTS_DIR/_ConversionsExt.so $PHYCAS_STAGE/conversions
    cp src/python/conversions/*.py $PHYCAS_STAGE/conversions
    cp $BOOST_ROOT/lib/libboost_python.dylib $PHYCAS_STAGE/conversions
    cp $NCL_ROOT/lib/ncl/libncl-2.1.18.dylib $PHYCAS_STAGE/conversions
    cp $PYTHON_ROOT/lib/libpython2.7.dylib $PHYCAS_STAGE/conversions
else
    if [[ ! -e "phycas/conversions" ]]; then
        echo $0 expected to find a directory named "phycas/conversions" inside `pwd`
        exit 1
    fi

    rm -f phycas/conversions/*.py
    cp src/python/conversions/*.py phycas/conversions
fi

#######################################
# Copy files for datamatrix extension #
#######################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/datamatrix
    cp $BUILT_PRODUCTS_DIR/_DataMatrixExt.so $PHYCAS_STAGE/datamatrix
    cp src/python/datamatrix/*.py $PHYCAS_STAGE/datamatrix
else
    if [[ ! -e "phycas/datamatrix" ]]; then
        echo $0 expected to find a directory named "phycas/datamatrix" inside `pwd`
        exit 1
    fi

    rm -f phycas/datamatrix/*.py
    cp src/python/datamatrix/*.py phycas/datamatrix
fi

#######################################
# Copy files for likelihood extension #
#######################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/likelihood
    cp $BUILT_PRODUCTS_DIR/_LikelihoodExt.so $PHYCAS_STAGE/likelihood
    cp src/python/likelihood/*.py $PHYCAS_STAGE/likelihood
else
    if [[ ! -e "phycas/likelihood" ]]; then
        echo $0 expected to find a directory named "phycas/likelihood" inside `pwd`
        exit 1
    fi

    rm -f phycas/likelihood/*.py
    cp src/python/likelihood/*.py phycas/likelihood
fi

######################################
# Copy files for phylogeny extension #
######################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/phylogeny
    cp $BUILT_PRODUCTS_DIR/_PhylogenyExt.so $PHYCAS_STAGE/phylogeny
    cp src/python/phylogeny/*.py $PHYCAS_STAGE/phylogeny
else
    if [[ ! -e "phycas/phylogeny" ]]; then
        echo $0 expected to find a directory named "phycas/phylogeny" inside `pwd`
        exit 1
    fi

    rm -f phycas/phylogeny/*.py
    cp src/python/phylogeny/*.py phycas/phylogeny
fi

#####################################
# Copy files for probdist extension #
#####################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/probdist
    cp $BUILT_PRODUCTS_DIR/_ProbDistExt.so $PHYCAS_STAGE/probdist
    cp src/python/probdist/*.py $PHYCAS_STAGE/probdist
else
    if [[ ! -e "phycas/probdist" ]]; then
        echo $0 expected to find a directory named "phycas/probdist" inside `pwd`
        exit 1
    fi

    rm -f phycas/probdist/*.py
    cp src/python/probdist/*.py phycas/probdist
fi

######################################
# Copy files for readnexus extension #
######################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/readnexus
    cp $BUILT_PRODUCTS_DIR/_ReadNexusExt.so $PHYCAS_STAGE/readnexus
    cp src/python/readnexus/*.py $PHYCAS_STAGE/readnexus
else
    if [[ ! -e "phycas/readnexus" ]]; then
        echo $0 expected to find a directory named "phycas/readnexus" inside `pwd`
        exit 1
    fi

    rm -f phycas/readnexus/*.py
    cp src/python/readnexus/*.py phycas/readnexus
fi

################################
# Copy files for pdfgen module #
################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/pdfgen
    cp src/python/pdfgen/*.py $PHYCAS_STAGE/pdfgen
    cp -R src/python/pdfgen/AFM $PHYCAS_STAGE/pdfgen
else
    if [[ ! -e "phycas/pdfgen" ]]; then
        cd phycas
        rm -rf pdfgen
        mkdir pdfgen
        if [ $? -ne 0 ]; then
            echo $0 could not create pdfgen directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp src/python/pdfgen/*.py phycas/pdfgen
    cp -R src/python/pdfgen/AFM phycas/pdfgen
fi

####################################
# Copy files for treeviewer module #
####################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/treeviewer
    cp src/python/treeviewer/*.py $PHYCAS_STAGE/treeviewer
else
    if [[ ! -e "phycas/treeviewer" ]]; then
        cd phycas
        rm -rf treeviewer
        mkdir treeviewer
        if [ $? -ne 0 ]; then
            echo $0 could not create treeviewer directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp src/python/treeviewer/*.py phycas/treeviewer
fi

###################################
# Copy files for utilities module #
###################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/utilities
    cp src/python/utilities/*.py $PHYCAS_STAGE/utilities
else
    if [[ ! -e "phycas/utilities" ]]; then
        cd phycas
        rm -rf utilities
        mkdir utilities
        if [ $? -ne 0 ]; then
            echo $0 could not create utilities directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp src/python/utilities/*.py phycas/utilities
fi

##################################
# Copy files for commands module #
##################################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/commands
    cp src/python/commands/*.py $PHYCAS_STAGE/commands
else
    if [[ ! -e "phycas/commands" ]]; then
        cd phycas
        rm -rf commands
        mkdir commands
        if [ $? -ne 0 ]; then
            echo $0 could not create commands directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp src/python/commands/*.py phycas/commands
fi

##############
# Copy tests #
##############
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/tests
    cp -R tests/* $PHYCAS_STAGE/tests
    cp scripts/wingdbstub.py $PHYCAS_STAGE/tests/gtrtest
else
    if [[ ! -e "phycas/tests" ]]; then
        cd phycas
        rm -rf tests
        mkdir tests
        if [ $? -ne 0 ]; then
            echo $0 could not create tests directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp -R tests/* phycas/tests
fi

#################
# Copy examples #
#################
if [[ $RUNNING_FROM_XCODE -eq 1 ]]; then
    mkdir $PHYCAS_STAGE/examples
    cp -R examples/* $PHYCAS_STAGE/examples
else
    if [[ ! -e "phycas/examples" ]]; then
        cd phycas
        rm -rf examples
        mkdir examples
        if [ $? -ne 0 ]; then
            echo $0 could not create examples directory inside `pwd`
            exit 1
        fi
        cd ..
    fi

    cp -R examples/* phycas/examples
fi

echo "$0 done."

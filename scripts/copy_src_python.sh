#!/bin/bash

################################
# Copy master __init__.py file #
################################

if [ ! -e "phycas" ]; then
    echo $0 expected to find a directory named "phycas" inside `pwd`
    exit 1
fi

cd phycas
cp ../src/python/__init__.py .

###############################################
# Copy python files for conversions extension #
###############################################

if [ ! -e "conversions" ]; then
    echo $0 expected to find a directory named "conversions" inside `pwd`
    exit 1
fi

cd conversions
rm -f *.py
cp ../../src/python/conversions/*.py .
cd ..

##############################################
# Copy python files for datamatrix extension #
##############################################

if [ ! -e "datamatrix" ]; then
    echo $0 expected to find a directory named "datamatrix" inside `pwd`
    exit 1
fi

cd datamatrix
rm -f *.py
cp ../../src/python/datamatrix/*.py .
cd ..

##############################################
# Copy python files for likelihood extension #
##############################################

if [ ! -e "likelihood" ]; then
    echo $0 expected to find a directory named "likelihood" inside `pwd`
    exit 1
fi

cd likelihood
rm -f *.py
cp ../../src/python/likelihood/*.py .
cd ..

##############################################
# Copy python files for phylogeny extension #
##############################################

if [ ! -e "phylogeny" ]; then
    echo $0 expected to find a directory named "phylogeny" inside `pwd`
    exit 1
fi

cd phylogeny
rm -f *.py
cp ../../src/python/phylogeny/*.py .
cd ..

############################################
# Copy python files for probdist extension #
############################################

if [ ! -e "probdist" ]; then
    echo $0 expected to find a directory named "probdist" inside `pwd`
    exit 1
fi

cd probdist
rm -f *.py
cp ../../src/python/probdist/*.py .
cd ..

#############################################
# Copy python files for readnexus extension #
#############################################

if [ ! -e "readnexus" ]; then
    echo $0 expected to find a directory named "readnexus" inside `pwd`
    exit 1
fi

cd readnexus
rm -f *.py
cp ../../src/python/readnexus/*.py .
cd ..

#######################################
# Copy python files for pdfgen module #
#######################################

if [ -e "pdfgen" ]; then
    rm -rf pdfgen
fi

mkdir pdfgen
if [ $? -ne 0 ]; then
	echo $0 could not create pdfgen directory inside `pwd`
    exit 1
fi

cd pdfgen
cp ../../src/python/pdfgen/*.py .
cp -r ../../src/python/pdfgen/AFM .
cd ..

###########################################
# Copy python files for treeviewer module #
###########################################

if [ -e "treeviewer" ]; then
    rm -rf treeviewer
fi

mkdir treeviewer
if [ $? -ne 0 ]; then
	echo $0 could not create treeviewer directory inside `pwd`
    exit 1
fi

cd treeviewer
cp ../../src/python/treeviewer/*.py .
cd ..

##########################################
# Copy python files for utilities module #
##########################################

if [ -e "utilities" ]; then
    rm -rf utilities
fi

mkdir utilities
if [ $? -ne 0 ]; then
	echo $0 could not create utilities directory inside `pwd`
    exit 1
fi

cd utilities
cp ../../src/python/utilities/*.py .
cd ..

#########################################
# Copy python files for commands module #
#########################################

if [ -e "commands" ]; then
    rm -rf commands
fi

mkdir commands
if [ $? -ne 0 ]; then
	echo $0 could not create commands directory inside `pwd`
    exit 1
fi

cd commands
cp ../../src/python/commands/*.py .
cd ..

##############
# Copy tests #
##############

if [ -e "tests" ]; then
    rm -rf tests
fi

mkdir tests
if [ $? -ne 0 ]; then
	echo $0 could not create tests directory inside `pwd`
    exit 1
fi

cd tests
cp -r ../../tests/* .
cd ..

#################
# Copy examples #
#################

if [ -e "examples" ]; then
    rm -rf examples
fi

mkdir examples
if [ $? -ne 0 ]; then
	echo $0 could not create examples directory inside `pwd`
    exit 1
fi

cd examples
cp -r ../../examples/* .
cd ..


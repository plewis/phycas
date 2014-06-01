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
# Copy python files for Conversions extension #
###############################################

if [ ! -e "Conversions" ]; then
    echo $0 expected to find a directory named "Conversions" inside `pwd`
    exit 1
fi

cd Conversions
rm -f *.py
cp ../../src/python/conversions/*.py .
cd ..

##############################################
# Copy python files for DataMatrix extension #
##############################################

if [ ! -e "DataMatrix" ]; then
    echo $0 expected to find a directory named "DataMatrix" inside `pwd`
    exit 1
fi

cd DataMatrix
rm -f *.py
cp ../../src/python/datamatrix/*.py .
cd ..

##############################################
# Copy python files for Likelihood extension #
##############################################

if [ ! -e "Likelihood" ]; then
    echo $0 expected to find a directory named "Likelihood" inside `pwd`
    exit 1
fi

cd Likelihood
rm -f *.py
cp ../../src/python/likelihood/*.py .
cd ..

##############################################
# Copy python files for Likelihood extension #
##############################################

if [ ! -e "Phylogeny" ]; then
    echo $0 expected to find a directory named "Phylogeny" inside `pwd`
    exit 1
fi

cd Phylogeny
rm -f *.py
cp ../../src/python/phylogeny/*.py .
cd ..

############################################
# Copy python files for ProbDist extension #
############################################

if [ ! -e "ProbDist" ]; then
    echo $0 expected to find a directory named "ProbDist" inside `pwd`
    exit 1
fi

cd ProbDist
rm -f *.py
cp ../../src/python/probdist/*.py .
cd ..

#############################################
# Copy python files for ReadNexus extension #
#############################################

if [ ! -e "ReadNexus" ]; then
    echo $0 expected to find a directory named "ReadNexus" inside `pwd`
    exit 1
fi

cd ReadNexus
rm -f *.py
cp ../../src/python/readnexus/*.py .
cd ..

#######################################
# Copy python files for PDFGen module #
#######################################

if [ -e "PDFGen" ]; then
    rm -rf PDFGen
fi

mkdir PDFGen
if [ $? -ne 0 ]; then
	echo $0 could not create PDFGen directory inside `pwd`
    exit 1
fi

cd PDFGen
cp ../../src/python/pdfgen/*.py .
cp -r ../../src/python/pdfgen/AFM .
cd ..

###########################################
# Copy python files for TreeViewer module #
###########################################

if [ -e "TreeViewer" ]; then
    rm -rf TreeViewer
fi

mkdir TreeViewer
if [ $? -ne 0 ]; then
	echo $0 could not create TreeViewer directory inside `pwd`
    exit 1
fi

cd TreeViewer
cp ../../src/python/treeviewer/*.py .
cd ..

##########################################
# Copy python files for Utilities module #
##########################################

if [ -e "Utilities" ]; then
    rm -rf Utilities
fi

mkdir Utilities
if [ $? -ne 0 ]; then
	echo $0 could not create Utilities directory inside `pwd`
    exit 1
fi

cd Utilities
cp ../../src/python/utilities/*.py .
cd ..

#########################################
# Copy python files for Commands module #
#########################################

if [ -e "Commands" ]; then
    rm -rf Commands
fi

mkdir Commands
if [ $? -ne 0 ]; then
	echo $0 could not create Commands directory inside `pwd`
    exit 1
fi

cd Commands
cp ../../src/python/commands/*.py .
cd ..

##############
# Copy tests #
##############

if [ -e "Tests" ]; then
    rm -rf Tests
fi

mkdir Tests
if [ $? -ne 0 ]; then
	echo $0 could not create Tests directory inside `pwd`
    exit 1
fi

cd Tests
cp -r ../../tests/* .
cd ..

#################
# Copy examples #
#################

if [ -e "Examples" ]; then
    rm -rf Examples
fi

mkdir Examples
if [ $? -ne 0 ]; then
	echo $0 could not create Examples directory inside `pwd`
    exit 1
fi

cd Examples
cp -r ../../examples/* .
cd ..


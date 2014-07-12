#!/bin/bash

# Call this script from same level as Jamroot file to create a source distribution.
# Specify desired name of phycas directory on the command line:
#
# scripts/create_src_distr.sh phycas-2.0.0
#
# This will create a dist directory and a directory phycas-2.0.0 within dist.
# The phycas-2.0.0 directory will be archived to create a phycas-2.0.0-src.tar.gz file
# inside dist.

#####################################################
# Sanity checks to ensure we are in the right place #
#####################################################

if [ ! -e "scripts" ]; then
    echo $0 expected to find a directory named "scripts" inside `pwd`
    exit 1
fi

if [ ! -e "src" ]; then
    echo $0 expected to find a directory named "src" inside `pwd`
    exit 1
fi

if [ ! -e "src/cpp" ]; then
    echo $0 expected to find a directory named "src/cpp" inside `pwd`
    exit 1
fi

if [ ! -e "src/python" ]; then
    echo $0 expected to find a directory named "src/python" inside `pwd`
    exit 1
fi

if [ $# -ne 1 ]; then
    echo $0 expected one command line argument \(name of phycas source distribution directory\)
    exit 1
fi

#########################
# Create dist directory #
#########################

if [ -e "dist" ]; then
    rm -rf dist
fi

mkdir dist
if [ $? -ne 0 ]; then
	echo $0 could not create dist directory inside `pwd`
    exit 1
fi

#######################################
# Create phycas directory inside dist #
#######################################

mkdir dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not create $1 directory inside `pwd`/dist
    exit 1
fi

##########################
# Copy scripts directory #
##########################

cp -R scripts dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy scripts directory to `pwd`/dist/$1
    exit 1
fi

######################
# Copy src directory #
######################

cp -R src dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy src directory to `pwd`/dist/$1
    exit 1
fi

########################
# Copy tests directory #
########################

rm -f tests/diffs.txt

rm -f tests/codontest/params.p
rm -f tests/codontest/trees.t

rm -f tests/exploreprior/nodata.nex.p
rm -f tests/exploreprior/nodata.nex.t
rm -f tests/exploreprior/nodata.nex.txt

rm -f tests/fixedparams/fixed.p
rm -f tests/fixedparams/fixed.t
rm -f tests/fixedparams/output.txt

rm -f tests/gtrtest/gtr_test.p
rm -f tests/gtrtest/gtr_test.t
rm -f tests/gtrtest/output.txt

rm -f tests/likelihoodtest/check.nex
rm -f tests/likelihoodtest/simulated.nex

rm -f tests/mtss/ss-multiss.p
rm -f tests/mtss/ss-multiss.t
rm -f tests/mtss/ss-multiss.txt
rm -f tests/mtss/ss-multisump.txt

rm -f tests/pdftree/test.pdf

rm -f tests/polytomies/mcmcoutput.txt
rm -f tests/polytomies/params.p
rm -f tests/polytomies/simHKY.nex
rm -f tests/polytomies/sump-log.txt
rm -f tests/polytomies/sumt-log.txt
rm -f tests/polytomies/sumt-splits.pdf
rm -f tests/polytomies/sumt-trees.pdf
rm -f tests/polytomies/sumt-trees.tre
rm -f tests/polytomies/trees.t

rm -f tests/randgen/first100.txt

rm -f tests/rzyprior/rzy-explore-prior.p
rm -f tests/rzyprior/rzy-explore-prior.t
rm -f tests/rzyprior/rzy-explore-prior.txt

rm -f tests/rzyss/rzy-ss-refdist.txt
rm -f tests/rzyss/rzy-ss-sump.txt
rm -f tests/rzyss/rzy-ss.p
rm -f tests/rzyss/rzy-ss.t
rm -f tests/rzyss/rzy-ss.txt

rm -f tests/simulator/simulated.nex

rm -f tests/splittest/out.txt

rm -f /tests/sump/logfile.txt
rm -f /tests/sump/ss-jc4-params.p

rm -f tests/sumt/logfile.txt
rm -f tests/sumt/splits.pdf
rm -f tests/sumt/trees.pdf
rm -f tests/sumt/trees.tre

rm -f tests/underflow/output.txt

cp -R tests dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy tests directory to `pwd`/dist/$1
    exit 1
fi

###########################
# Copy examples directory #
###########################

cp -R examples dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy examples directory to `pwd`/dist/$1
    exit 1
fi

##################
# Copy README.md #
##################

cp README.md dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy README.md file to `pwd`/dist/$1
    exit 1
fi

################
# Copy INSTALL #
################

cp INSTALL dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy INSTALL file to `pwd`/dist/$1
    exit 1
fi

################
# Copy CHANGES #
################

cp CHANGES dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy CHANGES file to `pwd`/dist/$1
    exit 1
fi

#############
# Copy BUGS #
#############

cp BUGS dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy BUGS file to `pwd`/dist/$1
    exit 1
fi

################
# Copy Jamroot #
################

cp Jamroot dist/$1
if [ $? -ne 0 ]; then
	echo $0 could not copy Jamroot file to `pwd`/dist/$1
    exit 1
fi

######################
# Create tar.gz file #
######################

cd dist
tar zcvf $1-src.tar.gz $1
cd ..



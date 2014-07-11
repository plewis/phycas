#!/bin/bash

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

##################################
# Copy scripts directory to dist #
##################################

cp -R scripts dist
if [ $? -ne 0 ]; then
	echo $0 could not copy scripts directory to `pwd`/dist
    exit 1
fi

##############################
# Copy src directory to dist #
##############################

cp -R src dist
if [ $? -ne 0 ]; then
	echo $0 could not copy src directory to `pwd`/dist
    exit 1
fi

##########################
# Copy README.md to dist #
##########################

cp README.md dist
if [ $? -ne 0 ]; then
	echo $0 could not copy README.md file to `pwd`/dist
    exit 1
fi

########################
# Copy INSTALL to dist #
########################

cp INSTALL dist
if [ $? -ne 0 ]; then
	echo $0 could not copy INSTALL file to `pwd`/dist
    exit 1
fi

########################
# Copy CHANGES to dist #
########################

cp CHANGES dist
if [ $? -ne 0 ]; then
	echo $0 could not copy CHANGES file to `pwd`/dist
    exit 1
fi

#####################
# Copy BUGS to dist #
#####################

cp CHANGES dist
if [ $? -ne 0 ]; then
	echo $0 could not copy BUGS file to `pwd`/dist
    exit 1
fi

########################
# Copy Jamroot to dist #
########################

cp Jamroot dist
if [ $? -ne 0 ]; then
	echo $0 could not copy Jamroot file to `pwd`/dist
    exit 1
fi

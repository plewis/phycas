#!/bin/bash

# define Phycas script name, data file name, and destination folder
PHYCAS_SCRIPT="ss.py"
DATA_FILE="marshall.nex"
PROJECT_NAME="out-ss"

#BSUB -J ss

# you should not need to change anything below here
# from one run to the next 

WORKING_DIR="/work/$USER/$PROJECT_NAME"

# If working directory does not exist, create it
# The -p means "create parent directories as needed"
if [ ! -d "$WORKING_DIR" ]; then
mkdir -p $WORKING_DIR
fi

# If destination directory does not exist, create it
# The -p in mkdir means "create parent directories as needed"
if [ ! -d "$HOME/$PROJECT_NAME" ]; then
mkdir -p $HOME/$PROJECT_NAME
fi

# navigate to the working directory
cd $WORKING_DIR  

# Get script and input data from home directory
# and put a copy in the working directory
cp $HOME/$PHYCAS_SCRIPT .
cp $HOME/$DATA_FILE .

# Run the program
export PYTHONPATH="/shared/pol02003"
export LD_LIBRARY_PATH="/shared/pol02003/phycas/Conversions:/shared/pol02003/ncl"
python $PHYCAS_SCRIPT $LSB_JOBID $LSB_JOBINDEX

# copy output files back to your home directory
cp * $HOME/$PROJECT_NAME

# clean up scratch directory
rm -rf $WORKING_DIR

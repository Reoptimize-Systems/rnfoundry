#!/bin/bash
#
## file: release.sh
#
# This script creates releases of the RenewNet Foundry code
#
# Run this script from the directory containing it

# stop on first error
set -e

VERSION="0_1"
WORKING_COPY_DIR=$(pwd)
#echo $WORKING_COPY_DIR

RELEASE_NAME="RenewNet_$VERSION"

# create release directory
RELEASE_DIR="$WORKING_COPY_DIR/RenewNet_$VERSION"
mkdir $RELEASE_DIR

echo "Creating release $RELEASE_NAME in directory $RELEASE_DIR"

# export from the working directory to the release directory
hg archive $RELEASE_DIR
# remove file created by mercurial
rm $RELEASE_DIR/.hg_archival.txt
#remove hgignore file
rm $RELEASE_DIR/.hgignore
# remove the release scripts
rm $RELEASE_DIR/release.sh
rm $RELEASE_DIR/test_release.sh

# tar up the result
cd $WORKING_COPY_DIR
zip -r ${RELEASE_NAME}.zip $RELEASE_DIR/


# test
#cd $WORKING_COPY_DIR
#./test_release.sh $RELEASE_DIR


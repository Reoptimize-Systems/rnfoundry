#!/bin/bash
#
## file: release.sh
#
# This script creates releases of the RenewNet Foundry code
#
# Run this script from the directory containing it

# stop on first error
set -e

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
version="0_4"
working_copy_dir=$(pwd)
run_tests=false
copy_win_libs=true

while getopts "h?v:t:w" opt; do
    case "$opt" in
    h|\?)
        echo "Release script for Renewnet Foundry"
        exit 0
        ;;
    v)  version=$OPTARG
        echo "Version string changed to: $version"
        ;;
    t)  run_tests=true
        echo "run_tests: $run_tests"
        ;;
    w)  copy_win_libs=false
        echo "copy_win_libs: $copy_win_libs"
        ;;
    esac
done

#echo $working_copy_dir

echo "Releasing with version string: $version"

release_name="RenewNet_${version}"

# create release directory
release_dir="${working_copy_dir}/${release_name}"
mkdir $release_dir

echo "Creating release $release_name in directory $release_dir"

# export from the working directory to the release directory
hg archive $release_dir
# remove file created by mercurial
rm $release_dir/.hg_archival.txt
#remove hgignore file
rm $release_dir/.hgignore
# remove the release scripts
rm $release_dir/release.sh
rm $release_dir/test_release.sh

if [ "$copy_win_libs" = true ]; then

    # we need to copy a bunch of files cross-compiled using MXE to the
    # release so it can be built on windows machines
    mkdir $release_dir/x86_64-w64-mingw32_static
    mkdir $release_dir/x86_64-w64-mingw32_static/include
    mkdir $release_dir/x86_64-w64-mingw32_static/lib

    # gsl
    # octave require standard unix lib names (.a)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgsl.a  $release_dir/x86_64-w64-mingw32_static/lib/
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a $release_dir/x86_64-w64-mingw32_static/lib/
    # matlab needs libraries to have a different name (.lib)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgsl.a  $release_dir/x86_64-w64-mingw32_static/lib/gsl.lib
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a $release_dir/x86_64-w64-mingw32_static/lib/gslcblas.lib
    cp -r /opt/mxe/usr/x86_64-w64-mingw32.static/include/gsl/ $release_dir/x86_64-w64-mingw32_static/include/

    # f2c
    # octave requires standard unix lib names (.a)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libf2c.a  $release_dir/x86_64-w64-mingw32_static/lib/
    # matlab needs libraries to have a different name (.lib)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libf2c.a  $release_dir/x86_64-w64-mingw32_static/lib/f2c.lib

else
  echo "Not copying gsl gslcblas and f2c libraries"
fi

# zip up the result
cd $working_copy_dir
zip -qr ${release_name}.zip $release_name/

if [ "$run_tests" = true ]; then
  # test
  cd $working_copy_dir
  ./test_release.sh $release_dir
else
  echo "Skipping Tests"
fi



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
version="dev"
working_copy_dir=$(pwd)
run_tests=false
copy_win_libs=true
skip_mex=false
make_zip=true

usage="$(basename "$0") [-h] [-v <version>] [-t] [-w] [-m] [-z]  -- creates Renewnet Foundry release

where:
    -h  show this help text
    -v  set the release version string (default: $version)
    -t  run tests
    -w  copy windows libraries etc
    -m  skip building mex files (requires matlab)
    -z  don't create zip file"

while getopts "h?:vtwmz" opt; do
    case "$opt" in
    h|\?)
        echo "$usage"
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
    m)  skip_mex=true
        echo "skip_mex: $skip_mex"
        ;;
    z)  make_zip=false
        echo "make_zip: $make_zip"
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
    mkdir $release_dir/x86_64-w64-mingw32
    mkdir $release_dir/x86_64-w64-mingw32/include
    mkdir $release_dir/x86_64-w64-mingw32/lib

    # gsl
    # octave requires standard unix lib names (.a)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgsl.a  $release_dir/x86_64-w64-mingw32/lib/
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a $release_dir/x86_64-w64-mingw32/lib/
    # matlab needs libraries to have a different name (.lib)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgsl.a  $release_dir/x86_64-w64-mingw32/lib/gsl.lib
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a $release_dir/x86_64-w64-mingw32/lib/gslcblas.lib
    cp -r /opt/mxe/usr/x86_64-w64-mingw32.static/include/gsl/ $release_dir/x86_64-w64-mingw32/include/

    # f2c
    # octave requires standard unix lib names (.a)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libf2c.a  $release_dir/x86_64-w64-mingw32/lib/
    # matlab needs libraries to have a different name (.lib)
    cp /opt/mxe/usr/x86_64-w64-mingw32.static/lib/libf2c.a  $release_dir/x86_64-w64-mingw32/lib/f2c.lib

    # mbdyn
    cp -r ~/build/mbdyn/windows/* $release_dir/x86_64-w64-mingw32/
    # matlab needs libraries to have a different name (.lib)
    cp ~/build/mbdyn/windows/lib/libmbc.a  $release_dir/x86_64-w64-mingw32/lib/mbc.lib

else
  echo "Not copying gsl gslcblas and f2c libraries"
fi

if [ "$make_zip" = true ]; then
  # zip up the result
  cd $working_copy_dir
  zip -qr ${release_name}.zip $release_name/
fi

if [ "$skip_mex" = false ]; then
  if ! [ -x "$(command -v matlab)" ]; then
    echo 'matlab is not installed, not building mex files using Matlab.' >&2
  else
    # buld the mex files using matlab
    matlab -nodesktop -r "restoredefaultpath; cd('$release_dir'); rnfoundry_setup('Runtests', false, 'PreventXFemmCheck', true); quit"
    #matlab -nodesktop -r "restoredefaultpath; cd('$release_dir'); rnfoundry_setup('Runtests', false, 'PreventXFemmCheck', true, 'CrossBuildW64', true); quit"
  fi
fi

if [ "$run_tests" = true ]; then
  # test
  cd $working_copy_dir
  ./test_release.sh $release_dir
else
  echo "Skipping Tests"
fi



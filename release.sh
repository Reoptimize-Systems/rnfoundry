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
make_zip=false
make_docs=true
matlab_cmds=""
verbose=false

usage="$(basename "$0") [-h] [-v <version>] [-t] [-w] [-m] [-z] [-c <matlab_commands>] [-o] -- creates Renewnet Foundry release

where:
    -h  show this help text
    -v  set the release version string (default: $version)
    -t  run tests
    -w  copy windows libraries etc
    -m  skip building mex files (requires matlab)
    -z  create zip file
    -c  additional matlab commands to run before running rnfoundry_release
    -o  verbose output (default: false)"

while getopts "h?v:twmzc:o" opt; do
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
    z)  make_zip=true
        echo "make_zip: $make_zip"
        ;;
    c)  matlab_cmds=$OPTARG
        echo "matlab_cmds : $matlab_cmds"
        ;;
    o)  verbose=true
        echo "verbose: $verbose"
        ;;
    esac
done

#echo ${working_copy_dir}

echo "Releasing with version string: ${version}"

release_name="RenewNet_${version}"

# create release directory
release_dir="/tmp/${release_name}"

mkdir -p ${release_dir}

rm -rf ${release_dir}/*

echo "Creating release ${release_name} in directory ${release_dir}"

# export from the working directory to the release directory
hg archive ${release_dir}
# remove file created by mercurial
rm ${release_dir}/.hg_archival.txt
#remove hgignore file
rm ${release_dir}/.hgignore
# remove the release scripts
rm ${release_dir}/release.sh
rm ${release_dir}/test_release.sh
# remove test scripts
rm ${release_dir}/test_hg.sh

# create version.txt
echo ${version} > ${release_dir}/version.txt
echo ${version} > ${release_dir}/wave/matlab-octave/wec-sim/version.txt

if [ "$copy_win_libs" = true ]; then

    mxedir=/opt/mxe

    cd ${mxedir}
    make gsl libf2c

    # we need to copy a bunch of files cross-compiled using MXE to the
    # release so it can be built on windows machines
    mkdir ${release_dir}/x86_64-w64-mingw32
    mkdir ${release_dir}/x86_64-w64-mingw32/include
    mkdir ${release_dir}/x86_64-w64-mingw32/lib

    # gsl
    # octave requires standard unix lib names (.a)
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libgsl.a  ${release_dir}/x86_64-w64-mingw32/lib/
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a ${release_dir}/x86_64-w64-mingw32/lib/
    # matlab needs libraries to have a different name (.lib)
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libgsl.a  ${release_dir}/x86_64-w64-mingw32/lib/gsl.lib
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libgslcblas.a ${release_dir}/x86_64-w64-mingw32/lib/gslcblas.lib
    cp -r ${mxedir}/usr/x86_64-w64-mingw32.static/include/gsl/ ${release_dir}/x86_64-w64-mingw32/include/

    # f2c
    # octave requires standard unix lib names (.a)
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libf2c.a  ${release_dir}/x86_64-w64-mingw32/lib/
    # matlab needs libraries to have a different name (.lib)
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libf2c.a  ${release_dir}/x86_64-w64-mingw32/lib/f2c.lib

    # mbdyn
    cp -r ~/build/mbdyn/x86_64-w64-mingw32_static/* ${release_dir}/x86_64-w64-mingw32/
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libws2_32.a ${release_dir}/x86_64-w64-mingw32/lib/
    # matlab needs libraries to have a different name (.lib)
    cp ~/build/mbdyn/x86_64-w64-mingw32_static/lib/libmbc.a  ${release_dir}/x86_64-w64-mingw32/lib/mbc.lib
    cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libws2_32.a ${release_dir}/x86_64-w64-mingw32/lib/ws2_32.lib

    # boost (for shared memory communication)
    #mkdir -p ${release_dir}/x86_64-w64-mingw32/include/boost
    #cp -r ${mxedir}/usr/x86_64-w64-mingw32.static/include/boost/interprocess ${release_dir}/x86_64-w64-mingw32/include/boost/interprocess
    #cp -r ${mxedir}/usr/x86_64-w64-mingw32.static/include/boost/date_time ${release_dir}/x86_64-w64-mingw32/include/boost/date_time

else
  echo "Not copying gsl gslcblas and f2c libraries, or mbdyn program"
fi

if [ "$skip_mex" = false ]; then
  if ! [ -x "$(command -v matlab)" ]; then
    echo 'matlab is not installed, not building mex files using Matlab.' >&2
  else
    # buld the mex files using (oldish) version of matlab
    /usr/local/MATLAB/R2016b/bin/matlab -nodesktop -r "restoredefaultpath; cd('${release_dir}'); ${matlab_cmds}; rnfoundry_release ('Throw', true, 'Verbose', ${verbose}); quit"
  fi
fi

rm ${release_dir}/rnfoundry_release.m

if [ "$make_zip" = true ]; then
  # zip up the result
  cd ${working_copy_dir}
  zip -qr ${release_name}.zip ${release_name}/
fi

#if [ "$make_docs" = true ]; then
#  # copy the mbdyn sphinx docs directory to /tmp
#  cp ${working_copy_dir}
#  zip -qr ${release_name}.zip $release_name/
#fi

if [ "$run_tests" = true ]; then
  # test
  cd ${working_copy_dir}
  ./test_release.sh ${release_dir}
else
  echo "Skipping Tests"
fi





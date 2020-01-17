#!/bin/bash
#
## file: release.sh
#
# This script creates releases of the RenewNet Foundry code
#
# Run this script from the directory containing it

# stop on first error
set -ex

# A POSIX variable
OPTIND=1  # Reset in case getopts has been used previously in the shell.

# Initialize our own variables:
version="dev"
working_copy_dir=$(pwd)
run_tests=false
copy_win_libs=true
skip_mex=false
make_win_mex=false
make_zip=false
make_docs=true
matlab_cmds=""
verbose=false
launch_matlab_cmd="matlab"
mbdyn_release_dir="${HOME}/build/mbdyn-official-git/x86_64-w64-mingw32_static"
mxedir="/opt/mxe"
matlab_windows_lib_dir="${HOME}/Nextcloud/Personal/work/matlab_windows_libs/r2016b/extern/lib/win64/mingw64"
mbdyn_src_dir=""
build_windows_mbdyn=false
output_dir="/tmp"

usage="$(basename "$0") [-h] [-v <version>] [-t] [-w] [-m] [-W] [-z] [-c <matlab_commands>] [-o] [b <matlab_cmd>] [-M <mbdyn_release_dir>] [-x <mxedir>] [-l <mat_lib_dir>] -- creates Renewnet Foundry release

where:
    -h  show this help text
    -v  set the release version string (default: $version)
    -t  run tests
    -w  don't build windows libraries using MXE or copy them
    -m  skip building any mex files (building mex requires matlab)
    -W  attempt to cross-build mex files for windows
    -z  create zip file
    -c  additional matlab commands to run before running rnfoundry_release (default: "${matlab_cmds}")
    -n  (noisy) verbose output (default: ${verbose})
    -b  command to run matlab (default: "${launch_matlab_cmd}")
    -M  mbdyn release dir (default: "${mbdyn_release_dir}")
    -x  MXE install dir (default: "${mxedir}")
    -l  directory containing the windows mingw64 mex libraries (default: "${matlab_windows_lib_dir}")
    -s  mbdyn source directory (default: "${mbdyn_src_dir}")
    -t  cross build mbdyn. (default: ${build_windows_mbdyn})
    -o  output directory where the directory containing the release will be created (default: "${output_dir}")"

while getopts "h?v:twmWzc:nM:x:l:b:s:to:" opt; do
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
    W)  make_win_mex=true
        echo "make_win_mex: $make_win_mex"
        ;;
    z)  make_zip=true
        echo "make_zip: $make_zip"
        ;;
    c)  matlab_cmds=$OPTARG
        echo "matlab_cmds : $matlab_cmds"
        ;;
    n)  verbose=true
        echo "verbose: $verbose"
        ;;
    b)  launch_matlab_cmd=$OPTARG
        echo "launch_matlab_cmd: $launch_matlab_cmd"
        ;;
    M)  mbdyn_release_dir=$OPTARG
        echo "mbdyn_release_dir: $mbdyn_release_dir"
        ;;
    s)  mbdyn_src_dir=$OPTARG
        echo "mbdyn_src_dir: $mbdyn_src_dir"
        ;;
    t)  build_windows_mbdyn=true
        echo "build_windows_mbdyn: $build_windows_mbdyn"
        ;;
    x)  mxedir=$OPTARG
        echo "mxedir: $mxedir"
        ;;
    l)  matlab_windows_lib_dir=$OPTARG
        echo "matlab_windows_lib_dir: $matlab_windows_lib_dir"
        ;;
    o)  output_dir=$OPTARG
        echo "output_dir: $output_dir"
        ;;
    esac
done

#echo ${working_copy_dir}

echo "Releasing with version string: ${version}"

release_name="RenewNet_${version}"

# create release directory
release_dir="${output_dir}/${release_name}"

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
if [ -f "${release_dir}/upload_release.sh" ]; then
  rm ${release_dir}/upload_release.sh
fi
rm ${release_dir}/test_release.sh
# remove test scripts
rm ${release_dir}/test_hg.sh

# create version.txt
echo ${version} > ${release_dir}/version.txt
echo ${version} > ${release_dir}/wave/matlab-octave/wec-sim/version.txt



if [ "$copy_win_libs" = true ]; then

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
  if [ "$build_windows_mbdyn" = true ]; then
    if [ -d "$mbdyn_src_dir" ]; then
      cd ${mbdyn_src_dir}/packaging/mswindows
      ./win_package.sh
    fi
  fi

  # mbdyn
  cp -r ${mbdyn_release_dir}/* ${release_dir}/x86_64-w64-mingw32/
  cp ${mxedir}/usr/x86_64-w64-mingw32.static/lib/libws2_32.a ${release_dir}/x86_64-w64-mingw32/lib/
  # matlab needs libraries to have a different name (.lib)
  cp ${mbdyn_release_dir}/lib/libmbc.a  ${release_dir}/x86_64-w64-mingw32/lib/mbc.lib
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
    ${launch_matlab_cmd} -nodesktop -r "restoredefaultpath; cd('${release_dir}'); ${matlab_cmds}; rnfoundry_release ('Throw', true, 'Verbose', ${verbose}, 'Version', '${version}', 'W64CrossBuild', ${make_win_mex}, 'W64CrossBuildMexLibsDir', '${matlab_windows_lib_dir}'); quit"
  fi
fi

rm ${release_dir}/rnfoundry_release.m

if [ "$make_zip" = true ]; then
  # zip up the resultto maximise profit.
  cd ${release_dir}/..
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





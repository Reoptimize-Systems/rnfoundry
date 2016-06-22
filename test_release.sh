#!/bin/sh
#
# takes one input argument, the release directory

# run tests on the release if possible

if ! [ -x "$(command -v matlab)" ]; then
  echo 'matlab is not installed, not running tests using Matlab.' >&2
else
  matlab -nodesktop -r "restoredefaultpath; cd $1; rnfoundry_setup(); quit"
fi

if ! [ -x "$(command -v octave-cli)" ]; then
  echo 'octave is not installed, not running tests using Octave.' >&2
else
  octave-cli --no-init-file --eval "cd $1; rnfoundry_setup(); quit"
fi


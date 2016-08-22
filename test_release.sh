#!/bin/sh
#
# takes one input argument, the release directory

# run tests on the release if possible

RED='\033[0;31m'
BRed='\033[1;31m'
NC='\033[0m' # No Color

printf "Testing RenewNet release in dir ${1}\n"

if ! [ -x "$(command -v matlab)" ]; then
  printf 'matlab is not installed, not running tests using Matlab.\n' >&2
else
  printf "${BRed}Testing release using Matlab.${NC}\n"
  matlab -nodesktop -r "restoredefaultpath; fprintf(1, 'changing dir to %s\n', strtrim('${1}')); cd(strtrim('${1}')); rnfoundry_setup('RunTests', true); quit"
fi

if ! [ -x "$(command -v octave-cli)" ]; then
  printf 'octave is not installed, not running tests using Octave.\n' >&2
else
  printf "${BRed}Testing release using Octave.${NC}\n"
  octave-cli --no-init-file --eval "pkg load odepkg; fprintf(1, 'changing dir to %s\n', strtrim('${1}')); cd(strtrim('${1}')); rehash(); rnfoundry_setup('RunTests', true); quit"
fi


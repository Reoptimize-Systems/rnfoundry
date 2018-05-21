#!/bin/sh

set -ex

repo_co_dir="/tmp/rnfoundry-hg-for-testing"

if [ ! -d "${repo_co_dir}" ]; then
  hg clone http://hg.code.sf.net/p/rnfoundry/hgrepo ${repo_co_dir}
fi

cd ${repo_co_dir}

hg purge --all
hg pull
hg update

octave --no-init-file --eval "dbstop if error; fprintf(1, 'changing dir to %s\n', strtrim('${repo_co_dir}')); cd(strtrim('${repo_co_dir}')); rehash(); rnfoundry_setup('RunTests', true); quit"

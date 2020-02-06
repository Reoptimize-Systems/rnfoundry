#!/bin/bash

set -ex

version=$1
release_dir=$2

release_name="RenewNet_${version}"

# create release directory
#release_dir="/tmp/${release_name}"

cd ${release_dir}/doc/sphinx
make html

cd ${release_dir}/..

zip -qr ${release_name}.zip ${release_name}/

# upload the releases
rsync -ravP -e ssh ${release_dir}/../RenewNet_${version}.zip ${release_dir}/{README.rst,CHANGELOG.md}  crobarcro@frs.sourceforge.net:"/home/frs/project/rnfoundry/Release\\ ${version}/"
#rsync -ravP -e ssh ${release_dir}/../README.rst  crobarcro@frs.sourceforge.net:"/home/frs/project/rnfoundry/Release\\ ${version}/README.rst"
#rsync -ravP -e ssh ${release_dir}/../CHANGELOG.md  crobarcro@frs.sourceforge.net:"/home/frs/project/rnfoundry/Release\\ ${version}/CHANGELOG.md"

# upload documentation
rsync -ravP -e ssh ${release_dir}/doc/sphinx/_build/html/  ${release_dir}/documentation  crobarcro@web.sourceforge.net:/home/project-web/rnfoundry/htdocs/


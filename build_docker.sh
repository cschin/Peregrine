#!/bin/bash
set -e o pipefail
pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
popd

pushd src
make clean
popd

tar czvf src.tgz src/ falcon/ py/ .git/
mv src.tgz docker/

pushd docker/
if [ $1 == 'tag' ] 
then
    tag=$(git describe --always --abbrev=0 --tags)
    tag=${tag:2}
else
    tag=latest
fi
echo current docker tag: ${tag}
docker build . --tag cschin/peregrine:${tag}
popd



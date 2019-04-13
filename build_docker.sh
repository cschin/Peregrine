#!/bin/bash

pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
popd

pushd src
make clean
popd

tar czvf src.tgz src/ falcon/ py/ .git/
mv src.tgz docker/

pushd docker/
tag=$(git describe --abbrev=0 --tags)
tag=${tag:2}
docker build . --tag cschin/peregrine:${tag}
popd

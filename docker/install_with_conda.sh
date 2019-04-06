#!/bin/bash
. /opt/conda/bin/activate
conda create -n peregrine -y python=3.7

conda activate peregrine
conda install -c conda-forge -y pypy3.6

pushd py
rm -rf .eggs/ dist/ build/ peregrine.egg-info/ peregrine_pypy.egg-info get-pip.py
python3 setup.py install
python3 setup.py clean --all
popd
git clone -b peregrine https://github.com/cschin/pypeFLOW.git
pushd pypeFLOW
python3 setup.py install
popd
pushd py
wget -q https://bootstrap.pypa.io/get-pip.py
wget -q https://bootstrap.pypa.io/get-pip.py
pypy3 get-pip.py
pypy3 setup_pypy.py install
popd

pushd src
make all
make install
popd

python3 -m pip install cffi==1.12.2

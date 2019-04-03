from setuptools import setup
import os
os.environ["peregrine_base"] = os.path.abspath(os.path.pardir)

setup(name='peregrine',
      version='0.1',
      packages=['peregrine'],
      package_dir = {'peregrine': 'peregrine'},
      scripts = ["scripts/path_to_contig.py", "scripts/cns_prototype.py"],
      setup_requires=["cffi>=1.12.0"],
      cffi_modules=["peregrine/build_shimmer4py.py:ffibuilder",
                    "peregrine/build_falcon4py.py:ffibuilder"],
      install_requires=["cffi>=1.12.0", "numpy>=1.16.2"])

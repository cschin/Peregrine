from setuptools import setup
import versioneer
import os
os.environ["peregrine_base"] = os.path.abspath(os.path.pardir)

setup(name='peregrine',
      version=versioneer.get_version(),
      cmdclass=versioneer.get_cmdclass(),
      packages=['peregrine'],
      package_dir = {'peregrine': 'peregrine'},
      scripts = ["scripts/path_to_contig.py",
                 "scripts/pg_asm_cns.py",
                 "scripts/pg_run.py",
                 "scripts/pg_run_dev.py"],
      setup_requires=["cffi>=1.12.0",
                      "versioneer==0.18"],
      cffi_modules=["peregrine/build_shimmer4py.py:ffibuilder",
                    "peregrine/build_falcon4py.py:ffibuilder"],
      install_requires=["cffi>=1.12.0",
                        "docopt>=0.6.2",
                        "numpy>=1.16.2"])

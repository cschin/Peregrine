import sys
from ._version import get_versions
__version__ = get_versions()['version']
sys.stderr.write(f"Peregrine Assembler & SHIMMER ASMKit({__version__})\n")
del get_versions

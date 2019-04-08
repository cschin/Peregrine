# Peregrine &  SHIMMER Genome Assembly Toolkit

Peregrine is a fast genome assembler for accurate long
reads (length > 10kb, accuraccy > 99%). It can assemble
a human genome from 30x reads within 20 cpu hours from
reads to polished consensus. It uses Sparse HIereachical
MimiMizER (SHIMMER) for fast read-to-read overlaping
without quadratic comparisions used in other OLC
assemblers.

This code base includes code that uses SHIMMER (Sparse 
HIerarchical MiniMimER) for genome assembly and other
related applications.

Currently, the assembly graph process is more or less
identical to the approaches used in the FALCON assembler
developed by Jason Chin and others in Pacific Biosciences,
Inc. There are a number of other possible ways to generate
contigs without a  string graph but it will need some
research work to make it happening. The FALCON graph
module is also not very efficient as python scripts
running single thread mode.


## Install

See the `docker/Dockerfile` and `docker/install_with_conda.sh`
as examples to install from scratch within a Conda
environemnt.


## LICENSE

### Peregrine &  SHIMMER Genome Assembly Toolkit

Peregrine Assembler and SHIMMER Genome Assembly Toolkit
Copyright (c) 2019- by Jason, Chen-Shan, Chin

Peregrine Assembler and  SHIMMER Genome Assembly Toolkit 
is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International 
License.

You should have received a copy of the license along with this
work. If not, see <http://creativecommons.org/licenses/by-nc-sa/4.0/>.


### Minimap2

SHIMMER genome assembly toolkit uses C library developed by
Heng Li for Minimap2.  See LICENSE.minimap2


### FALCON

See LICENSE.falcon for license for the code from FALCON 

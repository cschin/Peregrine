  <img src="misc/logo.png" alt="PeregrineLogo" width="120"/>

# Peregrine & SHIMMER Genome Assembly Toolkit

Peregrine is a fast genome assembler for accurate long reads (length > 10kb,
accuraccy > 99%). It can assemble a human genome from 30x reads within 20 cpu
hours from reads to polished consensus. It uses Sparse HIereachical MimiMizER
(SHIMMER) for fast read-to-read overlaping without quadratic comparisions used
in other OLC assemblers.

This code base includes code that uses SHIMMER (Sparse HIerarchical MiniMimER)
for genome assembly and other related applications.

Currently, the assembly graph process is more or less identical to the
approaches used in the FALCON assembler developed by Jason Chin and others in
Pacific Biosciences, Inc. There are a number of other possible ways to generate
contigs without a  string graph but it will need some research work to make it
happening. The FALCON graph module is also not very efficient as python scripts
running single thread mode.


## Install

We *do not* recommend that you install the software from the source code unless
you are comfortable handling the required dependences for your system
independently. Unless you have full control (e.g. root access) of the computer 
system you use to build Peregrine and you can install the proper GCC 
compiler/python/pypy/conda version, then you should try to learn to use [Docker 
images](https://hub.docker.com/r/cschin/peregrine/tags) that 
we provide, it will make your life and our life easier. 

As independent deverlopers with limit resource, we cannot provide free support for 
solving dependence problem of your specific system. Instead, we can provide 
docker image so you can run the executables and their dependency using Docker.  

If you want to build for your system without using Docker, please see the 
`docker/Dockerfile` and `docker/install_with_conda.sh` as examples to
install from scratch within a clean Conda environemnt.

## Run the assembler

Peregrine is designed to run on single compute node. It does not need a grid
computing job scheduling system. It uses Pypeflow to coordinate multiple
concurrent processes.  

After revsion 0.1.5.3, You can test a small assembly using simulated E. Coli 
reads with Docker:

```
# please substitue $PWD and $IMAGETAG with proper values
docker run -it --rm -v $PWD:/wd cschin/peregrine:$IMAGETAG test
```

The assembly results are in `$PWD/ecoli_test_results/`. The testing case will
download an E. Coli reference and generate simulated reads. After the assembly
is done, it also installs `nucmer` to run `dandiff` comparing the assembled 
contigs with the original E. coli reference. You can check the ouput by using 
`cat $PWD/ecoli_test_results/out.report` command.

Here is the general usage for `pg_run.py` which starts the workflow for 
assembling a genome from input `fasta`, `fastq`, `fasta.gz` or 
`fastq.gz` files. 

```
Usage:
  pg_run.py asm <reads.lst> <index_nchunk> <index_nproc>
                            <ovlp_nchunk> <ovlp_nproc>
                            <mapping_nchunk> <mapping_nproc>
                            <cns_nchunk> <cns_nproc>
                            <sort_nproc>
                            [--with-consensus]
                            [--with-L0-index]
                            [--output <output>]
                            [--shimmer-k <shimmer_k>]
                            [--shimmer-w <shimmer_w>]
                            [--shimmer-r <shimmer_r>]
                            [--shimmer-l <shimmer_l>]
                            [--best_n_ovlp <n_ovlp>]
                            [--mc_lower <mc_lower>]
                            [--mc_upper <mc_upper>]
                            [--aln_bw <aln_bw>]
                            [--ovlp_upper <ovlp_upper>]
  pg_run.py (-h | --help)
  pg_run.py --verison

Options:
  -h --help                   Show this help
  --version                   Show version
  --with-consensus            Generate consensus after getting the draft contigs
  --with-L0-index             Keep level-0 index
  --output <output>           Set output directory (will be created if not exist) [default: ./wd]
  --shimmer-k <shimmer_k>     Level 0 k-mer size [default: 16]
  --shimmer-w <shimmer_w>     Level 0 window size [default: 80]
  --shimmer-r <shimmer_r>     Reduction factore for high level SHIMMER [default: 6]
  --shimmer-l <shimmer_l>     number of level of shimmer used, the value should be 1 or 2 [default: 2]
  --best_n_ovlp <n_ovlp>      Find best n_ovlp overlap [default: 4]
  --mc_lower <mc_lower>       Does not cosider SHIMMER with count less than mc_low [default: 2]
  --mc_upper <mc_upper>       Does not cosider SHIMMER with count greater than mc_upper [default: 240]
  --aln_bw <aln_bw>           Max off-diagonal gap allow during overlap confirmation [default: 100]
  --ovlp_upper <ovlp_upper>   Ignore cluster with overlap count greater ovlp_upper [default: 120]
```

The first required option is `reads.lst`.  The `reads.list` should a
path to a file that contains the list of the paths of the input sequence files.

The rest required options specify how to partition the data for different part
of the pipeline and the number of the processors used for each of the step.

`<index_nchunk>`  and `<index_nproc>` control the number of "chunks" and the
number of cpu used concurrently for the initial SHIMMER index generation.

`<ovlp_nchunk>`  and `<ovlp_nproc>` control the number of "chunks" and the
number of cpu used concurrently for generating overlap inforrmation between
reads. This part typically use the most memory and the exact size of RAM used
concurrently depends on the size of input sequence data and the index file
size. 

You can use larger number of `<ovlp_nchunk>` and smaller number of
`<ovlp_nproc>` on a smaller memory mechine. For example, I was able to finish
this part using a machine with 32G RAM with `ovlp_nchunk=24` and
`ovlp_nproc=1`. 

If there is enough memory, for example, AWS bothe m5d.metal and r5d.12xlarge
have 384G RAM, they can support running 24 to 48 cpu cores at once. However,
the overlap step needs to do random access the sequence data through shared
memory mapped file, it will be great to reserve some RAM to cache the sequence
in memory in RAM. In our test, 48 cores does not provide significant speeding
comparing to use 24 cores. Also, if there is not enough memory, you may need
fast SSD or nvme drives and reduce the number or CPU core concurrently
accessing the sequence data.

`<mapping_nchunk>` and `<mapping_nproc>` control the partitioning and the
number of cores used for mapping the sequence reads to draft contigs for the
following consensus step.

`<sort_nproc>` controls the number of cpu cores used for sorting the reads to
contigs map.

`<cns_nchunk>` and  `<cns_nproc>` control the partitioning and the number of
cores used for generating the consensus from draft contigs.


## Runing Peregrine Using Docker

Here is an example running Peregrine with Docker for a Peregrin build 
of tag 0.1.5.0 using an AWS m5d.metal or r5d.12xlarge instance. (You will
need to configure the AWS instance to utilize the NVME drives and a 
docker environment.)

```
find /wd/chm13-fastq/ -name "*.fastq" | sort > chm13-seqdata.lst 

docker run -it -v /wd:/wd --user $(id -u):$(id -g) cschin/peregrine:0.1.5.0 asm \
    /wd/chm13-seqdata.lst 24 24 24 24 24 24 24 24 24 \ 
    --with-consensus --shimmer-r 3 --best_n_ovlp 8 \ 
    --output /wd/chm13-asm-r3-pg0.1.5.0 
```

Note that the paths in the `<reads.lst>` should be the full paths to the
sequuence files inside the docker container.


## LICENSE

### Peregrine & SHIMMER Genome Assembly Toolkit

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

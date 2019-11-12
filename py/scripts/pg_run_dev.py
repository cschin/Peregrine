import peregrine
from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, Dist)
from pypeflow.tasks import gen_task
from docopt import docopt
import logging
import os
import sys

__version__ = peregrine.__version__

__doc__ = """
Peregrine
===============================================================
Peregrine is a fast genome assembler for accurate long
reads (length > 10kb, accuraccy > 99%). It can assemble
a human genome from 30x reads within 20 cpu hours from
reads to polished consensus. It uses Sparse HIereachical
MimiMizER (SHIMMER) for fast read-to-read overlaps without
explicitly quadratic comparisions used in other OLC
assemblers.

Currently, the assembly graph process is more or less
identical to the approaches used in the FALCON assembler
developed by Jason Chin and others in Pacific Biosciences, Inc.

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
                            [--min_len <min_len>]
                            [--min_idt <min_idt>]
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
  --min_len <min_len>         Minimum overlap length for assembly graph construction [default: 4000]
  --min_idt <min_idt>         Minimum identity for considering two reads that are properly overlaps [default: 96]

Licenses:

Peregrine Assembler and SHIMMER Genome Assembly Toolkit
Copyright (c) 2019- by Jason, Chen-Shan, Chin

Peregrine Assembler and  SHIMMER Genome Assembly Toolkit
is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International
License.

You should have received a copy of the license along with
this work. If not, see
<http://creativecommons.org/licenses/by-nc-sa/4.0/>.

************************************************************
If you want to use it for any commericial purposes
(including promotion activity for a commerical product),
please contact Jason Chin for a commericial license.
************************************************************

This software uses the following libraray from Heng Li's
Minimap2 codebase under MIT License:

mm_sketch.c kvec.h kseq.h khash.h kalloc.h kalloc.c

The MIT License

Copyright (c) 2018-     Dana-Farber Cancer Institute
2017-2018 Broad Institute, Inc.

Permission is hereby granted, free of charge, to any person obtaining
a copy of this software and associated documentation files (the
"Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish,
distribute, sublicense, and/or sell copies of the Software, and to
permit persons to whom the Software is furnished to do so, subject to
the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


License for code from FALCON

#################################################################################$$
# Copyright (c) 2011-2015, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted (subject to the limitations in the
# disclaimer below) provided that the following conditions are met:
#
#  * Redistributions of source code must retain the above copyright
#  notice, this list of conditions and the following disclaimer.
#
#  * Redistributions in binary form must reproduce the above
#  copyright notice, this list of conditions and the following
#  disclaimer in the documentation and/or other materials provided
#  with the distribution.
#
#  * Neither the name of Pacific Biosciences nor the names of its
#  contributors may be used to endorse or promote products derived
#  from this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
# GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
# BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
# CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
# USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
# OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
# SUCH DAMAGE.
#################################################################################$$
"""


LOG = logging.getLogger(__name__)


def Task(script="", inputs={}, outputs={}, parameters=None, dist=None):
    if parameters is None:
        parameters = dict()
    if dist is None:
        dist = Dist()

    # Make paths relative to CWD. (But ok if caller does this.)
    def get_rel(maybe_abs):
        rel = dict()
        for (k, v) in maybe_abs.items():
            try:
                if os.path.isabs(v):
                    v = os.path.relpath(v)
                rel[k] = v
            except Exception:
                LOG.exception('Error for {!r}->{!r}'.format(k, v))
                raise
        return rel
    inputs = get_rel(inputs)
    outputs = get_rel(outputs)

    # All outputs must be in same directory.

    params = dict(parameters)

    pt = gen_task(script, inputs, outputs, params, dist)

    return pt


def run_build_db(wf, args, seq_dataset_lst_fn):

    build_db = """
/usr/bin/time shmr_mkseqdb \
    -p {params.read_db_prefix} \
    -d {input.seq_dataset}
"""

    read_db_dir = os.path.join(os.path.abspath(args["--output"]), "0-seqdb")
    read_db_prefix = "seq_dataset"
    seq_dataset = seq_dataset_lst_fn
    read_db = os.path.join(read_db_dir, f"{read_db_prefix}.seqdb")
    seqidx = os.path.join(read_db_dir, f"{read_db_prefix}.idx")
    read_db_abs_prefix = os.path.join(read_db_dir, read_db_prefix)
    outputs = {'read_db': read_db,
               'seqidx': seqidx}

    wf.addTask(Task(
        script=build_db,
        inputs={'seq_dataset': seq_dataset},
        outputs=outputs,
        parameters={
            'read_db_prefix': read_db_abs_prefix
        },
        dist=Dist(NPROC=1, local=True)
    ))

    wf.max_jobs = 1
    wf.refreshTargets()

    return read_db_abs_prefix, outputs


def run_build_idx(wf, args, read_db_abs_prefix):
    n_chunk = int(args["<index_nchunk>"])
    n_proc = int(args["<index_nproc>"])
    shimmer_l = int(args["--shimmer-l"])

    build_idx = """
/usr/bin/time shmr_index\
    -m {params.output_L0_index}\
    -p {params.read_db_prefix}\
    -k {params.shimmer_k}\
    -w {params.shimmer_w}\
    -r {params.shimmer_r}\
    -l {params.shimmer_l}\
    -t {params.n_chunk}\
    -c {params.my_chunk}\
    -o {params.index_prefix}
ln -s {params.index_prefix}* {params.index_dir}
"""
    index_dir = os.path.join(os.path.abspath(args["--output"]), "1-index")
    if shimmer_l == 2:
        index_abs_prefix = os.path.join(index_dir, "shmr-L2")
    elif shimmer_l == 1:
        index_abs_prefix = os.path.join(index_dir, "shmr-L1")
    else:
        sys.exit(1)

    outputs = {}
    for my_chunk in range(1, n_chunk+1):
        index_chunk_dir = os.path.join(index_dir, f"chunk-{my_chunk:02d}")
        index_chunk_abs_prefix = os.path.join(index_chunk_dir, "shmr")
        # index_L0_fn = f"{index_chunk_abs_prefix}-L0-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        # index_L0_MC_fn = f"{index_chunk_abs_prefix}-L0-MC-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        if shimmer_l == 2:
            index_shmr_fn = f"{index_chunk_abs_prefix}-L2-{my_chunk:02d}-of-{n_chunk:02d}.dat"
            index_shmr_MC_fn = f"{index_chunk_abs_prefix}-L2-MC-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        elif shimmer_l == 1:
            index_shmr_fn = f"{index_chunk_abs_prefix}-L1-{my_chunk:02d}-of-{n_chunk:02d}.dat"
            index_shmr_MC_fn = f"{index_chunk_abs_prefix}-L1-MC-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        else:
            sys.exit(1)

        outputs[f'index_shmr_{my_chunk:02d}'] = index_shmr_fn
        outputs[f'index_shmr_MC_{my_chunk:02d}'] = index_shmr_MC_fn
        wf.addTask(Task(
            script=build_idx,
            inputs={
                'read_db': f"{read_db_abs_prefix}.seqdb",
                'seqidx': f"{read_db_abs_prefix}.idx"
            },
            outputs={
                'index_shmr': index_shmr_fn,
                'index_shmr_MC': index_shmr_MC_fn
            },
            parameters={
                'read_db_prefix': read_db_abs_prefix,
                'index_prefix': index_chunk_abs_prefix,
                'index_dir': index_dir,
                'output_L0_index': 1 if args["--with-L0-index"] else 0,
                'shimmer_k': int(args["--shimmer-k"]),
                'shimmer_w': int(args["--shimmer-w"]),
                'shimmer_r': int(args["--shimmer-r"]),
                'shimmer_l': int(args["--shimmer-l"]),
                'n_chunk': n_chunk,
                'my_chunk': my_chunk
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = n_proc
    wf.refreshTargets()

    return index_abs_prefix, outputs


def run_overlapper(wf, args,
                   read_db_abs_prefix, index_abs_prefix, ovlp_in):
    n_chunk = int(args["<ovlp_nchunk>"])
    n_proc = int(args["<ovlp_nproc>"])
    shmr_ovlp_script = """
/usr/bin/time shmr_overlap\
    -b {params.best_n_ovlp}\
    -m {params.mc_lower}\
    -M {params.mc_upper}\
    -w {params.align_bandwidth}\
    -n {params.ovlp_upper}\
    -p {params.read_db_prefix}\
    -l {params.index_prefix}\
    -t {params.n_chunk}\
    -c {params.my_chunk}\
    -o {output.ovlp_out}
"""
    ovlp_dir = os.path.join(os.path.abspath(args["--output"]), "2-ovlp")
    outputs = {}
    for my_chunk in range(1, n_chunk+1):
        ovlp_chunk_dir = os.path.join(ovlp_dir, f"chunk-{my_chunk:02d}")
        ovlp_chunk_abs_prefix = os.path.join(ovlp_chunk_dir, "ovlp")
        ovlp_fn = f"{ovlp_chunk_abs_prefix}-{my_chunk:02d}.dat"
        outputs[f'ovlp_{my_chunk:02d}'] = ovlp_fn

        wf.addTask(Task(
            script=shmr_ovlp_script,
            inputs=ovlp_in,
            outputs={'ovlp_out': ovlp_fn},
            parameters={
                'read_db_prefix': read_db_abs_prefix,
                'index_prefix': index_abs_prefix,
                'n_chunk': n_chunk,
                'my_chunk': my_chunk,
                'best_n_ovlp': int(args["--best_n_ovlp"]),
                'mc_lower': int(args["--mc_lower"]),
                'mc_upper': int(args["--mc_upper"]),
                'align_bandwidth': int(args["--aln_bw"]),
                'ovlp_upper': int(args["--ovlp_upper"]),
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = n_proc
    wf.refreshTargets()

    return outputs


def run_ovlp_to_ctg(wf, args, read_db_abs_prefix, read_db, ovlps):
    asm_script = """\
cat {params.ovlps} | shmr_dedup > preads.ovl; echo "-" >> preads.ovl
/usr/bin/time ovlp_to_graph.py --min_len {params.min_len} \
    --min_idt {params.min_idt} >& asm.log
/usr/bin/time graph_to_path.py >& to_path.log
/usr/bin/time path_to_contig.py {params.read_db_prefix} \
    p_ctg_tiling_path > {output.p_ctg} 2> to_contig.log
"""
    asm_dir = os.path.join(os.path.abspath(args["--output"]), "3-asm")
    ovlps_list = " ".join(sorted([v for v in ovlps.values()]))
    inputs = {}
    inputs.update(ovlps)
    inputs.update(read_db)
    outputs = {}
    outputs["p_ctg"] = os.path.join(asm_dir,  "p_ctg.fa")
    wf.addTask(Task(
        script=asm_script,
        inputs=inputs,
        outputs=outputs,
        parameters={
            'read_db_prefix': read_db_abs_prefix,
            'ovlps': ovlps_list,
            'min_len': int(args["--min_len"]),
            'min_idt': int(args["--min_idt"])
        },
        dist=Dist(NPROC=1, local=True)
    ))
    wf.max_jobs = 1
    wf.refreshTargets()
    return outputs


def run_cns(wf, args, read_db_abs_prefix, read_db,
            index_abs_prefix, read_index, p_ctg):
    mapping_nchunk = int(args["<mapping_nchunk>"])
    mapping_nproc = int(args["<mapping_nproc>"])
    cns_nchunk = int(args["<cns_nchunk>"])
    cns_nproc = int(args["<cns_nproc>"])
    sort_nproc = int(args["<sort_nproc>"])
    shimmer_k = int(args["--shimmer-k"])
    shimmer_w = int(args["--shimmer-w"])
    shimmer_r = int(args["--shimmer-r"])
    shimmer_l = int(args["--shimmer-l"])
    build_index_script = """\
echo {input.p_ctg} > p_ctg.lst

/usr/bin/time shmr_mkseqdb -p p_ctg \
    -d p_ctg.lst 2> build_p_ctg_db.log
"""

    build_index_script += f"""
/usr/bin/time shmr_index \
    -p p_ctg -t 1 -c 1 \
    -k {shimmer_k}\
    -w {shimmer_w}\
    -r {shimmer_r}\
    -l {shimmer_l}\
    -o p_ctg 2> build_p_ctg_index.log
"""
    cns_dir = os.path.join(os.path.abspath(args["--output"]), "4-cns")
    inputs = {}
    inputs.update(p_ctg)
    output_dir = os.path.join(cns_dir, "p_ctg_index")
    outputs = {}
    outputs["p_ctg_db"] = os.path.join(output_dir, "p_ctg.seqdb")
    outputs["p_ctg_idx"] = os.path.join(output_dir, "p_ctg.idx")
    if shimmer_l == 2:
        outputs["p_ctg_shmr_idx"] = os.path.join(output_dir, "p_ctg-L2-01-of-01.dat")
        p_ctg_idx_abs_prefix = os.path.join(output_dir, "p_ctg-L2")
    elif shimmer_l == 1:
        outputs["p_ctg_shmr_idx"] = os.path.join(output_dir, "p_ctg-L1-01-of-01.dat")
        p_ctg_idx_abs_prefix = os.path.join(output_dir, "p_ctg-L1")

    p_ctg_db_abs_prefix = os.path.join(output_dir, "p_ctg")

    wf.addTask(Task(
        script=build_index_script,
        inputs=inputs,
        outputs=outputs,
        parameters={
            'read_db_prefix': read_db_abs_prefix,
            'index_prefix': index_abs_prefix
        },
        dist=Dist(NPROC=1, local=True)
    ))
    wf.max_jobs = 1
    wf.refreshTargets()

    mapping_script = """\
/usr/bin/time shmr_map \
    -r {params.p_ctg_db_prefix} \
    -m {params.p_ctg_idx_prefix} \
    -p {params.read_db_prefix} \
    -l {params.index_prefix} \
    -t {params.n_chunk} -c {params.my_chunk}  > {output.readmap}
"""
    inputs = {}
    inputs["p_ctg_db"] = outputs["p_ctg_db"]
    inputs["p_ctg_idx"] = outputs["p_ctg_idx"]
    inputs["p_ctg_shmr_idx"] = outputs["p_ctg_shmr_idx"]
    inputs.update(read_db)
    inputs.update(read_index)
    outputs = {}
    for my_chunk in range(1, mapping_nchunk+1):
        mapping_chunk_dir = os.path.join(cns_dir, f"map-{my_chunk:02d}")
        mapping_chunk_abs_prefix = os.path.join(mapping_chunk_dir, "reads2ref")
        map_fn = f"{mapping_chunk_abs_prefix}-{my_chunk:02d}.dat"
        outputs[f'readmap_{my_chunk:02d}'] = map_fn

        wf.addTask(Task(
            script=mapping_script,
            inputs=inputs,
            outputs={'readmap': map_fn},
            parameters={
                'read_db_prefix': read_db_abs_prefix,
                'index_prefix': index_abs_prefix,
                'p_ctg_db_prefix': p_ctg_db_abs_prefix,
                'p_ctg_idx_prefix': p_ctg_idx_abs_prefix,
                'n_chunk': mapping_nchunk,
                'my_chunk': my_chunk,
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = mapping_nproc
    wf.refreshTargets()
    mapping_chunk_outputs = outputs
    map_files = " ".join(sorted([v for v in mapping_chunk_outputs.values()]))

    mapping_merge_script = """
mkdir -p {params.tmp_dir}
cat {params.map_files} | \
    sort -T {params.tmp_dir} -S 8g --parallel {params.sort_nproc}\
        -k 1 -g -k 2 -g > {output.merged_mapping_file}
"""
    mapping_merge_dir = os.path.join(cns_dir, "map-merge")
    merged_mapping_fn = os.path.join(mapping_merge_dir, "reads2ref_all.out")

    wf.addTask(Task(
        script=mapping_merge_script,
        inputs=mapping_chunk_outputs,
        outputs={"merged_mapping_file": merged_mapping_fn},
        parameters={
            'tmp_dir': os.path.join(cns_dir, "tmp"),
            'sort_nproc': sort_nproc,
            'map_files': map_files
        },
        dist=Dist(NPROC=1, local=True)
    ))

    cns_script = """\
/usr/bin/time pg_asm_cns.py {params.read_db_prefix} \
    {params.p_ctg_db_prefix} {input.merged_mapping_file} \
    {params.n_chunk} {params.my_chunk} > {output.cns_file} 2> cns.log
"""
    inputs.update({"merged_mapping_file": merged_mapping_fn})
    outputs = {}
    for my_chunk in range(1, cns_nchunk+1):
        cns_chunk_dir = os.path.join(cns_dir, f"cns-{my_chunk:02d}")
        cnd_chunk_abs_prefix = os.path.join(cns_chunk_dir, "p_ctg_cns")
        cns_fn = f"{cnd_chunk_abs_prefix}-{my_chunk:02d}.fa"
        outputs[f'cns_{my_chunk:02d}'] = cns_fn
        wf.addTask(Task(
            script=cns_script,
            inputs=inputs,
            outputs={"cns_file": cns_fn},
            parameters={
                'read_db_prefix': read_db_abs_prefix,
                'p_ctg_db_prefix': p_ctg_db_abs_prefix,
                'n_chunk': cns_nchunk,
                'my_chunk': my_chunk
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = cns_nproc
    wf.refreshTargets()

    gather_cns_script = """\
cat {params.cns_chunk_files} > {output.cns_file}
ln -sf {params.cns_merge_dir}/{output.cns_file} {params.workdir}
"""
    cns_chunk_files = " ".join(sorted([v for v in outputs.values()]))
    inputs = outputs
    cns_merge_dir = os.path.join(cns_dir, f"cns-merge")
    cns_fn = os.path.join(cns_merge_dir, "p_ctg_cns.fa")
    outputs = {"cns_file": cns_fn}
    wf.addTask(Task(
        script=gather_cns_script,
        inputs=inputs,
        outputs=outputs,
        parameters={
            'cns_chunk_files': cns_chunk_files,
            'workdir': os.path.abspath(args["--output"]),
            'cns_merge_dir': "./4-cns/cns-merge"
        },
        dist=Dist(NPROC=1, local=True)
    ))
    wf.max_jobs = 1
    wf.refreshTargets()

    return outputs


# Simple local-only submit-string.
submit = 'bash -C ${CMD} >| ${STDOUT_FILE} 2>| ${STDERR_FILE}'


def main(args):

    job_defaults = dict(
        njobs=1,
        NPROC=1,
        MB=24000,
        submit=submit,
        job_type='local',
        pwatcher_type='blocking',
    )

    wf = PypeProcWatcherWorkflow(
            job_defaults=job_defaults,
    )

    seq_dataset_lst = os.path.abspath(args["<reads.lst>"])
    read_db_abs_prefix, read_db = run_build_db(wf, args, seq_dataset_lst)
    LOG.info('Finished: {}'.format(read_db_abs_prefix))

    index_abs_prefix, read_idx = run_build_idx(wf, args,
                                               read_db_abs_prefix)
    LOG.info('Finished: {}'.format(index_abs_prefix))

    ovlp_in = {}
    ovlp_in.update(read_db)
    ovlp_in.update(read_idx)
    ovlp_out = run_overlapper(wf, args,
                              read_db_abs_prefix,
                              index_abs_prefix,
                              ovlp_in)

    LOG.info('Finished: {}'.format(ovlp_out))

    ctg_out = run_ovlp_to_ctg(wf, args,
                              read_db_abs_prefix,
                              read_db,
                              ovlp_out)
    LOG.info('Finished: {}'.format(ctg_out))

    if args['--with-consensus']:
        cns_out = run_cns(wf, args,
                          read_db_abs_prefix,
                          read_db,
                          index_abs_prefix,
                          read_idx,
                          ctg_out)
        LOG.info('Finished: {}'.format(cns_out))


if __name__ == "__main__":
    import pkg_resources
    short_doc = """
Peregrine
=========

Peregrine is a fast genome assembler for accurate long
reads (length > 10kb, accuraccy > 99%). It can assemble
a human genome from 30x reads within 20 cpu hours from
reads to polished consensus. It uses Sparse HIereachical
MimiMizER (SHIMMER) for fast read-to-read overlaps without
explicitly quadratic comparisions used in other OLC
assemblers.

Peregrine Assembler and SHIMMER Genome Assembly Toolkit
Copyright (c) 2019- by Jason, Chen-Shan, Chin

Peregrine Assembler and  SHIMMER Genome Assembly Toolkit
is licensed under a Creative Commons
Attribution-NonCommercial-ShareAlike 4.0 International
License.

You should have received a copy of the license along with
this work. If not, see
<http://creativecommons.org/licenses/by-nc-sa/4.0/>.

************************************************************
If you want to use it for any commericial purposes
(including promotion activity for a commerical product),
please contact Jason Chin for a commericial license.
************************************************************

run `pg_run.py -h` for help and other license information


"""
    sys.stderr.write(short_doc)

    sys.stderr.write('using {}\n\n'.format(pkg_resources.get_distribution('pypeflow')))

    logging.basicConfig(level=logging.INFO)
    args = docopt(__doc__, version=__version__)
    print(f"Peregrine Assembler ({__version__}) has been started with the following option:\n", args, file=sys.stderr)
    if not int(args["--shimmer-l"]) in (1, 2):
        print(f"\nERROR: <shimmer-l> should be 1 or 2, stopping\n", file=sys.stderr)
        sys.exit(128)

    main(args)

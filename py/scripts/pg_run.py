from pypeflow.simple_pwatcher_bridge import (PypeProcWatcherWorkflow, Dist)
from pypeflow.tasks import gen_task
from docopt import docopt
import peregrine
import logging
import os
import sys

__version__ = peregrine.__version__

__doc__ = """Peregrine Long Accurate Read Fast Genome Assembler

Usage:
  pg_run.py asm <reads.lst> <index_nchunk> <index_nproc>
                            <ovlp_nchunk> <ovlp_nproc>
                            <mapping_nchunk> <mapping_nproc>
                            <cns_nchunk> <cns_nproc>
                            <sort_nproc> [--with-consensus]
                            [--output <output>]
  pg_run.py (-h | --help)
  pg_run.py --verison

Options:
  -h --help     Show this help
  --version     Show version
  --output <output>     set output directory (will be created if not exist) [default: ./wd]
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

    build_idx = """
/usr/bin/time shmr_index\
    -p {params.read_db_prefix}\
    -t {params.n_chunk}\
    -c {params.my_chunk}\
    -o {params.index_prefix}
ln -s {params.index_prefix}* {params.index_dir}
"""
    index_dir = os.path.join(os.path.abspath(args["--output"]), "1-index")
    index_abs_prefix = os.path.join(index_dir, "shmr-L2")
    outputs = {}
    for my_chunk in range(1, n_chunk+1):
        index_chunk_dir = os.path.join(index_dir, f"chunk-{my_chunk:02d}")
        index_chunk_abs_prefix = os.path.join(index_chunk_dir, "shmr")
        index_L0_fn = f"{index_chunk_abs_prefix}-L0-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        index_L0_MC_fn = f"{index_chunk_abs_prefix}-L0-MC-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        index_L2_fn = f"{index_chunk_abs_prefix}-L2-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        index_L2_MC_fn = f"{index_chunk_abs_prefix}-L2-MC-{my_chunk:02d}-of-{n_chunk:02d}.dat"
        outputs[f'index_L2_{my_chunk:02d}'] = index_L2_fn
        outputs[f'index_L2_MC_{my_chunk:02d}'] = index_L2_MC_fn
        wf.addTask(Task(
            script=build_idx,
            inputs={
                'read_db': f"{read_db_abs_prefix}.seqdb",
                'seqidx': f"{read_db_abs_prefix}.idx"
            },
            outputs={
                'index_L0': index_L0_fn,
                'index_L0_MC': index_L0_MC_fn,
                'index_L2': index_L2_fn,
                'index_L2_MC': index_L2_MC_fn
            },
            parameters={
                'read_db_prefix': read_db_abs_prefix,
                'index_prefix': index_chunk_abs_prefix,
                'index_dir': index_dir,
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
                'my_chunk': my_chunk
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = n_proc
    wf.refreshTargets()

    return outputs


def run_ovlp_to_ctg(wf, args, read_db_abs_prefix, read_db, ovlps):
    asm_script = """\
cat {params.ovlps} | shmr_dedup > preads.ovl; echo "-" >> preads.ovl
/usr/bin/time ovlp_to_graph.py >& asm.log
/usr/bin/time graph_to_path.py >& to_path.log
/usr/bin/time path_to_contig.py {params.read_db_prefix} \
    p_ctg_tiling_path > {output.p_ctg} 2> to_contig.log
"""
    asm_dir = os.path.join(os.path.abspath(args["--output"]), "3-asm")
    ovlps_list = " ".join([v for k, v in ovlps.items()])
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
            'ovlps': ovlps_list
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
    build_index_script = """\
echo {input.p_ctg} > p_ctg.lst

/usr/bin/time shmr_mkseqdb -p p_ctg \
    -d p_ctg.lst 2> build_p_ctg_db.log

/usr/bin/time shmr_index -p p_ctg -t 1 -c 1 \
    -o p_ctg 2> build_p_ctg_index.log

"""
    cns_dir = os.path.join(os.path.abspath(args["--output"]), "4-cns")
    inputs = {}
    inputs.update(p_ctg)
    output_dir = os.path.join(cns_dir, "p_ctg_index")
    outputs = {}
    outputs["p_ctg_db"] = os.path.join(output_dir, "p_ctg.seqdb")
    outputs["p_ctg_idx"] = os.path.join(output_dir, "p_ctg.idx")
    outputs["p_ctg_L2_idx"] = os.path.join(output_dir, "p_ctg-L2-01-of-01.dat")
    p_ctg_db_abs_prefix = os.path.join(output_dir, "p_ctg")
    p_ctg_idx_abs_prefix = os.path.join(output_dir, "p_ctg-L2")

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
    inputs["p_ctg_L2_idx"] = outputs["p_ctg_L2_idx"]
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

    cns_script = """\
mkdir -p {params.tmp_dir}
cat {params.map_files} | \
    sort -T {params.tmp_dir} -S 8g --parallel {params.sort_nproc}\
        -k 1 -g -k 2 -g > reads2ref_all.out
/usr/bin/time cns_prototype.py {params.read_db_prefix} \
    {params.p_ctg_db_prefix} reads2ref_all.out \
    {params.n_chunk} {params.my_chunk} > {output.cns_file} 2> cns.log
"""
    map_files = " ".join([v for k, v in outputs.items()])
    inputs.update(outputs)
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
                'my_chunk': my_chunk,
                'tmp_dir': 'wd/tmp',
                'sort_nproc': sort_nproc,
                'map_files': map_files
            },
            dist=Dist(NPROC=1, local=True)
        ))

    wf.max_jobs = cns_nproc
    wf.refreshTargets()

    gather_cns_script = """\
cat {params.cns_chunk_files} > {output.cns_file}
"""
    cns_chunk_files = " ".join([v for k, v in outputs.items()])
    inputs = outputs
    cns_merge_dir = os.path.join(cns_dir, f"cns-merge")
    cns_fn = os.path.join(cns_merge_dir, "p_ctg_cns.fa")
    outputs = {"cns_file": cns_fn}
    wf.addTask(Task(
        script=gather_cns_script,
        inputs=inputs,
        outputs=outputs,
        parameters={
            'cns_chunk_files': cns_chunk_files
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

    nchunk = int(args['<index_nchunk>'])
    nproc = int(args['<index_nproc>'])
    index_abs_prefix, read_idx = run_build_idx(wf, args,
                                               read_db_abs_prefix)
    LOG.info('Finished: {}'.format(index_abs_prefix))

    ovlp_in = {}
    ovlp_in.update(read_db)
    ovlp_in.update(read_idx)
    nchunk = int(args['<ovlp_nchunk>'])
    nproc = int(args['<ovlp_nproc>'])
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
        mapping_nchunk = int(args["<mapping_nchunk>"])
        mapping_nproc = int(args["<mapping_nproc>"])
        cns_nchunk = int(args["<cns_nchunk>"])
        cns_nproc = int(args["<cns_nproc>"])
        sort_nproc = int(args["<sort_nproc>"])
        cns_out = run_cns(wf, args,
                          read_db_abs_prefix,
                          read_db,
                          index_abs_prefix,
                          read_idx,
                          ctg_out)
        LOG.info('Finished: {}'.format(cns_out))


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    args = docopt(__doc__, version=__version__)
    print(args)
    main(args)

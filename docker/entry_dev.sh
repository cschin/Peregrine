#!/bin/bash
. /opt/conda/etc/profile.d/conda.sh
conda activate peregrine
pg_run_dev.py $@

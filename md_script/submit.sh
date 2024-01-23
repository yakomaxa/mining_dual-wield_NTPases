#!/bin/bash
#PJM -L rscgrp=cxgfs-share
#PJM -L elapse=168:00:00
#PJM -j
#PJM -S
export OMP_NUM_THREADS=10
bash ./nvt_npt_shortmd.sh INPUTST

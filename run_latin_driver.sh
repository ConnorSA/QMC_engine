#!/bin/bash
python3 init_params.py
#export OMP_NUM_THREADS=4
./latin_driver_exec
python3 plotting.py
#!/bin/bash
python3 init_params.py
./bond_driver_exec
python3 ./energy_plotting.py
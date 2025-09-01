#!/bin/bash

set -e

python generate_all.py

for f in bolsig_script_*.txt; do
    echo " " | bolsigminus "$f"
done

for f in bolsig_result_*.txt; do
    python convert_bolsig.py "$f" "../${f/bolsig_result_/TD_}"
done

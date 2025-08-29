#!/bin/bash

set -e

bolsig=~/opt/sw/bolsigminus/bolsigminus

python generate_all.py

for f in bolsig_script_*.txt; do
    echo " " | $bolsig "$f"
done

for f in bolsig_result_*.txt; do
    python convert_bolsig.py "$f" "../${f/bolsig_result_/TD_}"
done

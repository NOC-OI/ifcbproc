#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR
cd ..
bash scripts/get-test-data.sh
mkdir -p testout
python3 cli.py features ./testdata/D20191211T034109_IFCB010.roi -o ./testout
python3 cli.py ecotaxa ./testdata/D20191211T034109_IFCB010.roi -o ./testout/D20191211T034109_IFCB010.zip

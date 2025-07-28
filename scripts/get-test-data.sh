#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR
cd ..
mkdir -p testdata
cd testdata
SHAOUT=($(cat *.adc *.hdr *.roi *.csv | sha256sum -b))
SHAOUT=${SHAOUT[0]}
if [ "$SHAOUT" = "80d497e4a8616ab11abbcbb6fe3beb5fcd44f7a0bb23e4da2a563661cebccfd4" ]; then
    echo "Test data OK, skipping download"
else
    echo "Test data hash ($SHAOUT) did not match expected output, redownloading data"
    rm *.adc *.hdr *.roi *.csv
    ROOT_URI="https://ifcb-data.whoi.edu/mvco/D20191211T034109_IFCB010"
    wget ${ROOT_URI}.adc
    wget ${ROOT_URI}.hdr
    wget ${ROOT_URI}.roi
    wget ${ROOT_URI}_features.csv
# src = https://ifcb-data.whoi.edu/timeline?dataset=mvco&bin=D20191211T034109_IFCB010
fi

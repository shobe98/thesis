#!/bin/bash

echo "1"

./riboscript.sh /data/caitken/files/CA1_1_S1_L001_R1_001.fq > ../logs/single_run.ribo.out.log 2> ../logs/single_run.ribo.err.log

echo "2"

./mrnascript.sh /data/caitken/files/CA1_3_S1_L001_R1_001.fq > ../logs/single_run.mrna.out.log 2> ../logs/single_run.mrna.err.log

echo "3"

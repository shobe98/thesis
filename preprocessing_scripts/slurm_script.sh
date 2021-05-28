#!/bin/bash

echo "Running the full Ingolia pipeline on June 12."
echo "from initial files to alignments"

./script.sh > ../logs/full_pipeline.out.log 2> ../logs/full_pipeline.err.log 

echo "Done!"

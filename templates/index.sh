#!/bin/bash
set -euo pipefail
salmon index --no-clip                  \
             --threads     ${task.cpus} \
             --kmerLen     ${kmerLen}   \
             --transcripts ${refSeqs}   \
             --index       ${outDir}

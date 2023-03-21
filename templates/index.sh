#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Runs salmon index for a given set of transcript sequences and specified params

set -euo pipefail
salmon index --no-clip                  \
             --threads     ${task.cpus} \
             --kmerLen     ${kmerLen}   \
             --transcripts ${refSeqs}   \
             --index       ${outDir}

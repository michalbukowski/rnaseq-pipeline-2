#!/bin/bash
set -euo pipefail
salmon quant --threads ${task.cpus} \
             --libType ${libType}   \
             --index   ${index}     \
             --mates1  ${reads_1}   \
             --mates2  ${reads_2}   \
             --output  ${outDir}

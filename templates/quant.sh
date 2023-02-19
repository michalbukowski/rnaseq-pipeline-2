#!/bin/bash
set -euo pipefail
salmon quant --threads ${task.cpus} \
             --libType ${libType}   \
             --index   ${index}     \
             --mates1  ${reads1}    \
             --mates2  ${reads2}    \
             --output  ${outDir}

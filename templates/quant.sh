#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Runs salmon quant for a given index, paired-read archives and specified params

set -euo pipefail
salmon quant --threads ${task.cpus} \
             --libType ${libType}   \
             --index   ${index}     \
             --mates1  ${reads1}    \
             --mates2  ${reads2}    \
             --output  ${outDir}

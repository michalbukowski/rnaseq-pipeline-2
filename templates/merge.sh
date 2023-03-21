#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io) under GPL-3.0 license

# Runs salmon quantmerge for a given set of salmon quant results and specified params

set -euo pipefail
salmon quantmerge --quants ${fsorted}  \
                  --names  ${colnames} \
                  --column ${metric}   \
                  --output ${outFile}

#!/bin/bash
set -euo pipefail
salmon quantmerge --quants ${fsorted}  \
                  --names  ${colnames} \
                  --column ${metric}   \
                  --output ${outFile}

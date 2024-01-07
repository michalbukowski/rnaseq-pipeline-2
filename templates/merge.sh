#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license. For detail see the following publication:
#
# Bukowski M, Kosecka-Strojek M, Wladyka B. Disruption of saoB and saoC genes
# in Staphylococcus aureus alters transcription of genes involved in amino acid
# metabolism and virulence. [awaiting publication]

# Runs salmon quantmerge for a given set of salmon quant results and specified params

salmon quantmerge --quants ${fsorted}  \
                  --names  ${colnames} \
                  --column ${metric}   \
                  --output ${outFile}

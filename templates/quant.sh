#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license. For detail see the following publication:
#
# Bukowski M, Kosecka-Strojek M, Wladyka B. Disruption of saoB and saoC genes
# in Staphylococcus aureus alters transcription of genes involved in amino acid
# metabolism and virulence. [awaiting publication]

# Runs salmon quant for a given index, paired-read archives and specified params

salmon quant --threads ${task.cpus} \
             --libType ${libType}   \
             --index   ${index}     \
             --mates1  ${reads1}    \
             --mates2  ${reads2}    \
             --output  ${outDir}

#!/bin/bash
# Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
# under GPL-3.0 license. For detail see the following publication:
#
# Bukowski M, Kosecka-Strojek M, Wladyka B. Disruption of saoB and saoC genes
# in Staphylococcus aureus alters transcription of genes involved in amino acid
# metabolism and virulence. [awaiting publication]

# Runs salmon index for a given set of transcript sequences and specified params

salmon index --no-clip                  \
             --threads     ${task.cpus} \
             --kmerLen     ${kmerLen}   \
             --transcripts ${refSeqs}   \
             --index       ${outDir}

#!/usr/bin/env nextflow
/* Created by Michal Bukowski (michal.bukowski@tuta.io, m.bukowski@uj.edu.pl)
   under GPL-3.0 license. For detail see the following publication:
   
   Bukowski M, Kosecka-Strojek M, Wladyka B. Disruption of saoB and saoC genes
   in Staphylococcus aureus alters transcription of genes involved in amino acid
   metabolism and virulence. [awaiting publication]
   
   Differential gene expression pipeline utilising Salmon and Bioconductor DESeq2
   for paired-end Illumina RNA-Seq squencing results.
   
   FASTQ files with paired reads are named according to a pattern:
   {strain}_{group}_{replica}_{reads}.fastq, eg. rn_wt_lg_1_R2.fastq
   for RN4220 (rn), wild type in logarithmic growth phase (wt_lg), replica 1 of 3,
   1 of 2 files with reads. R1 and R2 are known, the remianing
   are inferred experiment descriptions (defined in the config file). The read files
   are located in reads/ directory.
   
   In conda/ directory YML files describe conda environments used by the pipeline.
   Miniconda/Anaconda installation is a prerequisite. The default environment that
   must be installed prior running the pipeplie is described in rnaseq-base.yml.
   It contains nexflow 22.10.6 sra-tools 3.0.3. Run the following commands from
   the pipeline directory to recreate and activate the environment:
   
   conda env create -f conda/rnaseq-base.yml
   conda activate rnaseq-base
   
   The processes run scripts located in templates/ directory. To start the pipeline
   run the command from the pipeline directory:
   
   nextflow main.nf
   
   Detailed comments on what each process do and how are provided at the end of
   this file in the workflow{} block.
*/

process salmonIndex {
    conda    'envs/rnaseq-salmon.yml'
    publishDir 'output/salmon', mode: 'link'
    cpus     4
    
    input:
        val  kmerLen
        each strain
    
    exec:
        refSeqs = file("refs/${strain}.fna", checkIfExists: true)
        outDir  = "index/${strain}_index"
    
    output:
        tuple val(strain), path(outDir)
    
    script:
        template 'index.sh'
}

process salmonQuant {
    conda    'envs/rnaseq-salmon.yml'
    publishDir 'output/salmon/quant', mode: 'link'
    cpus     4
    
    input:
        tuple val(strain), val(group), val(replica), path(index)
        each  libType
    
    exec:
        reads1 = file("reads/${strain}_${group}_${replica}_R1.fastq.gz",
                       checkIfExists: true)
        reads2 = file("reads/${strain}_${group}_${replica}_R2.fastq.gz",
                       checkIfExists: true)
        outDir = "${libType}/${strain}_${group}_${replica}"
    
    output:
        tuple val(strain), val(group), val(libType), path(outDir)
    
    script:
        template 'quant.sh'
}

process salmonQuantMerge {
    conda    'envs/rnaseq-salmon.yml'
    publishDir 'output/salmon/quantmerge', mode: 'link'
    
    input:
        tuple val(exp), val(strain), val(libType), val(groups), path(files)
        each  metric
    
    exec:
        outFile  = "${libType}/${metric}/${strain}_${groups.join('_')}_${libType}_${metric}.tsv"
        fsorted  = []
        colnames = []
        groups.each{ group ->
            sublist = []
            num     = 1
            files.each{ fpath ->
                if (fpath.getName().indexOf(group) != -1) {
                    sublist  << fpath
                    colnames << "${group}_${num}"
                    num++
                }
            }
            fsorted.addAll(sublist.sort())
        }
        fsorted  = fsorted.join(' ')
        colnames = colnames.join(' ')
    
    output:
        tuple val(exp), val(strain), val(libType), val(metric), val(groups), path(outFile)
    
    script:
        template 'merge.sh'
}

process calculateDGE {
    conda    'envs/rnaseq-data.yml'
    publishDir 'output/DGE', mode: 'link'
    
    input:
        tuple val(exp), val(strain), val(libType), val(metric), val(groups), path(file)
        each  alpha
    
    exec:
        outFile = "${libType}/${strain}_${groups.join('_')}_${libType}_${alpha}.tsv"
    
    output:
        tuple val(exp), val(strain), val(libType), val(groups), path(outFile)
    
    script:
        template 'dge.r'
}

workflow {
    /* STAGE 1 -- create indices for transcriptomes (salmon index)
    */
    /*
       From the list of dictionaries describing experiments (defined in the config file)
       fetch a set of strain names, for which salmon index must be created, and
       convert it into a channel.
    */
    strains = Channel.fromList(
        params.experiments
        .collect{ exp -> exp['strain'] }
        .unique()
    )
    /* Set kmerLen to a fixed value of 15.
    */
    kmerLen = 15
    /* Run salmonIndex process for given kmerLen and a the channel of strains.
       It returns a channel of tupels: (strain, outDir), where outDir is a path
       to the directory with index for a given strain.
    */
    salmonIdx = salmonIndex(kmerLen, strains)
    
    
    /* STAGE 2 -- calculate TPMs and NumReads for all replicas (salmon quant)
    */
    /* Based on the list of dictionaries describing experiments (defined
       in the config file) create a nested list describing each replica that
       samlon quant must be run for. Each replica is described by
       strain, group (experimental group), and a number {1..n}.
       These correspond to file names: strain_group_replica_RN.fastq.gz, where
       N is 1 or 2 and refers to reads 1. General readsParams list structure:
       [ [strain, group, replica],
         ...
       ]
    */
    readsParams = [] as Set
    params.experiments.each{ exp ->
        strain = exp['strain']
        exp['groups'].each{ group ->
            for(replica in 1..exp['repNo']) {
                readsParams << [strain, group, replica]
            }
        }
    }
    /* Turn the list of libTypes for salmon quant (defined in the config file) into
       a channel.
    */
    libTypes = Channel.fromList(params.libTypes)
    /* Turn readsParams nested list describing replicas into a channel and combine
       it with salmonIdx output channel by the first value, which is a strain name.
       It results in a channel of tupels: (strain, group, replica, outDir), where
       outDir is a path to the directory with index for a given strain.
    */
    allParams = Channel
        .fromList(readsParams)
        .combine(salmonIdx, by: 0)
    /* Run salmonQuant process for given (strain, group, replica, outDir) tuples,
       combined with all libTypes of interest. For each analysed replica,
       it returns a channel of tupels: (strain, group, libType, outDir), where
       outDir is a path to the directory with quant files for a given replica.
    */
    salmonQnt = salmonQuant(allParams, libTypes)
    
    
    /* STAGE 3 -- for each experiment merge all replica TPMs and NumReads (salmon quantmerge)
    */
    /* Based on the list of dictionaries describing experiments (defined
       in the config file) create a nested list describing each experiment
       to be able to assign replica quant files to experiments. Give each
       experiment a label (exp), which is the index of an experiment
       in the said list (i):
       [ [strain, group, exp],
         ...
       ]
    */
    expParams = []
    params.experiments.eachWithIndex{ exp, i ->
        strain = exp['strain']
        groups = exp['groups']
        groups.each{ group ->
            expParams << [strain, group, i]
        }
    }
    /* Turn the list of metrics (TPM and/or NumReads, defined in the config file)
       for salmon quantmerge into a channel.
    */
    metrics = Channel.fromList(params.metrics)
    /* Turn expParams nested list describing experiments into a channel and combine
       it with salmonQnt output channel by the first and the second value,
       which are a strain name and an experimental group. It results in
       a channel of tupels: (strain, group, exp, libType, outDir), where
       outDir is a path to the directory with quant files for a given replica.
       Then group channel items in respect to experiment (exp) and libType.
       The obtained structure should be:
       ( [str_1, str_1, ...], [gr_1, gr_1, ..., gr_2, gr_2...], exp, libType,
         [outDir_1, outDir_2, outDir_3, ...] )
       Then rearange it to obtain a channel of items of the following structure:
       (exp, strain, libType, [gr_1, gr_2], [outDir_1, outDir_2, outDir_3, ...]).
       List of groups within experiments are probed directly from config data
       to preserve the order of groups as the first is treated as a reference one.
       In this way, each channel item describes one experiment, in particular
       it gives a list of experimental groups and a list of paths to all directories
       with all files that are associated with all replicas assigned to the experiment.
    */
    expParams = Channel.fromList(expParams)
    .combine(salmonQnt, by: [0, 1])
    .groupTuple(by:[2,3])
    .map{
        tuple(it[2], it[0].first(), it[3], params.experiments[it[2]]['groups'], it[4])
    }
    /* Run salmonQuantMerge process for given (exp, strain, libType,
       [gr_1, gr_2], [outDir_1, outDir_2, outDir_3, ...]) tuples,
       combined with all metrics of interest. For each analysed experiment,
       it returns a channel of tupels: (exp, strain, libType, metric, groups, outFile),
       where outDir is a path to a TSV file with merged values of a given metrics for
       all transcripts across all replicas associated with an experiment.
       The columns are named by experimental group labels (group) and
       replica numbers (replica).
    */
    salmonQntM = salmonQuantMerge(expParams, metrics)
    
    
    /* STAGE 4 -- for each experiment perform DGE analysis (DESeq2)
    */
    /* Turn the list of alpha thresholds (defined in the config file) for DESeq2
       into a channel. In the salmonQntM output channel leave items related to NumReads
       metric, which DESeq2 DGE analysis is based on.
    */
    alpha      = Channel.fromList(params.alpha)
    salmonQntM = salmonQntM.filter{ it[3] == 'NumReads' }
    /* Run calculateDGE process for (exp, strain, libType, metric, groups, outFile),
       where outDir is a path to a TSV file with merged values of a given metrics for
       all transcripts across all replicas associated with an experiment, and the only
       value of metric is "NumReads". Run the analysis for each alpha threshold.
       In output file names, the group that is mentioned first is the reference one.
    */
    calculateDGE(salmonQntM, alpha)
}


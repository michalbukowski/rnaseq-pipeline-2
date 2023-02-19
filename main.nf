#!/usr/bin/env nextflow

process salmonIndex {
    conda 'conda/rnaseq-salmon.yml' //'/store/miniconda3/envs/rnaseq-salmon'//'bioconda::salmon=1.9.0'
    cpus 4
    storeDir 'output/salmon'
    
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
    conda 'conda/rnaseq-salmon.yml' //'/store/miniconda3/envs/rnaseq-salmon'
    cpus 4
    storeDir 'output/salmon/quant'
    
    input:
        tuple val(strain), val(group), val(replica), path(index)
        each  libType
    
    exec:
        reads1 = file("reads/${strain}_${group}_${replica}_R1.fastq.gz",
                       checkIfExists: true)
        reads2 = file("reads/${strain}_${group}_${replica}_R2.fastq.gz",
                       checkIfExists: true)
        outDir  = "${libType}/${strain}_${group}_${replica}"
    
    output:
        tuple val(strain), val(group), val(libType), path(outDir)
    
    script:
        template 'quant.sh'
}

process salmonQuantMerge {
    conda 'conda/rnaseq-salmon.yml'
    storeDir 'output/salmon/quantmerge'
    
    input:
        tuple val(exp), val(strain), val(libType), val(groups), path(files)
        each  metric
    
    exec:
        outFile  = "${libType}/${metric}/${strain}_${groups.join('_')}_${libType}_${metric}.tsv"
        fsorted  = []
        colnames = []
        groups.each{ group ->
            sublist = []
            count = 1
            files.each{ fpath ->
                if (fpath.getName().indexOf(group) != -1) {
                    sublist  << fpath
                    colnames << "${group}_${count}"
                    count++
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
    conda 'conda/rnaseq-data.yml'
    storeDir 'output/DGE'
    
    input:
        tuple val(exp), val(strain), val(libType), val(metric), val(groups), path(file)
        each  alpha
    
    exec:
        outFile  = "${libType}/${strain}_${groups.join('_')}_${libType}_${alpha}.tsv"
    
    output:
        tuple val(exp), val(strain), val(libType), val(groups), path(outFile)
    
    script:
        template 'dge.r'
}

workflow {
    strains = params.experiments.collect{ exp -> exp['strain'] }.unique()
    kmerLen   = 15
    salmonIdx = salmonIndex(kmerLen, strains)
    
    readsParams = [] as Set
    params.experiments.each{ exp ->
        strain = exp['strain']
        exp['groups'].each{ group ->
            for(replica in 1..exp['repNo']) {
                readsParams << [strain, group, replica]
            }
        }
    }
    libTypes  = Channel.fromList(params.libTypes)
    allParams = Channel
        .fromList(readsParams)
        .combine(salmonIdx, by: 0)
    salmonQnt = salmonQuant(allParams, libTypes)
    
    expParams = []
    params.experiments.eachWithIndex{ exp, i ->
        strain = exp['strain']
        groups = exp['groups']
        groups.each{ group ->
            expParams << [strain, group, "exp_${i}"]
        }
    }
    metrics = Channel.fromList(params.metrics)
    expParams = Channel.fromList(expParams)
    .combine(salmonQnt, by: [0, 1])
    .groupTuple(by:[2,3])
    .map{
        tuple(it[2], it[0].first(), it[3], it[1].unique(), it[4])
    }
    salmonQntM = salmonQuantMerge(expParams, metrics)
    
    alpha = Channel.fromList(params.alpha)
    salmonQntM = salmonQntM.filter{ it[3] == 'NumReads' }
    calculateDGE(salmonQntM, alpha)
}


## Simple RNA-Seq pipeline (rnaseq-pipeline-2)
Differential gene expression pipeline utilising Salmon and Bioconductor DESeq2 for Illumina RNA-Seq sequencing results.

### 1. Environment
Miniconda/Anaconda installation is a prerequisite for running the pipeline. The pipeline was created and tested using the following conda environments and crucial packages:
- rnaseq-base:
  - sra-tools (3.0.3), only for test data acquisition
  - nextflow (22.10.6)
- rnaseq-salmon:
  - salmon (1.9.0)
- rnaseq-data:
  - python (3.11.0)
  - r-base (4.2.2)
  - bioconductor-deseq2 (1.38.0)

The first environment is required for obtaining the test data and running the pipeline. Remaining environemnts (see YML files in `conda/` directory) are installed by the pipeline. Run the following commands from the pipeline directory to recreate and activate the `rnaseq-base` environment:
```
conda env create -f conda/rnaseq-base.yml
conda activate rnaseq-base
```

### 2. Test data acquisition
To download test data (read archives) from NCBI SRA run the `test_data.sh` script. This script will fetch 12 read archives to `reads/` directory (~3.7 GB). It may take a while. Any problems are usually related to the accesibility of NCBI servers. If something does not work, try again later.
```
./test_data.sh
```

### 3. Directory structure and pipeline files
Names of the FASTQ files with reads in [NCBI BioProject](https://www.ncbi.nlm.nih.gov/bioproject/) [`PRJNA798259`](https://www.ncbi.nlm.nih.gov/bioproject?term=PRJNA798259%5BProject%20Accession%5D) follow the pattern, which is required by the pipeline: `{strain}_{genotype}_{replica}_{reads}.fastq`, e.g. `rn_wt_lg_1_R1.fastq`. Regarding the strain, in the project data there are files only for `rn` (_Staphylococcus aureus_ RN4220). There are two genotypes: `wt` (wild type, the reference group), `mt` (mutant, _saoC_ gene mutant), sampled in two settings/conditions: `lg` (logarithmic growth phase) and `lt` (late growth phase). That two pairs of groups (`wt_lg`, `mt_lg` and `wt_lt`, `mt_lt`). For each there are 3 replicas (`1` -- `3`). Reads of each are written to two files: `R1` and `R2`. All in all, there are 24 files.

The pipeline utilises the following directory structure:
```
your_pipeline_location/
├── conda/
│   ├── rnaseq-base.yml
│   ├── rnaseq-data.yml
│   └── rnaseq-salmon.yml
├── refs/
│   └── rn.fna
├── templates/
│   ├── dge.r
│   ├── index.sh
│   ├── merge.sh
│   └── quant.sh
├── reads/
├── work/
├── output/
├── nextflow.config
└── main.nf
```
In the working directory you can find the `main.nf` file that describes the pipeline. The experimental design is described in `nextflow.config` file. Template scripts are in `template/`. In `input/refs/` you will have `rn.fna` multiple FASTA file with reference transcript sequences for Staphylococcus aureus RN4220. In `input/reads/` you will find 24 read archives downloaded with `test_data.sh` script. Directories `work/` and `output/` will be created automatically once the pipeline is run. The latter will contain final output from each stage of the analysis.

### 4. Pipeline architecture
The pipline described in the `main.nf` Nexflow file implements the following processes:
1. **salmonIndex** -- preparation of an index of reference transcript sequences with Salmon.
1. **salmonQuant** -- read mapping and counts with Salmon.
1. **salmonQuantMerge** -- gathering TPM and NumRead values for replicas assigned to one experiment.
1. **calculateDGE** -- differential gene expression with utilising DESeq2 library for each experiment.

More detailed description on how the pipeline works you will find in comments included in the `main.nf` and `nexflow.config` files.

### 5. Running the pipeline
To start the pipeline, once you have read data and experiment design description done as well as rnaseq-base conda environment activated, run the following command from the pipeline directory:
```
nextflow main.nf
```

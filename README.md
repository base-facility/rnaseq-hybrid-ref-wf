## Hybrid reference generation pipeline

Generate hybrid reference genome, transcriptome, and annotation files based on a primary reference (e.g., Human reference files) from a list of input FASTA files. Quantify RNA-seq data using Salmon. Align data to hybrid reference using STAR. 

**Input:**

- Contig sequence fasta file(s) in gzip format (`.fa.gz`)
- Primary reference genome (e.g., Human reference files) (`.fa.gz`)
- Reference annotation file (GTF) (`.gtf.gz`)
- Reference transcriptome file (`.fa.gz`)
- RNA-seq read 1 dataset (`.fq.gz`)
- RNA-seq read 2 dataset (`.fq.gz`)

**System requirements:**

- Snakemake v7.0 and above

### Usage

Populate your input file path in `config/config.yaml`, `config/reads.tsv`, `config/units.tsv`, and `config/samples.tsv` as described below:

- Modify `config/config.yaml` to add your reference files (i.e., genome, transcriptome, and annotation of sample's primary species). Genome reference: `ref_fa`, Annotation: `ref_gtf`, Transcript reference sequences: `ref_tx`.
- Modify `config/reads.tsv` to add your RNA-seq read 1 and 2 datasets. Note that `sample_id` column values will be utilised for naming output files.
- Modify `config/units.tsv` to add your synthetic contigs in fasta `.fa` format file path. Each contig sequence could be passed as a single fasta file with a given ID and name in `config/samples.tsv`.

If `snakemake` is not installed on your local system/server, you can install it using mamba (highly recommended) or conda:

```
mamba create -n snakemake -c conda-forge -c bioconda snakemake;
conda activate snakemake
```

Then in activated `snakemake` environment execute the following command from [pipeline_root](#pipeline-directory-structure):

```
snakemake -j 8 --use-conda --rerun-incomplete
```

This Snakemake pipeline takes advantage of isolated conda environments generated for each processing step and therefore, Snakemake is the only pre-installed required package. All tools will be downloaded and installed automatically by Snakemake and conda. `j` specifies allocated number of threads for running the pipeline and can be adjusted to control system workload.

#### Pipeline directory structure

User specific inputs should be provided via `config/` files. Input reference files (e.g., species specific annotation and genome fasta files) must be specified by user in `config/config.yaml`. Sample metadata including a `sample_id` must be supplied in `config/samples.tsv`. Sample input files (i.e., genome mapped sorted and indexed BAM) must be supplied in `config/units.tsv`. **Note that `sample_id` is used as primary key and therefore, must be the same value for a given sample in both `config/samples.tsv` and `config/units.tsv` files.**

Pipeline output files are stored in `results/` and each sample output is prefixed with the supplied `sample_id` value in `config/reads.tsv`.

```
pipeline_root
├── workflow
│   ├── rules 
|   │   ├── commons.smk
|   │   └── process.smk
│   ├── envs
│   ├── scripts
|   └── Snakefile
├── config
│   ├── config.yaml (modify config and provide required reference input files)
│   ├── reads.tsv (add your RNA-seq read 1 and 2 file path per sample) 
│   ├── samples.tsv (add your synthetic construct meta data)
|   └── units.tsv (add your synthetic construct input fasta file path)
├── results (pipeline output)
└── resources
```
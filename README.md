## workflow title

Generate hybrid reference genome, transcriptome, and annotation files based on a primary reference (e.g., Human reference files) from a list of input FASTA files.

**Input:**

- Contig sequence fasta file(s) in gzip format (`.fa.gz`)
- Primary reference genome (e.g., Human reference files) (`.fa.gz`)
- Reference annotation file (GTF) (`.gtf.gz`)
- Reference transcriptome file (`.fa.gz`)

**System requirements:**

- Snakemake v7.0 and above

### Usage

Populate your input file path in `config/config.yaml`, sample names and input files in `config/samples.yaml` and `config/units.yaml`.

If `snakemake` is not installed on your local system, simplest way to install `snakemake` is by creating a new conda environment with an isolated Snakemake installation:

```
mamba create -n snakemake -c conda-forge -c bioconda snakemake;
mamba activate snakemake
```

Then in activated `snakemake` environment execute the following command from [pipeline_root](#pipeline-directory-structure):

```
snakemake -j 16 --use-conda --rerun-incomplete
```

This Snakemake pipeline takes advantage of isolated conda environments generated for each processing step and therefore, Snakemake is the only pre-installed required package. All tools will be downloaded and installed automatically by Snakemake and conda. `j` specifies allocated number of threads for running the pipeline and can be adjusted to control system workload.

#### Pipeline directory structure

User specific inputs should be provided via `config/` files. Input reference files (e.g., species specific annotation and genome fasta files) must be specified by user in `config/config.yaml`. Sample metadata including a `sample_id` must be supplied in `config/samples.tsv`. Sample input files (i.e., genome mapped sorted and indexed BAM) must be supplied in `config/units.tsv`. **Note that `sample_id` is used as primary key and therefore, must be the same value for a given sample in both `config/samples.tsv` and `config/units.tsv` files.**

Pipeline output files are stored in `results/` and each sample output is prefixed with the supplied `sample_id` value.

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
│   └── samples.tsv (add your sample meta data)
|   └── units.tsv (add your sample input file path)
├── results (pipeline output)
└── resources
```
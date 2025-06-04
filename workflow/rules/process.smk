# Concatenate reference sequences
rule concatenate:
    input: all_files
    output: 
        genome = "results/cat/hybrid_ref.fa.gz",
        transcriptome = "results/cat/hybrid_ref_tx.fa.gz"
    log: "logs/cat.log"
    benchmark: "benchmarks/cat.benchmark"
    params:
        ref_fa = config['ref_fa'],
        ref_tx = config['ref_tx']
    shell:
        '''
        cat {params.ref_fa} {input} > {output.genome};
        cat {params.ref_tx} {input} > {output.transcriptome} 2> {log}
        '''

# Generate GTF annotation files for each construct
rule gen_gtf:
    input: all_files
    output: "results/gtf/hybrid.gtf.gz"
    log: "logs/gen_gtf.log"
    benchmark: "benchmarks/gen_gtf.benchmark"
    conda: "../envs/gffutils.yaml"
    params:
        ref_gtf = config['ref_gtf'],
        contigs_gtf = "results/gtf/contigs.gtf"
    shell:
        '''
        python workflow/scripts/genSingleFeatureGTF.py -o {params.contigs_gtf} {input};
        bgzip {params.contigs_gtf};
        cat {params.ref_gtf} {params.contigs_gtf}.gz > {output} 2> {log}
        '''

# Generate Salmon index
rule salmon_idx:
    input: "results/cat/hybrid_ref_tx.fa.gz"
    output: "results/salmon_idx/info.json"
    log: "logs/salmon_idx.log"
    benchmark: "benchmarks/salmon_idx.benchmark"
    conda: "../envs/salmon.yaml"
    threads: 16
    params:
        idx_dir = "results/salmon_idx"
    shell:
        '''
        salmon index \
        -t {input} \
        -i {params.idx_dir} \
        -p {threads} \
        --gencode
        '''

# Sequencing QC of full dataset
rule fastqc:
    input:
        r1 = lambda wildcards: read_1[wildcards.id],
        r2 = lambda wildcards: read_2[wildcards.id]
    output:
        done = "results/qc/{id}_fastqc.done",
    conda: "../envs/fastqc.yaml"
    log: "logs/fastqc_{id}.log"
    benchmark: "benchmarks/fastqc_{id}.benchmark"
    threads: 4
    shell:
        '''
        fastqc -t {threads} \
        -o results/qc/ \
        {input.r1} \
        {input.r2} 2> {log};
        touch {output.done}
        '''

# Subsample reads
rule subsample:
    input:
        r1 = lambda wildcards: read_1[wildcards.id],
        r2 = lambda wildcards: read_2[wildcards.id],
    output:
        r1 = "results/subsample/{id}_sub_R1.fastq.gz",
        r2 = "results/subsample/{id}_sub_R2.fastq.gz"
    threads: 8
    conda: "../envs/seqtk_samtools.yaml"
    log: "logs/subsample_{id}.log"
    benchmark: "benchmarks/subsample_{id}.benchmark"
    params:
        size = 20000000
    shell:
       '''
       seqtk sample -s80 {input.r1} {params.size} | bgzip -@ {threads} > {output.r1};
       seqtk sample -s80 {input.r2} {params.size} | bgzip -@ {threads} > {output.r2} 2> {log}
       '''

# Quality filter sub-sampled reads
rule fastp_trim:
    input:
        r1 = "results/subsample/{id}_sub_R1.fastq.gz",
        r2 = "results/subsample/{id}_sub_R2.fastq.gz"
    output:
        html = "results/fastp/{id}_fastp.html",
        json = "results/fastp/{id}_fastp.json",
        tr1 = "results/qt_reads/{id}_R1_qt.fastq.gz",
        tr2 = "results/qt_reads/{id}_R2_qt.fastq.gz"
    threads: 8
    conda: "../envs/fastp.yaml"
    log: "logs/{id}_fastp_trim.log"
    benchmark: "benchmarks/{id}_fastp_trim.benchmark"
    shell:
        '''
        fastp --thread {threads} \
        -i {input.r1} \
        -I {input.r2} \
        -h {output.html} \
        -j {output.json} \
        -o {output.tr1} \
        -O {output.tr2} 2> {log}
        '''

rule salmon:
    input:
        tr1 = "results/qt_reads/{id}_R1_qt.fastq.gz",
        tr2 = "results/qt_reads/{id}_R2_qt.fastq.gz",
        index = "results/salmon_idx/info.json"
    output:
        quant = "results/salmon/{id}/quant.sf"
    threads: 12
    conda: "../envs/salmon.yaml"
    log: "logs/salmon_{id}.log"
    benchmark: "benchmarks/salmon_{id}.benchmark"
    params:
        idx = "results/salmon_idx",
        outdir = "results/salmon/{id}/"
    shell:
        '''
        salmon quant \
        -i {params.idx} \
        -l A \
        -1 {input.tr1} \
        -2 {input.tr2} \
        --validateMappings \
        --seqBias \
        --gcBias \
        --posBias \
        --threads {threads} \
        -o {params.outdir} 2> {log}
        '''

# Generate STAR index
rule star_idx:
    input: 
        ref_genome = "results/cat/hybrid_ref.fa.gz",
        ref_gtf = "results/gtf/hybrid.gtf.gz"
    output: "results/star_idx/sjdbList.out.tab"
    log: "logs/star_idx.log"
    benchmark: "benchmarks/star_idx.benchmark"
    conda: "../envs/star.yaml"
    threads: 16
    params:
        idx_dir = "results/star_idx",
        unzip_fa = "results/cat/hybrid_ref.fa",
        unzip_gtf = "results/gtf/hybrid.gtf"
    shell:
        '''
        gunzip -k {input.ref_genome};
        gunzip -k {input.ref_gtf};
        STAR \
        --runThreadN {threads} \
        --runMode genomeGenerate \
        --genomeDir {params.idx_dir} \
        --genomeFastaFiles {params.unzip_fa} \
        --sjdbGTFfile {params.unzip_gtf} \
        --sjdbOverhang 100;
        rm {params.unzip_fa};
        rm {params.unzip_gtf} 2> {log}
        '''

# Align reads with STAR
rule star:
    input:
        tr1 = "results/qt_reads/{id}_R1_qt.fastq.gz",
        tr2 = "results/qt_reads/{id}_R2_qt.fastq.gz",
        star_idx = "results/star_idx/sjdbList.out.tab"
    output: "results/star/{id}.SJ.out.tab"
    threads: 16
    conda: "../envs/star.yaml"
    log: "logs/star_{id}.log"
    benchmark: "benchmarks/star_{id}.benchmark"
    params:
        star_idx = "results/star_idx",
        prefix = "results/star/{id}.",
        tmp = "tmp/{id}"
    shell:
        '''
        STAR \
        --runThreadN {threads} \
        --genomeDir {params.star_idx} \
        --readFilesIn {input.tr1} {input.tr2} \
        --readFilesType Fastx \
        --readFilesCommand gunzip -c \
        --outFilterType BySJout \
        --outFilterMultimapNmax 20 \
        --alignSJoverhangMin 8 \
        --alignSJDBoverhangMin 1 \
        --outFilterMismatchNmax 999 \
        --outFilterMismatchNoverReadLmax 0.04 \
        --alignIntronMin 20 \
        --alignIntronMax 1000000 \
        --alignMatesGapMax 1000000 \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM Unsorted 2> {log}
        '''

# Sort and index BAM files
rule sort_bam:
    input: "results/star/{id}.SJ.out.tab"
    output: "results/star/{id}.sorted.bam"
    threads: 8
    conda: "../envs/samtools.yaml"
    log: "logs/{id}_sort_bam.log"
    benchmark: "benchmarks/{id}_sort_bam.benchmark"
    params:
        bam = "results/star/{id}.Aligned.out.bam"
    shell:
        '''
        samtools sort \
        -@ {threads} \
        -O BAM \
        -o {output} \
        {params.bam};
        samtools index \
        -@ {threads} \
        {output};
        rm {params.bam} 2> {log}
        '''

# Generate synthetic construct coverage
rule synth_coverage:
    input: "results/star/{id}.sorted.bam"
    output: "results/samtools_depth/{id}_samtools.depth.tsv"
    conda: "../envs/samtools.yaml"
    log: "logs/{id}_synth_coverage.log"
    benchmark: "benchmarks/{id}_synth_coverage.benchmark"
    threads: 8
    params:
        region = " ".join(samples['name'].to_list())
    shell:
        '''
        for region in {params.region}; do
            samtools depth -@ {threads} -a -r $region {input} >> {output}
        done 2> {log}
        '''

# Generate report for synthetic contruct quantitation
rule report:
    input: 
        salmon = "results/salmon/{id}/quant.sf",
        coverage = "results/samtools_depth/{id}_samtools.depth.tsv"
    output: "results/report/{id}.html"
    conda: "../envs/r.yaml"
    log: "logs/{id}_report.log"
    benchmark: "benchmarks/{id}_report.benchmark"
    threads: 32
    params:
        sample_id = "{id}",
        synth = samples['name'].to_list(),
        outdir = "results/report",
        outfile = "{id}.html",
        knitroot = cwd,
        txid_mapping = "resources/txid_mapping.tsv"
    shell:
        '''
        Rscript workflow/scripts/report.R \
        --id {params.sample_id} \
        --quant {input.salmon} \
        --txid_mapping {params.txid_mapping} \
        --coverage {input.coverage} \
        --output_dir {params.outdir} \
        --output_file {params.outfile} \
        --knit_root {params.knitroot} \
        --synth {params.synth}
        '''
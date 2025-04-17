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

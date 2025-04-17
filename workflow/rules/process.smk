# Concatenate reference sequences
rule concatenate:
    input: all_files
    output: 
        genome = "results/cat/hybrid_ref.fa",
        transcriptome = "results/cat/hybrid_ref_tx.fa"
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

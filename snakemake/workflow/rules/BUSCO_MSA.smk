rule BUSCO_MSA:
    """
    Perform MSA on BUSCO OGs
    of each class
    """
    input:
        os.path.join(out_dir, 'all_species', 'BUSCO_sequences', '{class}/{BUSCO}.faa')
    output:
        os.path.join(out_dir, 'all_species', 'BUSCO_sequences', '{class}/{BUSCO}_MSA.faa')
    log:
        os.path.join(logs_dir, 'BUSCO_MSA', '{class}.{BUSCO}.BUSCO_MSA.log')
    conda:
        os.path.join(envs_dir, 'mafft.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        mafft --auto --thread 1 {input} > {output}
        """

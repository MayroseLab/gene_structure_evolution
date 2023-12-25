rule collect_class_BUSCO_sequences:
    """
    Collect protein sequences of
    BUSCO OGs from all class  species
    into fasta files
    """
    input:
        sp_list = os.path.join(out_dir, 'all_species', '{class_}.list'),
        stats = expand(os.path.join(out_dir, 'per_species', '{species}/BUSCO.stats'), species=species_list)
    output:
        expand(os.path.join(out_dir, 'all_species', 'BUSCO_sequences', '{{class_}}/{BUSCO}.faa'), BUSCO=busco_list)
    log:
        os.path.join(logs_dir, 'collect_class_BUSCO_sequences', '{class_}.collect_class_BUSCO_sequences.log')
    params:
        in_dir = os.path.join(out_dir, 'per_species'),
        out_dir = lambda wc: os.path.join(out_dir, 'all_species', 'BUSCO_sequences', wc.class_),
        collect_script = os.path.join(scripts_dir, 'collect_BUSCO_proteins.py')
    conda:
        os.path.join(envs_dir, 'python_analyze.yml')
    resources:
        mem_mb = 10000
    shell:
        """
        python {params.collect_script} {input.sp_list} {params.in_dir} {params.out_dir}
        """

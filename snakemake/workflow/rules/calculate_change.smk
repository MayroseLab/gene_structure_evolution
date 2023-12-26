rule calculate_change:
    """
    Calculate the ammount of
    sequence and intron content
    change for a given BUSCO in a
    class
    """
    input:
        config['species_tree'],
        os.path.join(out_dir, 'all_species', 'BUSCO_trees', '{class_}/{BUSCO}.treefile'),
        #os.path.join(out_dir, 'all_species', 'BUSCO_intron_frac.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_total_intron_len.stats')
    output:
        os.path.join(out_dir, 'all_species', 'BUSCO_change/{class_}/{BUSCO}_change.tsv')
    log:
        os.path.join(logs_dir, 'calculate_change', '{class_}.{BUSCO}.calculate_change.log')
    params:
        calculate_change_script = os.path.join(scripts_dir, 'calculate_change.R')
    conda:
        'gene_structure_evolution'
    resources:
        mem_mb = 32000
    shell:
        """
        Rscript {params.calculate_change_script} {input} {wildcards.BUSCO} {output}
        """

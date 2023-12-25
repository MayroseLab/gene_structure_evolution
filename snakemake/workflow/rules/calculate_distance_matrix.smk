rule calculate_distance_matrix:
    """
    calculate pairwise KS distances
    between species, based on
    several gene structure features
    """
    input:
        expand(os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3'), species=species_list)
    output:
        expand(os.path.join(out_dir, 'all_species', '{feature}_KS_dist.tsv'), feature=features)
    log:
        os.path.join(logs_dir, 'calculate_distance_matrix', 'all_species.calculate_distance_matrix.log')
    params:
        matrix_script = os.path.join(scripts_dir, 'calculate_distance_matrix.py'),
        out_dir = os.path.join(out_dir, 'all_species')
    conda:
        os.path.join(envs_dir, 'python_analyze.yml')
    resources:
        mem_mb = 20000
    shell:
        """
        python {params.matrix_script} {params.out_dir} {input}
        """

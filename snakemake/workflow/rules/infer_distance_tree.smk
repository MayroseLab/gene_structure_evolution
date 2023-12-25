rule distance_matrix_to_phylip:
    """
    Convert a distance matrix to PHYLIP format
    """
    input:
        os.path.join(out_dir, 'all_species', '{feature}_KS_dist.tsv')
    output:
        os.path.join(out_dir, 'all_species', '{feature}_KS_dist.phylip')
    log:
        os.path.join(logs_dir, 'infer_distance_tree', '{feature}.distance_matrix_to_phylip.log')
    resources:
        mem_mb = 20000
    run:
        def binom_species(s):
            genus, species = s.split('_')[:2]
            return f'{genus}_{species}'

        ks_dist_df = pd.read_csv(input[0], sep='\t', index_col=0)
        ks_dist_df.index = ks_dist_df.index.to_series().apply(binom_species)
        ks_dist_df.index.name = ks_dist_df.shape[0]
        ks_dist_df.columns = ['']*ks_dist_df.shape[1]
        ks_dist_df.to_csv(output[0], sep='\t')

rule infer_distance_tree:
    """
    Infer phylogenetic tree based
    on the distance matrix of a feature.
    Keep the original topology and
    only adjust brach lengths.
    """
    input:
        matrix = os.path.join(out_dir, 'all_species', '{feature}_KS_dist.phylip'),
        species_tree = config['species_tree']
    output:
        os.path.join(out_dir, 'all_species', '{feature}_KS_dist.nwk')
    log:
        os.path.join(logs_dir, 'infer_distance_tree', '{feature}.infer_distance_tree.log')
    resources:
        mem_mb = 20000
    conda:
        os.path.join(envs_dir, 'fastME.yml')
    shell:
        """
        fastme -i {input.matrix} -u {input.species_tree} -o {output}
        """

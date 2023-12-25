rule collect_intron_length_stats:
    """
    Collect intron lengths stats from
    all species into one table
    """
    input:
        intron_length = expand(os.path.join(out_dir, 'per_species', '{species}', 'intron_lengths.stats'), species=species_list),
        genome_size = expand(os.path.join(out_dir, 'per_species', '{species}', 'genome.size'), species=species_list)
    output:
        os.path.join(out_dir, 'all_species', 'intron_lengths.stats')
    log:
        os.path.join(logs_dir, 'collect_intron_length_stats', 'all_species.collect_intron_length_stats.log')
    run:
        len_df = pd.concat([pd.read_csv(tsv, sep='\t', index_col=0) for tsv in input['intron_length']])
        gs_df = pd.concat([pd.read_csv(tsv, sep='\t', index_col=0, names=['Genome_size']) for tsv in input['genome_size']])
        gs_df = pd.concat([species_df, gs_df], axis=1)
        len_df = len_df.join(gs_df, how='outer')
        len_df.to_csv(output[0], sep='\t', float_format=lambda x: round(x,2))

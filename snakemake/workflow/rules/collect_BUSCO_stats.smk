rule collect_BUSCO_stats:
    """
    Collect BUSCO stats for all species.
    Creates two matrices, one for total
    intron lengt and one for n introns
    """
    input:
        expand(os.path.join(out_dir, 'per_species', '{species}', 'BUSCO.stats'), species=species_list),
    output:
        os.path.join(out_dir, 'all_species', 'BUSCO_total_intron_len.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_n_introns.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_intron_frac.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_total_exon_len.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_intron_lens.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_protein_len.stats'),
        os.path.join(out_dir, 'all_species', 'BUSCO_transcript_id.stats')
    log:
        os.path.join(logs_dir, 'collect_BUSCO_stats', 'all_species.collect_BUSCO_stats.log')
    run:
        def read_species_busco_stats(p, species):
            df = pd.read_csv(p, sep='\t', index_col=0)
            df['species'] = species
            return df
        
        def create_busco_matrix(column):
            return all_data_df[['Busco_id', 'species', column]].pivot(index='Busco_id', columns='species', values=column)

        all_data = [read_species_busco_stats(p, sp) for p,sp in zip(input,species_list)]
        all_data_df = pd.concat(all_data)
        all_data_df['tot_exon_len'] = all_data_df['transcript_len'] - all_data_df['tot_intron_len']

        tot_intron_len_df = create_busco_matrix('tot_intron_len')
        tot_intron_len_df.to_csv(output[0], sep='\t', float_format=lambda x: round(x,2))
        tot_exon_len_df = create_busco_matrix('tot_exon_len')
        tot_exon_len_df.to_csv(output[3], sep='\t', float_format=lambda x: round(x,2))
        n_introns_df = create_busco_matrix('n_introns')
        n_introns_df.to_csv(output[1], sep='\t')
        intron_frac_df = create_busco_matrix('intron_frac')
        intron_frac_df.to_csv(output[2], sep='\t', float_format=lambda x: round(x,3))
        intron_lens_df = create_busco_matrix('intron_lens')
        intron_lens_df.to_csv(output[4], sep='\t')
        protein_len_df = create_busco_matrix('Protein_length')
        protein_len_df.to_csv(output[5], sep='\t')
        transcript_id_df = create_busco_matrix('canon_mRNA')
        transcript_id_df.to_csv(output[6], sep='\t')

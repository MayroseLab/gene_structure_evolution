rule create_class_species_list:
    """
    Create file with a list
    of species in each group
    """
    output:
        expand(os.path.join(out_dir, 'all_species', '{class_}.list'), class_=classes)
    log:
        os.path.join(logs_dir, 'create_class_species_list', 'all_classes.create_class_species_list.log')
    resources:
        mem_mb = 10000
    run:
        for cl in classes:
            cl_species = species_df.query('`class` == @cl').index.to_list()
            class_file = os.path.join(out_dir, 'all_species', f'{cl}.list')
            with open(class_file, 'w') as fo:
                print('\n'.join(cl_species), file=fo)

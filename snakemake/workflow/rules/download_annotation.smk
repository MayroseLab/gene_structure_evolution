ensembl_vert_ftp = "ftp.ensembl.org"
ensembl_genomes_ftp = "ftp.ensemblgenomes.ebi.ac.uk"

rule download_annotation:
    """
    Download genome annotation
    in GFF3 format and related
    proteins FASTA
    """
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.gff3'),
        os.path.join(out_dir, 'per_species', '{species}', 'prot.fasta')
    log:
        os.path.join(logs_dir, 'download_annotation', '{species}.download_annotation.log')
    run:
        species = wildcards.species
        group = species_df.loc[species]['group']
        if group == 'vertebrates':
            ftp = ftputil.FTPHost(ensembl_vert_ftp, "anonymous", "")
            gff_ftp = f'pub/current_gff3/{species}'
            pep_ftp = f'pub/current_fasta/{species}/pep'
        else:
            ftp = ftputil.FTPHost(ensembl_genomes_ftp, "anonymous", "")
            gff_ftp = f'pub/current/{group}/gff3/{species}'
            pep_ftp = f'pub/current/{group}/fasta/{species}/pep'

        # download gff3
        ftp.chdir(gff_ftp)
        all_files = ftp.listdir(ftp.curdir)
        full_gff = [f for f in all_files if f.endswith('.gff3.gz') and '.chr.' not in f and '.chromosome.' not in f and '.primary_assembly.' not in f and  '.abinitio.' not in f and 'group' not in f and '.chr_patch_hapl_scaff.' not in f and '.supercontig.' not in f and '.scaffold.' not in f]
        if len(full_gff) != 1:
            print(full_gff)
            sys.exit(f"Can't determine gff3 for species {species}")
        full_gff = full_gff[0]
        ftp.download(full_gff, output[0] + '.gz')
        os.system('gzip -d ' + output[0])

        # download proteins
        ftp.chdir('/')
        ftp.chdir(pep_ftp)
        all_files = ftp.listdir(ftp.curdir)
        full_fasta = [f for f in all_files if f.endswith('.pep.all.fa.gz')]
        if len(full_fasta) != 1:
            print(full_fasta)
            sys.exit(f"Can't determine prot fasta for species {species}")
        full_fasta = full_fasta[0]
        ftp.download(full_fasta, output[1] + '.gz')
        os.system('gzip -d ' + output[1])
        os.system(r"sed -i 's/>.*\(transcript:[^ ]*\) .*/>\1/' " + output[1])

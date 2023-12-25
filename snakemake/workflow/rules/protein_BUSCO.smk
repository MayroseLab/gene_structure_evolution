rule download_lineage:
    """
    Download the BUSCO lineage dataset
    """
    output:
        tgz=os.path.join(out_dir, 'all_species', 'BUSCO_lineage.tar.gz'),
        dir=directory(os.path.join(out_dir, 'all_species', 'busco_downloads', 'lineages', lineage))
    log:
        os.path.join(logs_dir, 'protein_BUSCO', 'download_lineage.log')
    params:
        lineage = lineage,
        out_dir = os.path.join(out_dir, 'all_species'),
        ftp = busco_lineages_ftp
    shell:
         """
         cd {params.out_dir}
         #mkdir busco_downloads/lineages
         # find the file name from list
         file_name=$(wget {params.ftp}"/file_versions.tsv" --no-check-certificate -O - | awk '$1 == "{params.lineage}" {{print $1"."$2".tar.gz"}}')
         file_url={params.ftp}"/lineages/"$file_name
         # download lineage dataset
         wget $file_url -O {output.tgz} --no-check-certificate
         # extract
         tar -zxvf {output.tgz} -C ./busco_downloads/lineages
         """

rule run_BUSCO:
    """
    Run BUSCO on proteins
    """
    input:
        prot=os.path.join(out_dir, 'per_species', '{species}', 'prot.canon.fasta'),
        lineage=os.path.join(out_dir, 'all_species', 'busco_downloads', 'lineages', lineage)
    output:
        os.path.join(out_dir, 'per_species', '{species}', f'BUSCO/run_{lineage}/full_table.tsv')
    log:
        os.path.join(logs_dir, 'protein_BUSCO', '{species}.run_BUSCO.log')
    conda:
        os.path.join(envs_dir, 'BUSCO.yml')
    resources:
        mem_mb = 4000
    params:
        lineage = lineage,
        lineages = os.path.join(out_dir, 'all_species', 'busco_downloads'),
        out_dir = os.path.join(out_dir, 'per_species', '{species}')
    shell:
         """
         cd {params.out_dir}
         busco -i {input.prot} -m prot --offline --download_path {params.lineages} --lineage {params.lineage} -o BUSCO -f
         """

rule parse_BUSCO_table:
    """
    Create a table with single
    enry per BUSCO and canonical
    transcript for this BUSCO
    (if any)
    """
    input:
        busco_tbl = os.path.join(out_dir, 'per_species', '{species}', f'BUSCO/run_{lineage}/full_table.tsv'),
        canon_tbl = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3.mRNA_to_canon'),
        gff = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'BUSCO.stats')
    log:
        os.path.join(logs_dir, 'protein_BUSCO', '{species}.parse_BUSCO_table.log')
    conda:
        os.path.join(envs_dir, 'gffutils.yml')
    params:
        parse_script = os.path.join(scripts_dir, 'parse_BUSCO_table.py')
    shell:
        """
        python {params.parse_script} {input.busco_tbl} {input.canon_tbl} {input.gff} {output}
        """

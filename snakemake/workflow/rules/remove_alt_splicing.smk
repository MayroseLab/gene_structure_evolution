rule remove_alt_splicing:
    """
    Remove alternative splice
    variants for all genes, leaving
    only the canonical transcript.
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.gff3')
    output:
        gff = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3'),
        tbl = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3.gene_to_mRNA'),
        tbl2 = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3.mRNA_to_canon')
    log:
        os.path.join(logs_dir, 'remove_alt_splicing', '{species}.remove_alt_splicing.log')
    conda:
        os.path.join(envs_dir, 'gffutils.yml')
    params:
        remove_alt_splicing_script = os.path.join(scripts_dir, 'remove_alt_splicing_from_gff.py')
    shell:
        """
        python {params.remove_alt_splicing_script} {input} {output.gff}
        """

rule extract_canonical_proteins:
    """
    Based on the list of canonical
    transcripts, extract matching
    proteins from fasta
    """
    input:
        tbl = os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3.gene_to_mRNA'),
        prot = os.path.join(out_dir, 'per_species', '{species}', 'prot.fasta')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'prot.canon.fasta')
    log:
        os.path.join(logs_dir, 'remove_alt_splicing', '{species}.extract_canonical_proteins.log')
    conda:
        os.path.join(envs_dir, 'biopython.yml')
    params:
        filter_fasta_script = os.path.join(scripts_dir, 'filter_fasta_by_id.py')
    shell:
        """
        cut -f2 {input.tbl} | python {params.filter_fasta_script} {input.prot} > {output}
        """

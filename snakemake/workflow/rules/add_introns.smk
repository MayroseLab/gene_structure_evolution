rule add_introns:
    """
    Add intron features to
    annotation gff
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.gff3') 
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3')
    log:
        os.path.join(logs_dir, 'add_introns', '{species}.add_introns.log')
    conda:
        os.path.join(envs_dir, 'genomeTools.yml')
    resources:
        mem_mb = 4000
    shell:
         """
         gt gff3 -tidy -addintrons -retainids yes {input} > {output}
         """

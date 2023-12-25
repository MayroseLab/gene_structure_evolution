rule analyze_intron_lengths:
    """
    Analyze and extract stats
    about intron length from gff
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.canon.introns.gff3'),
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'intron_lengths.stats')
    log:
        os.path.join(logs_dir, 'analyze_intron_lengths', '{species}.analyze_intron_lengths.log')
    params:
        analyze_script = os.path.join(scripts_dir, 'analyze_intron_lengths.py')
    conda:
        os.path.join(envs_dir, 'python_analyze.yml')
    resources:
        mem_mb = 4000
    shell:
        """
        python {params.analyze_script} {input} {output} {wildcards.species}
        """

rule collect_change_stats:
    """
    """
    input:
        expand(os.path.join(out_dir, 'all_species', 'BUSCO_change/{{class_}}/{BUSCO}_change.tsv'), BUSCO=busco_list)
    output:
        os.path.join(out_dir, 'all_species', '{class_}_change.tsv')
    log:
        os.path.join(logs_dir, 'collect_change_stats', '{class_}.collect_change_stats.log')
    resources:
        mem_mb = 8000
    shell:
        """
        header=$(head -1 {input[0]})
        echo $header | tr ' ' '\t' > {output}
        tail -q -n 1 {input} >> {output}
        """

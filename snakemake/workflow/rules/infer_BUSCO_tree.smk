rule infer_BUSCO_tree:
    """
    Tree inference based on
    BUSCO MSA
    """
    input:
        os.path.join(out_dir, 'all_species', 'BUSCO_sequences', '{class_}/{BUSCO}_MSA.faa')
    output:
        os.path.join(out_dir, 'all_species', 'BUSCO_trees', '{class_}/{BUSCO}.treefile')
    log:
        os.path.join(logs_dir, 'infer_BUSCO_tree', '{class_}.{BUSCO}.infer_BUSCO_tree.log')
    params:
        pref = lambda wc: os.path.join(out_dir, 'all_species', 'BUSCO_trees', f'{wc.class_}/{wc.BUSCO}')
    conda:
        os.path.join(envs_dir, 'iqtree.yml')
    resources:
        mem_mb = 2000
    shell:
        """
        if [ -s {input} ]
        then
            iqtree -nt 1 -m LG+I+G -s {input} --prefix {params.pref} 
        else
            touch {output}
        fi
        """

rule calculate_genome_size_from_gff3:
    """
    Calculate and write out
    genome (assembly) size
    """
    input:
        os.path.join(out_dir, 'per_species', '{species}', 'annotation.gff3')
    output:
        os.path.join(out_dir, 'per_species', '{species}', 'genome.size')
    log:
        os.path.join(logs_dir, 'calculate_genome_size', '{species}.calculate_genome_size.log')
    shell:
        """
        set +e
        seq=$(grep '##sequence-region' {input})
        if [ -z "$seq" ]
        then
            gs="NA"
        else
            gs=$(echo "${{seq}}" | awk '{{SUM+=$4}}END{{print SUM}}')
        fi
        echo -e "{wildcards.species}\t$gs" > {output}
        """

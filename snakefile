######
# Ensemble MAG annotator
# Antton Alberdi
# 2024/08/25
######

import os
import re

#### Fetch input genomes

def detect_extension(directory, valid_extensions):
    extensions = set() 
    pattern = re.compile(r'\.([a-zA-Z0-9]+(?:\.[a-zA-Z0-9]+)*)$')

    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isdir(file_path):
            continue
        match = pattern.search(filename)
        if match:
            ext = match.group(1)
            if ext in valid_extensions:
                extensions.add(ext)

    if len(extensions) == 1:
        return extensions.pop()  # Return the single common extension
    else:
        print("Error: Genome files have different extensions.")
        return None

valid_extensions = ["fa", "fna", "fasta", "fa.gz", "fna.gz", "fasta.gz"]
extension = detect_extension("resources/genomes/",valid_extensions)

# List genome and target wildcards
genomes, = glob_wildcards(f"resources/genomes/{{genome}}.{extension}")

# Expand target files
rule all:
    input:
        expand("results/output/{genome}.tsv", genome=genomes)

rule prepare_input:
    input:
        lambda wildcards: f"resources/genomes/{wildcards.genome}.{extension}"
    output:
        "results/input/{genome}.fna"
    params:
        jobname="{genome}.pi",
        extension=lambda wildcards: extension
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    shell:
        """
            if [[ "{params.extension}" == "gz" ]]; then
                gunzip -c {input} > {output}
            else
                cp {input} {output}
            fi
        """

rule prodigal:
    input:
        "results/input/{genome}.fna"
    output:
        fna="results/prodigal/{genome}.fna",
        faa="results/prodigal/{genome}.faa"
    params:
        jobname="{genome}.pr",
    threads:
        1
    resources:
        mem_gb=8,
        time=15
    shell:
        """
        module load prodigal/2.6.3
        prodigal -i {input} -d {output.fna} -a {output.faa}
        """

rule final:
    input:
        "results/prodigal/{genome}.faa"
    output:
        "results/output/{genome}.tsv"
    params:
        jobname="allgenomes.fi"
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    shell:
        """
        cat {input} > {output} 
        """

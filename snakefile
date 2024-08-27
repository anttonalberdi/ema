######
# Ensemble MAG annotator
# Antton Alberdi
# 2024/08/25
######

import os
import re
import requests

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

###########################
#### Prepare databases #### 
###########################

rule prepare_kofams:
    output:
        h3f="resources/databases/kofams/kofams.h3f",
        h3i="resources/databases/kofams/kofams.h3i",
        h3m="resources/databases/kofams/kofams.h3m",
        h3p="resources/databases/kofams/kofams.h3p"
    params:
        jobname="pr.kofams"
    threads:
        1
    resources:
        mem_gb=16,
        time=300
    shell:
        """
        # Create directory
        if [ ! -f resources/databases/kofams ]; then
            mkdir resources/databases/kofams
        fi

        # Download
        if [ ! -f resources/databases/kofams/profiles.tar.gz ]; then
            cd resources/databases/kofams/
            wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
        fi

        # Decompress
        if [ ! -f resources/databases/kofams/kofams ]; then
            tar -xvzf profiles.tar.gz
            cat profiles/*hmm > kofams
            rm -rf profiles
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            module load hmmer/3.3.2
            hmmpress -f kofams
        fi
        """

rule prepare_cazy:
    output:
        h3f="resources/databases/cazy/cazy.h3f",
        h3i="resources/databases/cazy/cazy.h3i",
        h3m="resources/databases/cazy/cazy.h3m",
        h3p="resources/databases/cazy/cazy.h3p"
    params:
        jobname="pr.cazy"
    threads:
        1
    resources:
        mem_gb=16,
        time=300
    shell:
        """
        # Create directory
        if [ ! -d resources/databases/cazy ]; then
            mkdir resources/databases/cazy
        fi

        # Download
        if [ ! -f resources/databases/cazy/cazy ]; then
            cd resources/databases/cazy/
            wget https://bcb.unl.edu/dbCAN2/download/dbCAN-HMMdb-V13.txt
            mv dbCAN-HMMdb-V13.txt cazy
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            module load hmmer/3.3.2
            hmmpress -f cazy
        fi
        """

rule prepare_pfam:
    output:
        h3f="resources/databases/pfam/pfam.h3f",
        h3i="resources/databases/pfam/pfam.h3i",
        h3m="resources/databases/pfam/pfam.h3m",
        h3p="resources/databases/pfam/pfam.h3p"
    params:
        jobname="pr.pfam"
    threads:
        1
    resources:
        mem_gb=16,
        time=300
    shell:
        """
        # Create directory
        if [ ! -d resources/databases/pfam ]; then
            mkdir resources/databases/pfam
        fi

        # Download
        if [ ! -f resources/databases/pfam/pfam ]; then
            cd resources/databases/pfam/
            wget https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
            gunzip Pfam-A.hmm.gz
            mv Pfam-A.hmm pfam
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            module load hmmer/3.3.2
            hmmpress -f pfam
        fi
        """

### Run pipeline

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

rule kofams:
    input:
        faa="results/prodigal/{genome}.faa",
        db="resources/databases/kofams/kofams.h3p"
    output:
        txt="results/kofams/{genome}.txt",
        tsv="results/kofams/{genome}.tsv"
    params:
        jobname="{genome}.kf",
        db="resources/databases/kofams/kofams"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    shell:
        """
        module load hmmer/3.3.2
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input.faa}
        """

rule cazy:
    input:
        faa="results/prodigal/{genome}.faa",
        db="resources/databases/cazy/cazy.h3p"
    output:
        txt="results/cazy/{genome}.txt",
        tsv="results/cazy/{genome}.tsv"
    params:
        jobname="{genome}.cz",
        db="resources/databases/cazy/cazy"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    shell:
        """
        module load hmmer/3.3.2
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input.faa}
        """

rule pfam:
    input:
        faa="results/prodigal/{genome}.faa",
        db="resources/databases/pfam/pfam.h3p"
    output:
        txt="results/pfam/{genome}.txt",
        tsv="results/pfam/{genome}.tsv"
    params:
        jobname="{genome}.pf",
        db="resources/databases/pfam/pfam"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    shell:
        """
        module load hmmer/3.3.2
        hmmscan -o {output.txt} --tblout {output.tsv} --noali {params.db} {input.faa}
        """

rule final:
    input:
        kofams="results/kofams/{genome}.txt",
        cazy="results/cazy/{genome}.txt",
        pfam="results/pfam/{genome}.txt"
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

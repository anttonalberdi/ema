######
# Ensemble MAG annotator
# Antton Alberdi
# 2024/08/25
######

import os
import re
import requests
import glob


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
        expand("results/final/{genome}.tsv", genome=genomes)
        
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
        h3p="resources/databases/kofams/kofams.h3p",
        json="resources/databases/kofams/kegg.json"
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
        if [ ! -d resources/databases/kofams ]; then
            mkdir resources/databases/kofams
        fi

        # Download
        if [ ! -f resources/databases/kofams/profiles.tar.gz ]; then
            wget -O resources/databases/kofams/profiles.tar.gz https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
        fi

        # Download KEGG hierarchy
        if [ ! -f {output.json} ]; then
            curl -o {output.json} "https://www.kegg.jp/kegg-bin/download_htext?htext=ko00001.keg&format=json&filedir="
        fi

        # Decompress
        if [ ! -f resources/databases/kofams/kofams ]; then
            tar -xvzf resources/databases/kofams/profiles.tar.gz -C resources/databases/kofams
            cat resources/databases/kofams/profiles/*hmm > resources/databases/kofams/kofams
            rm -rf resources/databases/kofams/profiles
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            module load hmmer/3.3.2
            hmmpress -f resources/databases/kofams/kofams
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

        # Change to working directory
        cd resources/databases/cazy/

        # Download
        if [ ! -f cazy ]; then
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
        h3p="resources/databases/pfam/pfam.h3p",
        ec="resources/databases/pfam/pfam_ec.tsv"
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

        # Download pfams
        if [ ! -f resources/databases/pfam/pfam ]; then
            wget -O resources/databases/pfam/Pfam-A.hmm.gz https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz
            gunzip resources/databases/pfam/Pfam-A.hmm.gz
            mv resources/databases/pfam/Pfam-A.hmm resources/databases/pfam/pfam
        fi

        # Download pfam-ec mapping
        if [ ! -f resources/databases/pfam/pfam_ec.tsv ]; then
            wget -O resources/databases/pfam/EC-Pfam_calculated_associations_Extended.csv https://ecdm.loria.fr/data/EC-Pfam_calculated_associations_Extended.csv
            mv resources/databases/pfam/EC-Pfam_calculated_associations_Extended.csv {output.ec}
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            cd resources/databases/pfam/
            module load hmmer/3.3.2
            hmmpress -f pfam
        fi
        """

rule prepare_vfdb:
    output:
        db="resources/databases/vfdb/vfdb.idx",
        mapping="resources/databases/vfdb/vfdb.tsv"
    params:
        jobname="pr.vfdb",
        db="resources/databases/vfdb/vfdb",
    threads:
        1
    resources:
        mem_gb=16,
        time=15
    conda:
        "workflow/envs/environment.yml"
    shell:
        """
        # Create directory
        if [ ! -d resources/databases/vfdb ]; then
            mkdir resources/databases/vfdb
        fi

        # Download
        if [ ! -f resources/databases/vfdb/VFDB_setB_pro.fas ]; then
            wget -O resources/databases/vfdb/VFDB_setB_pro.fas.gz http://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
            gunzip resources/databases/vfdb/VFDB_setB_pro.fas.gz
        fi

        #Generate mapping file
        if [ ! -f {output.mapping} ]; then
            python workflow/scripts/get_fvid.py resources/databases/vfdb/VFDB_setB_pro.fas {output.mapping}
        fi

        #Create mmseqs2 db
        if [ ! -f {output.db} ]; then
            module load mmseqs2/14.7e284
            mmseqs createdb resources/databases/vfdb/VFDB_setB_pro.fas {params.db}
            mmseqs createindex {params.db} tmp
        fi
        """

rule prepare_amr:
    output:
        h3f="resources/databases/amr/amr.h3f",
        h3i="resources/databases/amr/amr.h3i",
        h3m="resources/databases/amr/amr.h3m",
        h3p="resources/databases/amr/amr.h3p",
        tsv="resources/databases/amr/amr.tsv"
    params:
        jobname="pr.amr"
    threads:
        1
    resources:
        mem_gb=16,
        time=300
    shell:
        """
        # Create directory
        if [ ! -d resources/databases/amr ]; then
            mkdir resources/databases/amr
        fi

        # Download and decompress
        if [ ! -f resources/databases/amr/amr ]; then
            wget -O resources/databases/amr/NCBIfam-AMRFinder.HMM.tar.gz https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz
            tar -xvzf resources/databases/amr/NCBIfam-AMRFinder.HMM.tar.gz -C resources/databases/amr
            cat resources/databases/amr/HMM/*.HMM > resources/databases/amr/amr
            rm -rf HMM
        fi

        # Download metadata
        if [ ! -f {output.tsv} ]; then
            wget -O {output.tsv} https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.tsv
        fi

        # Build index
        if [ ! -f {output.h3p} ]; then
            module load hmmer/3.3.2
            hmmpress -f resources/databases/amr/amr
        fi
        """

##########################
####   Run pipeline   #### 
##########################

# Predict genes
rule prodigal:
    input:
        "results/input/{genome}.fna"
    output:
        fna="results/prodigal/{genome}.fna",
        faa="results/prodigal/{genome}.faa",
        gff="results/prodigal/{genome}.gff"
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
        prodigal -i {input} -d {output.fna} -a {output.faa} -o {output.gff} -f gff
        """

# Annotate KEGG orthologs
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

# Annotate carbohydrate active enzymes (CAZYs)
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

# Annotate protein families and assign enzyme commission codes (EC)
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

# Annotate virulence factors
rule vfdb:
    input:
        faa="results/prodigal/{genome}.faa",
        db="resources/databases/vfdb/vfdb.idx"
    output:
        "results/vfdb/{genome}.txt"
    params:
        jobname="{genome}.vf",
        db="resources/databases/vfdb/vfdb"
    threads:
        1
    resources:
        mem_gb=8,
        time=60
    shell:
        """
        module load mmseqs2/14.7e284
        mmseqs easy-search {input.faa} {params.db} {output} results/vfdb/{wildcards.genome}
        """

# Annotate antimicrobial resistance genes
rule amr:
    input:
        faa="results/prodigal/{genome}.faa",
        db="resources/databases/amr/amr.h3p"
    output:
        txt="results/amr/{genome}.txt",
        tsv="results/amr/{genome}.tsv"
    params:
        jobname="{genome}.ar",
        db="resources/databases/amr/amr"
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

# Annotate signal peptides
rule signalp:
    input:
        "results/prodigal/{genome}.faa"
    output:
        "results/signalp/{genome}.txt"
    params:
        jobname="{genome}.sp",
        outputdir="results/signalp/{genome}"
    threads:
        4
    resources:
        mem_gb=4,
        time=120
    shell:
        """
        module load signalp/6h
        signalp6 --fastafile {input} --output_dir {params.outputdir} --write_procs {threads}
        cat {params.outputdir}/output.gff3  | cut -f1,3,6 | awk -F' # |[ \t]+' '!/^#/ {{print $1, $6, $7}}' OFS='\t' > {output}
        """

# Filter and merge all annotations
rule final:
    input:
        gff="results/prodigal/{genome}.gff",
        kofamsdb="results/kofams/{genome}.tsv",
        kofams="resources/databases/kofams/kegg.json",
        cazy="results/cazy/{genome}.tsv",
        pfam="results/pfam/{genome}.tsv",
        ec="resources/databases/pfam/pfam_ec.tsv",
        vfdb="results/vfdb/{genome}.txt",
        vf="resources/databases/vfdb/vfdb.tsv",
        amrdb="results/amr/{genome}.tsv",
        amr="resources/databases/amr/amr.tsv",
        sp="results/signalp/{genome}.txt"
    output:
        "results/final/{genome}.tsv"
    params:
        jobname="merge_annotations"
    threads:
        1
    conda:
        "workflow/envs/environment.yml"
    resources:
        mem_gb=8,
        time=5
    shell:
        """
        python workflow/scripts/merge_annotations.py \
            -gff {input.gff} \
            -kofamsdb {input.kofamsdb} \
            -kofams {input.kofams} \
            -pfam {input.pfam} \
            -ec {input.ec} \
            -cazy {input.cazy} \
            -vfdb {input.vfdb} \
            -vf {input.vf} \
            -amrdb {input.amrdb} \
            -amr {input.amr} \
            -signalp {input.sp} \
            -o {output}
        """

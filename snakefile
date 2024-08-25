######
# Ensemble MAG annotator
# Antton Alberdi
# 2024/08/25
######

#### Fetch input genomes

def detect_extension(directory, valid_extensions):
    extensions = set() 
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isdir(file_path):
            continue
        _, ext = os.path.splitext(filename)
        ext = ext.lstrip('.')
        if ext:
            # Check if the extension is in the list of valid extensions
            if ext in valid_extensions:
                extensions.add(ext)
            elif any(ext.endswith(valid_ext) for valid_ext in valid_extensions):
                # For cases like .gz extensions, check the base extension
                base_ext = ext.split('.')[-2]  # Assumes .gz is last if present
                if base_ext in valid_extensions:
                    extensions.add(base_ext)
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
        "results/annotations.tsv"


rule prepare_input:
    input:
        "resources/genomes/{genome}.{extension}"
    output:
        "results/input/{genome}.fna"
    params:
        jobname="{genome}.pi",
        extension="{extension}"
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    shell:
        """
            if [[ "{params.ext}" == "gz" ]]; then
                echo "Decompressing {input} to {output}"
                gunzip -c {input} > {output}
            else
                echo "Processing {input} as {output}"
                cp {input} {output}
            fi
        """

rule final:
    input:
        expand("results/input/{genome}.fna", genome=genomes)
    output:
        "results/annotations.tsv"
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

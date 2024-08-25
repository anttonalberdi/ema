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
        expand("results/output/{genome}.tsv", genome=genomes)

rule prepare_input:
    input:
        lambda wildcards: f"resources/genomes/{wildcards.genome}.{extension}"
    output:
        "results/input/{genome}.fna"
    params:
        jobname="{genome}.pi",
        extension=lambda wildcards: extension  # Use lambda to pass the extension
    threads:
        1
    resources:
        mem_gb=8,
        time=5
    shell:
        """
            if [[ "{params.extension}" == "gz" ]]; then
                echo "Decompressing {input} to {output}"
                gunzip -c {input} > {output}
            else
                echo "Processing {input} as {output}"
                cp {input} {output}
            fi
        """

rule final:
    input:
        "results/input/{genome}.fna"
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

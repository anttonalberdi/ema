# Ensemble MAG annotator

Annotates metagenome-assembled genomes (MAGs) or complete bacterial genomes against the KEGG, PFAM, CAZY, VFDB and AMR databases using HMMER and MMSEQS2, and classifies annotations into ecologically meaningful traits in a straightforward way: you just need to place the genome fasta files in the `resources/genomes` folder.

## Steps

1. Download and prepare databases
2. Annotate genomes
3. Merge annotations

## Usage

1. Clone this repo.
2. Place genomes in the `resources/genomes` folder.
3. Create a screen session.
4. Launch the snakemake using the following code:
```
module purge && module load snakemake/7.20.0 mamba/1.3.1
snakemake -j 20 --cluster 'sbatch -o log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600 --keep-going
```

### Resistance type (resistance_type)

resistance_type == AMRFinderPlus subtype
from: https://github.com/ncbi/amr/wiki/Interpreting-results

AMR	Antimicrobial resistance gene
AMR-SUSCEPTIBLE	Not used in AMRFinderPlus output, but in the Reference Gene Catalog (see note below)
POINT	Known point mutation associated with antimicrobial resistance
VIRULENCE	Virulence gene
ANTIGEN	Gene codes for a known antigen (these are often used for typing)
ACID	Acid resistance gene
BIOCIDE	Biocide resistance gene
HEAT	Heat resistance gene
METAL	Metal resistance gene

## Structure

### 1. METABOLISM
Proteins involved in catalyzing metabolic reactions can be further divided based on their specific biochemical functions:

- **Carbohydrate metabolism**: Proteins involved in breaking down or synthesizing sugars (e.g., glycolysis, gluconeogenesis, fermentation).
- **Lipid metabolism**: Proteins involved in fatty acid synthesis, degradation, and lipid-based energy storage.
- **Amino acid metabolism**: Proteins involved in the synthesis and breakdown of amino acids.
- **Nucleotide metabolism**: Enzymes related to the synthesis and breakdown of nucleotides.
- **Cofactor/vitamin biosynthesis**: Proteins involved in the production of cofactors (e.g., NADH, FAD) and vitamins (e.g., biotin, B12).

### Resistance target (resistance_target)

resistance_target == AMRFinderPlus subclass
from: https://github.com/ncbi/amr/wiki/class-subclass

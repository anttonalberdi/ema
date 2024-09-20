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

## Structure

### 1. METABOLISM
Proteins involved in catalyzing metabolic reactions can be further divided based on their specific biochemical functions:

- **Carbohydrate metabolism**: Proteins involved in breaking down or synthesizing sugars (e.g., glycolysis, gluconeogenesis, fermentation).
- **Lipid metabolism**: Proteins involved in fatty acid synthesis, degradation, and lipid-based energy storage.
- **Amino acid metabolism**: Proteins involved in the synthesis and breakdown of amino acids.
- **Nucleotide metabolism**: Enzymes related to the synthesis and breakdown of nucleotides.
- **Cofactor/vitamin biosynthesis**: Proteins involved in the production of cofactors (e.g., NADH, FAD) and vitamins (e.g., biotin, B12).

### 2. TRANSPORTER
Transporter proteins can be divided based on the types of molecules they transport and the mechanism of transport:

- **Ion transporters**: Proteins that move ions such as Na⁺, K⁺, Cl⁻, Ca²⁺ across membranes.
- **Nutrient transporters**: Proteins that facilitate the uptake of essential nutrients, including sugars, amino acids, and vitamins.
- **Secondary transporters**: Proteins that couple the transport of one molecule to another, typically using ion gradients (e.g., symporters, antiporters).
- **ABC transporters**: ATP-powered transporters that move a wide variety of molecules, including lipids, ions, and organic molecules.

### 3. RECEPTOR
Receptor proteins can be further categorized by the type of signal they detect or the downstream signaling pathways they trigger:

- **Chemoreceptors**: Proteins that detect chemical gradients (e.g., chemotaxis receptors).
- **Photoreceptors**: Proteins that sense light (e.g., rhodopsins, photolyases).
- **Histidine kinases**: Part of two-component signaling systems that detect environmental stimuli and trigger a phosphorylation cascade.
- **Quorum-sensing receptors**: Proteins that detect signaling molecules used in bacterial communication for population-density sensing.
- **Toll-like receptors (TLRs)**: Involved in recognizing microbial-associated molecular patterns in immune-related contexts.

### 4. STRUCTURE
Structural proteins can be divided based on the cellular structure they contribute to:

- **Cytoskeletal proteins**: Proteins like MreB, FtsZ, and crescentin involved in cell shape and division.
- **Cell wall proteins**: Proteins involved in peptidoglycan synthesis and maintenance (e.g., penicillin-binding proteins).
- **Flagellar and pili proteins**: Proteins forming components of bacterial motility structures.
- **Capsule proteins**: Proteins involved in the synthesis and assembly of the capsule that protects the microbe from the host immune system.
- **Spore coat proteins**: Proteins involved in forming the outer coat of bacterial spores (in spore-forming bacteria like Bacillus).

### 5. SIGNAL
Signal-related proteins can be subdivided based on their function in cell-to-cell communication or intracellular signaling:

- **Quorum sensing molecules**: Autoinducers and other molecules involved in population-density-based signaling.
- **Second messengers**: Small molecules like cAMP, cGMP, or ppGpp that mediate internal signaling.
- **Phosphorelays**: Components of two-component systems that pass on signals through phosphorylation cascades.
- **Extracellular signaling peptides**: Proteins or peptides secreted to communicate with other microbes.
- **Autoinducers**: Small diffusible molecules like acyl-homoserine lactones (AHLs) that regulate quorum sensing.

### 6. REGULATOR
Regulatory proteins can be classified based on the processes they regulate:

- **Transcription factors**: Proteins that bind to DNA and regulate transcription.
- **Post-transcriptional regulators**: Proteins involved in RNA processing or degradation (e.g., RNA-binding proteins).
- **Sigma factors**: Specialized subunits of RNA polymerase that regulate sets of genes during specific conditions (e.g., stress response).
- **Transcriptional repressors**: Proteins that bind to operator sequences and block transcription.
- **Epigenetic regulators**: Proteins that modify DNA or histones (in eukaryotic microbes) to control gene expression.

### 7. MOBILITY
Mobile genetic elements can be subdivided based on the type of mobility they promote:

- **Integrases**: Enzymes that integrate or excise DNA elements (e.g., integrons).
- **Transposases**: Enzymes that move transposable elements (jumping genes).
- **Phage-related proteins**: Proteins associated with bacteriophages, including capsid proteins, tail fibers, and packaging enzymes.
- **Conjugation proteins**: Proteins involved in plasmid transfer between bacteria through conjugation (e.g., Tra proteins).
- **Recombinases**: Enzymes involved in site-specific recombination and gene shuffling.

### 8. DEFENSE
Defense-related proteins can be split based on the specific type of threat they protect against:

- **Antibiotic resistance proteins**: Proteins that degrade/inactivate antibiotics (e.g., beta-lactamases).
- **Efflux pumps**: Proteins that actively export toxic compounds or antibiotics out of the cell.
- **Toxin-antitoxin systems**: Pairs of proteins where one is a toxin and the other an antitoxin, often used in stress response.
- **Restriction-modification systems**: Proteins that protect against foreign DNA, such as restriction enzymes and DNA methyltransferases.
- **CRISPR-associated proteins**: Proteins involved in microbial adaptive immunity against viruses (e.g., Cas proteins).
- **Oxidative stress proteins**: Proteins that protect against reactive oxygen species (e.g., catalase, superoxide dismutase).

### 9. ENERGY
Energy-related proteins can be divided based on the energy-conversion process they are involved in:

- **Electron transport proteins**: Proteins involved in electron transfer, such as cytochromes, quinones, and iron-sulfur proteins.
- **ATP synthase component**s: Subunits of the ATP synthase complex involved in oxidative phosphorylation.
- **Photosynthetic proteins**: Proteins involved in capturing light energy (e.g., photosystem components).
- **Fermentation enzymes**: Proteins that enable energy production through fermentation.
- **Proton pumps**: Proteins that move protons across membranes to generate a proton motive force (e.g., bacteriorhodopsin).

### 10. CHAPERONE
Chaperones can be split based on their specific roles in protein maintenance:

- **Folding chaperones**: Proteins that assist in proper protein folding (e.g., GroEL, Hsp70, DnaK).
- **Unfolding chaperones**: Proteins involved in unfolding misfolded proteins for degradation.
- **Disaggregation chaperones**: Proteins that help recover proteins from aggregates.
- **Protein export chaperones**: Proteins that assist in the export or secretion of newly synthesized proteins.
- **Heat shock proteins**: Chaperones specifically upregulated during heat or stress conditions.

## Resistance details

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

### Resistance target (resistance_target)

resistance_target == AMRFinderPlus subclass
from: https://github.com/ncbi/amr/wiki/class-subclass

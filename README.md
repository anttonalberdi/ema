# Ensemble MAG annotator

Annotates metagenome-assembled genomes (MAGs) or complete bacterial genomes against the KEGG, PFAM, CAZY, VFDB and AMR databases, and classifies annotations into ecologically meaningful traits. 

## Usage

1. Clone this repo.
2. Place genomes in the `resources/genomes` folder.
3. Create a screen session.
4. Launch the snakemake using the following code:
```
module purge && module load snakemake/7.20.0 mamba/1.3.1
snakemake -j 20 --cluster 'sbatch -o log/{params.jobname}-slurm-%j.out --mem {resources.mem_gb}G --time {resources.time} -c {threads} --job-name={params.jobname} -v'   --use-conda --conda-frontend mamba --conda-prefix conda --latency-wait 600 --keep-going
```


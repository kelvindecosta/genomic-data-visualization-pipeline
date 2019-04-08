# Genomic Data Visualization Pipeline

This is my attempt at the selection tests for the [GSoC 2019](https://summerofcode.withgoogle.com/) project, "Human history and data visualization" under the [Canadian Center for Computational Genomics](http://computationalgenomics.ca/)

## Usage

```
usage: main.py [-h] [--log_dir] [--plot_dir] [--coords_dir] [--genedata]
               [--labels] [--pops] [--seed] [--downsample_n] [--prune_size]
               [--prune_step] [--prune_thresh] [--prune_iter] [--pca_n]
               [--pca_scaler] [--plot_limit_components] [--title TITLE]

Applying Dimensionality Reduction on Genomic Data

optional arguments:
  -h, --help            show this help message and exit
  --log_dir             log directory
  --plot_dir            plot directory
  --coords_dir          coordinates directory
  --genedata            genomic data file
  --labels              sample population mappings file
  --pops                population description file
  --seed                seed for prng
  --downsample_n        number of SNPs to randomly choose
  --prune_size          size for pruning
  --prune_step          step for pruning
  --prune_thresh        threshold for pruning
  --prune_iter          iterations for pruning
  --pca_n               number of principal components
  --pca_scaler          scaler for pca
  --plot_limit_components 
                        number of components to consider plotting, [None]
                        means plot all
  --title TITLE         process title
```

## Example

With the dataset in a folder ```data``` in the repository directory, run the following command :

```bash
python main.py --genedata=data/ALL.wgs.nhgri_coriell_affy_6.20140825.genotypes_has_ped.vcf --labels=data/affy_samples.20141118.panel --pops=data/20131219.populations.tsv --title="Chromosome 1 Assembly B37" --plot_limit_components=4
```

After running command, the project structure changes to :

```
.
├── coords
│   └── [1554763502304] Chromosome 1 Assembly B37
│       ├── [1554763502304] Chromosome 1 Assembly B37 - [PCA] coords.txt
│       ├── [1554763502304] Chromosome 1 Assembly B37 - [UMAP] coords.txt
│       └── [1554763502304] Chromosome 1 Assembly B37 - [UMAP_PCA] coords.txt
├── logs
│   └── [1554763502304] Chromosome 1 Assembly B37.log
├── main.py
├── plots
│   └── [1554763502304] Chromosome 1 Assembly B37 plots
│       ├── [1554763502304] Chromosome 1 Assembly B37 [PCA]
│       │   ├── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[1][2].png
│       │   ├── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[1][3].png
│       │   ├── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[1][4].png
│       │   ├── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[2][3].png
│       │   ├── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[2][4].png
│       │   └── [1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[3][4].png
│       ├── [1554763502304] Chromosome 1 Assembly B37 [UMAP]
│       │   └── [1554763502304] Chromosome 1 Assembly B37 [UMAP] - PC[1][2].png
│       └── [1554763502304] Chromosome 1 Assembly B37 [UMAP_PCA]
│           └── [1554763502304] Chromosome 1 Assembly B37 [UMAP_PCA] - PC[1][2].png
└── README.md
```

The [`coords`](coords/) directory contains the lower dimensional coordinates written into text files.

The [`logs`](logs/) directory contains the log reports.

The [`plots`](plots/) directory has the plots of the components with one another.

### PCA

<img src="plots/[1554763502304] Chromosome 1 Assembly B37 plots/[1554763502304] Chromosome 1 Assembly B37 [PCA]/[1554763502304] Chromosome 1 Assembly B37 [PCA] - PC[1][2].png" alt="alt text" title="Chromosome 1 Assembly B37 [PCA] - PC[1][2]">

### UMAP

<img src="plots/[1554763502304] Chromosome 1 Assembly B37 plots/[1554763502304] Chromosome 1 Assembly B37 [UMAP]/[1554763502304] Chromosome 1 Assembly B37 [UMAP] - PC[1][2].png" alt="alt text" title="[1554763502304] Chromosome 1 Assembly B37 [UMAP] - PC[1][2]">

### UMAP + PCA

<img src="plots/[1554763502304] Chromosome 1 Assembly B37 plots/[1554763502304] Chromosome 1 Assembly B37 [UMAP_PCA]/[1554763502304] Chromosome 1 Assembly B37 [UMAP_PCA] - PC[1][2].png" alt="alt text" title="[1554763502304] Chromosome 1 Assembly B37 [UMAP_PCA] - PC[1][2]">
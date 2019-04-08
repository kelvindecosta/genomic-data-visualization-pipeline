import allel
import argparse
import logging
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import pandas as pd
import random
import seaborn as sns
import sys
import time
import umap
from os.path import join, exists
from os import makedirs

def package_version():
    """
    Returns a dictionary with the package names and versions
    """
    packages = {}
    packages['scikit-allel'] = allel.__version__
    packages['matplotlib'] = matplotlib.__version__
    packages['numpy'] = np.__version__
    packages['pandas'] = pd.__version__
    packages['seaborn'] = sns.__version__
    packages['umap'] = umap.__version__

    return packages

def setup_log(filename):
    """
    @Params: filename: relative path to log file
    """
    logging.StreamHandler(sys.stdout)
    logging.basicConfig(format='[%(asctime)s] : %(message)s', level=logging.INFO, filename=filename)


def population_descriptions(filename):
    """
    @Params: filename: relative path to population description file
    """
    df = pd.read_csv(filename , sep = '\t')  
    
    label_name_dict = {} 
    for i in range(df.shape[0]):
        label_name_dict[df.loc[i]['Population Code']] = df.loc[i]['Population Description']
    
    logging.info(f"Constructed dictionary of population codes to population descriptions from {filename}")
    return label_name_dict


def samples_population_mapping(filename):
    """
    @Params: filename: relative path to samples to population mapping file
    """
    df = pd.read_csv(filename , sep = '\t') 
    
    sample_label_dict = {} 
    for i in range(df.shape[0]):
        sample_label_dict[df.loc[i]['sample']] = df.loc[i]['pop']

    logging.info(f"Constructed dictionary of samples to population codes from {filename}")
    return sample_label_dict


def genotype_array_from_vcf(filename):
    """
    @Params: filename: relative path to genotpye data file
    """
    g =  allel.GenotypeChunkedArray(allel.read_vcf(filename)['calldata/GT'])
    logging.info(f"Loaded Genotype Data from file {filename}")
    return g

def filter_genotype_array(g):
    """
    @Params: g: unfiltered GenotypeChunkedArray
    """
    ac = g.count_alleles()[:]
    flt = (ac.max_allele() == 1) & (ac[:, :2].min(axis=1) > 1)
    g = g.compress(flt, axis=0)
    logging.info(f"Filtered Genotype Data")
    return g



def input_data_from_genotype_array(g):
    """
    @Params: g: filtered GenotypeChunkedArray
    """
    return g.to_n_alt()


def downsample(g, n, seed):
    """
    @Params: g: input data format
             n: Number of SNPs to randomly choose
             seed: seed for prng
    """
    random.seed(seed)
    np.random.seed(seed)
    vidx = np.random.choice(g.shape[0], n, replace=False)
    vidx.sort()
    g = g.take(vidx, axis=0)
    logging.info(f"Applied downsampling to data")
    return g


def ld_prune(gn, size, step, threshold, n_iter):
    """
    Applies pruning; removing SNPs that are not correlated
    """
    for i in range(n_iter):
        loc_unlinked = allel.locate_unlinked(gn, size=size, step=step, threshold=threshold)
        n = np.count_nonzero(loc_unlinked)
        n_remove = gn.shape[0] - n
        logging.info(f"Pruning iteration {i+1} retaining  {n} removing {n_remove} variants")
        gn = gn.compress(loc_unlinked, axis=0)
    logging.info(f"Applied pruning with paramaters, size : {size}, step : {step}, threshold : {threshold}, n_iter : {n_iter}")
    return gn[:]


def apply_pca(g, outfile, seed,  n, s):
    """
    Applies PCA to data and saves low-dimensional coordinates in outfile
    @Params: g: input data format
             outfile: path to coordinate file
             seed: seed for prng
             n: Number of principal components
             s: scaler
    """
    coords, _ = allel.randomized_pca(g, n_components=n, scaler=s, random_state=seed)
    logging.info(f"Applied PCA with {n} components with scaler {s} and seed {seed}")
    np.savetxt(outfile, coords)
    logging.info(f"Saved coordinates to {outfile}")


def apply_umap(g, outfile):
    """
    Applies UMAP to data and saves low-dimensional coordinates in outfile
    @Params: g: input data format
             outfile: path to coordinate file
    """
    embedding = umap.UMAP(n_neighbors=5, min_dist=0.3, metric='correlation').fit_transform(g)
    logging.info(f"Applied UMAP")
    np.savetxt(outfile, embedding)
    logging.info(f"Saved coordinates to {outfile}")


def apply_umap_pca(pca_coords_file, outfile):
    """
    Applies UMAP to data from PCA reduced coordinates file and saves low-dimensional coordinates in outfile
    @Params: pca_coords_file: input data coordinates file
             outfile: path to coordinate file
    """
    coords = np.loadtxt(pca_coords_file)
    embedding = umap.UMAP(n_neighbors=5, min_dist=0.3, metric='correlation').fit_transform(coords)
    logging.info(f"Applied UMAP on PCA reduced space")
    np.savetxt(outfile, embedding)
    logging.info(f"Saved coordinates to {outfile}")


def generate_plots(coords_file, samples_mapping, populations_dict, title, current_type_dir, limit, figsize=(10, 10)):
    """
    Generates plots for coordinates based on labels and descriptions
    @Params:    coords_file: path to file with low-dimensional coordinates
                samples_mapping: path to labels
                populations_dict: path to descriptions
                title: plot title
                current_type_dir: directory to save plots
                limit: number of components to plot
                figsize: tuple describing plot size in inches
    """
    coords = np.loadtxt(coords_file)
    logging.info(f"Loaded coordinates from {coords_file}")

    n = coords.shape[1]
    if limit and limit < n:
        n = limit
    

    logging.info(f"Plotting for {n} components")
    colors = len(populations_dict)
    cm = plt.get_cmap('gist_rainbow')
    
    for i in range(n):
        for j in range(i, n):
            if i != j:
                fig = plt.subplots(figsize=figsize)
                ax = plt.subplot(111)
                ax.set_prop_cycle(color=[cm(float(x/colors)) for x in range(colors)])
                sns.despine(ax=ax, offset=5)
                x = coords[:, i]
                y = coords[:, j]
                for pop, label in populations_dict.items():
                    filtr = (samples_mapping == pop)
                    ax.plot(x[filtr], y[filtr], marker='o', linestyle=' ', label=label, markersize=6, mec='k', mew=.5)
                ax.set_xlabel(f"PC{i+1}")
                ax.set_ylabel(f"PC{j+1}")
                box = ax.get_position()
                ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
                plt.title(title)
                logging.info(f"Saving plot of PC{i+1} vs PC{j+1}")
                plt.savefig(join(current_type_dir, f"{title} - PC[{i+1}][{j+1}].png"), bbox_inches = 'tight')
                plt.clf()

def get_args():
    parser = argparse.ArgumentParser(description='Applying Dimensionality Reduction on Genomic Data')
    parser.add_argument('--log_dir' , metavar='', help = 'log directory', default="./logs/") 
    parser.add_argument('--plot_dir' , metavar='', help = 'plot directory', default="./plots/")
    parser.add_argument('--coords_dir' , metavar='', help = 'coordinates directory', default="./coords/") 
    parser.add_argument('--genedata' , metavar='', help = 'genomic data file')
    parser.add_argument('--labels' , metavar='', help = 'sample population mappings file')
    parser.add_argument('--pops' , metavar='', help= 'population description file')
    parser.add_argument('--seed' , metavar='', help= 'seed for prng', default=42)
    parser.add_argument('--downsample_n' , metavar='', help= 'number of SNPs to randomly choose', default=10000)
    
    parser.add_argument('--prune_size' , metavar='', help= 'size for pruning', default=500)
    parser.add_argument('--prune_step' , metavar='', help= 'step for pruning', default=200)
    parser.add_argument('--prune_thresh' , metavar='', help= 'threshold for pruning', default=0.1)
    parser.add_argument('--prune_iter' , metavar='', help= 'iterations for pruning', default=5)
    
    parser.add_argument('--pca_n' , metavar='', help= 'number of principal components', default=10)
    parser.add_argument('--pca_scaler' , metavar='', help= 'scaler for pca', default='patterson')

    parser.add_argument('--plot_limit_components' , type=int ,metavar='', help= 'number of components to consider plotting, [None] means plot all', default=None)


    parser.add_argument('--title', help= 'process title')
    return parser.parse_args() 


def dictstring(d, l):
    r = ""
    for k, v in d.items():
        r += f"\t\t{k:{l}}: {v}\n"
    return r[:-1]

# Seaborn
sns.set_style('white')
sns.set_style('ticks')

if __name__ == "__main__":
    start_time = time.time()
    args = get_args()
    for d in [args.log_dir, args.plot_dir, args.coords_dir]:
        if not exists(d):
            makedirs(d)
    
    pid = f"[{int(start_time * 1000)}] {args.title}"
    logfile = join(args.log_dir, f"{pid}.log")

    setup_log(logfile)
    logging.info(f"Package Versions\n{dictstring(package_version(), 25)}")
    logging.info(f"Parameters Passed\n{dictstring(vars(args), 25)}")
    g = genotype_array_from_vcf(args.genedata)
    g = filter_genotype_array(g)
    g = input_data_from_genotype_array(g)
    g = downsample(g, args.downsample_n, args.seed)
    g = ld_prune(g, args.prune_size, args.prune_step, args.prune_thresh, args.prune_iter)

    types = ["PCA", "UMAP", "UMAP_PCA"]
    times = {}

    coords_dir = join(args.coords_dir, f"{pid}")
    if not exists(coords_dir):
        makedirs(coords_dir) 
    coords_files = [join(coords_dir, f"{pid} - [{x}] coords.txt") for x in types]

    # PCA
    s = time.time()
    apply_pca(g, coords_files[0], args.seed, args.pca_n, args.pca_scaler)
    e = int((time.time() - s) * 1000)
    times[types[0]] = f"{e} ms"

    # UMAP
    s = time.time()
    apply_umap(np.transpose(g), coords_files[1])
    e = int((time.time() - s) * 1000)
    times[types[1]] = f"{e} ms"

    # UMAP + PCA
    s = time.time()
    apply_umap_pca(coords_files[0], coords_files[2])
    e = int((time.time() - s) * 1000)
    times[types[2]] = f"{e} ms"

    # Preparing data to plot
    # Loading population descriptions
    populations_dict = population_descriptions(args.pops)
    logging.info(f"Loaded descriptions from {args.pops}")

    # Loading samples population
    samples_mapping = np.array(list(samples_population_mapping(args.labels).values()))
    logging.info(f"Loaded samples mapping from {args.labels}")
    
    plotdir = join(args.plot_dir, f"{pid} plots")
    if not exists(plotdir):
        makedirs(plotdir)
    
    for i in range(len(types)):
        cdir = join(plotdir, f"{pid} [{types[i]}]")
        if not exists(cdir):
            makedirs(cdir)
        generate_plots(coords_files[i], samples_mapping, populations_dict, f"{pid} [{types[i]}]", cdir,  args.plot_limit_components)
    elapsed_time = int(time.time() - start_time)
    times["Total"] = f"{elapsed_time} s"

    logging.info(f"Elapsed Time\n{dictstring(times, 10)}")

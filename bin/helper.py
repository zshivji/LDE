import pandas as pd
import argparse
import glob
import os

import matplotlib.pyplot as plt
from matplotlib_venn import venn3

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
import Bio.Phylo.NewickIO

import gzip
import re
from functools import lru_cache

from ete3 import NCBITaxa

def load_hits(update_index, hits_file=None, archaea_only=False):
    if hits_file:
        nif = pd.read_csv(hits_file)
    elif archaea_only:
        nif = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
    else:
        nif_archaea = pd.read_csv('../results/hmmsearch_results/archaea_hits.csv')
        nif_bacteria = pd.read_csv('../results/hmmsearch_results/bacteria_hits.csv')
        nif = pd.concat([nif_archaea, nif_bacteria], ignore_index=True)

    # nif.reset_index(inplace = True)
    nif.set_index(update_index, inplace = True) # set index to selected columns

    return nif

def load_faa(dir: str = "../protein_faa_reps_232", genome=None):
    # load  protein fasta files from the specified directory
    # if using GTDB, pull the protein fasta file for the specified genome
    if dir.__contains__("protein_faa_reps_232"):
        candidates = glob.glob(f"{dir}/*/{genome}_protein.faa")
        if not candidates:
            raise FileNotFoundError(
                f"No protein FASTA files found under {dir}"
            )
        file = candidates[0]
    elif dir.endswith((".faa", ".fa", ".fasta", '.txt')): # if a file is passed (not dir)
        file = dir
    else: 
         # if running the API, the protein fasta files are organized into a single directory (accepts multiple file types)
        files = [
            f for f in glob.glob(f"{dir}/*")
            if f.endswith((".faa", ".fa", ".fasta", '.txt'))
            ]
        if not files:
            raise FileNotFoundError(
                f"No protein FASTA files found under {dir}"
            )
        file = files[0] if files else None

    return file

def get_seq(genome, hit, id=None, description=None, start=None, end=None, dir: str = "../protein_faa_reps_232"):
    # get the sequence of a hit from the protein fasta file
    record_id = genome if id is None else id
    record_description = hit if description is None else description

    file = load_faa(dir, genome)

    # get sequence from fasta file 
    for result in SeqIO.parse(file, "fasta"):
        if result.id == hit:
            # store seq
            seq = result.seq[start:end] # trim seq to match domain len
            # convert to seqrecord
            return SeqRecord(seq, id=record_id, description=record_description)
    raise KeyError(f"Sequence {hit!r} was not found in {file}")

def filter_groups_by_unique_counts(df, group_cols, requirements, exclude='nifB'):
    """
    Filter groups based on multiple unique-count requirements.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataframe
    group_cols : list[str]
        Columns defining groups
    requirements : dict
        Dictionary of {column: minimum_unique_count}

    Example
    -------
    requirements={
        "Gene": 3,
        "Hit": 3
    }
    """

    def has_excluded_protein(x):
        return (
            exclude is not None
            and 'Gene' in x.columns
            and exclude in x['Gene'].unique()
        )

    def meets_requirements(x):
        return all(
            x[col].nunique() >= threshold
            for col, threshold in requirements.items()
        )

    # Step 0: exit function if min_genes==1
    if requirements.get('gene', 0) <= 1 and requirements.get('hit', 0) <= 1:
        return df
    
    # Step 1: remove genomes that fewer than 3 unique nif genes.
    df = df.copy()
    df = df.groupby('GenomeID').filter(
        lambda x: meets_requirements(x)
    )

    # Step 2: remove groups (typically operons) that do not meet the requirements or contain the excluded protein.
        # an operon can be kept with fewer than (3) unique nif genes if it contains the excluded protein (e.g., nifB)
    return df.groupby(group_cols).filter(
        lambda x: has_excluded_protein(x) or meets_requirements(x)
    )

def tree_node_match_metadata(cluster_file,tree_file=None, tree_csv=None):

    clusters = pd.read_csv(cluster_file)

    if tree_file is not None: #rewrite newick tree with updated node names
        tree = Phylo.read(tree_file, "newick")

        for clade in tree.get_terminals():
            # for each node, update nodeID to match the metadata headers (except added outbranches)
            if clade.name.startswith('WP_'):
                new_name = clade.name
            else:
                protein = clade.name
                contig = '_'.join(protein.split('_')[:-1])
                position = protein.split('_')[-1]
                cluster = clusters[clusters["contig"] == contig]
                for _, row in cluster.iterrows():
                    positions = str(row['pos_num']).strip('[]').replace(',', ' ').split()
                    if position in {p.strip() for p in positions}:
                        organism = row['Organism']
                        clusterID = row['cluster']
                        genome = row['GenomeID']
                        new_name = f"{organism} | {clusterID} | {genome} | {protein}"
                        break

            clade.name = new_name

        Phylo.write(tree, tree_file, "newick")
    
    elif tree_csv is not None: #export tree node names to csv
        node_names = pd.read_csv(tree_csv, header=None, sep="\t", names=['current', 'new'])
        new_names = []

        for row in node_names.iterrows():
            current_name = row[1]['current']
            if current_name.startswith('AFDB'):
                current_name = current_name.replace("_",":",1)
            tmp = clusters[clusters['protein'] == current_name]
            domain = str(tmp['domain'].values[0]) + ' | ' if not tmp.empty else ''
            genome = str(tmp['GenomeID'].values[0]) + ' | ' if not tmp.empty else ''
            lde = str(tmp['LDE'].values[0]) + ' | ' if not tmp.empty else ''
            organism = str(tmp['Organism'].values[0]) + ' | ' if not tmp.empty else ''
            new_names.append(f"{domain}{organism}{lde}{genome}{current_name}")

        node_names['new'] = new_names
        node_names.to_csv(tree_csv, index=False, header=False, sep="\t")

def tree_comparison(cutoff=90):
    # for each comparison DB, create a list of nodes that have a match in diazoDB
    # matched are considered when the percent identity is above the cutoff
    #ref_tree = Phylo.read('../diazoDB-comparison/tree-comparison/nifD_anfD_vnfD_clustered.trim.treefile', 'newick')
    comparison_DBs = ['NSDB', 'Nif-Finder', 'NFixDB']
    #ref_nodes = {clade.name for clade in ref_tree.get_terminals()}
    node_matches = pd.read_csv('../diazoDB-comparison/tree-comparison/nifD_anfD_vnfD_clustered.fasta.tsv', sep='\t', names=['cluster', 'acc'], index_col='acc')
    node_matches[comparison_DBs] = False

    for comparison in comparison_DBs:
        blast = pd.read_csv(f'../diazoDB-comparison/venn/{comparison}.blast', 
                            sep='\t', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])

        blast = blast[blast['evalue'] < 1e-5]
        blast = blast.loc[blast.groupby('qseqid')['pident'].idxmax()]
        matches = blast.loc[blast['pident'] >= cutoff, 'sseqid']
        node_matches.loc[node_matches.index.isin(matches), comparison] = True        

    # collapse into cluster
    cluster_matches = node_matches.groupby('cluster', as_index=False)[comparison_DBs].any()

    # convert to boolean
    cluster_matches[comparison_DBs] = cluster_matches[comparison_DBs].astype(int)

    # export to csv
    cluster_matches.to_csv('../diazoDB-comparison/tree-comparison/nifD-tree-comparison.csv', index=False)

# def venn():
#     # create a venn diagram of the comparison DBs
#     comparison_DBs = ['NSDB', 'Nif-Finder', 'NFixDB']

#     # get the sets of clusters for each comparison DB


#     # create the venn diagram
#     plt.figure(figsize=(8, 8))
#     venn3(sets, set_labels=comparison_DBs)
#     plt.title('Venn Diagram of Comparison DBs')
#     plt.savefig('../diazoDB-comparison/tree-comparison/nifD-tree-comparison-venn.png')
#     plt.close()

def get_AFDB_domain(record, fasta_path = "../results/tmp/AFDB_hits.fasta"):
    """Yield one dictionary per FASTA record.

    domain is Bacteria, Archaea, Eukaryota, or None.
    Metagenome/environmental/unspecified sources return None,
    even if their lineage includes a domain.
    """
    ncbi = NCBITaxa()

    # Match stable lineage IDs rather than depending on rank names.
    domain_ids = {
        2: "Bacteria",
        2157: "Archaea",
        2759: "Eukaryota",
    }

    nonspecific = re.compile(
        r"\b(metagenom\w*|environmental|uncultured|unclassified|"
        r"unidentified|unspecified|synthetic|mixed)\b",
        re.IGNORECASE,
    )

    @lru_cache(maxsize=100_000)
    def lookup_domain(taxid):
        if taxid is None:
            return None

        try:
            lineage = ncbi.get_lineage(taxid)
        except ValueError:  # Unknown or deleted taxonomy ID
            return None

        names = ncbi.get_taxid_translator(lineage)

        # Also inspect the taxonomy itself, not just the FASTA OS field.
        if any(nonspecific.search(name) for name in names.values()):
            return None

        return next(
            (domain_ids[node] for node in lineage if node in domain_ids),
            None,
        )

    opener = gzip.open if str(fasta_path).endswith(".gz") else open

    with opener(fasta_path, "rt") as fasta:
        for line in fasta:
            if not line.startswith(">"):
                continue

            header = line[1:].strip()
            if header.split()[0] == record:
                
                ox = re.search(r"\bOX=(\d+)\b", header)
                os = re.search(
                    r"\bOS=(.*?)(?=\s+[A-Z][A-Z0-9]*=|$)", header
                )
                ua = re.search(r"\bUA=(\S+)", header)

                taxid = int(ox.group(1)) if ox else None
                organism = os.group(1).strip() if os else None # already included in AFDB fasta header

                domain = (
                    None
                    if organism and nonspecific.search(organism)
                    else lookup_domain(taxid)
                )

                return domain, organism

import json
import os
import subprocess
import re

def subset_results(
        prefix = "", 
        cutoff=20, 
        tmp_dir = "../results/tmp"):
    # subset hits with a higher passing_score
    prefix = prefix + "_" if prefix else ""
    gene = 'pqq'

    # 1. get subset
    subset = pd.read_csv(os.path.join(tmp_dir, f'{prefix}{gene}_residues.csv'))
    subset = subset[subset['score'] >= cutoff]

    ind = subset['hit'].values
    ind2 = subset['contig'].values

    # 2. subset final.csv, operons.csv
    result_dir = "../results/final"
    final = pd.read_csv(os.path.join(result_dir, f'{prefix}final.csv'), index_col=['hit'])
    final = final[final.index.isin(ind)]
    final.to_csv(os.path.join(result_dir, f'{prefix}final_{cutoff}.csv'), index=False)

    clusters = pd.read_csv(os.path.join(result_dir, f'{prefix}operons.csv'), index_col=['hit'])
    clusters = clusters[clusters.index.isin(ind)]
    clusters.to_csv(os.path.join(result_dir, f'{prefix}operons_{cutoff}.csv'), index=False)
    
    # 3. subset fastas
    with open("subset.txt", "w") as file:
        file.write("\n".join(ind2) + "\n")

    input = f'{result_dir}/fastas/{prefix}{gene}.fasta'
    output = f'{result_dir}/fastas/{prefix}{gene}_{cutoff}.fasta'

    subprocess.run(
        ["seqkit", "grep", "--quiet", "-f", "subset.txt", 
         input, "-o", output],
        check=True)

    # # 4. subset metadata
    # metadata_file = os.path.join(result_dir, "{prefix}metadata.json")
    # metadata_output = None

    # with open(metadata_file) as handle:
    #     metadata = json.load(handle)

    #     # Prefer the explicit genome/contig/operon fields, then fall back to
    #     # gene IDs for metadata produced by older versions of metadata.py.
    #     selected_locations = {
    #         (str(row.get("GenomeID", "")), str(row.get("contig", "")), str(row.get("operon", "")))
    #         for _, row in final_subset.iterrows()
    #     }
    #     filtered_metadata = {}
    #     for key, value in metadata.items():
    #         location = (
    #             str(value.get("genome", "")),
    #             str(value.get("contig", "")),
    #             str(value.get("operon", "")),
    #         )
    #         genes = value.get("operon", {}).get("genes", [])
    #         gene_ids = {str(gene.get("gene_id", "")) for gene in genes}
    #         if location in selected_locations or gene_ids & selected_ids:
    #             filtered_metadata[key] = value

    #     metadata_stem, metadata_extension = os.path.splitext(os.fspath(metadata_file))
    #     metadata_output = f"{metadata_stem}_{cutoff}{metadata_extension}"
    #     with open(metadata_output, "w") as handle:
    #         json.dump(filtered_metadata, handle, indent=2, allow_nan=False)

def subset_tree(clusters="../trees/config20/GTDB_AFDB_pqq_20_clustered_cluster.tsv"):
    subset = pd.read_csv(clusters, header=None, sep='\t').iloc[:,0].to_list()
    df2 = pd.read_csv('../results/final/GTDB_final_20.csv', index_col = 'protein')
    df3 = pd.read_csv('../results/final/AFDB_final_20.csv', index_col = 'protein')
    df = pd.concat([df2, df3])
    df = df[df.index.isin(subset)]

    df.to_csv('../results/final/tree20.csv')

def map_contig(file1, file2, col="contig"):
    # df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)
    # df3 = pd.read_csv('../results/hmmsearch_results/bacteria_hits.csv')
    # df = pd.concat([df1,df3], ignore_index=True)
    # df1 = df.drop_duplicates(subset=['hit'])
    # col2map = df.set_index('hit')[col]
    df2[col] = df2['protein']
    df2.to_csv(file2, index=False)


def main():
    parser = argparse.ArgumentParser(description="DiazoDB helper functions")
    subparsers = parser.add_subparsers(dest="command", required=True)

    tree_parser = subparsers.add_parser(
        "tree_node_match_metadata",
        help="Replace tree terminal IDs with metadata-matched IDs.",
    )
    tree_parser.add_argument("tree_file", help="Newick tree file to update in place.")

    args = parser.parse_args()
    if args.command == "tree_node_match_metadata":
        tree_node_match_metadata(args.tree_file)

    


if __name__ == "__main__":
    main()

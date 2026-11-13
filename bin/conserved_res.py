import pandas as pd
import numpy as np
import glob
import json
import os
import subprocess
import sys
from tqdm import tqdm
import argparse
import shutil

from Bio import AlignIO
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from helper import load_hits, filter_groups_by_unique_counts, get_seq
from test import test

NO_RESULTS_EXIT = 3

# get args (reload_fasta)
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make and align temporary fasta files based on filtered hits."
            "Check for conserved residues in alignments."
            "Update nif gene assignments based on conserved residue matching."
        )
    )
    parser.add_argument(
        "--reload_fasta",
        help="Whether or not to reload the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed).",
        action='store_true'
    )

    parser.add_argument(
        "--align_fasta",
        help="Whether or not to align the fasta files (slow but necessary if hmmsearch/parse_hmmp.py changed). Automatically defaults to True if --reload_fasta is specified.",
        action='store_true'
    )

    parser.add_argument(
        "--config_file",
        help="Path to the JSON configuration file containing reference sequences and important residues. Default bin/nif-config.json",
        default='LDE-config.json'
    )

    parser.add_argument(
        "--proteins_dir",
        help=(
            "Root containing <group>/<genome>_protein.faa files. "
            "The API runner supplies its uploaded proteins and packaged references."
        ),
        default='../protein_faa_reps_232',
    )

    parser.add_argument(
        "--min_genes",
        type=int,
        default=3,
        help="Minimum number of nif genes required to be considered a tophit cluster (default: 3)."
    )

    parser.add_argument(
        "--test",
        action='store_true',
        help="Whether or not to run the final results against the test file."
    )
    parser.add_argument(
        "--hits_file",
        help=(
            "Optional parsed HMM hits CSV. Defaults to the repository-wide "
            "archaea_hits.csv and bacteria_hits.csv files."
        ),
    )
    parser.add_argument(
        "--tmp",
        default="../results/tmp",
        help="Directory for temporary FASTA, alignment, and residue-check files.",
    )
    parser.add_argument(
        "--final_dir",
        default="../results/final",
        help="Directory for nif_final.csv, nif_clusters.csv, and final FASTAs.",
    )
    parser.add_argument(
        "--out_prefix",
        help="Prefix for output files (default: hits)."
    )
    parser.add_argument(
        "--skip_metadata",
        action="store_true",
        help="Do not read or join repository taxonomy metadata (used by API runner jobs).",
    )

    parser.add_argument(
        "--external",
        action="store_true",
        help="Export a cleaned-up nif_final.csv and nif_clusters.csv for external DiazoDB users.",
    )

    return parser.parse_args()

# Note: runtime data (nif table, config, per-gene frames) are loaded in `main`.
# Helper functions below accept explicit inputs and return results.

def write_fasta(
    gene_list: dict,
    config: dict,
    proteins_dir: str = "../protein_faa_reps_232",
    tmp: str = "../results/tmp",
    prefix: str = ""
):
    """Get fasta sequences for each gene & export to fasta of domains.
    """
    os.makedirs(tmp, exist_ok=True)
    for name, gene in gene_list.items():
        print('getting ' + name + ' fasta', flush=True)
        # filter out hits that don't meet length criteria (domain is too short, pre-filtering decreases runtime)
        minl = config[name]['min-length']
        gene = gene[abs(gene['Domain Start'] - gene['Domain End']) >= minl]
        gene_list[name] = gene

        records = []
        for genome, hit in gene.iterrows():
            seq = get_seq(
                genome,
                hit.hit[:-2],
                description=hit.hit,
                start=hit['Domain Start'],
                end=hit['Domain End'],
                dir=proteins_dir,
            )
            records.append(seq)

        print(len(records), "records found for " + name, flush=True)
        with open(os.path.join(tmp, f"{prefix}{name}.fasta"), "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

    gene_list = pd.concat(gene_list.values())
    gene_list.to_csv(os.path.join(tmp, f'{prefix}hits.csv'))

    return gene_list

def aln_fasta(
    genes: dict,
    config: dict,
    tmp_dir: str = "../results/tmp",
    prefix: str = ""
):
    # align fasta files
    splits_dir = os.path.join(tmp_dir, f'{prefix}fasta_splits')
    print(splits_dir, flush=True)
    shutil.rmtree(splits_dir) if os.path.exists(splits_dir) else None
    os.makedirs(splits_dir, exist_ok=True)
    for name in genes['gene'].unique():
        print("aligning "+ name, flush=True)
        tmp = genes[genes['gene']==name].copy()
        num = int(tmp.shape[0]/200)+1 # how many splits
        fasta_path = os.path.join(tmp_dir, f"{prefix}{name}.fasta")
        split_prefix = os.path.join(splits_dir, f"{name}_split")
        # split fasta file into smaller files for alignment
            # seqtk fails after 1024 so need to use seqkit
            # -f, force overwrite exisiting files
        subprocess.run(
            ["seqkit", "split2", "--quiet", "-f", "-p", str(num), "-O", splits_dir, "-P", f"{name}_split.", fasta_path],
            check=True
        )

        # manually zfill the split file names to ensure correct order for mafft
        for file in glob.glob(os.path.join(splits_dir, f"{name}_split.*.fasta")):
            base = os.path.basename(file)
            parts = base.split('.')
            if len(parts) == 3 and parts[1].isdigit():
                new_name = f"{parts[0]}.{str(int(parts[1])).zfill(5)}.fasta"
                os.rename(file, os.path.join(splits_dir, new_name))

        # for each split...
        for i in range(num):
            split_path = f"{split_prefix}.{str(i+1).zfill(5)}.fasta"
            alignment_path = f"{split_prefix}.{str(i+1).zfill(5)}.aln"
            
            # add reference sequences to split file
            ref = SeqRecord(Seq(config[name]['ref_seq']), id='reference', description=name)
            with open(split_path, "a") as f:
                SeqIO.write(ref, f, "fasta")

            # align split file
            with open(alignment_path, "wb") as alignment_file:
                subprocess.run(
                    ["mafft", "--auto", "--quiet", "--thread", "4", split_path],
                    stdout=alignment_file,
                    check=True,
                )


def check_gene(file, important_residues, residue_scores, passing_score, p=False):
    if len(important_residues) != len(residue_scores):
        raise ValueError(
            f"{file[-20:-16]}: residues and residue_scores must have the same length "
            f"({len(important_residues)} != {len(residue_scores)})"
        )

    residue_to_score = dict(zip(important_residues, residue_scores))
    alignment = AlignIO.read(f"{file}", "fasta")

    for record in alignment:
        if record.id == 'reference':
            aln = record.seq
            break

    # map important cols in aln to ref_seq
    residue_to_alignment = {}
    residue_idx = 0
    col2res = {}
    col2score = {}

    for aln_idx, char in enumerate(aln):
        if char != '-':
            residue_idx += 1
            if residue_idx in important_residues:
                residue_to_alignment[residue_idx] = aln_idx
                col2res[aln_idx] = f"{aln[aln_idx]}{residue_idx}"
                col2score[aln_idx] = residue_to_score[residue_idx]
                
            # early stop if we've found everything
            if len(residue_to_alignment) == len(important_residues):
                break

    # print the corresponding residues in the original sequence
    if p:
        for residue, aln_idx in residue_to_alignment.items():
            print(f"Residue {residue} corresponds to alignment index {aln_idx}: {aln[aln_idx]}", flush=True)

    # store in dataframe        
    acc = [result.id for result in alignment]
    seqs = [list(str(result.seq)) for result in alignment]
    hits = [result.description.split(' ')[-1] for result in alignment]

    pssm = pd.DataFrame(seqs, index = acc)

    pssm = pssm.iloc[:, list(residue_to_alignment.values())]
    pssm['hit'] = hits
    pssm['contig'] = pssm['hit'].str.split('_').str[:-1].str.join('_')


    # check if cols contain correct residue for function
    def check_res(row):
        score = 0
        for col in residue_to_alignment.values():
            if row[col] == aln[col]:
                score += 1 + col2score[col]
            elif row[col] != '-': # higher weight for correct C
                    score += 1
        if score >= passing_score:
            return score
        else:
            return np.nan

    pssm['score'] = pssm.apply(check_res, axis = 1)
    pssm.dropna(subset = ['score'], inplace = True)
    pssm.rename(columns=col2res, inplace=True)

    return pssm

def check_gene_wrapper(gene_list: dict, config: dict, tmp_dir: str = "../results/tmp", prefix: str = "") -> dict:
    print("Checking for conserved residues...", flush=True)

    all_genes_checked = {}

    for gene in gene_list['gene'].unique():
        print(f'checking {gene}', flush=True)

        # load json parameters
        important_residues = config[gene]['residues']
        residue_scores = config[gene].get('residue_scores', config[gene].get('reside_scores'))
        if residue_scores is None:
            residue_scores = [1] * len(important_residues)
        passing_score = config[gene]['passing_score']

        # loop through each split and check for conserved residues
            # store in list of dataframes
        checked = []
        for i, file in enumerate(glob.glob(os.path.join(tmp_dir, f'{prefix}fasta_splits', f"{gene}_split.*.aln"))):
            checked.append(check_gene(file, important_residues, residue_scores, passing_score, p=(i==1)))

        # concat list of df into df
        checked_df = pd.concat(checked) if checked else pd.DataFrame()
        checked_df.drop_duplicates(subset = ['hit'], inplace = True)
        checked_df.set_index(['hit'], append= True, inplace = True)
        checked_df = checked_df[checked_df.index != 'reference'] # drop reference row
        checked_df.to_csv(os.path.join(tmp_dir, f'{prefix}{gene}_residues.csv'))

        # check metal dependence (Asp320, LDE only***)
        checked_df['LDE'] = checked_df.apply(lambda row: 1 if row.get('D320') == 'D' else 0, axis=1)

        # only keep gene where residue matching passed
        tmp = gene_list[(gene_list['gene'] == gene) & (gene_list['hit'].isin(checked_df.index.get_level_values(1)))].copy()
        tmp['residue_match'] = gene
        # checked_df has a MultiIndex (sequence id, hit), while hit is a
        # column in tmp.  Reduce the residue results to a hit-indexed lookup
        # before assigning so values align on the matching hit.
        lde_by_hit = checked_df.reset_index().set_index('hit')['LDE']
        tmp['LDE'] = tmp['hit'].map(lde_by_hit) # apply metal dependence (Asp320, LDE only***)
        # drop duplicates (keep lowest e-val)
        tmp = tmp.sort_values('E-value').drop_duplicates('hit', keep='first')
        # tmp.to_csv(os.path.join(tmp_dir, f'{prefix}{gene}_rescheck_nofilt_config{passing_score}.csv'))

        all_genes_checked[gene] = tmp

        print(f"{len(tmp)} {gene} seqs \n", flush=True)

    # export
    os.makedirs(tmp_dir, exist_ok=True)
    all_genes_checked = pd.concat(all_genes_checked.values())
    # print(all_genes_checked, flush=True)
    all_genes_checked.set_index('GenomeID', inplace=True)
    all_genes_checked.to_csv(os.path.join(tmp_dir, f'{prefix}rescheck_nofilt_config{passing_score}.csv')) # index = GenomeID

    return all_genes_checked

# EXPORT
def export_results(
    # all_genes_checked,
    # min_genes: int,
    results_dir: str = "../results/final",
    proteins_dir: str = "../protein_faa_reps_232",
    include_metadata: bool = True,
    external: bool = False,
    prefix: str = ""
) -> None:
    os.makedirs(results_dir, exist_ok=True)
    tmp_dir = "../results/tmp"
    passing_score = 12
    hits = pd.read_csv(os.path.join(tmp_dir, f'{prefix}rescheck_nofilt_config{passing_score}.csv'))
    # hits = all_genes_checked.copy()

    # assign final annotation + group domains into proteins
    hits['protein'] = hits['hit'].str.rsplit('_', n=1).str[0] # remove domain suffix to get protein name
    genomes_to_keep = hits.copy()#pd.DataFrame(columns=hits.columns).set_index('protein')

    if prefix.__contains__('AFDB'):
        genomes_to_keep['Location_start'] = 0
        genomes_to_keep['Location_end'] = 0
        genomes_to_keep['pos_num'] = 0
    else:
        location_parts = genomes_to_keep['Location'].str.split('-', expand=True)
        genomes_to_keep['Location_start'] = location_parts[0].astype(int)
        genomes_to_keep['Location_end'] = location_parts[1].astype(int)
        genomes_to_keep['pos_num'] = genomes_to_keep['protein'].str.rsplit('_').str[-1].astype(int)

    # # get fasta sequences for each gene & export to fasta
    # os.makedirs(os.path.join(results_dir, "fastas"), exist_ok=True)
    # for gene in genomes_to_keep['gene'].unique().tolist():
    #     records = []

    #     for genome, hit in genomes_to_keep[genomes_to_keep.gene == gene].iterrows():
    #         record = get_seq(
    #             genome,
    #             hit.protein,
    #             id=hit.protein,
    #             description=f"{genome} {hit.gene}",
    #             dir=proteins_dir,
    #         )
    #         records.append(record)
                    
    #     # Write the records to a FASTA file
    #     with open(os.path.join(results_dir, f"fastas/{prefix}{gene}.fasta"), "w") as output_handle:
    #         SeqIO.write(records, output_handle, "fasta")

    #     print(f"{len(records)} final {gene} hits", flush=True)


    if include_metadata:
        if prefix.__contains__('GTDB'):
            # Database-build mode enriches representative genomes with taxonomy.
            gtdb_metadata = pd.read_csv(
                "GTDB_metadata.gz",
                sep="\t",
                usecols=[
                    "accession",
                    "gtdb_taxonomy",
                    "ncbi_taxonomy",
                    "ncbi_taxonomy_unfiltered",
                    "ncbi_country",
                    "ncbi_isolation_source",
                ],
            )

            genomes_to_keep = genomes_to_keep.join(
                gtdb_metadata.set_index("accession"),
                on="GenomeID",
                how="left",
            )

            genomes_to_keep["Organism"] = (
                genomes_to_keep["gtdb_taxonomy"].str.split(";").str[-1].str.split("__").str[-1]
            )
            genomes_to_keep["domain"] = genomes_to_keep["gtdb_taxonomy"].str.split(";").str[0].str.split("__").str[-1]
            # genomes_to_keep["Regulon"] = ""
            # genomes_to_keep["Operon"] = ""
            genomes_to_keep.rename(
                columns={
                    "gtdb_taxonomy": "GTDB Taxonomy",
                    "ncbi_taxonomy": "NCBI Taxonomy",
                    "ncbi_isolation_source": "Isolation Source",
                },
                inplace=True,
            )
        elif prefix.__contains__('AFDB'):
            from helper import get_AFDB_domain
            print('getting AFDB domains', flush=True)

            domains = {}
            organisms = {}

            for record in SeqIO.parse(f"{results_dir}/fastas/{prefix}pqq.fasta", "fasta"):
                domain, organism = get_AFDB_domain(record.id)
                domains[record.id] = domain
                organisms[record.id] = organism
            


            genomes_to_keep["domain"] = genomes_to_keep["protein"].map(domains)
            genomes_to_keep["Organism"] = genomes_to_keep["protein"].map(organisms)

            genomes_to_keep['GTDB Taxonomy'] = ""
            genomes_to_keep['NCBI Taxonomy'] = ""
            genomes_to_keep['Isolation Source'] = ""


    genomes_to_keep["Group"] = ""
    # export csv with each gene as individual row
    if external: # clean up export for external DiazoDB users
        external_export = genomes_to_keep.copy()
        external_export['Sequence Length'] = external_export['Domain End'] - external_export['Domain Start']
        external_export = external_export[[GenomeID','contig', 'operon', 'protein', 'pos_num', 'gene', 'Location', 
                                           'Sequence Length', 'Orientation', 
                                           'Location_start', 'Location_end', 'LDE', 'Group']]
        
        external_export.to_csv(os.path.join(results_dir, f'{prefix}final.csv'))
        operons = external_export.copy()
        operons.to_csv(os.path.join(results_dir, f'{prefix}operons.csv'))


    else:
        genomes_to_keep.to_csv(os.path.join(results_dir, f'{prefix}final.csv'))
        operons = genomes_to_keep.copy()
        operons.to_csv(os.path.join(results_dir, f'{prefix}operons.csv'))

    return genomes_to_keep


def main() -> int:
    args = parse_args()

    # load config and nif table at runtime
    config_file = args.config_file
    config = json.load(open(config_file, 'r'))

    prefix = args.out_prefix + '_' if args.out_prefix is not None else ''

    print('Running conserved-residue check with the following parameters:', flush=True)
    print(f"hmm hits: {args.hits_file}", flush=True)
    print('prefix:', prefix, flush=True)
    print(f"config passing score: {config['pqq']['passing_score']}", flush=True)

    # load hits dataframe
    hits = load_hits(update_index=['GenomeID'], hits_file=args.hits_file)
    
    if hits.empty:
        print("No qualifying HMM hits; stopping pipeline.", file=sys.stderr, flush=True)
        return NO_RESULTS_EXIT

    # store per-gene df in a dict
    gene_list = {gene: hits[hits.gene == gene].copy() for gene in config.keys()}
    print(f"checking conserved residues for: {list(gene_list)}", flush=True)

    # reload + align fasta files if specified
    if args.reload_fasta:
        gene_list = write_fasta(gene_list, config, proteins_dir=args.proteins_dir, tmp=args.tmp, prefix=prefix)
        args.align_fasta = True  # Automatically align if reloading fasta
    else: 
        gene_list = pd.read_csv(os.path.join(args.tmp, f'{prefix}hits.csv'))
    if args.align_fasta:
        aln_fasta(gene_list, config, tmp_dir=args.tmp, prefix=prefix)

    # check for conserved residues and save results
    checked = check_gene_wrapper(gene_list, config, tmp_dir=args.tmp, prefix=prefix)
    if checked.empty:
        print("No hits passed the conserved-residue check; stopping pipeline.", file=sys.stderr, flush=True)
        return NO_RESULTS_EXIT

    export_results(
        checked,
        args.min_genes,
        results_dir=args.final_dir,
        proteins_dir=args.proteins_dir,
        include_metadata=not args.skip_metadata,
        external=args.external,
        prefix=prefix
    )

    # check results against known nitrogenase
    if args.test:
        test(checked=pd.read_csv(f'../results/final/{prefix}final.csv'))

    return 0

if __name__ == "__main__":
    sys.exit(main())

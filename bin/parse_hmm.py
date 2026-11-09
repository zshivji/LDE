#!/usr/bin/env python3
"""Parse hmmsearch outputs and filter tophits results.

Inputs:
    --hits             Path to hmmsearch.domtblout results from the HMM search step.
    --outdir           Output directory for hits.feather and hits.csv.

Examples:
    python bin/parse_hmm.py \
        --hits ../results/hmmsearch_results/archaea/ \
        --outdir ../results/hmmsearch_results/archaea
"""

import argparse
import re
import warnings
import pandas as pd
import numpy as np
from pathlib import Path
from Bio import SearchIO
from scipy.cluster.hierarchy import linkage, fcluster
from helper import filter_groups_by_unique_counts

warnings.filterwarnings("ignore", category=FutureWarning)


def genome_id_from_path(file: Path) -> str | None:
    """Extract a recognized accession or runner input name from a domtblout path."""
    accession_match = re.search(r"([\w]+_[\w]+_[\d]+\.[\d])", str(file))
    if accession_match:
        return accession_match.group()
    if "__job_" in file.name:
        return file.name.split("__job_", 1)[0]
    if "AFDB" in file.name:
        return 'AFDB'
    return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse hmmsearch outputs and filter tophits. "
            "Writes hits.feather/hits.csv and tophits.feather/tophits.csv."
        )
    )
    parser.add_argument(
        "--hits",
        help="Path or glob to hmmsearch.out results from the HMM search step.",
    )
    parser.add_argument(
        "--outdir",
        help="Output directory for hmm hits that meet the following criteria (1) significant evalue, (2) have at least x nif genes within a y gene range.",
    )
    parser.add_argument(
        "--output_file",
        help=(
            "Optional exact CSV output path. If omitted, the historical "
            "<outdir>_hits.csv naming convention is used."
        ),
    )

    parser.add_argument(
        "--min_genes",
        type=int,
        default=3,
        help="Minimum number of nif genes required to be considered a tophit cluster (default: 3)."
    )

    parser.add_argument(
        "--gene_range",
        type=int,
        default=15,
        help="Maximum gene range (+/- from gene of interest) for nif genes to be considered a cluster (default: 15).")
    parser.add_argument(
        "--skip_taxonomy",
        action="store_true",
        help="Do not read or join repository taxonomy metadata (used by API runner jobs).",
    )
    parser.add_argument(
        "--refine",
        action="store_true",
        help="Refine hits by clustering into operons and maintaining operons with at --min_genes (3 for nitrogenase is default).",
    )
    
    return parser.parse_args()


# Get taxonomy from GTDB, https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/
# NCBI taxonomy from metadata files

def parse_hits(hits_path: Path) -> pd.DataFrame:

    # Create empty lists
    result_target = []
    query_id = []
    hit_id = []
    contig = [] 
    pos_num = []
    operon = []
    evalue = []
    bitscore = []
    bias = []
    dom_start = []
    dom_end = []
    location = []
    orientation = []
    flag = []

    ### NEED TO UPDATE TO ACCOMODATE ALL INPUT TYPES
    ind_map = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' # map domain index to letter for easier handling
    def append_hit(genome_id: str, gene: str, hit, hsp) -> None:
        # Save hit domains and relevant metadata
        result_target.append(genome_id)
        query_id.append(gene)
        hit_id.append(hit.id + "_" + ind_map[hsp.domain_index-1]) # hit in the form of contig_geneNo_domA/B/C
        parts = hit.id.split("_")
        if len(parts) > 1: # for AFDB (not contig, geneNo info)
            contig.append(str(("_").join(parts[:-1])))
            pos_num.append(int(parts[-1]))
        else: #for AFDB (not contig, geneNo info)
            contig.append(np.nan)
            pos_num.append(np.nan)
        operon.append(np.nan)
        evalue.append(hsp.evalue) # domain independent evalue
        bitscore.append(hsp.bitscore) # domain bitscore
        bias.append(hsp.bias) # domain bias
        dom_start.append(hsp.env_start) 
        dom_end.append(hsp.env_end)

        # grab description parameters
        if genome_id == 'AFDB':
            location.append('')
            orientation.append('')
        else:
            s = r"# ([0-9]+) # ([0-9]+) # (-?1) #"  ### NEED TO UPDATE TO ACCOMODATE ALL INPUT TYPES
            # location = full protein length
            location.append(
                re.match(s, hit.description).group(1)
                + "-"
                + re.match(s, hit.description).group(2) 
            )
            orientation.append(re.match(s, hit.description).group(3))
 

    # Parse through files in output directory
    for file in hits_path.glob("*.domtblout"):
        # Prefer a GTDB accession when present. Runner intermediates use
        # <input-stem>__job_<uuid>_nif.domtblout, so uploaded files without a
        # GTDB accession fall back to their sanitized input stem.
        genome_id = genome_id_from_path(file)
        if genome_id is None:
            print(f"Could not determine GenomeID from {file}", flush=True)
            continue

        # Parse file using SearchIO/HmmerIO w/ --domtblout option
        for result in SearchIO.parse(file, "hmmsearch3-domtab"):
            # grab gene name, ex. nifHDK, pchlide
            gene = 'pqq'#re.findall(r"([a-zA-Z]+)", result.id)[0]
            
            for hit in result.hits:
                # Check for positive bitscore and append the data to the corresponding lists
                if hit.bitscore <= 0 or hit.evalue >= 0.0001:
                    continue

                for hsp in hit.hsps:
                    if hsp.bitscore <= 0 or hsp.evalue >= 0.0001:
                        continue
                    append_hit(genome_id, gene, hit, hsp)

                    # check if bitscore >> bias (same order of magnitude) as bitscore
                    if hit.bias != 0 and hit.bitscore / hit.bias > 10:
                        flag.append(0)
                    elif hit.bias == 0:
                        flag.append(0)
                    else:
                        flag.append(1)

    # create and store dataframe
    hits = pd.DataFrame(
        {
            "GenomeID": result_target,
            "gene": query_id,
            "hit": hit_id,
            "contig": contig,
            "pos_num": pos_num,
            "operon": operon,
            "E-value": evalue,
            "Bit Score": bitscore,
            "Bias": bias,
            "Domain Start": dom_start,
            "Domain End": dom_end,
            "Location": location,
            "Orientation": orientation,
            "Flag_Bias": flag
        }
    )

    hits.drop_duplicates(inplace=True)
    return hits

def parse_tophits(
    hits: pd.DataFrame,
    min_genes: int,
    gene_range: int,
    include_taxonomy: bool = True,
) -> pd.DataFrame:

    # pre-filter to remove contigs strands with less than min_genes unique nif genes)
    hits = filter_groups_by_unique_counts(
        hits,
        group_cols=["GenomeID", "contig", "Orientation"],
        requirements={
        "gene": min_genes,
        "hit": min_genes
        },
        exclude='nifB'
    )

    # store genomes meeting criteria in new dataframe
    genomes_to_keep = pd.DataFrame(columns=hits.columns)

    # iterate through each genome and contig and split into operons
    for (genome, contig), tmp in hits.groupby(["GenomeID", "contig"]):

            # split contigs into operons when 1) strand changes or 2) distance between genes is greater than limit
            operon_counter = 1

            # sort hits by gene position and orientation
            # for each strand, cluster genes by distance and assign operon labels
            for _, strand in tmp.groupby("Orientation"):

                if len(strand) < 2: # skip strands with only nifB (linkage clustering breaks)
                    clusters = np.array([operon_counter])
                else: 
                    # use hierarchical clustering to group neighboring genes
                    strand = strand.sort_values("pos_num")
                    pos = strand["pos_num"].to_numpy().reshape(-1, 1)

                    # cluster genes by distance
                    Z = linkage(pos, method = 'single') # try method = 'ward' but singe prob better
                    clusters = fcluster(Z, t = gene_range, criterion = 'distance')
                    clusters = clusters.astype(int) + operon_counter -1 # convert to int for easier handling

                strand["operon"] = clusters

                # save clusters to tmp2
                operon_counter += clusters.max()
                genomes_to_keep = pd.concat([genomes_to_keep, strand])

    # refilter to only keep clusters with at least 3 unique genes, unless a gene is excluded from the requirement (ex. nifB)
    genomes_to_keep = filter_groups_by_unique_counts(
        genomes_to_keep,
        group_cols=["GenomeID", "contig", "operon"],
        requirements={
        "gene": min_genes,
        "hit": min_genes
        },
        exclude='nifB'
    )

    # record which gene has the highest bitscore for each protein
    genomes_to_keep["top_hit"] = (genomes_to_keep.sort_values("Bit Score", ascending=False).groupby(["GenomeID", "contig", "hit"])["gene"].transform("first"))

    if include_taxonomy:
        # Database-build mode enriches representative genomes with taxonomy.
        gtdb_taxonomy = pd.read_csv(
            "GTDB_metadata.gz",
            sep="\t",
            usecols=["accession", "gtdb_taxonomy"],
        ).rename(columns={"accession": "GenomeID", "gtdb_taxonomy": "GTDB"})
        genomes_to_keep = pd.merge(
            genomes_to_keep,
            gtdb_taxonomy,
            on="GenomeID",
            how="left",
        )
    else:
        # Keep the downstream schema stable without requiring repository metadata.
        genomes_to_keep["GTDB"] = ""

    return genomes_to_keep


def main() -> None:
    args = parse_args()

    if args.hits and args.outdir:
        print(f"Parsing HMM hits...\n\t\t from {args.hits}\n\t\t saving results to {args.outdir}_hits.csv", flush=True)

        hits = parse_hits(Path(args.hits))

        if args.refine: # only refine hits if --refine is specified (when we have useful info in description or looking for many genes in an operon)
            hits = parse_tophits(
                hits,
                args.min_genes,
                args.gene_range,
                include_taxonomy=not args.skip_taxonomy
            )

        # export
        Path(args.outdir).mkdir(parents=True, exist_ok=True)
        output_file = Path(args.output_file) if args.output_file else None
        output_path = output_file or Path(str(args.outdir) + "_hits.csv")
        output_path.parent.mkdir(parents=True, exist_ok=True)
        hits.to_csv(output_path, index=False)
        print(f"Saved {hits.shape[0]} hits to {output_path}", flush=True)

    elif args.hits or args.outdir:
        raise SystemExit("Both --hits and --outdir are required to parse HMM hits."
                         "Please provide both arguments and try again." )

if __name__ == "__main__":
    main()

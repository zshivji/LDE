import pandas as pd
import os
import glob
from tqdm import tqdm
import argparse
import json

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# get args (reload_fasta)
parser = argparse.ArgumentParser(description='Align MDH gene hits')
parser.add_argument('--reload_fasta', action='store_true', help='Whether to reload fasta sequences (True/False)')
args = parser.parse_args()

# grab both archaea + bacteria hits
MDH_archaea = pd.read_feather('../results/hmmsearch_results/MDH_archaea.feather')
MDH_bacteria = pd.read_feather('../results/hmmsearch_results/MDH_bacteria.feather')
MDH = pd.concat([MDH_archaea, MDH_bacteria])
#MDH = MDH_archaea.copy()

# read json
config = json.load(open('LDE-config.json', 'r'))

MDH.reset_index(inplace = True)
MDH.set_index(['GenomeID'], inplace = True)
MDH['Seq'] = ''

if args.reload_fasta: # should skip if fasta sequences have already been extracted
    
    # get fasta sequences for each gene & export to fasta
    for gene in config:
        print('getting ' + gene + ' fasta', flush=True)

#        for _ in tqdm(range(len(gene))): # progress bar
        records = []
        fasta = MDH[MDH.Gene == gene]
        for genome,hit in fasta.iterrows():
            hit = hit.Hit
            file = glob.glob(f"../../diazoDB-HPC/protein_faa_reps/*/{genome}_protein.faa")[0]
                
            for result in SeqIO.parse(file, "fasta"):
                if result.id == hit:
                    # store seq
                    fasta.loc[genome, 'Seq'] = str(result.seq)
                    # convert to seqrecord
                    record = SeqRecord(result.seq, id=genome, description=hit)
                    records.append(record)
                    # exit loop once sequence is found
                    break
                    
            # Write the records to a FASTA file ("w" means override)
        with open("../results/tmp/" + gene + ".fasta", "w") as output_handle:
            SeqIO.write(records, output_handle, "fasta")

#align fasta files
for gene in config:
    print("aligning "+ gene, flush=True)

    # remove exisiting fasta + aln files (we append to them below)
#    os.system("rm ../results/fasta_splits/*.fa*")
#    os.system("rm ../results/fasta_splits/*.aln")

    fasta = MDH[MDH.Gene == gene]
    num = int(fasta.shape[0]/200)+1 # how many splits
    print(num, flush=True)
    # seqtk fails after 1024 so need to use seqkit
#    os.system(f"seqkit split2 -p {num} -O ../results/fasta_splits/ -P {gene}_split. ../results/tmp/{gene}.fasta") # split fasta file
    for i in range(num):
        if i > 598:
            if i < 1000: # seqkit only fills 000 so need to count differently after 1000
                os.system(f"seqtk subseq ../results/tmp/{gene}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/{gene}_split.{str(i+1).zfill(3)}.fasta") # add reference sequences
                os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/{gene}_split.{str(i+1).zfill(3)}.fasta > ../results/fasta_splits/{gene}_split.{str(i+1).zfill(3)}.aln") # align files
            else:
                os.system(f"seqtk subseq ../results/tmp/{gene}.fasta ../results/ref_seq.ids >> ../results/fasta_splits/{gene}_split.{str(i+1)}.fasta") 
                os.system(f"mafft --auto --quiet --thread 4 ../results/fasta_splits/{gene}_split.{str(i+1)}.fasta > ../results/fasta_splits/{gene}_split.{str(i+1)}.aln")


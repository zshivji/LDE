# LDE

## Lanthanide Dependent Enzyme Annotation
Automated annotation and curation of LDEs genes using profile hidden Markov models and conserved residue matching.

## Install

## Scripts

+ **1_align_mafft.sh** - Clusters (mmseqs easy-cluster) HMM seed sequences at 90% AAI. Aligns (MAFFT) clustered seeds sequences. Builds HMM profiles and combines into single file. Stores alignmnets in ***alignments/*** and profiles in ***HMMs/***.

+ **2_hmmsearch.sh** - Runs pHMM search against GTDB R220 all_rep_proteins_aa database (bacteria and archaea). Stores results as individual files in ***results/{archaea | bacteria}*** (per GTDB accession number) 

+ **3_parse_hmm.sh** - Calls (1) **Parse_hmm_results.py** and (2) **Parse_hmm_results.py**.

+ **4_conserved_res.sh** - Calls (1) ***aln_nif_hits.py***, (2) ***conserved-res.py***, (3) ***final-fasta-export.py***, and (4) ***diazoDB-check.py***.

+ **5_make_trees.sh** - 

+ **6_operon_org.sh**

+ **7_SSN.sh**

+ **Parse_hmm_results.py** - Parse pHMM search results (see below criteria), combines with GTDB taxonomy data, and stores results in ***hits.feather***
    + positive bit score
    + full sequence evalue must be significant (<0.01)
    + best domain evalue should be significant (<0.01)
        + otherwise flagged for manual review to check if it is distant homolog or just short repeats

+ **Parse_tophits.py**

+ **aln_nif_hits.py** 

+ **cluster_pos.py**

+ **conserved-res.py**

+ **final-fasta-export.py**

+ **diazoDB-check.sh**

+ **get-operon.py**

+ **operon-org-plot.py**


 PDB residues: [76, 124, 125, 128, 130, 174, 189, 190, 191, 192, 258, 276, 318, 320, 345, 416, 494, 557, 558]
                E.   C.   C.   V.   R.   T.   S.   G.   G.   E.   W.   N.   D.   D.   R.   D.   W.   G.   W
CDD residues: [76, 124, 125, 128, 130, 174, 189, 190, 192, 256, 258, 276, 280, 320, 345, 416, 494, 558] # missing D-x-D
                E.   C.   C.   V.   R.   T.  S.   G.   E.   T.   W.   N.   W.   D.   R.   D.   W.   W
custom residues: [76, 124, 125, 128, 130, 174, 190, 192, 258, 276, 318, 320, 345, 494, 558]
                  E.   C.   C.   V.   R.   T.   G.   E.   W.   N.   D.   D.   R.   W.   W

W G W --> W not conserved in all MDH, 2nd XoxF in M. extorquens lacks LT962688.1_4238 but G is conserved

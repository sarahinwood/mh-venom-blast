library("data.table")
library("dplyr")
library("ggplot2")

venom_genes_blastx <- fread("output/full_seq_blast_annots/blastx.outfmt6")
setnames(venom_genes_blastx, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13"), new=c("gene_id", "Blastx_hit", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score", "Full_BlastX_annot"))
##order based off evalue and bitscore
setorder(venom_genes_blastx, gene_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
venom_min_evalues <- venom_genes_blastx[,.SD[which.min(evalue)], by=gene_id]
fwrite(venom_min_evalues, "output/full_seq_blast_annots/blastx_res.csv")

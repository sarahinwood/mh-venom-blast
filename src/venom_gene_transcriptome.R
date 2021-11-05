library(data.table)
library(ggplot2)
library(VennDiagram)

trinotate_report <- fread("data/trinotate_annotation_report.txt")
contig_blast_results <- fread("output/venom_gene_blast_transcriptome/blastn.outfmt6")

setnames(contig_blast_results, old=c("V2"), new=c("transcript_id"))

blast_w_annotations <- merge(contig_blast_results, trinotate_report, by.x = "transcript_id", by.y= "transcript_id")
setnames(blast_w_annotations, old=c("V1", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("contig_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

##best transcriptome hit per venom gene
##order based off evalue and bitscore
setorder(blast_w_annotations, contig_id, evalue, -bit_score)
##extract result with lowest evalue for each peptide - what if multiple rows with lowest min?
min_evalues <- blast_w_annotations[,.SD[which.min(evalue)], by=contig_id]
fwrite(min_evalues, "output/venom_gene_blast_transcriptome/blastn_trinotate_annotations.csv")

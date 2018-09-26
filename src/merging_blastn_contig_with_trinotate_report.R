library(data.table)
library(ggplot2)

trinotate_report <- fread("data/trinotate_annotation_report.txt")
contig_blast_results <- fread("output/contig_blast_transcriptome/blastno6.csv")

setnames(contig_blast_results, old=c("V2"), new=c("transcript_id"))

blast_w_annotations <- merge(contig_blast_results, trinotate_report, by.x = "transcript_id", by.y= "transcript_id")
setnames(blast_w_annotations, old=c("V1", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("contig_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

fwrite(blast_w_annotations, "output/contig_blast_transcriptome/blastn_w_trinotate_annotations.csv")

##Removed duplicated contigs (i.e. those with multiple annotations - only kept best)

##read back in annotated contigs & play with file
contigs_annots <- fread("output/contig_blast_transcriptome/blastn_w_trinotate_annotations_duplications_removed.csv")
contigs_annots_GO <- contigs_annots[,tstrsplit(gene_ontology_blast, "`", fixed = TRUE, keep=1)]
contigs_last_GO_term <- contigs_annots_GO[,tstrsplit(V1, "^", fixed=TRUE, keep=3)]
contigs_last_GO_term$V1[is.na(contigs_last_GO_term$V1)] <- "none"
GO_counts <- data.table(table(contigs_last_GO_term))
setnames(GO_counts, old=c("contigs_last_GO_term", "N"), new=c("Last_GO_term", "count"))
ordered_GO_counts <- GO_counts[order(-GO_counts$count),]

ggplot(ordered_GO_counts, aes(x = reorder(ordered_GO_counts$Last_GO_term, -count), y=ordered_GO_counts$count))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  geom_col()+xlab("GO Annotation Term")+ylab("No. of Contigs")

fwrite(ordered_GO_counts, "output/contig_blast_transcriptome/ordered_GO_counts.csv")


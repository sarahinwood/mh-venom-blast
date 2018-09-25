library(data.table)

trinotate_report <- fread("data/trinotate_annotation_report.txt")
contig_blast_results <- fread("output/contig_blast_transcriptome/blastno6.csv")

setnames(contig_blast_results, old=c("V2"), new=c("transcript_id"))

blast_w_annotations <- merge(contig_blast_results, trinotate_report, by.x = "transcript_id", by.y= "transcript_id")
setnames(blast_w_annotations, old=c("V1", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("contig_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

fwrite(blast_w_annotations, "output/contig_blast_transcriptome/blastn_w_trinotate_annotations.csv")

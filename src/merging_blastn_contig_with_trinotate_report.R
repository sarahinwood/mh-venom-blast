library(data.table)
library(ggplot2)
library(VennDiagram)

trinotate_report <- fread("data/trinotate_annotation_report.txt")
contig_blast_results <- fread("output/contig_blast_transcriptome/blastno6.csv")

setnames(contig_blast_results, old=c("V2"), new=c("transcript_id"))

blast_w_annotations <- merge(contig_blast_results, trinotate_report, by.x = "transcript_id", by.y= "transcript_id")
setnames(blast_w_annotations, old=c("V1", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("contig_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

fwrite(blast_w_annotations, "output/contig_blast_transcriptome/blastn_w_trinotate_annotations.csv")

##Removed duplicated contigs (i.e. those with multiple annotations - only kept best)

##read back in annotated contigs & play with file
contigs_annots <- fread("output/contig_blast_transcriptome/dedup_blastn_w_trinotate_annotations.csv", na.strings = ".")
contigs_annots_GO <- contigs_annots[,tstrsplit(gene_ontology_blast, "`", fixed = TRUE, keep=1)]
contigs_last_GO_term <- contigs_annots_GO[,tstrsplit(V1, "^", fixed=TRUE, keep=3)]
contigs_last_GO_term$V1[is.na(contigs_last_GO_term$V1)] <- "none"
GO_counts <- data.table(table(contigs_last_GO_term))
setnames(GO_counts, old=c("contigs_last_GO_term", "N"), new=c("Last_GO_term", "count"))
ordered_GO_counts <- GO_counts[order(-GO_counts$count),]

ggplot(ordered_GO_counts, aes(x = reorder(ordered_GO_counts$Last_GO_term, -count), y=ordered_GO_counts$count))+
  theme(axis.text.x = element_text(angle = 65, hjust = 1))+
  geom_col()+xlab("GO Annotation")+ylab("No. of Contigs")

fwrite(ordered_GO_counts, "output/contig_blast_transcriptome/ordered_GO_counts.csv")

##Get info for annotation venn diagram
pfam <- contigs_annots[!is.na(Pfam), unique(`contig_id`)]
blastx <- contigs_annots[!is.na(sprot_Top_BLASTX_hit), unique(`contig_id`)]
blastp <- contigs_annots[!is.na(sprot_Top_BLASTP_hit), unique(`contig_id`)]
kegg <- contigs_annots[!is.na(Kegg), unique(`contig_id`)]
number.contigs <- contigs_annots[!is.na(`contig_id`), length(unique(`contig_id`))]

#Draw Venn Diagram
Set1 <- RColorBrewer::brewer.pal(3, "Set1")
vd <- venn.diagram(x = list("Pfam"=pfam, "BlastX"=blastx, "Kegg"=kegg), filename=NULL, fill=Set1, alpha=0.5, cex = 1, cat.cex=1, lwd=1, main=paste("Total Number of Annotated Contigs = ", number.contigs))
grid.newpage()
grid.draw(vd)

#Sum of genes with any annotation
long.annotationreport <- melt(contigs_annots,id.vars = "contig_id", measure.vars = c("sprot_Top_BLASTX_hit", "sprot_Top_BLASTP_hit", "Pfam",  "eggnog", "Kegg"))
any.annotations <- long.annotationreport[,.(any_annotations = any(!is.na(value))),by=`contig_id`]
any.annotations[,length(unique(`contig_id`))]
any.annotations[,sum(any_annotations)]

#Sum of genes with any TransRate predicted protein
transrate.report <- melt(contigs_annots, id.vars="contig_id", measure.vars = c("prot_id"))
any.transrate <- transrate.report[,.(transrate_annotation = any(!is.na(value))),by=`contig_id`]
any.transrate[,length(unique(`contig_id`))]
any.transrate[,sum(transrate_annotation)]
sum(any.transrate$transrate_annotation==FALSE)

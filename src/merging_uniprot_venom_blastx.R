library(data.table)
library(Biostrings)

##read in contig annotations from uniprot venom db and rename columns
uniprot_venom_results <- fread("output/uniprot_venom_contigs/blastx.outfmt6")
setnames(uniprot_venom_results, old=c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12"), new=c("contig_id", "uniprot_venom_id", "%_identical_matches", "alignment_length", "no_mismatches", "no_gap_openings", "query_start", "query_end", "subject_start", "subject_end", "evalue", "bit_score"))

##Read in uniprotdb, keep only rownames and format as needed for merge
uniprot_db <- data.frame(readDNAStringSet("data/uniprot_toxin_venom_db.fasta"))
uniprot_db <- setDT(uniprot_db, keep.rownames = TRUE)
uniprot_db_names <- data.table(uniprot_db$rn)
uniprot_db_names_split <- uniprot_db_names[,tstrsplit(V1, "|", fixed = TRUE)]
uniprot_db_names_split$V1 <- NULL
setnames(uniprot_db_names_split, old=c("V2", "V3"), new=c("uniprot_venom_id", "uniprot_annotation"))

contig_uniprot_blastx <- merge(uniprot_venom_results, uniprot_db_names_split, by.x= "uniprot_venom_id", by.y = "uniprot_venom_id")
##Check sum of unique contig IDs
contig_uniprot_blastx[, unique(contig_id)]

fwrite(contig_uniprot_blastx, "output/uniprot_venom_contigs/contig_uniprot_blastx_annotations.csv")

##then manually sort and keep only lowest E-value annotation for each contig
dedup_annots <- fread("output/uniprot_venom_contigs/dedup_contig_uniprot_blastx_annotations.csv")
##check no. contigs = no. unique contig ids from earlier
unique <- dedup_annots[, unique(contig_id)]
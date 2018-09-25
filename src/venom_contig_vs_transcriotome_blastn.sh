#!/usr/bin/env bash

set -eu

blastn_results=output/contig_blast_transcriptome/blastn.outfmt6

outdir="$(dirname "${blastn_results}")"
if [[ ! -e "${outdir}" ]]; then
	mkdir -p "${outdir}"
fi

blastn \
	-query data/mh_venom_sequencing_contigs.fasta \
	-db output/mh_length_fil_db/mh_length_fil_db \
	-num_threads 50 \
	-outfmt 6 > output/contig_blast_transcriptome/blastn.outfmt6
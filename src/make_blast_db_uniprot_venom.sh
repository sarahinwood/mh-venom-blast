#!/usr/bin/env bash

set -eu

blast_db=output/mh_length_fil_db/mh_length_fil.db

outdir="$(dirname "${blast_db}")"
if [[ ! -e "${outdir}" ]]; then
	mkdir -p "${outdir}"
fi

makeblastdb \
	-in data/mh_transcriptome_lengthfil.fasta \
	-dbtype nucl \
	-title mh_length_filtered_transcriptome \
	-out output/mh_length_fil_db/mh_length_fil_db \
	-parse_seqids
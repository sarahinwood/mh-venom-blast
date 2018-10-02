#!/usr/bin/env bash

set -eu

blast_db=output/uniprot_venom_db/uniprot_venom.db

outdir="$(dirname "${blast_db}")"
if [[ ! -e "${outdir}" ]]; then
	mkdir -p "${outdir}"
fi

makeblastdb \
	-in data/uniprot_toxin_venom_db.fasta \
	-dbtype prot \
	-title uniprot_venom_db \
	-out output/uniprot_venom_db/uniprot_venom_db \
	-parse_seqids
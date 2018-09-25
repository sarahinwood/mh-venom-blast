#!/usr/bin/env bash

set -eu

blastx_results=output/full_seq_blast_annots/blastx.outfmt3

outdir="$(dirname "${blastx_results}")"
if [[ ! -e "${outdir}" ]]; then
	mkdir -p "${outdir}"
fi

blastx \
	-query data/mh_venom_nt.fasta \
	-db bin/db/uniprot_sprot.pep \
	-num_threads 50 \
	-outfmt 3 > output/full_seq_blast_annots/blastx.outfmt3

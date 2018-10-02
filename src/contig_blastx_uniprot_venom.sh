#!/usr/bin/env bash

set -eu

blastx_results=output/uniprot_venom_contigs/blastx.outfmt6

outdir="$(dirname "${blastx_results}")"
if [[ ! -e "${outdir}" ]]; then
	mkdir -p "${outdir}"
fi


blastx \
	-query data/mh_venom_sequencing_contigs.fasta \
	-db output/uniprot_venom_db/uniprot_venom_db \
	-num_threads 50 \
	-outfmt 6 > output/uniprot_venom_contigs/blastx.outfmt6
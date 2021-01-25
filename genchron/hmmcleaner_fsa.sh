#!/bin/bash

set -euo pipefail

help="USAGE: hmmcleaner_fsa.sh <input prot fasta> [<input cDNA>]

If prot fasta is named:      'inputdir/realign/inputbase_protfsa.fa',
unless given, use cDNA from: 'inputdir/inputbase_fsa.fa'
and output to:               'inputdir/inputbase_fsa_hmm.fa'
"
# Modified from hmmcleaner_codons to reuse the protein alignment already existing.

suffix='fsa'

OPTIND=1
while getopts 's:' opt; do
    case "$opt" in
        s) suffix="$OPTARG" ;;
        *) echo "$help" && exit 2 ;;
    esac
done
shift $(( OPTIND-1 ))

[[ -n "${1:-}" ]] || { echo "$help" >&2; exit 2; }

protinput="$1"

# Remove extension
if [[ "$protinput" = *.* ]]; then
    inputroot="${protinput%.*}"
else
    inputroot="$protinput"
fi
# Find the dirname and basename
if [[ "$inputroot" = */* ]]; then
    inputdir="${inputroot%/*}"
    inputbase="${inputroot##*/}"
else
    inputdir='.'
    inputbase="$inputroot"
fi

# Output directory is the parent folder of "realign/"
outputroot="${inputdir%/realign}/${inputbase}"

# The second script argument ($2) gives the source file for cDNA, otherwise use a default
input_cdna="${2:-${outputroot%_prot$suffix}_${suffix}.fa}"

SCRIPTS="$HOME/scripts"

[[ -d "$SCRIPTS" ]] || { echo "Script directory does not exist." >&2 && exit 1; }

# Run HmmCleaner on the proteins
HmmCleaner.pl --log_only "${protinput}" > "${inputdir}/${inputbase}_cleanedcounts.txt"

total=0

while read -r filename seqname count; do
    (( total+=count )) || true
done < "${inputdir}/${inputbase}_cleanedcounts.txt"

outfile="${outputroot%_prot$suffix}_${suffix}_hmm.fa"

if (( total )); then
# Find out which of the codon position correspond to the trimmed aa.
$SCRIPTS/seqtools/fillpositions.py -c -o "$outfile" \
         "${inputdir}/${inputbase}_hmm.log" "${input_cdna}"

elif [[ -f "$input_cdna" ]]; then
    ln -s "$input_cdna" "$outfile"
else
    cat "$input_cdna" > "$outfile"
fi


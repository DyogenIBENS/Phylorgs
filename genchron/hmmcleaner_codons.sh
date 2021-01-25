#!/bin/bash

set -euo pipefail

usage="USAGE: hmmcleaner_codons.sh [-t] [-p <prot fasta name>] <input cDNA fasta>"

help="$usage

-t: do translate the input cDNA
-p: the protein alignment (if none, the cDNA will be translated).
-s: suffix to strip from the basename to construct the output (protein file)

Output the filtered cDNA alignment
"

protinput=""
dotranslate=0
suffix='_prot'

OPTIND=1
while getopts 'htp:s:' opt; do
    case "$opt" in
        t) dotranslate=1 ;;
        p) protinput="$OPTARG" ;;
        s) suffix="$OPTARG" ;;
        h) echo "$help" && exit ;;
    esac
done
shift $(( OPTIND-1 ))

inputal="${1:-}"

[[ -n "${1:-}" ]] || { echo -e "ERROR: Missing arg: -p=${protinput} inputal=${inputal}.\n$usage" >&2 && exit 2; }


# Remove extension
if [[ "$inputal" = *.* ]]; then
    inputbase="${1%.*}"
else
    inputbase="$inputal"
fi


if [[ -z "$protinput" ]]; then
    protinputbase="${inputbase}_prot"
    protinput="${inputbase}_prot.fa"
else
    # the output name should be based on the protein input
    protinputbase="${protinput%.fa*}"
    inputbase="${protinputbase%$suffix}"
fi


SCRIPTS="$HOME/scripts"

[[ -d "$SCRIPTS" ]] || { echo "Script directory does not exist." >&2 && exit 1; }

if (( dotranslate )); then

# Ungap and Convert to protein alignment
#$SCRIPTS/seqtools/ungap.py "$inputal" |\
$SCRIPTS/seqtools/seq_translate.py "${inputal}" "${protinput}"

fi

# Run HmmCleaner on the proteins
HmmCleaner.pl --log_only "${protinput}" >"${protinputbase}_cleanedcounts.txt"

total=0

while read -r filename seqname count; do
    (( total+=count )) || true
done < "${protinputbase}_cleanedcounts.txt"

if (( total )); then
    # Find out which of the codon position correspond to the trimmed aa.
    $SCRIPTS/seqtools/fillpositions.py -c -o "${inputbase}_hmm.fa" \
            "${protinputbase}_hmm.log" "${inputal}"
elif [[ -f "$inputal" ]]; then
    ln -s "$inputal" "${inputbase}_hmm.fa"
else
    cat "$inputal" > "${inputbase}_hmm.fa"
fi

# Back-translate
#treebest backtrans "${inputal}_prot_hmm.fa" "${}" > ""

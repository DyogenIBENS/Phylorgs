#!/bin/bash

set -euo pipefail

usage="USAGE: hmmcleaner_codonzilla.sh -p <input prot fastas> [ -n <n lines for 1 dataset in cdna> ] [ -o <substitution-to-name-output> | -c <input cDNA fastas> ]"

help="$usage

Output the filtered cDNA alignment
"

SCRIPTS="$HOME/scripts"
[[ -d "$SCRIPTS" ]] || { echo "Script directory does not exist." >&2 && exit 1; }

dryrun=0
input_cdnas=()
input_prots=()
subst=''
nlinesdata=0

(( "$#" > 0 )) || { echo "$help" ; exit 2; }

add_prot=0
add_cdna=0


while (( $# )) ; do
    if [[ "$1" = "-o" ]]; then
        shift
        [[ -n "${1:-}" ]] || { echo "ERROR: value needed for -s" >&2; exit 2; }
        subst="$1"
        add_prot=0
        add_cdna=0
    elif [[ "$1" = "-n" ]]; then
        shift
        nlinesdata=$1
        [[ "${1:-}" =~ ^[0-9]+$ ]] || { echo "ERROR: -n should be an integer" >&2 ; exit 2; }
        add_cdna=0
        add_prot=0
    elif [[ "$1" = "-d" ]]; then
        dryrun=1
        add_cdna=0
        add_prot=0
    elif [[ "$1" = "-p" ]]; then
        add_prot=1
        add_cdna=0
    elif [[ "$1" = "-c" ]]; then
        add_cdna=1
        add_prot=0
    elif [[ "${1:0:1}" = "-" ]]; then
        echo "ERROR: Invalid option: '$1'" >&2
        exit 2
    elif (( !add_prot )) && (( !add_cdna )); then
        echo "ERROR: Unexpected argument: '$1'" >&2
        exit 2
    elif (( add_prot )); then
        input_prots+=("$1")
    else
        input_cdnas+=("$1")
    fi
    shift
done

if [[ -n "$subst" ]] && [[ ${#input_cdnas[@]} -eq 0 ]]; then
    for (( i=0; i<${#input_prots[@]}; i++ )); do
        input_cdnas[$i]=$(echo -n "${input_prots[$i]}" | sed -r "$subst")
    done
fi

(( dryrun )) && echo -n "[DRYRUN]" >&2
echo -e "DEBUG: cdnas=${input_cdnas[@]}\n     prots=${input_prots[@]}\n"\
  "    outsubst='$subst'\n     nlinesdata=${nlinesdata}" >&2

if (( !dryrun )); then
    HmmCleaner.pl --log_only ${input_prots[@]}
fi

#i=0
#ith_cleanedcounts=""


#if (( ${#input_cdnas[@]} < ${#input_prots[@]} )); then
#    for i in ${!input_prots[@]}; do
#        sed -n "$((i*nlinesdata)),$(( (i+1)*nlinesdata - 1 ))p" ${input_cdnas}
#    done
#fi

#while IFS= read -r line; do

    #total=0

    #while read -r filename seqname count; do
    #    (( total+=count )) || true
    #done < "_cleanedcounts.txt"

for i in ${!input_prots[@]}; do
    protinput=${input_prots[$i]}
    # Find out which of the codon position correspond to the trimmed aa.
    inputal=$(echo -n "$protinput" | sed -r "$subst")
    if (( nlinesdata )); then
        echo "DEBUG: input #$i: $protinput, seqconv ${input_cdnas:-}:$i -> ${inputal%.*}_hmm.fa" >&2
        if (( !dryrun )); then
            $SCRIPTS/seqtools/fillpositions.py -c -o "${inputal%.*}_hmm.fa" \
                "${protinput%.*}_hmm.log" <(seqconv.py -t fasta "${input_cdnas}:$i")
            #<(sed -n "$((i*nlinesdata)),$(( (i+1)*nlinesdata - 1 ))p" ${input_cdnas})
        fi
    else
        echo "DEBUG: input #$i : $protinput, ${input_cdna[i]} -> ${inputal%.*}_hmm.fa" >&2
        if (( !dryrun )); then
            $SCRIPTS/seqtools/fillpositions.py -c -o "${inputal%.*}_hmm.fa" \
                    "${protinput%.*}_hmm.log" "${input_cdnas[$i]}"
        fi
    #elif [[ -f "$inputal" ]]; then
    #    ln -s "$inputal" "${inputbase}_hmm.fa"
    #else
    #    cat "$inputal" > "${inputbase}_hmm.fa"
    fi
done

if (( dryrun )); then echo "DRYRUN finished." >&2; fi
#done <<< _cleanedcounts.txt

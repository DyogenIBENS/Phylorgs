#!/usr/bin/env bash

set -euo pipefail

help="
Convert a newick tree with PhylTree notations (pipes and spaces in labels, multiline, internal node names) to a 'standard' newick format

    ${0} [-s space_repl] [-b] inputtree

-s  character to replace spaces in labels [.]
-b  even more basic: no internal node labels, no root length.

NOTE: the input must be indented: one label per line. See indent_nwk.py
"

space_repl='.'
basic=0

while getopts 's:b' opt; do
    case "$opt" in
        s) space_repl="$OPTARG" ;;
        b) basic=1 ;;
        *) echo "$help"; exit 2 ;;
    esac
done
shift $((OPTIND-1))

inputtree=${1:-}

[ -n "$inputtree" ] || { echo "$help"; exit 2; }


if [ "$basic" -eq 0 ]; then

sed -r 's/([A-Za-z0-9]+) /\1'"${space_repl}"'/g; s/\|[A-Za-z0-9_.|-]+:/:/' \
    "${inputtree}" \
    | sed ':a; N;s/\n\s*//; ta'

else
# Additionally, remove internal node labels, and the root branch length.
sed -r -e 's/([A-Za-z0-9]+) /\1'"${space_repl}"'/g  #Replace spaces in labels' \
    -e 's/\|[A-Za-z0-9 _.|-]+:/:/  #Drop pipe-separated alternative names, preceeding colon' \
    -e 's/\)[A-Za-z0-9 _.-]+([:;])/)\1/  #Remove internal node labels' \
    -e 's/:[0-9.-];/;/  #Remove root length' \
    "${inputtree}" \
    | sed ':a; N;s/\n\s*//; ta  #Read the whole file and join lines'
fi

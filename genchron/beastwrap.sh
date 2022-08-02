#!/usr/bin/env bash

#echo -e "PATH=$PATH\nUSER=$USER\nTMPDIR=$TMPDIR\nLD_LIBRARY_PATH=$LD_LIBRARY_PATH\nhostname:$(hostname -d)"
#exit

# HARDCODE (Condor cluster won't let me use the environment value)
[ -d "/localtmp" ] && export TMPDIR="/localtmp/$USER"
[ -d "$TMPDIR" ] || mkdir "$TMPDIR"


set -euo pipefail


usage="    beastwrap.sh [-d <temp subdir>|-a|-t] -b <beast> -p -- <beast arguments> [treeannotator <treeannotator options>]

    -d <temp subdir>: output in this directory (in \$TEMPDIR, or /tmp).
    -a (automatic): use a directory in \$TEMPDIR named with the dirname of the xml file
    -t (mktemp): create a directory using mktemp.
    -p : print stdout directly to screen instead of saving.
    <beast arguments>: the xml file must be the last argument

    optional: treeannotator <treeannotator options>: run treeannotator with the given options
    output file will have the suffix '.annot.tree'

beast is run with '-working', producing output files in <temp subdir>.
stdout is redirected as well.
loganalyser is run with burnin 0 (output to .summary).
Output traces are bzipped and transfered back to the template file directory
(overwriting existing bz2 outputs).


# xml preparation:
    - Use fileName="'"'"\$(filebase).log"'"'" in tracelog, and treelog (.trees)
    - Disable state saving in the '<run><state>' element with 'storeEvery="'"'"-1"'"'"'
    - optionally, write less often to stdout (screenlog)
"


# First argument: the temporary directory for running Beast
currenttmpdir="${TMPDIR:-/tmp}"
beast='beast'
printstdout=0

# First check if there is enough space in the tmp dir (>100M to be safe)
available="$(df "$currenttmpdir" | awk 'NR==2 {print $4}')"

(( "$available" >= 102400 )) || { echo "Available space < 100M in $currenttmpdir: $available. Abort." >&2; exit 3; }

#OPTIND=1
while getopts 'b:d:ph' opt; do
    case "$opt" in
        h) echo "$usage" && exit ;;
        d) subdir="$OPTARG" ;;
        b) [[ "$OPTARG" =~ ^eagle ]] && break || beast="$OPTARG" ;;
        p) printstdout=1 ;;
        ?) break ;;
    esac
done

# Special case: if we encountered a -beagle (beast) option, keep it for beast args: OPTIND-2.
# Otherwise it is recommended to use '--' before beast args.
if [ -n "${OPTARG-}" ]; then shift $(( OPTIND - 2 )); else shift $(( OPTIND - 1 )); fi

[ "$#" -eq 0 ] && { echo "$usage"; exit 2; }
[ -z "${subdir:-}" ] && { echo 'Give a temporary directory to use inside $TMPDIR' >&2; exit 2; }
# suggestion:
# 1. automatically use dirname of the xmlfile as subdir.
# 2. or make a new safe temp name with mktemp -d

beast_args=()
while [ "${1:-treeannotator}" != "treeannotator" ]; do
    beast_args+=("$1")
    [[ "$1" =~ .+\.xml$ ]] && beast_xml="$1"
    shift
done

[ -n "$beast_xml" ] || { echo "xml file argument not given!" >&2; exit 2; }

srcbase="$(basename "$beast_xml" .xml)"
finaldir="$(dirname "$beast_xml")"
workdir="$currenttmpdir/$subdir"
tmproot="$workdir/$srcbase"

echo -e "DEBUG:
    LD_LIBRARY_PATH=$LD_LIBRARY_PATH
    beast=${beast}
    beast_args=${beast_args[@]}
    workdir=$workdir
    srcbase=$srcbase
    disk space left=$available" >&2

if [ "$#" -gt 0 ]; then
    shift
    treeannot_args="$@ $tmproot.trees $tmproot.annot.tree"
    echo "treeannot_args=$treeannot_args" >&2
else
    treeannot_args=
fi

[ -r "$beast_xml" ] || { echo "xml file not found or not readable!" >&2; exit 2; }
# Check that '$(filebase) is in use in the xml'
n_matches="$(egrep -c '<logger.*tr(ace|ee)log.*file[Nn]ame="\$\(filebase(=[A-Za-z0-9/_.-]+)?\)\.(log|trees)"' "$beast_xml" || true)"
[ "$n_matches" -lt 2 ] && { echo "xml file $beast_xml must use \$(filebase).(log|trees) for tracelog and treelog" >&2; exit 4; }

# Using mkdir and not mktemp for each single beast run, because I don't want to
# create 1000 tmp dirs, just one per batch
mkdir -p "$workdir"

# If output files already exist, and -resume is given,
# move .log.bz2, .trees.bz2 and .state to the tmp dir, and uncompress.
#TODO: check if uncompressed traces also exist.
if [[ "${beast_args[@]}" =~ (^|[[:space:]])-resume([[:space:]]|$) ]]; then
    if [ -e "$finaldir/$srcbase.trees.bz2" ] && [ -e "$finaldir/$srcbase.log.bz2" ]; then
        if [ -e "$finaldir/$srcbase.trees" ] || [ -e "$finaldir/$srcbase.log" ]; then
            echo 'WARNING: Resuming from .bz2 files, but uncompressed .trees or .log exist!' >&2;
        fi
        cp -a "$finaldir/$srcbase".{trees.bz2,log.bz2,state} "$workdir/"
        bunzip2 "$workdir/$srcbase".{trees.bz2,log.bz2}
    else
        cp -a "$finaldir/$srcbase".{trees,log,state} "$workdir/"
    fi
fi

echo 'Starting beast' >&2

cleanup() {
    mv -v -t "$finaldir/" "$tmproot".* 1>&2
}

trap cleanup INT ERR

if [ $printstdout == 0 ]; then
    $beast -working -prefix "$workdir/" -statefile "$tmproot.state" ${beast_args[@]} > "${tmproot}.stdout"
    # stdout should not be needed since Condor writes anyway in a remote scratch directory
else
    $beast -working -prefix "$workdir/" -statefile "$tmproot.state" ${beast_args[@]}
fi

# We don't need the complete xml file outputted again.
sed -i '/^#/d' "$tmproot.log"

loganalyser -quiet -b 0 "$tmproot.log" > "$finaldir/$srcbase.summary"
if [ -n "$treeannot_args" ]; then
    treeannotator $treeannot_args
fi

#ionice -c2 -n5
bzip2 "$tmproot"{.log,.trees}

mv -t "$finaldir/" "$tmproot".*


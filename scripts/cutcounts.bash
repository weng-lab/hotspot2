#!/bin/bash
# This is a modified version of the cutcounts script used in stampipes
set -x -e -o pipefail

if [[ $# -lt 2 ]] ; then
  echo "Usage: $0 in.bam mappabilityFile.bed" >&2
  exit 2
fi

bam=$1
mappabilityFile=$2

awk_cut="$(dirname $0)/cutfragments.awk"

name=$(basename $bam .bam)
outputdir="$(dirname $bam)"

CUTS_BED=$outputdir/$name.cuts.sorted.bed.starch
CUTCOUNTS=$outputdir/$name.cutcounts.sorted.bed.starch
FRAGMENTS=$outputdir/$name.fragments.sorted.bed.starch

if [[ -z "$TMPDIR" ]] ;then
  TMPDIR=$(mktemp -d)
fi

# temp files
COUNTFILETMP="$TMPDIR/base-count.$name.bed"
ALLBASE="$TMPDIR/all-perBase.$name.temp"
BEDTMP="$TMPDIR/$name.count.uniques.bed"
CUTSTMP="$TMPDIR/cuts.bed"
FRAGMENTSTMP="$TMPDIR/fragments.bed"

# Create cut counts and fragments if they don't exist
if [[  ! -s "$CUTS_BED" || ! -s "$FRAGMENTS" ]]; then

  time bam2bed --do-not-sort < "$bam" \
    | awk -v "cutfile=$CUTSTMP" -v "fragmentfile=$FRAGMENTSTMP" -f "$awk_cut"
  sort-bed --max-mem 1G "$FRAGMENTSTMP" | starch - > "$FRAGMENTS"
  sort-bed --max-mem 1G "$CUTSTMP"      | starch - > "$CUTS_BED"

fi

if [[ ! -s "$CUTCOUNTS" ]]; then

  time unstarch "$CUTS_BED" \
    | cut -f1-3 \
    | bedops -m - \
    | awk '{ for(i = $2; i < $3; i += 1) { print $1"\t"i"\t"i + 1 }}' \
    > "$ALLBASE"

  time unstarch "$CUTS_BED" | \
    bedmap --echo --count --delim "\t" "$ALLBASE" - \
    | awk '{print $1"\t"$2"\t"$3"\tid-"NR"\t"$4}' \
    > "$COUNTFILETMP"

  time starch "$COUNTFILETMP" > "$CUTCOUNTS"

fi

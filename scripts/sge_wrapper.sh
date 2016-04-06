#!/bin/bash

usage(){
  "$hotspot_script" -h
  exit
}

# Alter this to whatever works for you
submit(){
  name=$1
  shift
  holds=$(sed 's/\s\+/,/g' <<< "$@")
  qsub -cwd -V -q all.q -hold_jid ".fake,$holds" -N "$name" -S /bin/bash
}

starch_merge(){
  ext=$1
  shift
  submit "$@" <<__EOF__
  starchcat $outdir.*/$base.$ext > $outdir/$base.$ext
__EOF__

}

density_bigwig(){
  in=$outdir/$base.density.starch
  out=$outdir/$base.density.bw
  submit "$@" <<__EOF__
  TMPFRAGS="\$(mktemp -t fragsXXXXX)"
  unstarch "$in" | cut -f1,2,3,5 > "\$TMPFRAGS"
  bedGraphToBigWig \
    "\$TMPFRAGS" \
    <(cut -f1,3 "$CHROM_SIZES") \
    "$out"
  rm -f "\$TMPFRAGS"
__EOF__
}

SPOT_score(){
  submit "$@" <<'__EOF__'
  num_cleaves=$(samtools idxstats "$bam" | awk '{s+=$3} END{print s}')
  cleaves_in_hotspots=$(bedops --ec -e -1 "$outdir/$base.cutcounts.starch" "$outdir/$base.hotspots.fdr$FDR.starch" | awk '{s=0} {s+=$5} END {print s}')
  echo "scale=4; $cleaves_in_hotspots / $num_cleaves" \
    | bc \
    > "$outdir/$base.SPOT.txt"
__EOF__
}

hotspot_script="$(dirname "$0")/hotspot2.sh"

# Treat -h, -c, and -f specially
# All arguments are passed to hotspot2.sh
otherargs=()
FDR=0.05
while getopts ':hc:e:f:F:m:n:p:s:w:O' opt ; do
  otherargs+=(-$opt $OPTARG)
  case "$opt" in
    h)
      usage
      ;;
    c)
      CHROM_SIZES=$OPTARG
      ;;
    f)
      FDR=$OPTARG
      ;;
  esac
done
shift $((OPTIND-1))

bam=$1
outdir=$2

if [[ -z "$outdir" ]] ; then
  usage
fi

if [[ ! -s "$bam.bai" ]] ; then
  samtools index "$bam"
fi

chroms=($(samtools idxstats "$bam" | cut -f1 | sed '/^\*$/d'))

base="$(basename "$bam" .bam)"

joblist=()

mkdir -p "$outdir"

for c in "${chroms[@]}" ; do
  jobname=".htspt$c.$base.$$"
  joblist+=($jobname)
  submit "$jobname" <<__EOF__
  samtools view -h -1 "$bam" $c > "\$TMPDIR/$base.bam"
  "$hotspot_script" ${otherargs[@]} "\$TMPDIR/$base.bam" "$outdir.$c"
  if [[ ! -s $outdir.$c/$base.density.starch ]] ; then
    rm -f $outdir.$c/$base.density.starch  # Hack to avoid 0-size chrM
  fi
__EOF__
done


#Command       Input                     Job Name            Dependencies
starch_merge   density.starch            ".htsptdn.$base.$$" "${joblist[@]}"
starch_merge   cutcounts.starch          ".htsptcc.$base.$$" "${joblist[@]}"
starch_merge   allcalls.starch           ".htsptac.$base.$$" "${joblist[@]}"
starch_merge   fragments.sorted.starch   ".htsptfg.$base.$$" "${joblist[@]}"
starch_merge   "hotspots.fdr$FDR.starch" ".htspths.$base.$$" "${joblist[@]}"
starch_merge   peaks.starch              ".htsptpk.$base.$$" "${joblist[@]}"
starch_merge   peaks.narrowpeaks.starch  ".htsptpk.$base.$$" "${joblist[@]}"

density_bigwig                           ".htsptdb.$base.$$" ".htsptdn.$base.$$"
SPOT_score                               ".htsptss.$base.$$" ".htsptcc.$base.$$" ".htspths.$base.$$"

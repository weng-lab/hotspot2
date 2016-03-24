#!/bin/bash

# proof-of-concept
# Usage: $0 in.bam outdir [regular options to hotspot2.sh]

# Alter this to whatever works for you
submit(){
  name=$1
  shift
  holds=$(sed 's/\s\+/,/g' <<< $@)
  qsub -cwd -V -q all.q -hold_jid ".fake,$holds" -N "$name" -S /bin/bash
}

starch_merge(){
  ext=$1
  shift
  submit "$@" <<__EOF__
  starchcat $outdir.*/$base.$ext > $outdir/$base.$ext
__EOF__

}

bam=$1
outdir=$2
shift
shift
otherargs=$@

hotspot_script=$(dirname $0)/hotspot2.sh

if [[ ! -s "$bam.bai" ]] ; then
  samtools index "$bam"
fi

chroms=($(samtools idxstats "$bam" | cut -f1 | sed '/^\*$/d'))

base="$(basename "$bam" .bam)"

joblist=""

mkdir -p "$outdir"

for c in "${chroms[@]}" ; do
  jobname=".htspt$c.$base.$$"
  joblist="$jobname $joblist"
  submit "$jobname" <<__EOF__
  samtools view -h -1 "$bam" $c > "\$TMPDIR/$base.bam"
  "$hotspot_script" ${otherargs[@]} "\$TMPDIR/$base.bam" "$outdir.$c"
  if [[ ! -s $outdir.$c/$base.density.starch ]] ;
    rm -f $outdir.$c/$base.density.starch  # Hack to avoid 0-size chrM
  fi
__EOF__
done


starch_merge allcalls.starch         ".htsptac.$base.$$" $joblist
starch_merge cutcounts.starch        ".htsptcc.$base.$$" $joblist
starch_merge fragments.sorted.starch ".htsptfg.$base.$$" $joblist
starch_merge hotspots.fdr0.05.starch ".htspths.$base.$$" $joblist
starch_merge peaks.starch            ".htsptpk.$base.$$" $joblist
starch_merge density.starch          ".htsptdn.$base.$$" $joblist

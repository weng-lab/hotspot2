#!/bin/bash
set -x -e -o pipefail

if [[ $# != 6 ]] ; then
  echo "Usage: $0 tmpdir wavelets-binary <tags.starch> <chrom-sizes.starch> <density-out.starch> <peaks-out.starch>" >&2
  exit 2
fi

tmpdir=$1
wavelets=$2
tags=$3
chrfile=$4
density=$5
pk=$6

if [ ! -d $tmpdir ] then
  mkdir -p $tmpdir
fi

## density params
bins=150
step=20
halfbin=$((bins/2))
rangepad=$((bins/2-step/2))

## wavelet peakfinding params
wavletlvl=3
filter_type=Haar
boundary_type=reflected

echo "density and peak-finding..."
pkouts=""
densouts=""
for chr in $(bedops --ec -u $chrfile | cut -f1)
do
  ## Tag density, 150bp window, sliding every 20bp, used for peak-finding.
  bedops --ec -u --chrom $chr $chrfile \
    | awk -v b=$bins -v s=$step \
       'BEGIN {OFS="\t"; hs=s/2; hb=b/2} ; { \
         for ( start = $2+hb-hs; start < $3-hb-hs; start+=s) { \
           print $1, start, start+s, "."; \
         } \
       }' \
    | bedmap --faster --chrom $chr --range $rangepad --delim "\t" --echo --count - $tags \
    | starch - \
   > $tmpdir/.dens.$chr.starch

  densouts="$densouts $tmpdir/.dens.$chr.starch"

  unstarch $tmpdir/.dens.$chr.starch \
    | cut -f5 \
    | wavelets --level $waveletlvl --to-stdout --boundary $boundary_type --filter $filter_type - \
    > $tmpdir/.waves

  unstarch $tmpdir/.dens.$chr.starch \
    | paste - $tmpdir/.waves \
    | cut -f1-4,6 \
    | awk 'BEGIN{incr=0} ; { \
            if ( NR > 0 ) { \
              if ( incr == 1 ) { \
                if ( $5-lastv < 0 ) { \
                  print lastl; incr=0; \
                } \
              } else { \
                if ( $5-lastv > 0 ) incr=1; \
              } \
            } \
            lastv=$5; lastl=$0; \
          }' \
    > $tmpdir/.wave.$chr

  pkouts="$pkouts $tmpdir/.wave.$chr"

  rm -f $tmpdir/.waves
done

cat $pkouts \
  | awk -v h=$halfbin '{m=($2+$3)/2; left=m-h; if(left < 0) left=0; print $1"\t"left"\t"m+h"\t"$4"\t"$5}' - \
  | starch - \
  > $pk

starchcat $densouts \
  > $density

#need to add pph

exit 0

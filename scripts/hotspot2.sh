#! /bin/bash

# IMPORTANT NOTE:
# REQUIRES BADSPOTS, UNMAPPABLE REGIONS, ETC.
# TO BE MERGED INTO A SINGLE FILE
# WHOSE PATH IS THEN ASSIGNED TO THE VARIABLE
# "EXCLUDE_THESE_REGIONS".
EXCLUDE_THESE_REGIONS=""

# ALSO, THE PATH TO A FILE OF CHROMOSOME NAMES AND SIZES
# CAN BE ASSIGNED TO THE VARIABLE
# "CHROM_SIZES".
# IMPORTANT NOTE:
# IF THIS IS DONE, AND chrM DOES NOT APPEAR IN "CHROM_SIZES",
# THEN DATA FOR chrM WILL BE INTERPRETED AS ERRORS.
CHROM_SIZES=""

COUNTING_EXE=tallyCountsInSmallWindows
HOTSPOT_EXE=hotspot2

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  "$0" infileCuts.starch outfileFromWhichHotspotsAreCalled.starch [localHalfWindowSize] [backgroundHalfWindowSize] [overlapping] [pvalDistnSize]"
    echo -e "\twhere an optional 3rd argument, the local neighborhood on both sides of each bp, is in bp (default = 100 = 201bp window),"
    echo -e "\tan optional 4th argument, the background region on both sides of each local site, is in bp (default = 25000 = 50001bp window),"
    echo -e "\tan optional 5th argument, the word \"overlapping\" or (without hyphen) \"nonoverlapping\", specifies 1bp-sliding overlapping windows or not,"
    echo -e "\tand an optional 5th argument specifies the number of P-values to use for FDR estimates (default = 1000000)."
    exit
fi

INFILE=$1
OUTFILE=$2

if [ "$3" != "" ]; then
    SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=$3
else
    SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE=100 # i.e., 201bp regions
fi

if [ "$4" != "" ]; then
    BACKGROUND_WINDOW_SIZE=`echo $4 | awk '{print 2*$1 + 1}'`
else
    BACKGROUND_WINDOW_SIZE=50001 # i.e., +/-25kb around each position
fi

if [ "$5" != "" ]; then
    OVERLAPPING_OR_NOT=$5
else
    OVERLAPPING_OR_NOT="overlapping"
fi

if [ "$6" != "" ]; then
    PVAL_DISTN_SIZE=$4
else
    PVAL_DISTN_SIZE=1000000
fi


unstarch $INFILE \
    | $COUNTING_EXE $SITE_NEIGHBORHOOD_HALF_WINDOW_SIZE $OVERLAPPING_OR_NOT "reportEachUnit" $CHROM_SIZES \
    | bedops -n - 1 - $EXCLUDE_THESE_REGIONS \
    | $HOTSPOT_EXE $BACKGROUND_WINDOW_SIZE $PVAL_DISTN_SIZE \
    | starch - \
    > $OUTFILE

# Sample code to create, e.g., FDR 5% hotspots.
#
# FDRthresh=0.05
# HOTSPOT_OUTFILE=hotspots_fdr${FDRthresh}.bed5
#
# # P-values of 0 will exist, and we don't want to do log(0).
# # Roughly 1e-308 is the smallest nonzero P usually seen,
# # so we can cap everything at that, or use a different tiny value.
# # The constant c below converts from natural logarithm to log10.
#
# unstarch $OUTFILE \
#    | awk -v threshold=$FDRthresh '{if($6 <= threshold){print}}' \
#    | bedops -m - \
#    | bedmap --delim "\t" --echo --min - $OUTFILE \
#    | awk 'BEGIN{OFS="\t";c=-0.4342944819}{if($4>1e-308){print $1, $2, $3, "id-"NR, c*log($4)}else{print $1,$2,$3,"id-"NR,"308"}' \
#    | starch - \
#    > $HOTSPOT_OUTFILE

exit

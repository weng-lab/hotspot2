hotspot2 is a program developed at the Altius Institute for Biomedical Sciences in Seattle, U.S.A.
for identifying genomic regions with statistically significant "hotspots," or enrichments,
of cleavage activity in DNase-seq experiments.  It is designed to run in a UNIX environment,
and to work with alignment files in BAM format.

hotspot2 requires each of the following programs to be installed and accessible from the user's path
before it can be used:
* g++ (or an equivalent C++ compiler; replace "g++" with its name in Makefile in this case)
* bedops
* samtools
* bedGraphToBigWig
* modwt

hotspot2 tallies cleavages within a small region ("window") around each site.  It slides the window
across the genome, and statistically evaluates cleavage tallies within their local context, i.e.,
within a larger "background window" that itself slides across the genome.  To enable hotspot2
to distinguish between "no data" (e.g., an unmappable region or the end of a chromosome) and
"no cleavages" (a valid region in which no cleavages were observed in the DNase-seq experiment),
a file of chromosome sizes is required.  Ideally, a file containing all the mappable regions
in the genome (e.g., regions uniquely mappable by 36mers) should also be supplied by the user.
(If a "blacklist" of, e.g., mappable but problematic microsatellite regions is available
for the genome, it should be subtracted from the file of mappable regions before being supplied
to hotspot2.)

Before hotspots can be identified, the set of viable positions that can serve as centers of
sliding windows must be determined.  The script extractCenterSites.sh in the scripts subdirectory
must be run to determine these positions.  This script requires a file of chromosome sizes
and the name of the file to which the positions should be written.  A file of mappable regions
can optionally be provided as well, and this is highly recommended.  The default radius around
each position is 100 bp, yielding a sliding window of width 201 bp; a different value can be
supplied if desired.  To see the usage information for this script, type

$ scripts/extractCenterSites.sh -h

Note:  extractCenterSites.sh only needs to be run once per genome.  (If analyses with different
window sizes are desired, extractCenterSites.sh needs to be run once per window size per genome.)

Hotspots are called for an input alignment file in BAM format via a set of scripts that are
each executed by the hotspot2.sh script.  This script also executes a program that needs to
be compiled or "made" on the computer where hotspot2.sh will run.  To make the hotspot2 program,
simply type the command "make" from within the hotspot2 directory.  (This will create the program
hotspot2 in the subdirectory src.)

Once the hotspot2 program has been compiled and the center sites file has been created
by running extractCenterSites.sh, hotspot2.sh will be ready to run.  To run hotspot2.sh
with default values for its various parameters, type

$ scripts/hotspot2.sh yourData.bam yourOutputDirectory

(hotspot2.sh will create the output directory "yourOutputDirectory" if it does not already exist.)

To see the various parameters for hotspot2.sh and their default settings, type

$ scripts/hotspot2.sh -h



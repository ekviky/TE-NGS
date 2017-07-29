# TE-NGS: TE calls from targeted TE-enriched NGS data
# subdirectories for distribution

1. scripts 
TE calling algorithm consists of two principle steps: 
(i) clustering on genomic coordinates 
(ii) annotation via comparison to public and local TE databases

# eg to run script 
# tested on rescomp and florence/dougal/dylan 
# run no arguments to see usage
# switch to run in test mode (NA12878) and output files for comparison to expectation; sample mode, eg unknown sample
$ cd scripts
$ module load R/3.1.2 
$ R --no-restore --no-save --no-readline --quiet < R_get_TE_calls_v0.1.r --args ../test/ ../test/seq_opt_table.txt 201 test 

2. annotations 
# contains local copies all input TE annotations, reference annotations, primer sequences, etc 
gw_TE_targets_window3flank600bp.interval
hg19_RM_allLINESINE_window3flank600bp.interval
hg19_access_custom2.bed
hg19_centroTelo_sorted.txt
phase3_20130502_NA12878_present_window3flank600bp.interval
polyTEdb_window3flank600bp.interval
P3_primers.txt

3. clusters ie working dir for results

4. test 
# contains files necessary for implementation and unit testing
seq_opt_table.txt # input file describing experimental design when running on multiple samples
input bams NA12878 # to do: eg chr1 or subset
cluster output files for NA12878 to compare output
annotated cluster files for NA12878 for comparison/unit testing

5. examples
# contains example output files
NA12878_TE-NGS_calls_annotated_v0.1.bed # all TE calls made for NA12878 and results described in manuscript

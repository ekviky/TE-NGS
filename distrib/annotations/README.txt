This directory contains flat files used by the TE detection script.

The following are interval/bed format files of TE insertions annotated in the reference genome (build GRCh37/hg19). Stard/end coordinate positions are extended by a 600bp flanking window downstream of the 3 prime TE insertion site wrt reference strand (ie., end+600 if inserted on forward strand; start-600 if inserted on reverse strand).

1. polyTEdb_window3flank600bp.interval.gz
A locally compiled database of polymorhic germline TE insertion calls obtained from public databases and TE-targeting protocols. TE calls are considered derived from the same insertion if coordinates are within a maximum 100bp distance between calls. Each TE is cross-referenced giving the number of sources, subfamily annotation in each source (where provided), and primary source providing support for the call. 

TE calls were obtained from the following sources:
	MEI, Mobile Element Insertions (Alu and LINE) described in Stewart et al 2011
	Witherspoon, Alu identified by ME-scan as described in Witherspoon et al 2013
	phase3, 1000 Genomes Phase 3 Alu and LINE calls in Sudmant et al 2015
	dbRIP, database of retrotransposon insertion polymorphisms (Alu and LINE) as described in Wang et al 2006
	Ewing, LINE-1 insertions described in Ewing and Kazazian 2011

2. phase3_20130502_NA12878_present_window3flank600bp.interval.gz
1000 Genomes Phase 3 (release 20130502) mobile element insertions observed in NA12878
NB, only used in test mode when run on NA12878 files provided in /test.

3. hg19_RM_allLINESINE_window3flank600bp.interval.gz
RepeatMasker SINE/LINE annotations for hg19

4. hg19_TE_targets_window3flank600bp.interval.gz
Targetable loci present in hg19 defined as loci with matches to both TE-target and TE-nested primers, allowing at most 1 mismatch excluding the ultimate 3 prime primer nucleotide, positioned a maximum distance of 200bp apart, and in proper order and orientation. 

The following are hg19 coordinates used to define regions accessible by short-read alignment
 - hg19_centroTelo_sorted.txt.gz, centromeres and telomeres obtained from UCSC Table Browser
 - hg19_NA12878_access_custom2.bed.gz, regions excluding centromeres/telomeres, reference gaps, low complexity regions obtained from UCSC Table Browser annotations and NA12878-specific CNV obtained from DGV. NB, only used in test mode when run on NA12878 files provided in /test.

A fasta file containing nested primer sequences used to generate TE-NGS libraries (Table S2 in manuscript) is also required as input.
 - primers.fa.gz

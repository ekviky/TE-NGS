This directory contains flat files used by the TE detection script.


The following are interval/bed format files of TE insertions annotated in the reference genome (build gr37/hg19). Stard/end coordinate positions are extended by a 600bp buffer downsream of the 3' end of the TE insertion wrt reference strand (ie., end+600 if inserted on forward strand; start-600 if inserted on reverse strand).

1. polyTEdb_window3flank600bp.interval.gz
A locally compiled database of polymorhic TE insertions obtained from public databases and TE-targeting protocols. Each TE is cross-referenced giving the number of sources, subfamily annotation in each source (where provided), and primary source reference. 
where nonreference annotations were obtained from the following sources:
	MEI, Stewart et al 2011
	Witherspoon, Witherspoon et al 2013
	phase3, Studmant et al 2015
	dbRIP, Wang 
	Ewing, Ewing et al 2011

2. phase3_20130502_NA12878_present_window3flank600bp.interval.gz
1000 Genomes Phase 3 MEI observed in NA12878

3. hg19_RM_allLINESINE_window3flank600bp.interval.gz
RepeatMasker SINE/LINE annotations

4. gw_TE_targets_window3flank600bp.interval.gz


The following are hg19 coordinates used to define regions accessible by short-read alignment
 - hg19_centroTelo_sorted.txt.gz
 - hg19_NA12878_access_custom2 # only used in test mode, when run on NA12878 files provided

Fasta file containing nested primer sequences used to generate TE-NGS libraries (Table SX) 
 - P3_primers.txt.gz

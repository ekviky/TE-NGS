# R script to analyze TE-enriched sequence data
# first step algorithm: identify read clusters by genome coordinates
# second step algorithm: intersection with TE annotations

options(echo=FALSE);

## arguments
args<-commandArgs(TRUE)
if(length(args)<4){
		stop("usage: R_get_TE_calls_v0.1.r --args [bam dir] [in:seq_opt_table] [in: minimum (bp) cluster inter-distance] [mode] 
Identify and annotate targeted TE sequencing clusters 
where:
	- bam dir is location bam alignment files
	- seq_opt_table corresponds to input file containing info: experiment, TE, primers used, and unique sequencing index	
	- cluster inter-distance, minimum distance required between unique clusters in bp [201]
	- run mode, one of test or sample; test mode using NA12878 for unit testing, sample mode to analyse sample bams
	- output is txt file containing summary level statistics per sample/unique index'\n",call.=FALSE)
}

print(paste("# Started",Sys.time(),sep=" "))

print(args)
bamDir<-file.path(args[1])
miseq_opt_table<-args[2]
T<-args[3]
T<-as.numeric(args[3])
mode<-args[4]
print(paste("# mode=",mode,sep=" "))

## environment
#library(Rsamtools)
#library(seqinr)
#library(ShortRead)
#library(data.table)
#library(GenomicRanges)

suppressMessages(require("Rsamtools",quietly=TRUE,warn.conflicts = FALSE));
suppressMessages(require("seqinr",quietly=TRUE,warn.conflicts = FALSE));
suppressMessages(require("ShortRead",quietly=TRUE,warn.conflicts = FALSE));
suppressMessages(require("data.table",quietly=TRUE,warn.conflicts = FALSE));
suppressMessages(require("GenomicRanges",quietly=TRUE,warn.conflicts = FALSE));

sessionInfo()

## files and locations
# input files
miseq_table<-read.delim(file.path(miseq_opt_table),sep="\t",header=TRUE);
nrow(miseq_table);

# annotations
# to do: loop to nicely read all gz files for github file size limit
annot_dir<-file.path("../annotations/");
primersP3<-read.fasta(file.path(annot_dir,"P3_primers.txt"),as.string=TRUE,strip.desc=TRUE) # object containing primer names and DNAString
targets<-fread(file.path(annot_dir,"gw_TE_targets_window3flank600bp.interval"),sep="\t",header=FALSE)
phase3<-fread(file.path(annot_dir,"phase3_20130502_NA12878_present_window3flank600bp.interval"),sep="\t",header=FALSE)
allTE<-fread(file.path(annot_dir, "hg19_RM_allLINESINE_window3flank600bp.interval") ,sep="\t",header=FALSE)
nonrefall<-fread(file.path(annot_dir, "polyTEdb_window3flank600bp.interval") ,sep="\t",header=FALSE)
centroTelo<-fread(file.path(annot_dir, "hg19_centroTelo_sorted.txt"),sep="\t",header=FALSE)
NA12878_access<-fread(file.path(annot_dir, "hg19_NA12878_access_custom2.bed"),sep="\t",header=FALSE)

# bam files
bamFls<-list.files(bamDir,"bam$",full=TRUE) # list all bam files in dir
names(bamFls)<-basename(bamFls)
length(bamFls)

## variables
# parameters for filtering bam records
flag1<-scanBamFlag(isUnmappedQuery=FALSE, hasUnmappedMate=FALSE)
flag2<-scanBamFlag(isUnmappedQuery=FALSE, hasUnmappedMate=FALSE,isProperPair=TRUE,isSecondMateRead=TRUE)
param1<-ScanBamParam(flag1,what=scanBamWhat()) 
param0<-ScanBamParam(what=scanBamWhat())
param2<-ScanBamParam(flag2,what=scanBamWhat())

# threshold for distance between unique clusters
print(paste("# Threshold cluster-D",T,sep=" "))

# list indices corresponding to each element
element_A<-as.character(miseq_table[which(miseq_table$category == "element_A"),"index"])
element_B<-as.character(miseq_table[which(miseq_table$category == "element_B"),"index"])
element_C<-as.character(miseq_table[which(miseq_table$category == "element_C"),"index"])

# targets and their respective primer sequence names
target_A="L1HS"
target_B="AluYa[58]"
target_C="AluYb[89]"

nestedprimerA="L1HS_nested"
nestedprimerB="AluYa58_nested"
nestedprimerC="AluYb89_nested"

# define sequence patterns thresholds
# TE consensus sequence 3'end exclude primer
TEcons_fwd_A<-"AGTATAAT"
TEcons_fwd_B<-TEcons_fwd_C<-"GCCTGGGCGACAGAGCGAGACTCCGTCTC"
# threshold for read mismatches to TE consensus
TE_T_A=3 
TE_T_B=TE_T_C=10

## functions
# convert data.frame object to genomic ranges
convert_GR<-function(x){
        a.gr=makeGRangesFromDataFrame(x,
        seqnames.field=c("V1"),
        start.field=c("V2"),
        end.field=c("V3"))
        return(a.gr)
}

# intersect genomic ranges objects
get_intersection_clusters<-function(clust,clust.gr,annot,annot.gr){
        dist<-data.frame(distanceToNearest(clust.gr, annot.gr,ignore.strand=TRUE))
        out<-data.frame(clust[dist[which(dist$distance<10),"queryHits"],], annot[dist[which(dist$distance<10),"subjectHits"],],dist[which(dist$distance<10),"distance"])
        return(out)
	rm(out)
}

# count number of reads per cluster
# binary search implemented in data.table using keys for fast subsetting 
# data.table v1.9.4 has known bug issues => fixed dev version 1.9.5, release to CRAN 1.9.6 tba
# 1.9.6 DOES NOT work as of 06/04/2016 => to debug
# setkey(DT[J(x)], id)[J(seq_along(vv)), list(lapply(list(pos),length))]$V1 # should work => setkey and subset by all possible id at same time creates NAs for missing observations
# R 3.1.2 and data.table 1.8.8 work on florence/dylan and rescomp1
count_read2<-function(x,y){
	# where x is vector of read2
	# y is list length == number of clusters, which reads per cluster 
	# create vector each value corresponds to length of list (number of reads) per item (cluster) in y
	vv<-vapply(y,length,0L)
	# create data.table rows == each item in 
	DT<-data.table(V1=unlist(y),id=rep(seq_along(y),vv),pos=sequence(vv))
	# setkey and subset by all possible id at same time creates NAs for missing observations
	# want this for calcuation no. of reads each cluster
	setkey(DT,V1)
	result<-setkey(DT[J(x)],id)
	result2<-(result[J(seq_along(vv)),list(pos)][,length(na.exclude(pos)),by=id]) 
	return(result2$V1);

}

# annotate observed reference RM 
get_ref_annot<-function(clust,clust.gr){
        tmp<-get_intersection_clusters(clust,clust.gr,targets,targets.gr)

        out_targets_annot<-tmp[,c(1:12,20,19)]
        out_targets_annot$ref_DB<-"RM"
        out_targets_annot$annot<-"ref"

        header_ref<-c("chrom","start","end","length","strand","reads","read2_count","read1_count","read1_mapq3","read1_mapq10","experiment","index","Ref_TE_family","Ref_subfamily","Ref_Db","Ref_annotation")
        names(out_targets_annot)<-header_ref
        return(out_targets_annot)
        rm(tmp)
}

# annotate observed NA12878 phase3 
get_phase3_annot<-function(clust,clust.gr){
        tmp<-get_intersection_clusters(clust,clust.gr,phase3,phase3.gr)

        out_phase3_annot<-tmp[,c(1:12,17,18,21)]
        out_phase3_annot$annot<-"phase3_poly"
        header_phase3NA12878<-c("chrom","start","end","length","strand","reads","read2_count","read1_count","read1_mapq3","read1_mapq10","experiment","index","phase3_TE_family","phase3_subfamily","phase3_Db","phase3_annotation")
        names(out_phase3_annot)<-header_phase3NA12878
        return(out_phase3_annot)
        rm(tmp)
}

# annotate observed known polymorphic 
get_polyTEdb_annot<-function(clust,clust.gr){
        tmp<-get_intersection_clusters(clust,clust.gr,nonrefall,nonrefall.gr)

        out_poly_annot<-tmp[,c(1:12,17:19)]
        out_poly_annot$annot<-"nonref_poly"
        header_poly<-c("chrom","start","end","length","strand","reads","read2_count","read1_count","read1_mapq3","read1_mapq10","experiment","index","Poly_TE_family","Poly_subfamily","Poly_Db","Poly_annotation")
        names(out_poly_annot)<-header_poly
        return(out_poly_annot)
        rm(tmp)
}

# annotate novel ie previously unreported
get_novel_annot<-function(clust,clust.gr){
        # lack annotation/fail overlap
        a<-which(mcols(distanceToNearest(clust.gr,targets.gr,ignore.strand=TRUE))$distance>1)
        b<-which(mcols(distanceToNearest(clust.gr,phase3.gr,ignore.strand=TRUE))$distance>1)
        c<-which(mcols(distanceToNearest(clust.gr,allTE.gr,ignore.strand=TRUE))$distance>1)
        d<-which(mcols(distanceToNearest(clust.gr,nonrefall.gr,ignore.strand=TRUE))$distance>1)

        tmp<-((clust)[Reduce(intersect, list(a,b,c,d)),])

        test_out_novel_annot<-tmp
        test_out_novel_annot$annot<-"novel"
        names(test_out_novel_annot)<-c("chrom","start","end","length","strand","reads","read2_count","read1_count","read1_mapq3","read1_mapq10","experiment","index","annotation")
        return(test_out_novel_annot)
        rm(a,b,c,d,tmp)
}

# get clusters and some prelim stats per cluster
get_clusters<-function(bf,p,T,primer,TE,TE_T,run_mode){
	# read bam into object in R
	sub<-readGAlignmentsFromBam(bf,param=p,use.names=TRUE)
	print(paste("# alignments read",Sys.time(),sep=" "))
	print(length(sub))

	# filter reads pre-clustering
	## Filter 1: primer matches
	### match primer 
	print(primer)
	print(primersP3[[match(primer,names(primersP3))]][1])
	Lp<-nchar(primersP3[[match(primer,names(primersP3))]][1])
       	mindex<-vmatchPattern(primersP3[[match(primer,names(primersP3))]][1], sub@elementMetadata$seq) # match primer on fwd strand
	mindex2<-vmatchPattern(reverseComplement(DNAString(primersP3[[match(primer,names(primersP3))]][1])), sub@elementMetadata$seq) # match primer on rev strand

	### match primer in read 2 and keep pairs with read 2 matching primer
	matches<-c(names(sub[which(countIndex(mindex)>0)]),names(sub[which(countIndex(mindex2)>0)]))
	print(length(unique(matches)))
	filtered<-sub[which(names(sub) %in% matches)]
	print(paste("# primer matches found and filtered",Sys.time(),sep=" "))
	print(length(filtered))

	## Filter 2: match TE consensus sequence
	### match TE consensus sequence on read 2 and keep pairs with read 2 matching consensus with max T mismatches
	print(TE);
	print(get(TE))
	print(TE_T)
	print(get(TE_T))
	mindex3<-vmatchPattern(get(TE), filtered@elementMetadata$seq,max.mismatch=get(TE_T)) # match seq on fwd strand
	mindex4<-vmatchPattern(reverseComplement(DNAString(get(TE))), filtered@elementMetadata$seq,max.mismatch=get(TE_T)) # match seq on rev strand

	### match TE consensus seq in read 2 and keep pairs with read 2 matching with max T mismatches
	matches_TEcons<-c(names(filtered[which(countIndex(mindex3)>0 ) ]),names(filtered[which(countIndex(mindex4)>0 )]))

	### restrict match to substring==TE just immediately flanking primer
	med_start_fwd<-median(unlist(startIndex(mindex3)));
	med_start_rev<-median(unlist(startIndex(mindex4)));

	polyA_hit<-c(names(filtered[which(sapply(startIndex(mindex3),min)<55)]), names(filtered[which(sapply(startIndex(mindex4),max)>151-55)])) # for universal?
	primer_and_polyA<-polyA_hit[which(polyA_hit %in% matches_TEcons)]

	filtered2<-filtered[which(names(filtered) %in% primer_and_polyA)]
	print(paste("# TE consensus matches found and filtered",Sys.time(),sep=" "))
	print(length(filtered2))

	## Filter 3: Minimum mapping quality
	print(paste("reads mapq==0", length(which(filtered2@elementMetadata$mapq==0)),sep=" " ))
	#filt_map<-which(filtered2@elementMetadata$mapq!=0) 
	#print(length(filt_map))
	#filtered3<-filtered2[filt_map]
	filtered3<-filtered2[which(filtered2@elementMetadata$mapq!=0)]
	print(length(filtered3))

	# get clusters
	## get genomic ranges object representing bams
	filt_r<-granges(filtered3) # primer and TEseq and filter min MAPQ>0
	print(length(filt_r))

	## cluster on genomic position
	clusters<-reduce(filt_r, min.gapwidth=T,ignore.strand=TRUE,with.revmap=TRUE) # get defined clustered regions as ranges, metadata "revmap" == indices to original ranges for each cluster
	print(length(clusters))
	print(paste("# clusters identified",Sys.time(),sep=" "))

	# convert ranges object to data.frame object
	out<-as.data.frame(clusters)
	print(nrow(out))

	# get some stats per cluster
	## count number reads
	out$reads<-as.vector(unlist(lapply(mcols(clusters)$revmap,length))) # get number of reads per cluster
	print(paste("# reads per cluster calculated",Sys.time(),sep=" "))
        
	## count read2 per cluster
	read2<-which(filtered3@elementMetadata$flag>128)
	out$read2<-count_read2(read2,mcols(clusters)$revmap);
	print(paste("# read 2 per cluster calculated",Sys.time(),sep=" "))
	
	## count read1 per cluster
	read1<-which(filtered3@elementMetadata$flag<128)
	out$read1<-count_read2(read1,mcols(clusters)$revmap);
        print(paste("# read 1 per cluster calculated",Sys.time(),sep=" "))
	
        ## count read1 with min mapq filter
	### min T=3
        filt2<- which(filtered3@elementMetadata$flag<128 & filtered3@elementMetadata$mapq>=3)
	out$read1_minMapq3<-count_read2(filt2,mcols(clusters)$revmap);
        print(paste("# read 1 min MAPQ 3 per cluster calculated",Sys.time(),sep=" "))

	### min T=10
        filt3<- which(filtered3@elementMetadata$flag<128 & filtered3@elementMetadata$mapq>=10)
        out$read1_minMapq10<-count_read2(filt3,mcols(clusters)$revmap);
        print(paste("# read 1 min MAPQ 10 per cluster calculated",Sys.time(),sep=" "))
	print(nrow(out))

	# convert to output in useful format
	bed<-out
	print(nrow(bed))
	bed$start<-bed$start-1 # transform coords to zero-based, half-open i.e. bed
	bed$strand<-sub("\\*","+",bed$strand) # convert funny to + strand
	bed<-bed[,-c(which(names(bed)=="revmap"))] # remove column list indices, NB replace 
	print(nrow(bed))

	# remove clusters on decoy non-canonical chromosomes & chrY
	if(run_mode!="test"){
		bed<-bed[-grep("_|M|Y",bed$seqnames),]
	}
	print(nrow(bed))	
	# filter clusters size/width minimum 100bp
	bed<-bed[-which(bed$width<100),]
	return(bed)
 }

# annotate clusters
get_annotations<-function(clust,clust.gr,target,run_mode){
	# get various annotations
        ## reference target calls
        test_out_targets_annot<-get_ref_annot(clust,clust.gr)
        print(paste("# reference targets annotated",Sys.time(),sep=" "))
        print(nrow(test_out_targets_annot))

        ## polyTEdb calls
        test_out_poly_annot<-get_polyTEdb_annot(clust,clust.gr)
        print(paste("# polyTEdb annotated",Sys.time(),sep=" "))
	print(nrow(test_out_poly_annot))

        ## novel calls
        novel<-get_novel_annot(clust,clust.gr)
        test_out_novel_annot<-novel
        print(paste("# novel calls annotated",Sys.time(),sep=" "))
	print(nrow(test_out_novel_annot))


        ## 1KG phase 3 MEI calls in NA12878
        if(run_mode=="test"){
                test_out_phase3_annot<-get_phase3_annot(clust,clust.gr)
	        print(paste("# 1KG phase 3 NA12878 annotated",Sys.time(),sep=" "))
                print(nrow(test_out_phase3_annot))
                annots<-c("test_out_targets_annot","test_out_phase3_annot","test_out_poly_annot", "test_out_novel_annot")
                cols<-c("chrom","start","end","Ref_TE_family","Ref_subfamily","Ref_Db","Ref_annotation","Poly_TE_family","Poly_subfamily","Poly_Db","Poly_annotation","phase3_TE_family","phase3_subfamily","phase3_Db","phase3_annotation","annotation")

        }else{
                annots<-c("test_out_targets_annot","test_out_poly_annot", "test_out_novel_annot")
                cols<-c("chrom","start","end","Ref_TE_family","Ref_subfamily","Ref_Db","Ref_annotation","Poly_TE_family","Poly_subfamily","Poly_Db","Poly_annotation","annotation")

        }

        # merge annotations 
	print(annots)
	print(exists(annots[1]))
        test_out_all<-Reduce(function(...) merge(...,all=TRUE,by=c("chrom","start","end")), mget(annots))

        # get relevant columns, merge back on filt clusters to get full details (including read info) and write out calls where min mapqT=2
        test_out_all_nice<-merge(clusters_filt,test_out_all[,cols],by=c("chrom","start","end"))
        calls<-test_out_all_nice[which(test_out_all_nice$read1_mapq3>=2),]

        # add one summary annotation column
        if(run_mode=="test"){
                p<-which(calls$phase3_annotation!="NA")
                calls[p,"annotation"]<-calls$phase3_annotation[p]
        }
        poly<-which(calls$Poly_annotation!="NA" & is.na(calls$annotation)==TRUE)
        calls[poly,"annotation"]<-calls$Poly_annotation[poly]
        ref<-which(calls$Ref_annotation!="NA" & is.na(calls$annotation)==TRUE)
        calls[ref,"annotation"]<-calls$Ref_annotation[ref]

        print(paste("# annotations completed",Sys.time(),sep=" "))
	# return annotations data.frame object
	return(calls)
	
}


# make genomic ranges objects 
targets.gr<-convert_GR(targets)
phase3.gr<-convert_GR(phase3)
allTE.gr<-convert_GR(allTE)
nonrefall.gr<-convert_GR(nonrefall)
centroTelo.gr<-convert_GR(centroTelo)
NA12878_access.gr<-convert_GR(NA12878_access)

# execute
# run on all samples
for (i in c(miseq_table$index)){
	print(paste("#",i,sep=" "));

	# infile
	bf<-grep(paste("_",i,sep=""),bamFls,value=TRUE)
        print(bf)

	# get experiment-level details including nested primer
        exp_i<-miseq_table[match(i,miseq_table$index),"experiment"];
        print(exp_i);

        primer_i<-miseq_table[match(i,miseq_table$index),"P3"];
        print(primer_i);

	category_i<-miseq_table[match(i,miseq_table$index),"category"];
	TEcons_i<-gsub("element_","TEcons_fwd_",category_i)
	print(TEcons_i);

	TE_cons_Ti<-gsub("element_","TE_T_",category_i)
	print(TE_cons_Ti);

	te_i<-gsub("element_","target_",category_i)
	print(te_i)

	# get read clusters
	clusters_i<-get_clusters(bf,param0,T,primer_i,TEcons_i,TE_cons_Ti,mode)
	print(nrow(clusters_i))

	# add info such as sample ID and index
	clusters_i$experiment<-rep(exp_i,nrow(clusters_i))
	clusters_i$index<-rep(i,nrow(clusters_i))

	# update clusters header 
	names(clusters_i)[1]<-paste("#",names(clusters_i)[1],sep="")
	print(names(clusters_i))
	print(nrow(clusters_i))

	# if mode == test, write all clusters and add unit test compare to expected out
        if(mode=="test"){
		print("# test mode write all clusters")
		outfile<-paste(sub(".*/","",bf),"clusters_v0.1.bed",sep="_")
		print(paste(sub(".*/","",bf),"clusters_v0.1.bed",sep="_"))
		write.table(clusters_i,file=outfile,quote=FALSE,sep="\t",row.names=FALSE)
		print(paste("# outfile written",Sys.time(),sep=" "))
	}
	
	# get GR object of clusters from bed file returned
        clusters_i.gr<-GRanges(seqnames= clusters_i[,1],
               ranges=IRanges(start= clusters_i[,2],
                              end= clusters_i[,3]),
               strand= clusters_i[,5])
        print(clusters_i.gr)

	# filter clusters by distance centro/telomeres
        filt1<-which(mcols(distanceToNearest(clusters_i.gr, centroTelo.gr,ignore.strand=TRUE))$distance>1000000)
        print(length(filt1))

        # filter NA12878 by NA12878_accessible genome
        if(mode=="test"){
                filt2<-which(mcols(distanceToNearest(clusters_i.gr, NA12878_access.gr,ignore.strand=TRUE))$distance<10)
                print(length(filt2))
                clusters_filt<-(clusters_i[Reduce(intersect, list(filt1,filt2)),])
        }else {
                clusters_filt<-(clusters_i[filt1,])
        }
        print(nrow(clusters_filt))
       	clusters_filt_i.gr<-GRanges(seqnames= clusters_filt[,1],
               ranges=IRanges(start= clusters_filt[,2],
                              end= clusters_filt[,3]),
               strand= clusters_filt[,5])

        print(clusters_filt_i.gr)

	header<-c("chrom","start","end","length","strand","reads","read2_count","read1_count","read1_mapq3","read1_mapq10","experiment","index")
        names(clusters_filt)<-header

        # annotate clusters
	TE_calls<-get_annotations(clusters_filt,clusters_filt_i.gr,te_i,mode)
	print(nrow(TE_calls))
	outfile2<-paste(sub(".*/","",bf),"TEcalls_annotated_v0.1.bed",sep="_")
	write.table(TE_calls,file=outfile2,quote=FALSE,sep="\t",row.names=FALSE)
	print(paste("# outfile written",Sys.time(),sep=" "))
}
	
print(paste("# Finished",Sys.time(),sep=" "))

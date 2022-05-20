#locus plots
#libraries
library(GenomicFeatures)
library(rtracklayer)
library(data.table)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)

#load data
#annotations
#anno_original_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.gtf")
anno_gtf <- makeTxDbFromGFF("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_sub.gtf")
anno <- readRDS("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.annotated.sorted_repeat_anno.rds")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/210726_shortRead_processing_deNovo_custom4/gffCompare.mergedTranscripts.gtf.tmap"))
anno_original <-  readRDS("/omics/groups/OE0219/internal/tinat/raw_data_repo/references/gencode.v29lift37.annotation.repeat_anno.rds")
repeats <- fread("/omics/groups/OE0219/internal/repguide/data/repeats/hg19_repeats.txt")

#peptides 
output.dir <- "/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression"
peptides_new <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/integration/peptidomics/comparison_gene_expression", "peptides_list_new.rds"))
peptide_coord <- readRDS(file.path(output.dir,"2of3DACSB_unique", "peptide_coordinates.rds") )
peptide_coord_unlist <- unlist(GRangesList(peptide_coord))
peptide_coord_sub <- peptide_coord_unlist[grep("_",as.character(seqnames(peptide_coord_unlist)), invert=TRUE),]
peptide_coord_sub <- keepStandardChromosomes(peptide_coord_sub, pruning.mode="coarse")

#rnaseq -> merged tracks
DMSO_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DMSO.normalized.bigWig")
SB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/SB939.normalized.bigWig")
DAC_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DAC.normalized.bigWig")
DACSB_rnaseq <- import.bw("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/bigwig/DACSB939.normalized.bigWig")


#cage
DMSO_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DMSO.bigWig")
SB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/SB939.bigWig")
DAC_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DAC.bigWig")
DACSB_cage <- import.bw("/omics/groups/OE0219/internal/tinat/raw_data_repo/CAGE/DACSB.bigWig")

#nanopore
DACSB_nanopore <- import.bw("/omics/groups/OE0219/internal/tinat/210709_nanopore_ashish_stringtie_assembly/pipeline-nanopore-ref-isoforms/bigwig/reads_aln_sorted.bw")

#methylation
DMSO_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150388_H2_DMSO_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
SB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150389_H2_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DAC_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150386_H2_DAC_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")
DACSB_meth <- fread("/omics/groups/OE0219/internal/tinat/raw_data_repo/ChIP/online/GSM2150387_H2_DAC_plus_SB939_2lanes_merged.CG.ALL.call.gz.BSmooth.csv.gz")


#prepare data
#create granges from repeats file
repeats <- makeGRangesFromDataFrame(repeats, seqnames.field="genoName",
    start.field="genoStart", end.field="genoEnd", keep.extra.columns= TRUE)
repeats$repName <- as.character(repeats$repName)

#subset annos for symbol parsing
# anno_original_transcript_symbol <- data.frame(row.names=anno_original[anno_original$type =="transcript",]$transcript_id,
# SYMBOL=  anno_original[anno_original$type =="transcript",]$transcript_name)
anno_transcript_symbol <- data.frame(row.names=anno[anno$type =="transcript",]$transcript_id,
    SYMBOL=  anno[anno$type =="transcript",]$gene_name)
#get class code information
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_transcript_class <- data.frame(row.names=anno_classi$qry_id,
    Class_code=  anno_classi$class_code_simple)
#replace with easy name for gviz function
anno_transcript_class$Class_code <- gsub("non-chimeric (novel)","NonChimeric",anno_transcript_class$Class_code,fixed=TRUE)
anno_transcript_class$Class_code <- gsub("chimeric (novel)","Chimeric",anno_transcript_class$Class_code,fixed=TRUE)
anno_transcript_class$Class_code <- gsub("known","Known",anno_transcript_class$Class_code,fixed=TRUE)

#make granges of methylation
DMSO_meth <- makeGRangesFromDataFrame(DMSO_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DMSO_meth)<- DMSO_meth$Smoothed_Methylation_Level_H2_DMSO
score(DMSO_meth)<- DMSO_meth$Smoothed_Methylation_Level_H2_DMSO
SB_meth <- makeGRangesFromDataFrame(SB_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(SB_meth)<- SB_meth$Smoothed_Methylation_Level_H2_SB939
score(SB_meth)<- SB_meth$Smoothed_Methylation_Level_H2_SB939
DAC_meth <- makeGRangesFromDataFrame(DAC_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DAC_meth)<- DAC_meth$Smoothed_Methylation_Level_H2_DAC
score(DAC_meth)<- DAC_meth$Smoothed_Methylation_Level_H2_DAC
DACSB_meth <- makeGRangesFromDataFrame(DACSB_meth, seqnames.field="Chromosome", start.field="Start", end.field="Start",keep.extra.columns=TRUE)
mcols(DACSB_meth)<- DACSB_meth$Smoothed_Methylation_Level_H2_DAC_plus_SB939 
score(DACSB_meth)<- DACSB_meth$Smoothed_Methylation_Level_H2_DAC_plus_SB939 

#parameters
#region
ext <- 0 
fontSize <- 10
#get color code
library(RColorBrewer)
col <- brewer.pal(12, "Paired")
treat_col <- col[c(4,2,8,6)]
names(treat_col)<-c("DMSO", "DAC", "SB939",  "DACandSB939")
class_col <- c("gray", col[c(9,10)])
names(class_col)<- c("known",  "chimeric (novel)", "non-chimeric (novel)")
ere_col <- col[c(1,12,5,7,3)]
names(ere_col)<-  c("LINE", "LTR", "no ERE", "other", "SINE")
exon_col <- c("darkgray", "whitesmoke")
names(exon_col)<- c("multi-exonic", "mono-exonic")



#with nanopore data or alignment data
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub")
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withNano")
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withAlign")
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withAlign_Peptide")
dir.create("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_Peptide_only")

#select region
temp <- "^GAPDH$"
roi <- anno_original[grep(temp, anno_original$gene_name),]
 roi_new <-GRanges(
     seqnames = unique(seqnames(roi)),
     ranges = IRanges(start = min(start(roi)),
                    end =  max(end(roi))))
roi_new$transcript_id <- unique(roi$transcript_id)[2]
anno_transcript <- anno[anno$type=="transcript",]

#subset peptide list based on our ORFs == Universe
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]

#get transcripts that give rise to uniquely presented peptides
peptides_new_ORF_noDMSO <-  peptides_new_ORF[peptides_new_ORF$DAC_SB>60 & peptides_new_ORF$DMSO ==0,]

#peptides_new_ORF_noDMSO_split <- split(peptides_new_ORF_noDMSO, peptides_new_ORF_noDMSO$DAC_SB)
#lapply(peptides_new_ORF_noDMSO_split, length)
#names(peptides_new_ORF_noDMSO_split)<- c("1/3 DAC + SB939 unique", "2/3 DAC + SB939 unique", "3/3 DAC + SB939 unique")
#roi <- anno_transcript[anno_transcript$transcript_id %in%  peptides_new_ORF_noDMSO_split$"3/3 DAC + SB939 unique"$transcript_id,]
roi <- anno_transcript[anno_transcript$transcript_id %in%  peptides_new_ORF_noDMSO$transcript_id,]
roi_new <- roi[grep("gl", as.character(seqnames(roi)), invert=TRUE),]

#specific selection
roi_new <- anno_transcript
roi_new <- roi_new[roi_new$transcript_id %in% c("MSTRG.16894.3","MSTRG.28868.1", "MSTRG.29191.1"),]

#peptide_coord
roi_new <- peptide_coord_sub[peptide_coord_sub$transcript_id %in% c("MSTRG.16894.3","MSTRG.28868.1", "MSTRG.29191.1"),]

#roi based on coordinates 
roi_new <- GRanges(
    seqnames = c("chr8", "chr8", "chr16", "chr10", "chr20","chr22","chr22", "chr9", "chr9","chr12", "chr1", "chr8", "chr19" ),
    ranges = IRanges(start = c(101267295, 101267295,20942433,114128500,55922500,30991718,30991718,91973702, 91973702, 7144000,119713854, 77347011,8443894),
                   end =c(101350432,101330400, 21172762,114190134,55956258,31039779,31027974,92112571,92116000,7179000,119727577,77361001,8454545)))
roi_new$transcript_id <- c("RNF19A_1","RNF19A_2", "DNAH3", "ACSL5" , "RAE1","TCN2", "TCN2","SEMA4D_1","SEMA4D_2" ,"C1S"  ,"MSTRG.1804","MSTRG.32198","MSTRG.16351"   )
roi_new$transcript_id <- paste0(roi_new$transcript_id ,"_selected")

#plot regions
for (i in 2:length(roi_new)){
    tryCatch({
    #Define Region
    lim <- c(start(roi_new[i,]), end(roi_new[i,]))
    Chr<- as.character(seqnames(roi_new[i,]))
    range<- GRanges(
    seqnames = Chr,
    ranges = IRanges(start = lim[1]-ext,
                   end = lim[2]+ext))

    #get ideogramm tracks
    itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr,fontcolor="black")

    #genome axis track
    getrack <- GenomeAxisTrack(fontcolor="black")

    #deNOvo anno
    grtrack <- GeneRegionTrack(anno_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, # fill="darkgray",
        geneSymbol=TRUE, name="deNovo",collapseTrack=TRUE,color="black", cex.feature = 0.75,
        fontsize=fontSize,fontcolor.title="black", fontcolor.group="black", 
        Known=class_col["known"], Chimeric=class_col["chimeric (novel)"], NonChimeric=class_col["non-chimeric (novel)"])
    symbol(grtrack) <- as.character(anno_transcript_symbol[transcript(grtrack),])
    feature(grtrack) <- anno_transcript_class[transcript(grtrack),]

    #get annotation tracks of LTR12 repeats
    anno_repeats <- AnnotationTrack(repeats[grep("LTR12", repeats$repName),], name="LTR12", 
        shape = "box",fill="darkgray",color="darkgray", fontsize=fontSize,fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="darkgray", featureAnnotation="feature", cex.feature = 0.5)
    feature(anno_repeats)<- repeats[grep("LTR12", repeats$repName),]$repName

     #get annotation tracks of peptides
    anno_peptides <- AnnotationTrack(reduce(peptide_coord_sub), name="Peptides", 
        shape = "box",fill="black",color="black", fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="black", cex.feature = 0.5)


    #Data Tracks
    #rnaseq
    #get ylim
    ylim_rnaseq <- c(-(max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,]))))/20), 
        max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,])))))    
    yTicksAt_rnaseq <- c(0, plyr::round_any(max(ylim_rnaseq), accuracy=10, f = floor))
    #prepare tracks
    DMSO_rnaseq_track <- DataTrack(range = DMSO_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_rnaseq_track <- DataTrack(range = SB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nSB939", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_rnaseq_track <- DataTrack(range = DAC_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_rnaseq_track <- DataTrack(range = DACSB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])
   
    #nanopore
    nano_track <- DataTrack(range = DACSB_nanopore, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",#ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "Nano\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])

    #bam track for dacsb
    DAC_SB_alignment <- AlignmentsTrack("/omics/groups/OE0219/internal/tinat/210712_shortRead_processing_knownRef/HISAT2/aligned_sorted/AS-641779-LR-57123.sorted.bam",
        isPaired = TRUE, type=c("sashimi"),col.sashimi=treat_col["DACandSB939"],sashimiHeight=.5,max.height=.5,
        fontsize=fontSize,fontcolor.title="black",col.axis="black")#,  sashimiFilter =  intronicParts(anno_gtf) ,sashimiFilterTolerance = 10000L
        
    #cage
    #get ylim
    ylim_cage <- c(-(max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,]))))/20), max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,])))))
    yTicksAt_cage <- c(0, plyr::round_any(max(ylim_cage), accuracy=10, f = floor))
    #prepare tracks
    DMSO_cage_track <- DataTrack(range = DMSO_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_cage_track <- DataTrack(range = SB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nSB", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_cage_track <- DataTrack(range = DAC_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage,  yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_cage_track <- DataTrack(range = DACSB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])
    
    #Methylation
    DMSO_meth_track <- DataTrack(range = DMSO_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05), yTicksAt=c(0,1),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_meth_track <- DataTrack(range = SB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nSB", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_meth_track <- DataTrack(range = DAC_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_meth_track <- DataTrack(range = DACSB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])


    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withNano",
        paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    plotTracks(list(itrack ,getrack,anno_repeats,grtrack,# grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,nano_track,
        DMSO_cage_track,SB_cage_track,DAC_cage_track,DACSB_cage_track,
        DMSO_meth_track,SB_meth_track,DAC_meth_track,DACSB_meth_track
        ),fontcolor.title="black",sizes=c(0.5,1,0.5,2,1,1,1,1,1,1,1,1,1,1,1,1,1),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()

    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withAlign",
        paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    plotTracks(list(itrack ,getrack,anno_repeats,grtrack, # grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,DAC_SB_alignment,
        DMSO_cage_track,SB_cage_track,DAC_cage_track,DACSB_cage_track,
        DMSO_meth_track,SB_meth_track,DAC_meth_track,DACSB_meth_track
        ),fontcolor.title="black",sizes=c(0.5,1,0.5,2,1,1,1,1,2.5,1,1,1,1,1,1,1,1),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()
    print(i)
    
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_withAlign_Peptide",
    paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    plotTracks(list(itrack ,getrack,anno_repeats,anno_peptides,grtrack, # grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,DAC_SB_alignment,
        DMSO_cage_track,SB_cage_track,DAC_cage_track,DACSB_cage_track,
        DMSO_meth_track,SB_meth_track,DAC_meth_track,DACSB_meth_track
        ),fontcolor.title="black",sizes=c(0.5,1,0.5,0.75,2,1,1,1,1,2.5,1,1,1,1,1,1,1,1),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()
    print(i)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}



#only peptide sequence
ext <- 50
#peptide_coord
roi_new <- peptide_coord_sub[peptide_coord_sub$transcript_id %in% c("MSTRG.16894.3","MSTRG.28868.1", "MSTRG.29191.1"),]

#plot regions
for (i in 1:length(roi_new)){
    tryCatch({
    #Define Region
    lim <- c(start(roi_new[i,]), end(roi_new[i,]))
    Chr<- as.character(seqnames(roi_new[i,]))
    range<- GRanges(
    seqnames = Chr,
    ranges = IRanges(start = lim[1]-ext,
                   end = lim[2]+ext))

    #get ideogramm tracks
    itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr,fontcolor="black")

    #genome axis track
    getrack <- GenomeAxisTrack(fontcolor="black")

    #sequence track
    sTrack <- SequenceTrack(BSgenome.Hsapiens.UCSC.hg19, chromosome =Chr, from =lim[1]-ext, to = lim[2]+ext,cex=0.5)

    #deNOvo anno
    grtrack <- GeneRegionTrack(anno_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, # fill="darkgray",
        geneSymbol=TRUE, name="deNovo",collapseTrack=TRUE,color="black", cex.feature = 0.75,
        fontsize=fontSize,fontcolor.title="black", fontcolor.group="black", 
        Known=class_col["known"], Chimeric=class_col["chimeric (novel)"], NonChimeric=class_col["non-chimeric (novel)"])
    symbol(grtrack) <- as.character(anno_transcript_symbol[transcript(grtrack),])
    feature(grtrack) <- anno_transcript_class[transcript(grtrack),]

    #get annotation tracks of LTR12 repeats
    anno_repeats <- AnnotationTrack(repeats[grep("LTR12", repeats$repName),], name="LTR12", 
        shape = "box",fill="darkgray",color="darkgray", fontsize=fontSize,fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="darkgray", featureAnnotation="feature", cex.feature = 0.5)
    feature(anno_repeats)<- repeats[grep("LTR12", repeats$repName),]$repName

     #get annotation tracks of peptides
    anno_peptides <- AnnotationTrack(reduce(peptide_coord_sub), name="Peptides", 
        shape = "box",fill="black",color="black", fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="black", cex.feature = 0.5)


    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub_Peptide_only",
    paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 3)
    plotTracks(list(itrack,sTrack, getrack,anno_repeats,anno_peptides,grtrack
        ),fontcolor.title="black",sizes=c(0.5,0.5,1,0.5,0.75,2),
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()
    print(i)
    },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
























































#different regions selection criteria
#select region
temp <- "^FBP2$"
roi <- anno_original[grep(temp, anno_original$gene_name),]
 roi_new <-GRanges(
     seqnames = unique(seqnames(roi)),
     ranges = IRanges(start = min(start(roi)),
                    end =  max(end(roi))))
roi_new$transcript_id <- unique(roi$transcript_id)[2]
anno_transcript <- anno[anno$type=="transcript",]

#subset peptide list based on our ORFs == Universe
peptides_new_ORF <- peptides_new[peptides_new$Species =="ORFs",]
#get transcripts that give rise to uniquely presented peptides
peptides_new_ORF_noDMSO <- peptides_new_ORF[which(peptides_new_ORF$DMSO ==0),]
peptides_new_ORF_noDMSO_split <- split(peptides_new_ORF_noDMSO, peptides_new_ORF_noDMSO$DAC_SB)
lapply(peptides_new_ORF_noDMSO_split, length)
names(peptides_new_ORF_noDMSO_split)<- c("1/3 DAC + SB939 unique", "2/3 DAC + SB939 unique", "3/3 DAC + SB939 unique")
roi <- anno_transcript[anno_transcript$transcript_id %in%  peptides_new_ORF_noDMSO_split$"3/3 DAC + SB939 unique"$transcript_id,]
roi_new <- roi

#gene of interest
roi_new <- roi_new[roi_new$transcript_id =="MSTRG.28868.1",]
#plot regions
for (i in 1:length(roi_new)){
    #Define Region
    lim <- c(start(roi_new[i,]), end(roi_new[i,]))
    Chr<- as.character(seqnames(roi_new[i,]))
    range<- GRanges(
    seqnames = Chr,
    ranges = IRanges(start = lim[1]-ext,
                   end = lim[2]+ext))

    #get ideogramm tracks
    itrack <- IdeogramTrack(genome = "hg19", chromosome =Chr,fontcolor="black")

    #genome axis track
    getrack <- GenomeAxisTrack(fontcolor="black")

    #deNOvo anno
    grtrack <- GeneRegionTrack(anno_gtf, chromosome=Chr, transcriptAnnotation="symbol", showId=TRUE, # fill="darkgray",
        geneSymbol=TRUE, name="deNovo",collapseTrack=TRUE,color="black", cex.feature = 0.75,
        fontsize=fontSize,fontcolor.title="black", fontcolor.group="black", 
        Known=class_col["known"], Chimeric=class_col["chimeric (novel)"], NonChimeric=class_col["non-chimeric (novel)"])
    symbol(grtrack) <- as.character(anno_transcript_symbol[transcript(grtrack),])
    feature(grtrack) <- anno_transcript_class[transcript(grtrack),]

    #get annotation tracks of LTR12 repeats
    anno_repeats <- AnnotationTrack(repeats[grep("LTR12", repeats$repName),], name="LTR12", 
        shape = "box",fill="darkgray",color="darkgray", fontsize=fontSize,fontcolor.title="black",fontcolor.feature ="black",
        chromosome=Chr,genome = "hg19", fill="darkgray", featureAnnotation="feature", cex.feature = 0.5)
    feature(anno_repeats)<- repeats[grep("LTR12", repeats$repName),]$repName

    #Data Tracks
    #rnaseq
    #get ylim
    ylim_rnaseq <- c(-(max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,]))))/20), 
        max(c(max(score(DMSO_rnaseq[DMSO_rnaseq %over%range,])), max(score(SB_rnaseq[SB_rnaseq %over%range,])), 
        max(score(DAC_rnaseq[DAC_rnaseq %over%range,])), max(score(DACSB_rnaseq[DACSB_rnaseq %over%range,])))))    
    yTicksAt_rnaseq <- c(0, plyr::round_any(max(ylim_rnaseq), accuracy=10, f = floor))
    #prepare tracks
    DMSO_rnaseq_track <- DataTrack(range = DMSO_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black", ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_rnaseq_track <- DataTrack(range = SB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nSB939", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_rnaseq_track <- DataTrack(range = DAC_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_rnaseq_track <- DataTrack(range = DACSB_rnaseq, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_rnaseq, yTicksAt=yTicksAt_rnaseq,
        type = c("histogram"), chromosome = Chr, name = "RNAseq\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])
   
    #cage
    #get ylim
    ylim_cage <- c(-(max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,]))))/20), max(c(max(score(DMSO_cage[DMSO_cage %over%range,])), max(score(SB_cage[SB_cage %over%range,])), 
        max(score(DAC_cage[DAC_cage %over%range,])), max(score(DACSB_cage[DACSB_cage %over%range,])))))
    yTicksAt_cage <- c(0, plyr::round_any(max(ylim_cage), accuracy=10, f = floor))
    #prepare tracks
    DMSO_cage_track <- DataTrack(range = DMSO_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_cage_track <- DataTrack(range = SB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nSB", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_cage_track <- DataTrack(range = DAC_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage,  yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_cage_track <- DataTrack(range = DACSB_cage, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=ylim_cage, yTicksAt=yTicksAt_cage,
        type = c("histogram"), chromosome = Chr, name = "CAGE\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])
    
    #Methylation
    DMSO_meth_track <- DataTrack(range = DMSO_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05), yTicksAt=c(0,1),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDMSO", fill.histogram=treat_col["DMSO"],col.histogram=treat_col["DMSO"])
    SB_meth_track <- DataTrack(range = SB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nSB", fill.histogram=treat_col["SB939"],col.histogram=treat_col["SB939"])
    DAC_meth_track <- DataTrack(range = DAC_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC", fill.histogram=treat_col["DAC"],col.histogram=treat_col["DAC"])
    DACSB_meth_track <- DataTrack(range = DACSB_meth, genome = "hg19",
        fontsize=fontSize,fontcolor.title="black",col.axis="black",ylim=c(-.05,1.05),
        type = c("histogram"), chromosome = Chr, name = "WGBS\nDAC+SB", fill.histogram=treat_col["DACandSB939"],col.histogram=treat_col["DACandSB939"])


    #Plot track
    pdf(file.path("/omics/groups/OE0219/internal/tinat/210727_shortRead_processing_deNovo_custom4_quantification_analysis/locus_plots/sub",paste0(roi_new[i,]$transcript_id, "_locus_","ext_",ext, "_v1.pdf")), width = 12, height = 8)
    plotTracks(list(itrack ,getrack,anno_repeats,grtrack,  # grtrack_original,
        DMSO_rnaseq_track,SB_rnaseq_track,DAC_rnaseq_track,DACSB_rnaseq_track,
        DMSO_cage_track,SB_cage_track,DAC_cage_track,DACSB_cage_track,
        DMSO_meth_track,SB_meth_track,DAC_meth_track,DACSB_meth_track
        ),fontcolor.title="black",
        from =lim[1]-ext, to = lim[2]+ext, chromosome=Chr)
    dev.off()
 
    print(i)
}
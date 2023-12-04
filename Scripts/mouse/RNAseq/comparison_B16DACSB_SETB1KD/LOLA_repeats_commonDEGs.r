#Perform differential epression analysis of different contrasts
#libraries
library(DESeq2)
library(ggpubr)
library(rafalib)
library(pheatmap)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(biomaRt)
library(knitr)
library(Glimma)
library(limma)
library(edgeR)
library(randomcoloR)
library(clusterProfiler)
library(readxl)
library(rtracklayer)
library(LOLA)
library(ChIPpeakAnno)
library(rtracklayer)
library(ggpubr)
library(LOLA)
library(GenomicFeatures)
#function 
label_func <- function(x){
    breaks <- x
    #breaks[breaks>=500] <-  ">=500"
    #breaks <- droplevels(breaks)
    breaks

}

# 设置相关目录

#Directories
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAseq_deNovo_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_DACSB<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))

#setb1
base.dir<- "/omics/groups/OE0219/internal/tinat/mouse_project/220811_RNAseqSETB1KD_deNovoB16DACSB_analysis"
base_results.dir <- file.path(base.dir, "results")
results.dir<- file.path(base_results.dir , "tables")
PreDE.dir <- file.path(base_results.dir,"PreDE")
PostDE.dir <- file.path(base_results.dir,"PostDE")
#Read in Data
DEG_results_list_SETB1<- readRDS(file.path(PostDE.dir, "DEG_results_group_list.rds"))


# 提取上下调基因

#subset upregulated genes
DEG_results_list_DACSB_up_id <- lapply(DEG_results_list_DACSB, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})
DEG_results_list_SETB1_up_id <- lapply(DEG_results_list_SETB1, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange>lfc),]
    x$transcript_id
})

#subset upregulated genes
DEG_results_list_DACSB_down_id <- lapply(DEG_results_list_DACSB, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange<(-lfc)),]
    x$transcript_id
})
DEG_results_list_SETB1_down_id <- lapply(DEG_results_list_SETB1, function(x){
    x<- x[which(x$padj < alpha & x$log2FoldChange<(-lfc)),]
    x$transcript_id
})

#Read in Dataanno <- readRDS("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.annotated.sorted_repeat_anno.rds")
anno_original <-  import.gff2("/omics/groups/OE0219/internal/genomes/Mmusculus/mm10/rnaseq/gencode.vM19.annotation.gtf")
anno_classi <- as.data.frame(data.table::fread("/omics/groups/OE0219/internal/tinat/mouse_project/220809_RNAse_processing_deNovo/gffCompare.mergedTranscripts.gtf.tmap"))
repeats <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/","data","repeats_mm10.rds"))

## 上下调基因的阈值
# 
# 在这个代码片段中,alpha和lfc代表两个重要的参数:

# alpha:

# 它表示多重测试校正后的p值阈值,通常设置为0.05或0.01
# 这里设置为0.01,代表选择FDR(假发现率)小于0.01来判断差异表达的统计显著性
# lfc:

# 它是log2(Fold change)的缩写,代表基因表达折变的对数值
# 如果一个基因的上调和下调超过此log2FC阈值,则判定其为显著差异表达
# 这里设置为1,表示选择两个样本间基因表达变化2倍(对数值为1)以上或以下的基因

#Take a look at design and annotation of samples
design(dds)
#Set specifications
alpha <- 0.01 #set FDR cutoff
lfc <- 2##set logfold2 cutoff

# 从这几个变量的命名可以看出,它们与基因组中重复元素(repeats)的注释信息相关:

# repName:重复元素的具体名字,比如特定家族的LTR。
# repClass:重复元素的类型,包括DNA转座子、LINE、LTR、SINE等分类。
# repFamily:重复元素分类的家族,更为广义的类别,如ERV家族。
# 所以这三个数据库分别储存了基因组重复元素的不同级别的注释信息:

# repName_LTR_mm10 - 详细到具体LTR元素的名称
# repClass_mm10 - 中级分类的transposon类
# repFamily_mm10 - 更宏观的transposon家族分类

#load database
regionDB_repFamily_mm10 <- loadRegionDB(file.path(datasets.dir,"repFamily_mm10"))
regionDB_repClass_mm10 <- loadRegionDB(file.path(datasets.dir,"repClass_mm10"))
regionDB_repName_LTR_mm10 <- loadRegionDB(file.path(datasets.dir,"repName_LTR_mm10"))

# 运行LOLA 通过预定义好的基因组注释区域数据库,快速判断输入的基因列表是否与这些感兴趣区域存在重叠富集
# 在使用LOLA包进行基因组区域富集分析时,userSet和userUniverse是两个重要的概念:

# userSet:
# 用户定义的基因/区域集合,通常是差异表达基因或者GWAS中的峰值基因等。
# 需要判断这些基因与某类基因组注释区域是否存在重叠或关联。
# userSet放入富集分析的“检验集”。
# userUniverse:
# 背景基因组范围,用于统计检验。
# 定义在此基因组范围内,任选基因组成的集合是用户的“对照组”。
# 一般将某整个注释类别的所有基因作为userUniverse。


#run enrichment analysis
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript[anno_transcript$transcript_id %in% c(DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO,DEG_results_list_SETB1_up_id$SETDB1_vs_control ) ],1)
results_repFamily_mm10 <- list()
results_repClass_mm10 <- list()
results_repName_LTR_mm10 <-list()



# 这段代码的主要目的是找到两个条件下共同的差异基因,也即差异表达共性,然后对这些共性基因进行区域富集分析。主要步骤如下:

# 首先定义了list_oi,其中comparison列表包含:

# common_up: 两个条件下共同上调的差异基因及相关信息
# common_down: 两个条件下共同下调的差异基因及相关信息
# 然后对list_oi的每个元素,也就是共性基因进行循环分析:

# 构建UserSet, 包含共同上调和共同下调两个基因集
# 对其分别在3个重复元素相关的区域数据库中进行富集分析
# 存储不同区域数据库的富集结果到不同列表中
# 这样,我们可以得到两个条件下的差异表达共性基因,在不同的重复元素区域数据库中的富集结果。

# 最后,这些共性基因的区域富集信息,可以帮助我们判断差异表达的共性机制,比如是否与重复元素的激活相关。


list_oi <- list(comparison=list(common_up=DEG_results_list_DACSB$DACandSB939_vs_DMSO[DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO[DEG_results_list_DACSB_up_id$DACandSB939_vs_DMSO %in%  DEG_results_list_SETB1_up_id$SETDB1_vs_control],],
common_down=DEG_results_list_DACSB$DACandSB939_vs_DMSO[DEG_results_list_DACSB_down_id$DACandSB939_vs_DMSO[DEG_results_list_DACSB_down_id$DACandSB939_vs_DMSO %in%  DEG_results_list_SETB1_down_id$SETDB1_vs_control],]))
for(i in names(list_oi)){
    #get start position(strand specific)
    UserSet <- list(common_up=resize(anno_transcript[anno_transcript$transcript_id %in% list_oi[[i]][["common_up"]]$transcript_id,],1),
        common_down=resize(anno_transcript[anno_transcript$transcript_id %in% list_oi[[i]][["common_down"]]$transcript_id,],1))
    results_repFamily_mm10[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repFamily_mm10, cores=4)
    results_repClass_mm10[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repClass_mm10, cores=4)
    results_repName_LTR_mm10[[i]]= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_mm10, cores=4)
    print(i)
}
# 保存结果
dir.create(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA"))
saveRDS(results_repFamily_mm10, file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repFamily_mm10.rds"))
saveRDS(results_repClass_mm10, file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repClass_mm10.rds"))
saveRDS(results_repName_LTR_mm10, file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repName_LTR_mm10.rds"))

#results_repFamily_mm10 <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repFamily_mm10.rds"))
#results_repClass_mm10 <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repClass_mm10.rds"))
#results_repName_LTR_mm10 <- readRDS(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", "results_repName_LTR_mm10.rds"))

for(i in names(results_repName_LTR_mm10)){
    temp <- as.data.frame(results_repName_LTR_mm10[[i]])
    write.table(temp,file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD","LOLA", paste0("results_repName_LTR_mm10_", i,".txt")), col.names=TRUE, row.names=FALSE, quote=FALSE)
}

# 可视化 
#Plot genomic regions in bubble plot
#subset result of interest
for(i in names(list_oi)){
dir.create(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i,"LOLA"), recursive=TRUE)

results <- results_repFamily_mm10[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]

    
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
    
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i,"LOLA", "results_repFamily_LTR_mm10_top20.pdf"), height=10)
print(g)
dev.off()


    
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repFamily_LTR_mm10_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
if(length(idx)>0){
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repFamily_LTR_mm10_Sign.pdf"), height=20)
print(g)
dev.off()
}

#Plot genomic regions in bubble plot
#subset result of interest
results <- results_repClass_mm10[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repClass_LTR_mm10_top20.pdf"), height=10)
print(g)
dev.off()
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repClass_LTR_mm10_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
if(length(idx)>0){
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repClass_LTR_mm10_Sign.pdf"), height=20)
print(g)
dev.off()
}
    
#Plot genomic regions in bubble plot
#subset result of interest
results <- results_repName_LTR_mm10[[i]]
#plotting
idx <- head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA", "results_repName_LTR_mm10_top20.pdf"), height=10)
print(g)
dev.off()

idx <- head(unique(results$filename),5)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA","results_repName_LTR_mm10_top5.pdf"), height=10)
print(g)
dev.off()
    
#all
#plotting
idx <- head(unique(results$filename),Inf)
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=12, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=12, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=12, family="sans"), 
        legend.title=element_text(size=12, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA","results_repName_LTR_mm10_topInf.pdf"), height=20)
print(g)
dev.off()

#plot sign ones
idx <- unique(results[results$qValue < 0.05, ]$filename)
if(length(idx)>0){
results_sub <- results[results$filename %in% idx,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
g <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_colour_manual(values=c("grey", "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")
pdf(file.path("/omics/groups/OE0219/internal/tinat/mouse_project/comparison_DACSB_SETB1KD",i, "LOLA","results_repName_LTR_mm10_Sign.pdf"), height=20)
print(g)
dev.off()
}
print(i)
}





#combined plot
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript,1)
#get degs
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
UserSet <- lapply(DEG_results_list_sub, function(x){
#run enrichment analysis
resize(anno_transcript[anno_transcript$transcript_id %in%rownames(x[which(x$padj < alpha & x$log2FoldChange>lfc),]),],1)
})

results_repFamily_mm10= runLOLA(UserSet, userUnisverse, regionDB_repFamily_mm10, cores=4)
results_repClass_mm10= runLOLA(UserSet, userUnisverse, regionDB_repClass_mm10, cores=4)
results_repName_LTR_mm10= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_mm10, cores=4)

#prepare plotting
#results_repFamily_mm10
results <- results_repFamily_mm10
#idx_results_repFamily_mm10<- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repFamily_mm10<-head(unique(results$filename),1)
head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx_results_repFamily_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repFamily_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

#results_repClass_mm10
results <- results_repClass_mm10
#idx_results_repClass_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repClass_mm10 <-head(unique(results$filename),1)
results_sub <- results[results$filename %in% idx_results_repClass_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repClass_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")


#results_repName_LTR_mm10
results <- results_repName_LTR_mm10
#idx_results_repName_LTR_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
#idx_results_repName_LTR_mm10 <-head(unique(results$filename),5)
idx_results_repName_LTR_mm10 <-unique(results$filename[grep("LTR12",results$filename)])

results_sub <- results[results$filename %in% idx_results_repName_LTR_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repName_LTR_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+
scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

g <-cowplot::plot_grid(g_results_repFamily_mm10+coord_flip()+rremove("legend")+rremove("xlab"),
        g_results_repClass_mm10+coord_flip()+rremove("y.text")+rremove("legend")+rremove("xlab"),
        g_results_repName_LTR_mm10+coord_flip()+rremove("y.text")+rremove("xlab"),
        ncol=3,align = c("h"),#labels=c("Family", "Class", "Name"),
        rel_widths=c(length(idx_results_repFamily_mm10)+4, length(idx_results_repClass_mm10), length(idx_results_repName_LTR_mm10)))

ggsave(plot=g,file.path(PostDE.dir,"LOLA", "repeats_combined_DEG_up.pdf"),height = 6, width = 10, useDingbats = FALSE)





#same for only novel genes
#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id 
#combined plot
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript,1)
#get degs
DEG_results_list_sub <- DEG_results_list[c( "DAC_vs_DMSO","SB939_vs_DMSO","DACandSB939_vs_DMSO")]
UserSet <- lapply(DEG_results_list_sub, function(x){
        x <- dplyr::left_join(x, anno_classi, by="transcript_id")
        x <- x[x$class_code_simple!="known",]
        resize(anno_transcript[anno_transcript$transcript_id %in%rownames(x[which(x$padj < alpha & x$log2FoldChange>lfc),]),],1)
})

#run enrichment analysis
results_repFamily_mm10= runLOLA(UserSet, userUnisverse, regionDB_repFamily_mm10, cores=4)
results_repClass_mm10= runLOLA(UserSet, userUnisverse, regionDB_repClass_mm10, cores=4)
results_repName_LTR_mm10= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_mm10, cores=4)

#prepare plotting
#results_repFamily_mm10
results <- results_repFamily_mm10
#idx_results_repFamily_mm10<- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repFamily_mm10<-head(unique(results$filename),1)
head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx_results_repFamily_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repFamily_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

#results_repClass_mm10
results <- results_repClass_mm10
#idx_results_repClass_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repClass_mm10 <-head(unique(results$filename),1)
results_sub <- results[results$filename %in% idx_results_repClass_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repClass_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")


#results_repName_LTR_mm10
results <- results_repName_LTR_mm10
#idx_results_repName_LTR_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
#idx_results_repName_LTR_mm10 <-head(unique(results$filename),5)
idx_results_repName_LTR_mm10 <-unique(results$filename[grep("LTR12",results$filename)])

results_sub <- results[results$filename %in% idx_results_repName_LTR_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c("DACandSB939_vs_DMSO", "DAC_vs_DMSO", "SB939_vs_DMSO"))
#prepare plotting
g_results_repName_LTR_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+
scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

g <-cowplot::plot_grid(g_results_repFamily_mm10+coord_flip()+rremove("legend")+rremove("xlab"),
        g_results_repClass_mm10+coord_flip()+rremove("y.text")+rremove("legend")+rremove("xlab"),
        g_results_repName_LTR_mm10+coord_flip()+rremove("y.text")+rremove("xlab"),
        ncol=3,align = c("h"),#labels=c("Family", "Class", "Name"),
        rel_widths=c(length(idx_results_repFamily_mm10)+4, length(idx_results_repClass_mm10), length(idx_results_repName_LTR_mm10)))

ggsave(plot=g,file.path(PostDE.dir,"LOLA", "repeats_combined_DEG_up_onlyNovel.pdf"),height = 5, width = 10, useDingbats = FALSE)




#same for only dac sb genes strat
#same for only novel genes
#get class code information of de novo assembly, joint with number of exons
#rename class codes
anno_classi$class_code_simple  <- NA
anno_classi$class_code_simple <- ifelse(anno_classi$class_code == "=", "known", NA)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("s", "x", "i", "y", "p", "u"), "non-chimeric (novel)", anno_classi$class_code_simple)
anno_classi$class_code_simple <- ifelse(anno_classi$class_code %in% c("c", "k", "m", "n", "j", "e", "o"), "chimeric (novel)", anno_classi$class_code_simple)
anno_classi$transcript_id <- anno_classi$qry_id 
#combined plot
anno_transcript <- anno[anno$type == "transcript",]
userUnisverse <- resize(anno_transcript,1)
#get degs
DEG_results_list_sub <- lapply( DEG_results_list[c( "DACandSB939_vs_DMSO")], function(x){
        x <- dplyr::left_join(x, anno_classi, by="transcript_id")
})
DEG_split <-split(DEG_results_list_sub[[c( "DACandSB939_vs_DMSO")]], DEG_results_list_sub[[c( "DACandSB939_vs_DMSO")]]$class_code_simple)
UserSet <- lapply(DEG_split, function(x){
        resize(anno_transcript[anno_transcript$transcript_id %in% x[which(x$padj < alpha & x$log2FoldChange>lfc),]$transcript_id,],1)
})
#run enrichment analysis
results_repFamily_mm10= runLOLA(UserSet, userUnisverse, regionDB_repFamily_mm10, cores=4)
results_repClass_mm10= runLOLA(UserSet, userUnisverse, regionDB_repClass_mm10, cores=4)
results_repName_LTR_mm10= runLOLA(UserSet, userUnisverse, regionDB_repName_LTR_mm10, cores=4)

#prepare plotting
#results_repFamily_mm10
results <- results_repFamily_mm10
#idx_results_repFamily_mm10<- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repFamily_mm10<-head(unique(results$filename),1)
head(unique(results$filename),20)
results_sub <- results[results$filename %in% idx_results_repFamily_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c( "known" , "non-chimeric (novel)", "chimeric (novel)"))
#prepare plotting
g_results_repFamily_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

#results_repClass_mm10
results <- results_repClass_mm10
#idx_results_repClass_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
idx_results_repClass_mm10 <-head(unique(results$filename),1)
results_sub <- results[results$filename %in% idx_results_repClass_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c( "known" , "non-chimeric (novel)", "chimeric (novel)"))
#prepare plotting
g_results_repClass_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")


#results_repName_LTR_mm10
results <- results_repName_LTR_mm10
idx_results_repName_LTR_mm10 <- unique(results[results$qValue < 0.05, ]$filename)
#idx_results_repName_LTR_mm10 <-head(unique(results$filename),5)i 
idx_results_repName_LTR_mm10 <-unique(results$filename[grep("LTR12",results$filename)])

results_sub <- results[results$filename %in% idx_results_repName_LTR_mm10,]
#data preparation
combined_data <- results_sub[,c("userSet","dbSet", "pValueLog", "oddsRatio", "cellType" ,"filename")]#[userSet=="closed",]
combined_data$significant<- ifelse(combined_data$pValueLog< -log10(0.05), "No", "Yes" )
combined_data$cellType<- c(rep("AM", nrow(combined_data)))
#change infinite values
combined_data$pValueLog[is.infinite(combined_data$pValueLog)] <- 500
combined_data$filename<-sapply(strsplit(combined_data$filename,".",fixed=TRUE),`[`, 1)
combined_data$userSet <-factor(combined_data$userSet, levels = c( "known" , "non-chimeric (novel)", "chimeric (novel)"))
#prepare plotting
g_results_repName_LTR_mm10 <- ggplot(data = combined_data, aes(y=filename, x=userSet))+coord_fixed()+
geom_point(aes(size=pValueLog, fill=oddsRatio, color=significant), pch=21)+
#scale_fill_gradient2( midpoint = 1, low="#273871", high="#6D0026", name = "Odds Ratio")+
scale_fill_gradientn(colours=rev(hcl.colors(20,"Reds")), name = "Odds Ratio")+
scale_colour_manual(values=c(No="grey",Yes= "black"), name="Significant", drop=FALSE)+
scale_size(name="P-value\n(-log10)", labels = label_func) +
scale_y_discrete(limits=rev(levels(as.factor(combined_data$filename))))+
theme(text =element_text(size=14, color="black", family = "sans"),
        axis.ticks = element_blank(), axis.line = element_blank(), 
        axis.text.x=element_text(size=18, angle = 90, vjust = 0, color="black", family="sans"),
        axis.text.y=element_text(size=18, family="sans", color="black"))+
scale_x_discrete(name=NULL)+
theme(legend.text=element_text(size=18, family="sans"), 
        legend.title=element_text(size=18, family= "sans"),
        legend.background = element_rect(fill="white", color="white"),
        panel.background =  element_rect(fill="white"), panel.grid.major = element_line(color="lightgrey"),
        legend.key = element_rect(fill="white"))+rremove("ylab")

g <-cowplot::plot_grid(g_results_repFamily_mm10+coord_flip()+rremove("legend")+rremove("xlab"),
        g_results_repClass_mm10+coord_flip()+rremove("y.text")+rremove("legend")+rremove("xlab"),
        g_results_repName_LTR_mm10+coord_flip()+rremove("y.text")+rremove("xlab"),
        ncol=3,align = c("h"),#labels=c("Family", "Class", "Name"),
        rel_widths=c(length(idx_results_repFamily_mm10)+4, length(idx_results_repClass_mm10), length(idx_results_repName_LTR_mm10)))

ggsave(plot=g,file.path(PostDE.dir,"LOLA", "repeats_combined_DEG_up_ClassCodeStrat.pdf"),height = 5, width = 10, useDingbats = FALSE)



g <-cowplot::plot_grid(g_results_repFamily_mm10+coord_flip()+rremove("legend")+rremove("xlab")+rremove("y.text"),
        g_results_repClass_mm10+coord_flip()+rremove("y.text")+rremove("legend")+rremove("xlab"),
        g_results_repName_LTR_mm10+coord_flip()+rremove("y.text")+rremove("xlab"),
        ncol=3,align = c("h"),#labels=c("Family", "Class", "Name"),
        rel_widths=c(length(idx_results_repFamily_mm10), length(idx_results_repClass_mm10), length(idx_results_repName_LTR_mm10)))

ggsave(plot=g,file.path(PostDE.dir,"LOLA", "repeats_combined_DEG_up_ClassCodeStrat_v2.pdf"),height = 5, width = 12, useDingbats = FALSE)

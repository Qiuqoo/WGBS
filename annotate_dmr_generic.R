# 加载必要包
library(GenomicRanges)

# 设置工作目录
setwd("/data4/users/Qiuqoo/DRM2_WGBS/makeBSseq")

# 获取所有 _DMRs.RData 文件路径
dmr_files <- list.files(pattern = "_DMRs\\.RData$", full.names = TRUE)

# 创建列表存储 DMRs 和 GRanges
dmr_list <- list()
granges_list <- list()

# 遍历文件并转换为 GRanges
for (file in dmr_files) {
  cat("Loading:", file, "\n")
  
  load(file)  # 加载后的对象名为 dmrs
  
  # 创建合法名称作为列表键（替换 "_" 和 "-" 为 "_"）
  base_name <- gsub("_DMRs\\.RData$", "", basename(file))
  clean_name <- gsub("-", "_", base_name)  # 将 "-" 替换为 "_"
  
  # 存储原始 DMRs 对象
  dmr_list[[clean_name]] <- dmrs
  
  # 检查是否为 data.frame 格式，且包含 chr/start/end 等字段
  if (is.data.frame(dmrs) &&
      all(c("chr", "start", "end") %in% colnames(dmrs))) {
    
    # 创建 GRanges 对象
    gr <- GRanges(
      seqnames = dmrs$chr,
      ranges = IRanges(start = dmrs$start, end = dmrs$end),
      strand = "*"
    )
    
    # 添加其余 metadata
    mcols(gr) <- dmrs[, setdiff(colnames(dmrs), c("chr", "start", "end"))]
    
    # 存入 GRanges 列表
    granges_list[[clean_name]] <- gr
    
  } else {
    warning(paste("File", file, "does not contain expected DMR format."))
  }
}

# 示例：查看某个组的 GRanges 对象
granges_list$CG_9930_0h_vs_9930_72h

##############################################################################################################

# gene            mRNA            exon            CDS             five_prime_UTR  three_prime_UTR
ChineseLong_v3_Genes = import("/data3/users/Qiuqoo/CsV3_genome/ChineseLong_v3.gff3")
gene = ChineseLong_v3_Genes[ ChineseLong_v3_Genes$type == "gene", c( "ID", "Name" ) ]
unique(ChineseLong_v3_Genes$type)

seqlevels(gene, pruning.mode = "coarse") <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7")

#gene_UpDown_2k <- GeneUpDownStream(gene, upstream = 2000, downstream = 2000)
gene_Up_2k <- GeneUpDownStream(gene, upstream = 2000)

#gene_Up_2k ≈ promoter
seqlevels(gene_Up_2k, pruning.mode = "coarse") <- c("chr1", "chr2", "chr3",
                                                    "chr4", "chr5", "chr6", "chr7")

promoters <- gene_Up_2k



##############################################################################################################

#########TE
filename = "/data3/users/Qiuqoo/CsV3_genome/genome.all.repeat.sorted.gff"
TE <- import(filename, format = "gff", genome = "Cs",
             colnames = c( "type") )
seqlevels(TE, pruning.mode = "coarse") <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7")

TE_Transposon <- TE[TE$type == "Transposon"]
TE_TEprotein <- TE[TE$type == "TEprotein"]


##########################################################
# 加载必要包
library(GenomicRanges)

# 读取 miRNA 位置信息文件
mirna_df <- read.csv("/data4/users/Qiuqoo/DRM2_WGBS/match/micoRNA_cluster_location.csv")

# 检查数据格式
head(mirna_df)

# 转换为 GRanges 对象
mirna_gr <- GRanges(
  seqnames = mirna_df$Chromosome,
  ranges = IRanges(start = mirna_df$Start, end = mirna_df$End),
  strand = mirna_df$Strand,
  miRNA_Name = mirna_df$miRNA.Name  # 注意列名中的点
)

# 查看结果
mirna_gr
TE_Transposon
promoters

##############################################################
ry(GenomicRanges)

# 通用注释函数（前面定义好的）
annotate_dmr_generic <- function(dmr, grange) {
  hits <- findOverlaps(dmr, grange, ignore.strand = TRUE)
  dmr_matched <- dmr[queryHits(hits)]
  gr_matched <- grange[subjectHits(hits)]
  
  gr_df <- as.data.frame(gr_matched)
  dmr_df <- as.data.frame(dmr_matched)
  
  # 判断类型并提取对应字段
  if ("miRNA_Name" %in% colnames(mcols(grange))) {
    result <- cbind(dmr_df, miRNA_Name = mcols(gr_matched)$miRNA_Name)
  } else if ("type" %in% colnames(mcols(grange))) {
    result <- cbind(dmr_df, TE_Type = mcols(gr_matched)$type)
  } else if (all(c("ID", "Name") %in% colnames(mcols(grange)))) {
    result <- cbind(dmr_df, Gene_ID = mcols(gr_matched)$ID, Gene_Name = mcols(gr_matched)$Name)
  } else {
    warning("Unknown grange type: expected miRNA_Name, type, or ID/Name fields.")
    result <- dmr_df
  }
  
  return(result)
}

# 对 granges_list 中每个 DMR，批量注释三种功能区
dmr_annotations_all <- lapply(granges_list, function(dmr) {
  list(
    miRNA = annotate_dmr_generic(dmr, mirna_gr),
    TE = annotate_dmr_generic(dmr, TE_Transposon),
    promoter = annotate_dmr_generic(dmr, promoters)
  )
})

#dmr_annotations_all$CG_9930_0h_vs_9930_72h$miRNA


# 创建一个输出文件夹（可选）
setwd("/data4/users/Qiuqoo/DRM2_WGBS/match/")

dir.create("dmr_annotations_csv", showWarnings = FALSE)

# 遍历保存
for (dmr_name in names(dmr_annotations_all)) {
  result <- dmr_annotations_all[[dmr_name]]
  
  # 保存 miRNA 注释
  write.csv(result$miRNA, file = file.path("dmr_annotations_csv", paste0(dmr_name, "_miRNA.csv")), row.names = FALSE)
  
  # 保存 TE 注释
  write.csv(result$TE, file = file.path("dmr_annotations_csv", paste0(dmr_name, "_TE.csv")), row.names = FALSE)
  
  # 保存 promoter 注释
  write.csv(result$promoter, file = file.path("dmr_annotations_csv", paste0(dmr_name, "_promoter.csv")), row.names = FALSE)
}

#########################################
# 创建一个统计表
dmr_annotation_summary <- data.frame(
  DMR_Name = character(),
  miRNA_Count = integer(),
  TE_Count = integer(),
  Promoter_Count = integer(),
  stringsAsFactors = FALSE
)

# 遍历每个注释结果，统计数目
for (dmr_name in names(dmr_annotations_all)) {
  result <- dmr_annotations_all[[dmr_name]]
  
  dmr_annotation_summary <- rbind(
    dmr_annotation_summary,
    data.frame(
      DMR_Name = dmr_name,
      miRNA_Count = nrow(result$miRNA),
      TE_Count = nrow(result$TE),
      Promoter_Count = nrow(result$promoter),
      stringsAsFactors = FALSE
    )
  )
}

# 保存为 CSV 文件
write.csv(dmr_annotation_summary, file = "dmr_annotation_summary.csv", row.names = FALSE)


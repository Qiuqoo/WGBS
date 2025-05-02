setwd("/data4/users/Qiuqoo/DRM2_WGBS/methylation_extractor/")

library(DSS)
library(dplyr)
library(parallel)
library(purrr)

# 样品定义
samples <- list(
  "DRM2OE-0h" = c("DRM2OE-0h-1", "DRM2OE-0h-2", "DRM2OE-0h-3"),
  "DRM2OE-72h" = c("DRM2OE-72h-1", "DRM2OE-72h-2", "DRM2OE-72h-3"),
  "DRM2CR-0h" = c("DRM2CR-0h-1", "DRM2CR-0h-2", "DRM2CR-0h-3"),
  "DRM2CR-72h" = c("DRM2CR-72h-1", "DRM2CR-72h-2", "DRM2CR-72h-3"),
  "9930-0h"    = c("9930-0h-1", "9930-0h-2", "9930-0h-3"),
  "9930-72h"   = c("9930-72h-1", "9930-72h-2", "9930-72h-3")
)

# 所有要比较的组
comparisons <- list(
  c("9930-0h", "9930-72h"),
  c("DRM2OE-0h", "DRM2OE-72h"),
  c("DRM2CR-0h", "DRM2CR-72h"),
  c("9930-72h", "DRM2OE-72h"),
  c("9930-72h", "DRM2CR-72h"),
  c("DRM2OE-72h", "DRM2CR-72h")
)

# 染色体过滤
chr_filter <- paste0("chr", 1:7)

# 读取并预处理单个样本
read_and_prepare <- function(sample, context) {
  file_path <- paste0(sample, "_R1_val_1_bismark_bt2_pe.", context, "_report.txt")
  if (!file.exists(file_path)) {
    stop(paste("Missing file:", file_path))
  }
  df <- read.delim(file_path, header = FALSE)
  colnames(df) <- c("chr", "pos", "strand", "X", "uC", "Context", "Seq")
  df <- df %>%
    mutate(N = X + uC) %>%
    select(chr, pos, N, X) %>%
    filter(chr %in% chr_filter)
  return(df)
}

# 执行差异甲基化分析
process_context <- function(context, group1_name, group2_name) {
  cat("Processing", group1_name, "vs", group2_name, "in context", context, "\n")

  group1_samples <- samples[[group1_name]]
  group2_samples <- samples[[group2_name]]
  all_samples <- c(group1_samples, group2_samples)

  # 读取数据
  all_data <- map(all_samples, read_and_prepare, context = context)

  sample_names <- paste0(c(rep(group1_name, 3), rep(group2_name, 3)),
                         "_", rep(1:3, 2))

  # 构建 BSseq 对象
  BSobj <- makeBSseqData(all_data, sample_names)

  # 输出文件名前缀
  prefix <- paste0("../makeBSseq/", context, "_", group1_name, "_vs_", group2_name)

  # 保存 BSseq 对象
  save(BSobj, file = paste0(prefix, "_BSobj.RData"))

  # 差异位点检测
  dml_test <- DMLtest(BSobj,
                      group1 = sample_names[1:3],
                      group2 = sample_names[4:6],
                      smoothing = TRUE,
                      smoothing.span = 100)

  save(dml_test, file = paste0(prefix, "_DMLtest.RData"))

  # 调用 DMR
  dmrs <- callDMR(dml_test, p.threshold = 0.001, minlen = 50, minCG = 3, dis.merge = 100)
  save(dmrs, file = paste0(prefix, "_DMRs.RData"))
}

# 并行运行所有比较（每个比较可并行 3 个 context）
mclapply(comparisons, function(pair) {
  mclapply(c("CG", "CHG", "CHH"), function(ct) {
    process_context(ct, pair[1], pair[2])
  }, mc.cores = 3)
}, mc.cores = 6)  # 外层并行处理比较组

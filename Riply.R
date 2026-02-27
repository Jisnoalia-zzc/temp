library(spatstat)
library(tidyverse)
library(ggplot2)
library(patchwork)
library(cowplot)

# 根据图片定义映射关系
mimer_mapping <- data.frame(
  MIMER = c(
    "CD4_C7_0X40", "CD4_C7_0X40",
    "CD8_C6_CD39", "CD8_C6_CD39", "CD8_C6_CD39",
    "Endo_C3_RGCC", "Endo_C3_RGCC",
    "FB_C3_COL1A1", "FB_C3_COL1A1",
    "Mac_C2_SPP1", "Mac_C2_SPP1"
  ),
  CODEX = c(
    "CD4T-aTreg", "Tcells.cc",
    "PDCD1+Tex", "LAG3+Tex", "PDCD1+Tex2",
    "PLVAP+Endo", "CollagenIV+Endo",
    "COL1A1+CAF", "POSTN+CAF",
    "ISG15+Mac", "SPP1+Mac"
  )
)

ct = data.frame(
  L4_C = c('Ki-67+Bcells','activateBcells','Bcells','ImmatureBcells','ISG15+Bcells','matureBcells',
           'PLVAP+Endo','Vimentin+Endo','PDPN+Endo','CollagenIV+Endo','Ki-67+Endo',
           'non-canonical_epi1','non-canonical_epi2','non-canonical_epi3','GP100+tumor','Ki-67+cancer','differentiating_cancer','basal-like_epi','differentiating_epi',
           'MyoFibroblast','COL1A1+CAF','ISG15+CAF','POSTN+CAF',
           'ISG15+Mac','M1like-Mac','M2like-Mac','SPP1+Mac','Ki-67+Mac',
           'CD4T-Treg','CD4T-aTreg','CD4Th','ISG15_CD4T','CD4TRTE','Tcells.cc',
           'CD8T','LAG3+Tex','PDCD1+Tex','PDCD1+Tex2','CD8Tm',
           'NKT'),
  L3_C = c('Bcells','Bcells','Bcells','Bcells','Bcells','Bcells',
           'Endo','Endo','Endo','Endo','Endo',
           'Epithelium','Epithelium','Epithelium','Epithelium','Epithelium','Epithelium','Epithelium','Epithelium',
           'Fibroblast','Fibroblast','Fibroblast','Fibroblast',
           'Myeloids','Myeloids','Myeloids','Myeloids','Myeloids',
           'CD4T','CD4T','CD4T','CD4T','CD4T','CD4T',
           'CD8T','CD8T','CD8T','CD8T','CD8T',
           'NKT')
)


# 从映射表中提取唯一的MIMER类型
mimer_types <- unique(mimer_mapping$MIMER)

load('~/workspace/proj_ESCC_STW_ZWM_2022_01/codex/codex_final.Rdata')
codex2 = subset(codex2,subCelltype=='Epithelium',invert=T)
Idents(codex2) = codex2$subCelltype
library(qs)
qsave(codex2,file = '~/workspace/proj_ESCC_STW_ZWM_2022_01/zzc/final_analysis/final_codex.qs')

head(ct)
setdiff(ct$L4_C,unique(codex2$subCelltype))
setdiff(unique(codex2$subCelltype),ct$L4_C)
# 设置分析参数
r_values <- seq(0, 1000, by = 50)
target_r <- 366
n_perm <- 1000

# 假设您的CODEX数据已经加载为一个名为codex_data的数据框
# 并且包含以下列：
# - L4_C: 细胞类型（精细分类）
# - L3_C: 细胞谱系（大类）
# - coord_x: x坐标
# - coord_y: y坐标

analyze_spatial_clustering <- function(codex_data, sample_name, mimer_mapping, r_values, target_r, n_perm) {
  
  # 使用提供的CODEX数据
  cell_trek_df <- codex_data
  
  # 使用逻辑条件排除特定行
  exclude_types <- c("Cancer", "Epithelium")
  cell_trek_df <- cell_trek_df[!cell_trek_df$L4_C %in% exclude_types, ]
  
  if (nrow(cell_trek_df) == 0) {
    stop("After excluding Cancer and Epithelium, no data remains.")
  }
  
  # 验证步骤1: 检查L4_C列的细胞类型组成
  cat("验证L4_C列的细胞类型组成:\n")
  
  # 获取所有唯一的L4_C细胞类型
  all_cell_types <- unique(cell_trek_df$L4_C)
  cat("总共有", length(all_cell_types), "种不同的L4_C细胞类型\n")
  
  # 获取MIMER对应的细胞类型
  mimer_cell_types <- unique(mimer_mapping$CODEX)
  cat("其中", length(mimer_cell_types), "种是MIMER对应的细胞类型\n")
  
  # 获取非MIMER细胞类型
  non_mimer_cell_types <- setdiff(all_cell_types, mimer_cell_types)
  cat("剩下的", length(non_mimer_cell_types), "种是非MIMER细胞类型\n")
  
  if (length(non_mimer_cell_types) > 0) {
    cat("前10个非MIMER细胞类型:", paste(head(non_mimer_cell_types, 10), collapse = ", "), "\n")
    
    # 统计非MIMER细胞类型的数量
    non_mimer_counts <- table(cell_trek_df$L4_C[cell_trek_df$L4_C %in% non_mimer_cell_types])
    cat("非MIMER细胞总数:", sum(non_mimer_counts), "\n")
  }
  
  # 使用映射表创建MIMER标记
  cell_trek_df$is_mimer <- ifelse(
    cell_trek_df$L4_C %in% mimer_mapping$CODEX, 
    "Mimer", 
    "Other"
  )
  
  # 检查MIMER和非MIMER细胞的比例
  mimer_count <- sum(cell_trek_df$is_mimer == "Mimer")
  other_count <- sum(cell_trek_df$is_mimer == "Other")
  total_count <- nrow(cell_trek_df)
  
  cat("\n细胞类型统计:\n")
  cat("MIMER细胞:", mimer_count, "(", round(mimer_count/total_count*100, 1), "%)\n")
  cat("非MIMER细胞:", other_count, "(", round(other_count/total_count*100, 1), "%)\n")
  cat("总细胞数:", total_count, "\n")
  
  # 检查每个MIMER类型对应的细胞类型是否存在
  cat("\n验证每个MIMER类型的覆盖情况:\n")
  validation_result <- validate_mimer_coverage(cell_trek_df, mimer_mapping)
  print(validation_result)
  
  # 创建研究区域窗口
  x_range <- range(cell_trek_df$coord_x)
  y_range <- range(cell_trek_df$coord_y)
  study_win <- owin(x_range, y_range)
  
  # 创建观测点模式
  obs_ppp <- ppp(x = cell_trek_df$coord_x, 
                 y = cell_trek_df$coord_y,
                 window = study_win,
                 marks = factor(cell_trek_df$is_mimer))
  
  # 计算观测L函数
  L_obs <- Lcross(obs_ppp, i = "Mimer", j = "Mimer",
                  r = r_values, correction = "iso")
  
  L_obs_df <- data.frame(r = L_obs$r, L_r = L_obs$iso - L_obs$r, Type = "Observed")
  obs_value <- L_obs$iso[which.min(abs(L_obs$r - target_r))] - target_r
  
  # 策略1：在谱系内部L3_C层面置换
  cat("\n开始策略1置换 (谱系内部置换)...\n")
  perm_curves_strat1 <- matrix(NA, nrow = n_perm, ncol = length(r_values))
  perm_results_strat1 <- numeric(n_perm)
  
  # 策略2：在所有细胞层面完全随机打乱标签
  cat("开始策略2置换 (完全随机打乱)...\n")
  perm_curves_strat2 <- matrix(NA, nrow = n_perm, ncol = length(r_values))
  perm_results_strat2 <- numeric(n_perm)
  
  for (perm in 1:n_perm) {
    
    # 策略1：在谱系内部L3_C层面置换
    cell_trek_df_perm1 <- cell_trek_df
    
    unique_l3 <- unique(cell_trek_df$L3_C)
    
    for (l3 in unique_l3) {
      l3_indices <- which(cell_trek_df$L3_C == l3)
      if (length(l3_indices) > 1) {
        # 在谱系内部随机打乱细胞类型标签
        original_labels <- cell_trek_df$L4_C[l3_indices]
        shuffled_labels <- sample(original_labels)
        cell_trek_df_perm1$L4_C[l3_indices] <- shuffled_labels
      }
    }
    
    # 重新标记MIMER状态
    cell_trek_df_perm1$is_mimer_perm <- ifelse(
      cell_trek_df_perm1$L4_C %in% mimer_mapping$CODEX,
      "Mimer", "Other"
    )
    
    # 计算置换后的L函数
    perm_ppp1 <- ppp(x = cell_trek_df_perm1$coord_x,
                     y = cell_trek_df_perm1$coord_y,
                     window = study_win,
                     marks = factor(cell_trek_df_perm1$is_mimer_perm))
    
    L_perm1 <- Lcross(perm_ppp1, i = "Mimer", j = "Mimer",
                      r = r_values, correction = "iso")
    
    perm_curves_strat1[perm, ] <- L_perm1$iso - L_perm1$r
    perm_value1 <- L_perm1$iso[which.min(abs(L_perm1$r - target_r))] - target_r
    perm_results_strat1[perm] <- perm_value1
    
    # 策略2：在所有细胞层面完全随机打乱标签
    cell_trek_df_perm2 <- cell_trek_df
    
    # 在所有细胞中完全随机打乱L4_C标签
    all_labels <- cell_trek_df$L4_C
    cell_trek_df_perm2$L4_C <- sample(all_labels)
    
    # 重新标记MIMER状态
    cell_trek_df_perm2$is_mimer_perm <- ifelse(
      cell_trek_df_perm2$L4_C %in% mimer_mapping$CODEX,
      "Mimer", "Other"
    )
    
    # 计算置换后的L函数
    perm_ppp2 <- ppp(x = cell_trek_df_perm2$coord_x,
                     y = cell_trek_df_perm2$coord_y,
                     window = study_win,
                     marks = factor(cell_trek_df_perm2$is_mimer_perm))
    
    L_perm2 <- Lcross(perm_ppp2, i = "Mimer", j = "Mimer",
                      r = r_values, correction = "iso")
    
    perm_curves_strat2[perm, ] <- L_perm2$iso - L_perm2$r
    perm_value2 <- L_perm2$iso[which.min(abs(L_perm2$r - target_r))] - target_r
    perm_results_strat2[perm] <- perm_value2
    
    # 显示进度
    if (perm %% 100 == 0) {
      cat(sprintf("  完成 %d/%d 次置换...\n", perm, n_perm))
    }
  }
  
  # 计算p值
  p_value_strat1 <- (sum(perm_results_strat1 >= obs_value) + 1) / (n_perm + 1)
  p_value_strat2 <- (sum(perm_results_strat2 >= obs_value) + 1) / (n_perm + 1)
  
  # 准备置换曲线数据
  prepare_perm_data <- function(perm_curves, strategy_name) {
    perm_curves_df <- as.data.frame(perm_curves)
    colnames(perm_curves_df) <- r_values
    perm_curves_df$perm <- 1:n_perm
    
    perm_curves_long <- perm_curves_df %>%
      pivot_longer(cols = -perm, names_to = "r", values_to = "L_value") %>%
      mutate(r = as.numeric(r), Strategy = strategy_name)
    
    return(perm_curves_long)
  }
  
  perm_curves_long_strat1 <- prepare_perm_data(perm_curves_strat1, "Within Lineage (L3_C)")
  perm_curves_long_strat2 <- prepare_perm_data(perm_curves_strat2, "Global Shuffle (All Cells)")
  
  # 计算均值曲线
  mean_curve_strat1 <- perm_curves_long_strat1 %>%
    group_by(r) %>%
    summarise(mean_L = mean(L_value), .groups = 'drop') %>%
    mutate(Strategy = "Within Lineage (L3_C)")
  
  mean_curve_strat2 <- perm_curves_long_strat2 %>%
    group_by(r) %>%
    summarise(mean_L = mean(L_value), .groups = 'drop') %>%
    mutate(Strategy = "Global Shuffle (All Cells)")
  
  # 创建曲线图
  curve_plot <- ggplot() +
    # 策略1的个体置换曲线
    geom_line(data = perm_curves_long_strat1,
              aes(x = r, y = L_value, group = perm),
              color = "lightgray", alpha = 0.3) +
    # 策略2的个体置换曲线
    geom_line(data = perm_curves_long_strat2,
              aes(x = r, y = L_value, group = perm),
              color = "lightpink", alpha = 0.2) +
    # 策略1的均值曲线
    geom_line(data = mean_curve_strat1, 
              aes(x = r, y = mean_L), 
              color = "blue", linetype = "dashed", linewidth = 0.8) +
    # 策略2的均值曲线
    geom_line(data = mean_curve_strat2, 
              aes(x = r, y = mean_L), 
              color = "purple", linetype = "dashed", linewidth = 0.8) +
    # 观测曲线
    geom_line(data = L_obs_df, 
              aes(x = r, y = L_r), 
              color = "red", linewidth = 1.2) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_vline(xintercept = target_r, linetype = "dotted", color = "darkgreen") +
    labs(title = paste("L-function:", sample_name),
         subtitle = paste("P-value (Within Lineage) =", round(p_value_strat1, 4), 
                          "| P-value (Global) =", round(p_value_strat2, 4)),
         x = "Distance r", 
         y = "L(r) - r") +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          plot.subtitle = element_text(size = 8, hjust = 0.5),
          axis.title = element_text(size = 8)) +
    scale_color_manual(name = "Strategy",
                       values = c("Within Lineage (L3_C)" = "blue", 
                                  "Global Shuffle (All Cells)" = "purple",
                                  "Observed" = "red")) +
    guides(color = guide_legend(override.aes = list(linetype = c("dashed", "dashed", "solid"),
                                                    linewidth = c(0.8, 0.8, 1.2))))
  
  # 创建直方图
  hist_data <- data.frame(
    perm_results = c(perm_results_strat1, perm_results_strat2),
    Strategy = rep(c("Within Lineage (L3_C)", "Global Shuffle (All Cells)"), each = n_perm)
  )
  
  hist_plot <- ggplot(hist_data, aes(x = perm_results, fill = Strategy)) +
    geom_histogram(bins = 20, alpha = 0.6, position = "identity", color = "black") +
    geom_vline(xintercept = obs_value, color = "red", linetype = "dashed", linewidth = 1) +
    geom_vline(data = hist_data %>% group_by(Strategy) %>% summarise(mean_val = mean(perm_results)),
               aes(xintercept = mean_val, color = Strategy), 
               linetype = "solid", linewidth = 0.8) +
    labs(title = paste("Permutation at r =", target_r),
         x = "L(r) - r",
         y = "Frequency") +
    annotate("text", x = Inf, y = Inf,
             label = paste("P (Within) =", round(p_value_strat1, 4), "\n",
                           "P (Global) =", round(p_value_strat2, 4), "\n",
                           "Obs =", round(obs_value, 3)),
             hjust = 1.1, vjust = 1.1, size = 2.5) +
    scale_fill_manual(values = c("Within Lineage (L3_C)" = "lightblue", 
                                 "Global Shuffle (All Cells)" = "pink")) +
    scale_color_manual(values = c("Within Lineage (L3_C)" = "blue", 
                                  "Global Shuffle (All Cells)" = "purple")) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, hjust = 0.5),
          axis.title = element_text(size = 8),
          legend.position = "bottom")
  
  return(list(
    curve_plot = curve_plot,
    hist_plot = hist_plot,
    p_value_strat1 = p_value_strat1,
    p_value_strat2 = p_value_strat2,
    obs_value = obs_value,
    sample_name = sample_name,
    validation = validation_result,
    cell_type_summary = data.frame(
      Total_cell_types = length(all_cell_types),
      MIMER_cell_types = length(mimer_cell_types),
      Non_MIMER_cell_types = length(non_mimer_cell_types),
      MIMER_cell_count = mimer_count,
      Non_MIMER_cell_count = other_count,
      Total_cell_count = total_count
    )
  ))
}

# 验证函数
validate_mimer_coverage <- function(cell_trek_df, mimer_mapping) {
  results <- data.frame()
  mimer_types <- unique(mimer_mapping$MIMER)
  
  for (mimer in mimer_types) {
    # 获取这个MIMER类型对应的所有CODEX细胞类型
    codex_types <- mimer_mapping$CODEX[mimer_mapping$MIMER == mimer]
    
    # 检查这些细胞类型在数据中是否存在
    present_types <- intersect(codex_types, unique(cell_trek_df$L4_C))
    
    # 统计细胞数量
    cell_counts <- sapply(codex_types, function(ct) {
      sum(cell_trek_df$L4_C == ct, na.rm = TRUE)
    })
    
    result <- data.frame(
      MIMER_Type = mimer,
      Expected_CODEX_Types = paste(codex_types, collapse = ", "),
      Present_CODEX_Types = paste(present_types, collapse = ", "),
      Missing_CODEX_Types = paste(setdiff(codex_types, present_types), collapse = ", "),
      N_Expected_Types = length(codex_types),
      N_Present_Types = length(present_types),
      Cell_Counts = paste(cell_counts, collapse = ", "),
      Total_Cells = sum(cell_counts)
    )
    
    results <- rbind(results, result)
  }
  
  return(results)
}

# 主分析流程
# 假设您已经有一个或多个CODEX数据集存储在列表中
# codex_datasets 应该是一个命名列表，每个元素是一个数据框
# 每个数据框包含 L4_C, L3_C, coord_x, coord_y 列

# 示例: 如果您有多个样本
# codex_datasets <- list(
#   sample1 = codex_data_sample1,
#   sample2 = codex_data_sample2,
#   ...
# )

codex2_df = codex2@meta.data
head(codex2_df)


codex2_df$Sample <- chartr("ABCDEFGHI", "BCDEFGHIJ", codex2_df$Studylevel2)
table(codex2_df$Celltype.x,codex2_df$subCelltype)

del = setdiff(unique(codex2$subCelltype),ct$L4_C)
codex2_df = codex2_df %>%
  filter(Celltype.x != 'Epithelium') %>%
  filter(!subCelltype %in% del)

codex2_df %>%
  mutate(L4_C = as.character(subCelltype),
         L3_C = as.character(Celltype.x),
         coord_x = x,
         coord_y = y) %>%
  dplyr::select(L4_C,L3_C,Sample,coord_x,coord_y) -> codex2_df_final
head(codex2_df_final)

codex2_list = split(codex2_df_final,codex2_df_final$Sample)

codex_datasets = codex2_list
save(codex_datasets,file = 'codex_datasets.Rdata')
# 循环处理所有CODEX数据集
all_plots <- list()
results_summary <- data.frame()

save(all_plots,results_summary,file = 'all_results_20251226.Rdata')
for (i in seq_along(codex_datasets)) {
  sample_name <- names(codex_datasets)[i]
  codex_data <- codex_datasets[[i]]
  
  cat("\n=========================================\n")
  cat("处理样本", i, "/", length(codex_datasets), ":", sample_name, "\n")
  cat("=========================================\n")
  
  tryCatch({
    # 执行分析
    result <- analyze_spatial_clustering(codex_data, sample_name, mimer_mapping, r_values, target_r, n_perm)
    
    # 存储图形
    all_plots[[paste0("curve_", i)]] <- result$curve_plot
    all_plots[[paste0("hist_", i)]] <- result$hist_plot
    
    # 存储结果摘要
    results_summary <- rbind(results_summary, 
                             data.frame(Sample = result$sample_name,
                                        P_value_Within_Lineage = result$p_value_strat1,
                                        P_value_Global_Shuffle = result$p_value_strat2,
                                        Observed_value = result$obs_value,
                                        Total_cell_types = result$cell_type_summary$Total_cell_types,
                                        MIMER_cell_types = result$cell_type_summary$MIMER_cell_types,
                                        Non_MIMER_cell_types = result$cell_type_summary$Non_MIMER_cell_types,
                                        MIMER_cell_count = result$cell_type_summary$MIMER_cell_count,
                                        Non_MIMER_cell_count = result$cell_type_summary$Non_MIMER_cell_count,
                                        Total_cell_count = result$cell_type_summary$Total_cell_count))
    
  }, error = function(e) {
    cat("处理样本时出错", sample_name, ":", e$message, "\n")
  })
}

# 构建拼接图形
plot_list <- list()
plot_counter <- 1

for (i in 1:length(codex_datasets)) {
  curve_name <- paste0("curve_", i)
  hist_name  <- paste0("hist_",  i)
  
  if (curve_name %in% names(all_plots) && hist_name %in% names(all_plots)) {
    # 每样本：curve 在上，hist 在下
    combined <- all_plots[[curve_name]] / all_plots[[hist_name]] + plot_layout(heights = c(2, 1))
    plot_list[[plot_counter]] <- combined
    plot_counter <- plot_counter + 1
  } else {
    warning(paste("样本", i, "的图形缺失:", names(codex_datasets)[i]))
  }
}

# 组合所有样本
if (length(plot_list) > 0) {
  n_cols <- min(8, length(plot_list))
  n_rows <- ceiling(length(plot_list) / n_cols)
  base_height_per_row <- 4
  
  final_plot <- wrap_plots(plot_list, ncol = n_cols)
  
  ggsave("spatial_clustering_summary_codex.png", 
         plot = final_plot, 
         width = 20, 
         height = n_rows * base_height_per_row,
         limitsize = FALSE,
         units = "in", dpi = 300)
  
  cat("\n图像已保存为 'spatial_clustering_summary_codex.png'\n")
} else {
  stop("No plots were successfully generated.")
}
save(all_plots,results_summary,file = 'all_results_20251226.Rdata')
# 保存结果摘要表格
write.csv(results_summary, "spatial_clustering_pvalues_codex.csv", row.names = FALSE)
cat("汇总结果已保存为 'spatial_clustering_pvalues_codex.csv'\n")

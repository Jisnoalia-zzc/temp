library(qs)
library(dplyr)
codex = qread(file = '~/workspace/proj_ESCC_STW_ZWM_2022_01/zzc/final_analysis/final_codex.qs')

codex@meta.data %>% 
  dplyr::select(subCelltype,Celltype.x) %>%
  distinct(subCelltype,.keep_all = T) %>% 
  arrange(desc(Celltype.x))
load('~/workspace/proj_ESCC_STW_ZWM_2022_01/liuliqiu/final_analysis/codex_datasets.Rdata')

library(tidyverse)
library(RANN)
library(ggplot2)
library(patchwork)

# 定义5个MIMER细胞类型
mimer_types <- c(
  "CD4_C7_0X40",
  "CD8_C6_CD39", 
  "Endo_C3_RGCC",
  "FB_C3_COL1A1",
  "Mac_C2_SPP1"
)

# 定义MIMER到CODEX的映射
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

# 1. 核心函数：计算每个位置的S_r(i)统计量
calculate_Sr_per_cell <- function(coords, mimer_labels, r, mimer_types) {
  n_cells <- nrow(coords)
  S_r_vector <- rep(0, n_cells)
  
  # 使用RANN包中的nn2函数进行高效近邻搜索
  nearest_neighbors <- RANN::nn2(
    data = coords, 
    query = coords, 
    k = n_cells,  # 查找所有可能的邻居
    radius = r,
    searchtype = "radius"
  )
  
  for (i in 1:n_cells) {
    # 获取当前细胞在半径r内的邻居索引（不包括自身）
    neighbor_indices <- nearest_neighbors$nn.idx[i, ]
    neighbor_indices <- neighbor_indices[neighbor_indices > 0]  # 去除填充的0
    neighbor_indices <- neighbor_indices[neighbor_indices != i]  # 排除自身
    
    if (length(neighbor_indices) > 0) {
      # 获取邻居的MIMER类型
      neighbor_types <- mimer_labels[neighbor_indices]
      neighbor_types <- neighbor_types[neighbor_types != "Non_MIMER"]
      
      # 检查是否5类MIMER细胞都存在
      if (all(mimer_types %in% neighbor_types)) {
        S_r_vector[i] <- 1
      }
    }
  }
  
  return(S_r_vector)
}

# 2. 计算平均S_bar(r) = (1/N) * ΣS_r(i)
calculate_S_bar <- function(coords, mimer_labels, r, mimer_types) {
  S_r <- calculate_Sr_per_cell(coords, mimer_labels, r, mimer_types)
  S_bar <- mean(S_r)
  return(list(S_bar = S_bar, S_r = S_r))
}

# 3. 准备数据函数
prepare_codrex_data <- function(codex_data, mimer_mapping, mimer_types) {
  required_cols <- c("L4_C", "L3_C", "coord_x", "coord_y")
  if (!all(required_cols %in% colnames(codex_data))) {
    stop("数据必须包含以下列: ", paste(required_cols, collapse = ", "))
  }
  
  # 标记MIMER细胞类型
  codex_data <- codex_data %>%
    mutate(
      mimer_type = case_when(
        L4_C %in% mimer_mapping$CODEX ~ 
          mimer_mapping$MIMER[match(L4_C, mimer_mapping$CODEX)],
        TRUE ~ "Non_MIMER"
      ),
      is_mimer_cell = L4_C %in% mimer_mapping$CODEX
    )
  
  return(codex_data)
}

# 4. 选择与MIMER细胞数量比例相近的其他细胞类型
select_similar_proportion_cells <- function(cell_data, mimer_mapping) {
  # 统计各种细胞类型的数量
  cell_type_counts <- table(cell_data$L4_C)
  
  # 计算MIMER细胞类型的总数量
  mimer_cell_types <- mimer_mapping$CODEX
  mimer_cell_counts <- cell_type_counts[names(cell_type_counts) %in% mimer_cell_types]
  mimer_total_count <- sum(mimer_cell_counts)
  
  # 找出与MIMER细胞总数量比例相近的其他细胞类型
  # 我们选择数量在0.5-2倍范围内的细胞类型
  all_cell_types <- names(cell_type_counts)
  non_mimer_types <- setdiff(all_cell_types, mimer_cell_types)
  
  # 计算每个非MIMER细胞类型的数量
  similar_types <- character(0)
  for (cell_type in non_mimer_types) {
    count <- cell_type_counts[cell_type]
    ratio <- count / mimer_total_count
    
    # 如果数量比例在0.2-5倍范围内，则认为是"比例相近"
    if (ratio >= 0.2 && ratio <= 5) {
      similar_types <- c(similar_types, cell_type)
    }
  }
  
  # 返回所有需要置换的细胞类型
  all_perm_types <- c(mimer_cell_types, similar_types)
  
  return(list(
    all_perm_types = all_perm_types,
    mimer_total_count = mimer_total_count,
    similar_types = similar_types,
    similar_type_counts = cell_type_counts[similar_types]
  ))
}

# 5. 新的置换策略：对MIMER细胞+比例相近的其他细胞进行置换
# 策略1：在谱系内部置换
permute_labels_within_lineage_proportional <- function(cell_data, mimer_mapping) {
  # 选择需要置换的细胞类型
  type_selection <- select_similar_proportion_cells(cell_data, mimer_mapping)
  perm_types <- type_selection$all_perm_types
  
  # 标记需要置换的细胞
  cell_data <- cell_data %>%
    mutate(
      to_permute = L4_C %in% perm_types
    )
  
  permuted_data <- cell_data
  
  for (lineage in unique(cell_data$L3_C)) {
    lineage_indices <- which(cell_data$L3_C == lineage)
    
    if (length(lineage_indices) > 1) {
      # 获取当前谱系中需要置换的细胞索引
      perm_indices <- lineage_indices[cell_data$to_permute[lineage_indices]]
      
      if (length(perm_indices) > 1) {
        # 只对这些细胞的L4_C标签进行随机打乱
        original_labels <- cell_data$L4_C[perm_indices]
        shuffled_labels <- sample(original_labels)
        permuted_data$L4_C[perm_indices] <- shuffled_labels
      }
    }
  }
  
  # 重新标记MIMER类型
  permuted_data <- permuted_data %>%
    mutate(
      mimer_type = case_when(
        L4_C %in% mimer_mapping$CODEX ~ 
          mimer_mapping$MIMER[match(L4_C, mimer_mapping$CODEX)],
        TRUE ~ "Non_MIMER"
      )
    )
  
  return(list(
    permuted_data = permuted_data,
    selected_types = type_selection
  ))
}

# 策略2：全局置换
permute_labels_global_proportional <- function(cell_data, mimer_mapping) {
  # 选择需要置换的细胞类型
  type_selection <- select_similar_proportion_cells(cell_data, mimer_mapping)
  perm_types <- type_selection$all_perm_types
  
  # 标记需要置换的细胞
  cell_data <- cell_data %>%
    mutate(
      to_permute = L4_C %in% perm_types
    )
  
  permuted_data <- cell_data
  
  # 获取所有需要置换的细胞索引
  perm_indices <- which(cell_data$to_permute)
  
  if (length(perm_indices) > 1) {
    # 只对这些细胞的L4_C标签进行随机打乱
    original_labels <- cell_data$L4_C[perm_indices]
    shuffled_labels <- sample(original_labels)
    permuted_data$L4_C[perm_indices] <- shuffled_labels
  }
  
  # 重新标记MIMER类型
  permuted_data <- permuted_data %>%
    mutate(
      mimer_type = case_when(
        L4_C %in% mimer_mapping$CODEX ~ 
          mimer_mapping$MIMER[match(L4_C, mimer_mapping$CODEX)],
        TRUE ~ "Non_MIMER"
      )
    )
  
  return(list(
    permuted_data = permuted_data,
    selected_types = type_selection
  ))
}

# 6. 置换检验函数
perform_permutation_test_single_r <- function(cell_data, r, mimer_types, mimer_mapping, 
                                              n_perm = 1000, method = "within_lineage") {
  
  # 计算观测值
  obs_coords <- as.matrix(cell_data[, c("coord_x", "coord_y")])
  obs_labels <- cell_data$mimer_type
  
  obs_result <- calculate_S_bar(obs_coords, obs_labels, r, mimer_types)
  obs_S_bar <- obs_result$S_bar
  
  # 执行置换
  perm_S_bar <- numeric(n_perm)
  selected_types_info <- NULL
  
  for (i in 1:n_perm) {
    if (method == "within_lineage") {
      perm_result <- permute_labels_within_lineage_proportional(cell_data, mimer_mapping)
      permuted_data <- perm_result$permuted_data
      if (i == 1) selected_types_info <- perm_result$selected_types
    } else if (method == "global") {
      perm_result <- permute_labels_global_proportional(cell_data, mimer_mapping)
      permuted_data <- perm_result$permuted_data
      if (i == 1) selected_types_info <- perm_result$selected_types
    } else {
      stop("置换方法必须是 'within_lineage' 或 'global'")
    }
    
    # 计算置换后的S_bar
    perm_coords <- as.matrix(permuted_data[, c("coord_x", "coord_y")])
    perm_labels <- permuted_data$mimer_type
    
    perm_result_val <- calculate_S_bar(perm_coords, perm_labels, r, mimer_types)
    perm_S_bar[i] <- perm_result_val$S_bar
    
    # 显示进度
    if (i %% 100 == 0) {
      cat(sprintf("  完成 %d/%d 次置换...\n", i, n_perm))
    }
  }
  
  # 计算p值
  p_value <- (sum(perm_S_bar >= obs_S_bar) + 1) / (n_perm + 1)
  
  return(list(
    obs_S_bar = obs_S_bar,
    perm_S_bar = perm_S_bar,
    p_value = p_value,
    obs_S_r = obs_result$S_r,
    selected_types_info = selected_types_info
  ))
}

# 7. 多尺度分析函数
analyze_multi_scale_with_envelope <- function(cell_data, mimer_types, mimer_mapping, 
                                              r_values = seq(10, 200, by = 10), 
                                              n_perm = 100, method = "within_lineage") {
  
  cat("进行多尺度分析 (r =", min(r_values), "到", max(r_values), "μm)...\n")
  
  # 计算观测曲线
  obs_S_bar_values <- numeric(length(r_values))
  obs_S_r_list <- list()
  
  cat("计算观测曲线...\n")
  for (j in 1:length(r_values)) {
    r <- r_values[j]
    
    obs_coords <- as.matrix(cell_data[, c("coord_x", "coord_y")])
    obs_labels <- cell_data$mimer_type
    
    obs_result <- calculate_S_bar(obs_coords, obs_labels, r, mimer_types)
    obs_S_bar_values[j] <- obs_result$S_bar
    obs_S_r_list[[j]] <- obs_result$S_r
    
    if (j %% 5 == 0) {
      cat(sprintf("  完成 r=%dμm\n", r))
    }
  }
  
  # 计算置换包络
  cat("计算置换包络...\n")
  perm_matrix <- matrix(NA, nrow = n_perm, ncol = length(r_values))
  selected_types_info <- NULL
  
  for (i in 1:n_perm) {
    if (i %% 10 == 0) cat(sprintf("  置换 %d/%d\n", i, n_perm))
    
    # 创建置换数据
    if (method == "within_lineage") {
      perm_result <- permute_labels_within_lineage_proportional(cell_data, mimer_mapping)
      permuted_data <- perm_result$permuted_data
      if (i == 1) selected_types_info <- perm_result$selected_types
    } else if (method == "global") {
      perm_result <- permute_labels_global_proportional(cell_data, mimer_mapping)
      permuted_data <- perm_result$permuted_data
      if (i == 1) selected_types_info <- perm_result$selected_types
    }
    
    # 对每个半径计算置换后的S_bar
    perm_coords <- as.matrix(permuted_data[, c("coord_x", "coord_y")])
    perm_labels <- permuted_data$mimer_type
    
    for (j in 1:length(r_values)) {
      r <- r_values[j]
      perm_result_val <- calculate_S_bar(perm_coords, perm_labels, r, mimer_types)
      perm_matrix[i, j] <- perm_result_val$S_bar
    }
  }
  
  # 计算包络统计量
  envelope_stats <- data.frame(
    r = r_values,
    obs_S_bar = obs_S_bar_values,
    perm_mean = apply(perm_matrix, 2, mean, na.rm = TRUE),
    perm_upper = apply(perm_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE)),
    perm_lower = apply(perm_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
  )
  
  return(list(
    envelope_stats = envelope_stats,
    obs_S_bar_values = obs_S_bar_values,
    perm_matrix = perm_matrix,
    obs_S_r_list = obs_S_r_list,
    selected_types_info = selected_types_info
  ))
}

# 8. 主分析函数
analyze_co_neighborhood_single_sample <- function(cell_data, sample_name, mimer_types, mimer_mapping, 
                                                  target_r = 50, n_perm = 1000, 
                                                  multi_scale = TRUE, r_values = seq(10, 200, by = 10)) {
  
  cat("=========================================\n")
  cat("分析样本:", sample_name, "\n")
  cat("=========================================\n")
  
  # 检查必需的5类MIMER细胞是否存在
  mimer_cells <- cell_data %>% 
    filter(mimer_type %in% mimer_types)
  
  cat("MIMER细胞类型统计:\n")
  mimer_counts <- table(mimer_cells$mimer_type)
  for (type in mimer_types) {
    count <- ifelse(type %in% names(mimer_counts), mimer_counts[type], 0)
    cat(sprintf("  %s: %d 个细胞\n", type, count))
  }
  
  # 检查是否所有5类都存在
  missing_types <- setdiff(mimer_types, unique(mimer_cells$mimer_type))
  if (length(missing_types) > 0) {
    warning(sprintf("样本 %s 缺少以下MIMER类型: %s\n", 
                    sample_name, paste(missing_types, collapse = ", ")))
    return(NULL)
  }
  
  # 1. 单个半径的置换检验
  cat("\n执行单个半径置换检验 (r =", target_r, "μm)...\n")
  perm_within <- perform_permutation_test_single_r(
    cell_data, target_r, mimer_types, mimer_mapping, 
    n_perm = n_perm, method = "within_lineage"
  )
  
  perm_global <- perform_permutation_test_single_r(
    cell_data, target_r, mimer_types, mimer_mapping, 
    n_perm = n_perm, method = "global"
  )
  
  # 显示选择的细胞类型信息
  cat("\n置换过程中选择的细胞类型:\n")
  cat("MIMER细胞类型总数:", perm_within$selected_types_info$mimer_total_count, "\n")
  cat("选中的非MIMER细胞类型数量:", length(perm_within$selected_types_info$similar_types), "\n")
  cat("选中的非MIMER细胞类型:", 
      paste(perm_within$selected_types_info$similar_types, collapse = ", "), "\n")
  
  cat("\n单个半径结果汇总 (r =", target_r, "μm):\n")
  cat(sprintf("观测 S_bar(r) = %.4f\n", perm_within$obs_S_bar))
  cat(sprintf("谱系内部置换 p值 = %.4f\n", perm_within$p_value))
  cat(sprintf("全局置换 p值 = %.4f\n", perm_global$p_value))
  
  # 2. 多尺度分析
  multi_scale_within <- NULL
  multi_scale_global <- NULL
  
  if (multi_scale) {
    cat("\n执行多尺度分析...\n")
    multi_scale_within <- analyze_multi_scale_with_envelope(
      cell_data, mimer_types, mimer_mapping, 
      r_values = r_values, 
      n_perm = 100, method = "within_lineage"
    )
    
    multi_scale_global <- analyze_multi_scale_with_envelope(
      cell_data, mimer_types, mimer_mapping, 
      r_values = r_values, 
      n_perm = 100, method = "global"
    )
  }
  
  # 3. 绘制图形
  
  # 3.1 单个半径置换分布图
  plot_permutation_distribution <- function(obs_value, perm_values, r_value, p_value, title) {
    perm_df <- data.frame(perm_S_bar = perm_values)
    
    p <- ggplot(perm_df, aes(x = perm_S_bar)) +
      geom_histogram(bins = 30, fill = "lightblue", color = "black", alpha = 0.7) +
      geom_vline(xintercept = obs_value, color = "red", 
                 linetype = "dashed", linewidth = 1.2) +
      geom_vline(xintercept = mean(perm_values), color = "blue", 
                 linetype = "solid", linewidth = 0.8) +
      annotate("text", x = Inf, y = Inf,
               label = paste("obs =", round(obs_value, 4), "\n",
                             "pvalue =", round(p_value, 4)),
               hjust = 1.1, vjust = 1.1, size = 4) +
      labs(
        title = paste(title, " (r =", r_value, "μm)"),
        x = expression(bar(S)(r)),
        y = "Freq"
      ) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    
    return(p)
  }
  
  p_within <- plot_permutation_distribution(
    perm_within$obs_S_bar, perm_within$perm_S_bar, target_r, 
    perm_within$p_value, "谱系内部置换"
  )
  
  p_global <- plot_permutation_distribution(
    perm_global$obs_S_bar, perm_global$perm_S_bar, target_r, 
    perm_global$p_value, "全局置换"
  )
  
  # 3.2 多尺度曲线图
  plot_multi_scale_curve <- function(envelope_stats, title) {
    p <- ggplot(envelope_stats, aes(x = r)) +
      # 置信区间
      geom_ribbon(aes(ymin = perm_lower, ymax = perm_upper), 
                  fill = "lightblue", alpha = 0.3) +
      # 置换均值
      geom_line(aes(y = perm_mean, color = "置换均值"), 
                linetype = "dashed", linewidth = 0.8) +
      # 观测曲线
      geom_line(aes(y = obs_S_bar, color = "观测值"), 
                linewidth = 1.2) +
      # 标记目标半径
      geom_vline(xintercept = target_r, color = "red", 
                 linetype = "dashed", alpha = 0.5) +
      labs(
        title = title,
        x = "半径 r (μm)",
        y = expression(bar(S)(r) ~ "比例"),
        color = ""
      ) +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.position = "bottom"
      ) +
      scale_color_manual(values = c("观测值" = "red", "置换均值" = "blue"))
    
    return(p)
  }
  
  p_multi_within <- NULL
  p_multi_global <- NULL
  
  if (multi_scale) {
    p_multi_within <- plot_multi_scale_curve(
      multi_scale_within$envelope_stats, 
      paste(sample_name, "- 谱系内部置换 (多尺度)")
    )
    
    p_multi_global <- plot_multi_scale_curve(
      multi_scale_global$envelope_stats,
      paste(sample_name, "- 全局置换 (多尺度)")
    )
  }
  
  # 返回结果
  result <- list(
    sample_name = sample_name,
    target_r = target_r,
    obs_S_bar_single = perm_within$obs_S_bar,
    p_value_within = perm_within$p_value,
    p_value_global = perm_global$p_value,
    selected_types_info = perm_within$selected_types_info,
    cell_data_with_Sr = cell_data %>% mutate(S_r = perm_within$obs_S_r),
    permutation_plot_within = p_within,
    permutation_plot_global = p_global
  )
  
  if (multi_scale) {
    result$multi_scale_within <- multi_scale_within
    result$multi_scale_global <- multi_scale_global
    result$multi_scale_plot_within <- p_multi_within
    result$multi_scale_plot_global <- p_multi_global
  }
  
  return(result)
}

# 9. 批量分析函数
batch_analyze_co_neighborhood <- function(codex_datasets, mimer_mapping, mimer_types, 
                                          target_r = 50, n_perm = 1000, 
                                          multi_scale = TRUE, r_values = seq(10, 200, by = 10)) {
  
  all_results <- list()
  summary_df <- data.frame()
  type_selection_info <- list()
  
  for (sample_name in names(codex_datasets)) {
    cat("\n", strrep("=", 50), "\n", sep = "")
    cat("处理样本:", sample_name, "\n")
    cat(strrep("=", 50), "\n", sep = "")
    
    cell_data <- codex_datasets[[sample_name]]
    
    # 预处理数据
    cell_data <- prepare_codrex_data(cell_data, mimer_mapping, mimer_types)
    
    # 分析
    result <- analyze_co_neighborhood_single_sample(
      cell_data, sample_name, mimer_types, mimer_mapping,
      target_r, n_perm, multi_scale, r_values
    )
    
    if (!is.null(result)) {
      all_results[[sample_name]] <- result
      type_selection_info[[sample_name]] <- result$selected_types_info
      
      # 添加到摘要
      summary_df <- rbind(summary_df, data.frame(
        Sample = sample_name,
        Target_r_um = target_r,
        Obs_S_bar_single = result$obs_S_bar_single,
        P_value_Within = result$p_value_within,
        P_value_Global = result$p_value_global,
        Sig_Within = ifelse(result$p_value_within < 0.05, "Yes", "No"),
        Sig_Global = ifelse(result$p_value_global < 0.05, "Yes", "No"),
        Total_Cells = nrow(cell_data),
        MIMER_Cells = sum(cell_data$mimer_type %in% mimer_types),
        Non_MIMER_Cells = sum(cell_data$mimer_type == "Non_MIMER"),
        Selected_NonMIMER_Types = length(result$selected_types_info$similar_types)
      ))
    }
  }
  
  # 生成汇总结果
  if (nrow(summary_df) > 0) {
    # 保存结果
    write.csv(summary_df, "co_neighborhood_summary.csv", row.names = FALSE)
    cat("\n结果已保存到: co_neighborhood_summary.csv\n")
    
    # 保存细胞类型选择信息
    type_selection_df <- do.call(rbind, lapply(names(type_selection_info), function(sample_name) {
      info <- type_selection_info[[sample_name]]
      data.frame(
        Sample = sample_name,
        MIMER_Total_Count = info$mimer_total_count,
        Selected_NonMIMER_Types = paste(info$similar_types, collapse = "; "),
        Selected_NonMIMER_Counts = paste(info$similar_type_counts, collapse = "; ")
      )
    }))
    write.csv(type_selection_df, "cell_type_selection_info.csv", row.names = FALSE)
    
    # 创建汇总条形图
    p_summary <- ggplot(summary_df, aes(x = Sample, y = Obs_S_bar_single, fill = Sig_Within)) +
      geom_bar(stat = "identity") +
      geom_hline(yintercept = 0, linetype = "solid") +
      geom_text(aes(label = sprintf("%.4f\n(p=%.3f)", Obs_S_bar_single, P_value_Within)), 
                vjust = -0.5, size = 3) +
      labs(
        title = "五类同时邻近统计量汇总",
        subtitle = paste("目标半径 r =", target_r, "μm"),
        x = "样本",
        y = expression(bar(S)(r) ~ "比例"),
        fill = "显著 (p<0.05)"
      ) +
      theme_bw() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5)
      )
    
    ggsave("co_neighborhood_summary.png", p_summary, 
           width = 8, height = 6, dpi = 300)
    
    # 组合置换检验图
    plot_list <- list()
    for (i in seq_along(all_results)) {
      sample_name <- names(all_results)[i]
      p_combined <- (all_results[[sample_name]]$permutation_plot_within + 
                       all_results[[sample_name]]$permutation_plot_global) +
        plot_annotation(title = paste("样本:", sample_name))
      plot_list[[i]] <- p_combined
    }
    
    n_cols <- 2
    n_rows <- ceiling(length(plot_list) / n_cols)
    combined_plots <- wrap_plots(plot_list, ncol = n_cols, nrow = n_rows)
    
    ggsave("co_neighborhood_permutation_plots.png", combined_plots,
           width = 10 * n_cols, height = 6 * n_rows, dpi = 300, limitsize = FALSE)
    
    # 保存多尺度图
    if (multi_scale) {
      multi_plot_list <- list()
      for (i in seq_along(all_results)) {
        sample_name <- names(all_results)[i]
        if (!is.null(all_results[[sample_name]]$multi_scale_plot_within)) {
          p_multi <- (all_results[[sample_name]]$multi_scale_plot_within + 
                        all_results[[sample_name]]$multi_scale_plot_global) +
            plot_annotation(title = paste("样本:", sample_name, "- 多尺度分析"))
          multi_plot_list[[i]] <- p_multi
        }
      }
      
      if (length(multi_plot_list) > 0) {
        combined_multi_plots <- wrap_plots(multi_plot_list, ncol = 1, nrow = length(multi_plot_list))
        
        ggsave("co_neighborhood_multi_scale_plots.png", combined_multi_plots,
               width = 12, height = 6 * length(multi_plot_list), dpi = 300, limitsize = FALSE)
      }
    }
    
    cat("\n图形已保存:\n")
    cat("  - co_neighborhood_summary.png: 汇总条形图\n")
    cat("  - co_neighborhood_permutation_plots.png: 置换检验图\n")
    cat("  - cell_type_selection_info.csv: 细胞类型选择信息\n")
    if (multi_scale) cat("  - co_neighborhood_multi_scale_plots.png: 多尺度曲线图\n")
  }
  
  return(list(
    all_results = all_results,
    summary = summary_df,
    type_selection_info = type_selection_info
  ))
}

# 使用示例
# 假设您的数据格式如下：
# tmp <- list(
#   B01 = data.frame(L4_C = ..., L3_C = ..., coord_x = ..., coord_y = ...),
#   B02 = data.frame(L4_C = ..., L3_C = ..., coord_x = ..., coord_y = ...),
#   ...
# )
tmp = codex_datasets[1]
# 运行分析
results <- batch_analyze_co_neighborhood(
  tmp,
  mimer_mapping,
  mimer_types,
  target_r = 50,      # 您关注的单个半径
  n_perm = 1000,      # 置换次数
  multi_scale = TRUE, # 进行多尺度分析
  r_values = seq(10, 200, by = 10)  # 多尺度分析的范围
)

library(tidyverse)
library(ez)
library(ggplot2)
library(reshape2)
library(RColorBrewer)

data_p200 <- read.csv("/Users/yuanxiaoyu/Desktop/project/dataset/data_P200.csv")
data_n400 <- read.csv("/Users/yuanxiaoyu/Desktop/project/dataset/data_N400.csv")
data_p600 <- read.csv("/Users/yuanxiaoyu/Desktop/project/dataset/data_P600.csv")

# 定义每个ERP成分对应的脑区功能分区
define_brain_regions_by_erp <- function() {
  brain_regions <- list(
    "P200" = list(
      "Frontal" = c("F3", "F4", "FZ", "F5", "F6"),
      "Frontocentral" = c("FC3", "FC4", "FCZ", "FC5", "FC6"),
      "Frontotemporal" = c("FT7", "FT8"),
      "Central" = c("C3", "C4", "CZ"),
      "Centroparietal" = c("CP3", "CP4", "CPZ"),
      "Occipital" = c("O1", "O2")
    ),
    "N400" = list(
      "Frontal" = c("F3", "F4", "FZ", "F5", "F7"),
      "Frontocentral" = c("FC3", "FC4", "FCZ", "FC5"),
      "Central" = c("C3", "C4", "CZ"),
      "Centroparietal" = c("CP3", "CP4", "CPZ"),
      "Parietal" = c("P3", "P4", "Pz"),
      "Occipital" = c("O1", "O2")
    ),
    "P600" = list(
      "Frontal" = c("F3", "F4", "FZ"),
      "Frontopolar" = c("FP1", "FP2"),
      "Frontocentral" = c("FC3", "FC4", "FCZ"),
      "Central" = c("C3", "C4", "CZ"),
      "Centroparietal" = c("CP3", "CP4", "CPZ"),
      "Parietal" = c("P3", "P4", "Pz", "P5", "P6"),
      "Temporal" = c("T7", "T8")
    )
  )
  return(brain_regions)
}

# 将电极位置映射到脑区
map_electrodes_to_regions <- function(data, erp_name) {
  brain_regions <- define_brain_regions_by_erp()
  erp_regions <- brain_regions[[erp_name]]
  
  # 创建位置到脑区的映射
  position_to_region <- data.frame(
    position = character(),
    brain_region = character(),
    stringsAsFactors = FALSE
  )
  
  for (region_name in names(erp_regions)) {
    positions <- erp_regions[[region_name]]
    temp_df <- data.frame(
      position = positions,
      brain_region = region_name,
      stringsAsFactors = FALSE
    )
    position_to_region <- rbind(position_to_region, temp_df)
  }
  
  # 合并数据
  data_with_regions <- data %>%
    left_join(position_to_region, by = "position") %>%
    filter(!is.na(brain_region))  # 移除未分类的电极
  
  return(data_with_regions)
}

# 对脑区进行平均，减少条件数量
aggregate_by_brain_region <- function(data, measure_col, erp_name) {
  data_with_regions <- map_electrodes_to_regions(data, erp_name)
  
  # 按脑区平均
  if (measure_col == "amplitude") {
    aggregated_data <- data_with_regions %>%
      group_by(sample_name, word_class, brain_region) %>%
      summarise(
        avg_amplitude = mean(amplitude, na.rm = TRUE),
        n_electrodes = n(),
        .groups = "drop"
      )
  } else {
    aggregated_data <- data_with_regions %>%
      group_by(sample_name, word_class, brain_region) %>%
      summarise(
        avg_latency = mean(latency, na.rm = TRUE),
        n_electrodes = n(),
        .groups = "drop"
      )
  }
  
  return(aggregated_data)
}

# 运行ERP分析
run_erp_analysis <- function(data, measure, erp_name, positions) {
  cat("\n", rep("=", 60), "\n")
  cat("分析", erp_name, "成分的", ifelse(measure == "amplitude", "峰值", "延迟"), "\n")
  cat(rep("=", 60), "\n")
  
  # 过滤指定位置的数据
  data_filtered <- data %>%
    filter(position %in% positions) %>%
    filter(!is.na(if(measure == "amplitude") amplitude else latency))
  
  cat("原始数据点数:", nrow(data_filtered), "\n")
  cat("被试数:", length(unique(data_filtered$sample_name)), "\n")
  cat("电极位置数:", length(unique(data_filtered$position)), "\n")
  cat("词类数:", length(unique(data_filtered$word_class)), "\n")
  
  # 聚合到脑区级别
  aggregated_data <- aggregate_by_brain_region(data_filtered, measure, erp_name)
  
  cat("\n聚合后的数据:\n")
  cat("数据点数:", nrow(aggregated_data), "\n")
  cat("脑区:", paste(unique(aggregated_data$brain_region), collapse = ", "), "\n")
  cat("词类:", paste(unique(aggregated_data$word_class), collapse = ", "), "\n")
  
  # 转换为因子
  aggregated_data$sample_name <- as.factor(aggregated_data$sample_name)
  aggregated_data$word_class <- as.factor(aggregated_data$word_class)
  aggregated_data$brain_region <- as.factor(aggregated_data$brain_region)
  
  # 检查数据完整性
  n_regions <- length(unique(aggregated_data$brain_region))
  n_word_classes <- length(unique(aggregated_data$word_class))
  expected_per_subject <- n_regions * n_word_classes
  
  complete_check <- aggregated_data %>%
    group_by(sample_name) %>%
    summarise(
      n_obs = n(),
      expected = expected_per_subject,
      complete = (n_obs == expected),
      .groups = "drop"
      
    )
  
  complete_subjects <- complete_check$sample_name[complete_check$complete]
  cat("完整数据的被试数:", length(complete_subjects), "\n")
  
  # 如果完整数据的被试太少，尝试放宽条件
  if (length(complete_subjects) < 15) {
    cat("完整数据的被试较少，尝试使用至少有80%数据的被试...\n")
    partial_complete_subjects <- complete_check$sample_name[complete_check$n_obs >= 0.8 * expected_per_subject]
    cat("至少80%数据的被试数:", length(partial_complete_subjects), "\n")
    
    if (length(partial_complete_subjects) >= 15) {
      complete_subjects <- partial_complete_subjects
    }
  }
  
  # 过滤数据
  final_data <- aggregated_data %>%
    filter(sample_name %in% complete_subjects) %>%
    droplevels()
  
  cat("最终分析的被试数:", length(unique(final_data$sample_name)), "\n")
  cat("最终数据点数:", nrow(final_data), "\n")
  
  # 计算描述性统计
  if (measure == "amplitude") {
    descriptive_stats <- final_data %>%
      group_by(word_class, brain_region) %>%
      summarise(
        mean_amp = mean(avg_amplitude, na.rm = TRUE),
        sd_amp = sd(avg_amplitude, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    cat("\n描述性统计 (振幅 μV):\n")
    # 修改：显示所有行
    print(descriptive_stats, n = Inf)
  } else {
    descriptive_stats <- final_data %>%
      group_by(word_class, brain_region) %>%
      summarise(
        mean_lat = mean(avg_latency, na.rm = TRUE),
        sd_lat = sd(avg_latency, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      )
    cat("\n描述性统计 (潜伏期 ms):\n")
    # 修改：显示所有行
    print(descriptive_stats, n = Inf)
  }
  
  # 运行ANOVA
  tryCatch({
    cat("\n正在运行重复测量ANOVA...\n")
    
    if (measure == "amplitude") {
      anova_result <- ezANOVA(
        data = final_data,
        dv = avg_amplitude,
        wid = sample_name,
        within = .(word_class, brain_region),
        detailed = TRUE,
        type = 2
      )
    } else {
      anova_result <- ezANOVA(
        data = final_data,
        dv = avg_latency,
        wid = sample_name,
        within = .(word_class, brain_region),
        detailed = TRUE,
        type = 2
      )
    }
    
    cat("\nANOVA结果:\n")
    anova_table <- anova_result$ANOVA
    print(anova_table)
    
    # 解释显著性结果
    cat("\n显著性结果解释:\n")
    significant_effects <- anova_table[anova_table$p < 0.05, ]
    if (nrow(significant_effects) > 0) {
      for (i in 1:nrow(significant_effects)) {
        effect <- significant_effects$Effect[i]
        p_val <- significant_effects$p[i]
        f_val <- significant_effects$F[i]
        cat(sprintf("- %s: F = %.3f, p = %.3f (显著)\n", effect, f_val, p_val))
      }
    } else {
      cat("- 没有发现显著效应 (p < 0.05)\n")
    }
    
    # 球形性检验
    if (!is.null(anova_result$`Mauchly's Test for Sphericity`)) {
      cat("\n球形性检验:\n")
      mauchly_table <- anova_result$`Mauchly's Test for Sphericity`
      print(mauchly_table)
      
      sphericity_violations <- mauchly_table[mauchly_table$p < 0.05, ]
      if (nrow(sphericity_violations) > 0) {
        cat("\n球形性校正 (Greenhouse-Geisser):\n")
        print(anova_result$`Sphericity Corrections`)
      }
    }
    
    return(list(
      data = final_data,
      anova = anova_result,
      descriptive = descriptive_stats,
      erp_name = erp_name,
      measure = measure
    ))
    
  }, error = function(e) {
    cat("ANOVA失败:", e$message, "\n")
    cat("尝试简化分析...\n")
    
    return(run_simplified_analysis(final_data, measure, erp_name))
  })
}

# 简化分析（当完整ANOVA失败时）
run_simplified_analysis <- function(data, measure, erp_name) {
  cat("\n运行简化分析 - 主效应检验:\n")
  
  # 词类主效应
  if (measure == "amplitude") {
    word_class_data <- data %>%
      group_by(sample_name, word_class) %>%
      summarise(mean_val = mean(avg_amplitude, na.rm = TRUE), .groups = "drop")
    
    brain_region_data <- data %>%
      group_by(sample_name, brain_region) %>%
      summarise(mean_val = mean(avg_amplitude, na.rm = TRUE), .groups = "drop")
  } else {
    word_class_data <- data %>%
      group_by(sample_name, word_class) %>%
      summarise(mean_val = mean(avg_latency, na.rm = TRUE), .groups = "drop")
    
    brain_region_data <- data %>%
      group_by(sample_name, brain_region) %>%
      summarise(mean_val = mean(avg_latency, na.rm = TRUE), .groups = "drop")
  }
  
  # 词类效应
  cat("\n--- 词类主效应 ---\n")
  tryCatch({
    word_anova <- ezANOVA(
      data = word_class_data,
      dv = mean_val,
      wid = sample_name,
      within = word_class,
      detailed = TRUE
    )
    print(word_anova$ANOVA)
  }, error = function(e) {
    cat("词类效应分析失败:", e$message, "\n")
  })
  
  # 脑区效应
  cat("\n--- 脑区主效应 ---\n")
  tryCatch({
    region_anova <- ezANOVA(
      data = brain_region_data,
      dv = mean_val,
      wid = sample_name,
      within = brain_region,
      detailed = TRUE
    )
    print(region_anova$ANOVA)
  }, error = function(e) {
    cat("脑区效应分析失败:", e$message, "\n")
  })
  
  return(list(
    data = data,
    erp_name = erp_name,
    measure = measure,
    simplified = TRUE
  ))
}

# 可视化函数
create_erp_plots <- function(analysis_result) {
  if (is.null(analysis_result$data)) return(NULL)
  
  data <- analysis_result$data
  erp_name <- analysis_result$erp_name
  measure <- analysis_result$measure
  
  # 创建热图
  if (measure == "amplitude") {
    plot_data <- data %>%
      group_by(word_class, brain_region) %>%
      summarise(mean_val = mean(avg_amplitude, na.rm = TRUE), .groups = "drop")
    
    p1 <- ggplot(plot_data, aes(x = word_class, y = brain_region, fill = mean_val)) +
      geom_tile() +
      scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                           name = "amplitude(μV)") +
      theme_minimal() +
      labs(title = paste(erp_name, "amplitude"),
           x = "word", y = "region") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    plot_data <- data %>%
      group_by(word_class, brain_region) %>%
      summarise(mean_val = mean(avg_latency, na.rm = TRUE), .groups = "drop")
    
    p1 <- ggplot(plot_data, aes(x = word_class, y = brain_region, fill = mean_val)) +
      geom_tile() +
      scale_fill_gradient(low = "lightblue", high = "darkblue", 
                          name = "latent(ms)") +
      theme_minimal() +
      labs(title = paste(erp_name, "latent"),
           x = "word", y = "region") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
  }
  
  return(p1)
}

# 主分析函数
run_complete_erp_analysis <- function() {
  # 存储所有结果
  results <- list()
  
  # P200分析
  cat("开始P200分析...\n")
  results$p200_amplitude <- run_erp_analysis(data_p200_uv, "amplitude", "P200", p200_positions)
  results$p200_latency <- run_erp_analysis(data_p200_ms, "latency", "P200", p200_positions)
  
  # N400分析
  cat("\n开始N400分析...\n")
  results$n400_amplitude <- run_erp_analysis(data_n400_uv, "amplitude", "N400", n400_positions)
  results$n400_latency <- run_erp_analysis(data_n400_ms, "latency", "N400", n400_positions)
  
  # P600分析
  cat("\n开始P600分析...\n")
  results$p600_amplitude <- run_erp_analysis(data_p600_uv, "amplitude", "P600", p600_positions)
  results$p600_latency <- run_erp_analysis(data_p600_ms, "latency", "P600", p600_positions)
  
  return(results)
}

# 运行完整分析
cat("开始ERP脑区功能分区与词类关系分析\n")
cat(rep("=", 80), sep = "", "\n")
all_results <- run_complete_erp_analysis()

# 生成可视化
cat("\n生成可视化图表...\n")
plots <- list()
for (result_name in names(all_results)) {
  if (!is.null(all_results[[result_name]])) {
    plots[[result_name]] <- create_erp_plots(all_results[[result_name]])
  }
}

# 打印图表
for (plot_name in names(plots)) {
  if (!is.null(plots[[plot_name]])) {
    print(plots[[plot_name]])
  }
}


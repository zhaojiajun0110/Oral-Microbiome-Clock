# 载入所需库
library(openxlsx)
library(dplyr)
library(tidyverse) # 包含 ggplot2, tidyr, readr 等
library(VennDiagram) # 用于绘制Venn图

# --- 1. 数据导入与初始化 ---
message("--- 1. 数据导入与初始化 ---")
# 导入新测序数据
serum_metabolite_new <- read.xlsx("/Volumes/dot_office/朋友文件/代谢组试测检验/代谢试测/HCS_血液/FD20250700048_analysis-report/report/1总结报告/代谢物鉴定定量列表.xlsx",sheet = "全鉴定（raw）")%>%
  select(1:19) # 选取相关列

# 导入旧测序数据
serum_metabolite_old <- read.xlsx("/Volumes/dot_office/数据库/Webirth/Webirth/MWY-23-1037-02-a2_2023-10-24-16-20-30/1.Data_Assess/all_group/ALL_sample_data.xlsx")%>%
  select(1:12) # 选取相关列

message("--- 比较两个测序代谢物的差异 ---")

# --- 预处理：清理和标准化分子式 ---
# 确保分子式列名一致，并清理格式 (去除空白，转为大写)
serum_metabolite_new <- serum_metabolite_new %>%
  mutate(formula_clean = toupper(trimws(as.character(formula))))

serum_metabolite_old <- serum_metabolite_old %>%
  mutate(formula_clean = toupper(trimws(as.character(Formula))))

# 移除NA或空字符串的分子式
formulas_new_all <- serum_metabolite_new %>%
  filter(!is.na(formula_clean) & formula_clean != "") %>%
  pull(formula_clean) %>% unique()

formulas_old_all <- serum_metabolite_old %>%
  filter(!is.na(formula_clean) & formula_clean != "") %>% # 使用清理后的列
  pull(formula_clean) %>% unique()


# --- 2. 两次测序代谢物的overlap有多少 ---
message("\n--- 2. 代谢物重叠比较 (基于分子式) ---")

n_new <- length(formulas_new_all)
n_old <- length(formulas_old_all)
n_overlap <- length(intersect(formulas_new_all, formulas_old_all))
n_only_new <- n_new - n_overlap
n_only_old <- n_old - n_overlap

message(paste0("新测序总共鉴定到的独特分子式数量: ", n_new))
message(paste0("旧测序总共鉴定到的独特分子式数量: ", n_old))
message(paste0("两者共同鉴定的独特分子式数量 (重叠): ", n_overlap))
message(paste0("仅新测序鉴定的独特分子式数量: ", n_only_new))
message(paste0("仅旧测序鉴定的独特分子式数量: ", n_only_old))

# 可视化：Venn Diagram
# 设置Venn图颜色
my_colors <- c("#4daf4a", "#377eb8") # 新测序绿色，旧测序蓝色

# 绘制Venn图并保存
venn.diagram(
  x = list(formulas_new_all, formulas_old_all),
  category.names = c("New Sequencing", "Old Sequencing"),
  filename = "/Volumes/dot_office/朋友文件/代谢组试测检验/Metabolite_Overlap_Venn.png", # 保存到文件
  output = TRUE,
  
  # 图形输出设置
  imagetype = "png",
  height = 4500,
  width = 4500,
  resolution = 600,
  compression = "lzw",
  
  # 图表样式
  col = my_colors,
  fill = my_colors,
  alpha = 0.5,
  cex = 1.5,
  fontfamily = "sans",
  
  cat.cex = 1.5,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27), # 调整标签位置
  cat.dist = c(0.055, 0.055), # 调整标签与圆圈距离
  cat.fontfamily = "sans"
)
message("Venn图已保存为 Metabolite_Overlap_Venn.png")


# --- 3. level级别的差异 (针对重叠代谢物) ---
message("\n--- 3. 针对重叠代谢物的鉴定Level级别差异分析 ---")

# 3.1 提取重叠代谢物的级别信息
# 对于每个重复的分子式，选择其最高置信度的鉴定级别（即数字最小的Level）
overlap_new_levels_processed <- serum_metabolite_new %>%
  filter(formula_clean %in% formulas_old_all) %>%
  mutate(msi_numeric = as.numeric(gsub("MSI_level", "", MSI_levels))) %>% # 将"MSI_levelX"转换为数字X
  group_by(formula_clean) %>%
  summarise(New_MSI_Level = min(msi_numeric, na.rm = TRUE), .groups = 'drop') # 取最小的Level

overlap_old_levels_processed <- serum_metabolite_old %>%
  filter(formula_clean %in% formulas_new_all) %>%
  group_by(formula_clean) %>%
  summarise(Old_Level = min(as.numeric(Level), na.rm = TRUE), .groups = 'drop') # 取最小的Level

# 合并重叠代谢物的级别信息
overlap_comparison_df <- inner_join(
  overlap_new_levels_processed,
  overlap_old_levels_processed,
  by = "formula_clean"
)

# 3.2 统计Level差异并进行可视化
# 定义Level标签，方便绘图显示
level_labels_new_all <- c("MSI_level1", "MSI_level2")
level_labels_old_all <- c("Level 1", "Level 2", "Level 3")

# 统计旧Level和新Level的组合数量
level_diff_counts <- overlap_comparison_df %>%
  mutate(
    New_Level_Factor = factor(paste0("MSI_level", New_MSI_Level), levels = level_labels_new_all),
    Old_Level_Factor = factor(paste0("Level ", Old_Level), levels = level_labels_old_all)
  ) %>%
  count(Old_Level_Factor, New_Level_Factor) %>%
  # 填充缺失的组合为0，确保热图完整
  complete(Old_Level_Factor, New_Level_Factor, fill = list(n = 0))

message("\n重叠代谢物的鉴定Level交叉统计：")
print(level_diff_counts)

# 可视化：Level差异的 Heatmap
p_level_heatmap <- ggplot(level_diff_counts, aes(x = Old_Level_Factor, y = New_Level_Factor, fill = n)) +
  geom_tile(color = "white") + # 绘制方块
  geom_text(aes(label = n), color = "black", size = 4) + # 在方块上显示数量
  scale_fill_gradient(low = "white", high = "#377eb8", name = "代谢物数量") + # 颜色梯度
  labs(
    title = "重叠代谢物鉴定Level对比 (新测序 vs. 旧测序)",
    x = "旧测序鉴定级别 (Old Level)",
    y = "新测序鉴定级别 (New MSI Level)"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_level_heatmap)

# 评估Level提升情况
level_improvement_summary <- overlap_comparison_df %>%
  mutate(
    Improvement = case_when(
      New_MSI_Level < Old_Level ~ "鉴定提升 (Improved)", # 新Level数字更小，置信度更高
      New_MSI_Level == Old_Level ~ "鉴定不变 (Same)",
      New_MSI_Level > Old_Level ~ "鉴定下降 (Worsened)", # 新Level数字更大，置信度更低
      TRUE ~ "未知" # 以防万一
    )
  ) %>%
  count(Improvement) %>%
  mutate(Percentage = n / sum(n) * 100) # 计算百分比

message("\n重叠代谢物鉴定Level提升情况：")
print(level_improvement_summary)

# 可视化：Level提升情况的条形图
p_level_improvement <- ggplot(level_improvement_summary, aes(x = Improvement, y = n, fill = Improvement)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = paste0(n, " (", round(Percentage, 1), "%)")), vjust = -0.5, size = 3) +
  labs(
    title = "重叠代谢物鉴定Level变化",
    x = "鉴定级别变化",
    y = "代谢物数量"
  ) +
  scale_fill_manual(values = c("鉴定提升 (Improved)" = "#4daf4a", "鉴定不变 (Same)" = "#ff7f00", "鉴定下降 (Worsened)" = "#e41a1c")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 15, hjust = 0.5))

print(p_level_improvement)


# --- 4. 其他需要做的分析 ---
message("\n--- 4. 其他分析 ---")

# 4.1 整体鉴定Level分布 (包含非重叠的全部代谢物)
message("\n--- 4.1 整体鉴定Level分布 ---")

# 新测序数据的Level分布
new_overall_level_dist <- serum_metabolite_new %>%
  mutate(msi_numeric = as.numeric(gsub("MSI_level", "", MSI_levels))) %>%
  filter(!is.na(msi_numeric) & !is.na(formula_clean) & formula_clean != "") %>%
  group_by(formula_clean) %>%
  summarise(Level_Numeric = min(msi_numeric, na.rm = TRUE), .groups = 'drop') %>%
  count(Level_Numeric, name = "Count") %>%
  mutate(Source = "新测序 (New Sequencing)")

# 旧测序数据的Level分布
old_overall_level_dist <- serum_metabolite_old %>%
  filter(!is.na(Level) & !is.na(formula_clean) & formula_clean != "") %>%
  group_by(formula_clean) %>%
  summarise(Level_Numeric = min(as.numeric(Level), na.rm = TRUE), .groups = 'drop') %>%
  count(Level_Numeric, name = "Count") %>%
  mutate(Source = "旧测序 (Old Sequencing)")

# 合并并计算百分比
combined_overall_level_dist <- bind_rows(new_overall_level_dist, old_overall_level_dist) %>%
  group_by(Source) %>%
  mutate(Percentage = Count / sum(Count) * 100) %>%
  ungroup()

message("\n整体鉴定Level分布统计：")
print(combined_overall_level_dist)

# 可视化：整体鉴定Level分布
# 确保X轴显示所有Level (1-4)，并提供有意义的标签
all_possible_levels <- c(1, 2, 3, 4)
combined_overall_level_dist_full <- combined_overall_level_dist %>%
  complete(Source, Level_Numeric = all_possible_levels, fill = list(Count = 0, Percentage = 0)) %>%
  mutate(
    Level_Label = case_when(
      Level_Numeric == 1 ~ "最高 (1)",
      Level_Numeric == 2 ~ "中等 (2)",
      Level_Numeric == 3 ~ "较低 (3)",
      Level_Numeric == 4 ~ "最低/未知 (4)",
      TRUE ~ as.character(Level_Numeric)
    ),
    Level_Factor = factor(Level_Numeric, levels = all_possible_levels) # 确保排序
  )

p_overall_level_dist <- ggplot(combined_overall_level_dist_full, aes(x = Level_Factor, y = Percentage, fill = Source)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8)) +
  geom_text(aes(label = paste0(round(Percentage, 1), "%")), 
            position = position_dodge(width = 0.8), vjust = -0.5, size = 3) +
  labs(
    title = "不同测序平台整体鉴定级别分布",
    x = "鉴定级别 (数字越小置信度越高)",
    y = "代谢物百分比 (%)",
    fill = "数据来源"
  ) +
  scale_x_discrete(labels = setNames(combined_overall_level_dist_full$Level_Label, combined_overall_level_dist_full$Level_Factor)) +
  scale_fill_manual(values = c("新测序 (New Sequencing)" = "#4daf4a", "旧测序 (Old Sequencing)" = "#377eb8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

print(p_overall_level_dist)


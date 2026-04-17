# Loading libraries
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(rstatix)
library(moments)
library(AICcmodavg)
library(lubridate)
library(survival)
library(survminer)
library(contsurvplot)
library(riskRegression)
library(ggrepel)
library(gt)
library(greekLetters)

# Clearing environment
rm(list = ls())

# Loading data --------------------------------------------------------


# #list.files(path = default working directory, pattern = regular expression to match, ...)
# aging_files <- list.files(pattern = "*intens.csv",
#                          full.names = TRUE,
#                          ignore.case = TRUE)
# 
# rt_files <- list.files(pattern = "*rt.csv",
#                        full.names = TRUE,
#                        ignore.case = TRUE)
# 
# 
# #lapply applies read_csv over all elements of the list
# aging_data <- lapply(aging_files, read_csv)
# 
# rt_data <- lapply(rt_files, read_csv)
# 
# 
# # Merging Data
# merged_data <- aging_data[[1]]
# for (i in 2:length(aging_data)) {
#   merged_data <- left_join(merged_data, aging_data[[i]], by = c("Name", "Formula", "rt", "mz"))
# }
# 
# merged_rt_data <- rt_data[[1]]
# for (i in 2:length(rt_data)) {
#   merged_rt_data <- left_join(merged_rt_data, rt_data[[i]], by = c("Name", "Formula", "rt", "mz"))
# }


# Peak Picking --------------------------------------------------------


# # apply(data, 1 indicates that we are iterating over rows, function)
# # Calculation of tic median per metabolite (per row)
# merged_data$tic_median = apply(merged_data[, 5:ncol(merged_data)] |> select(!matches("[A-Za-z]")),
#                                1, median, na.rm=FALSE)
# 
# 
# 
# 
# # Calculation of blank to sample ratio per metabolite (per row)
# merged_data$blank_median = apply(merged_data |> select(matches("[Bb]...k.*")),
#                                      1, median, na.rm = FALSE)
# merged_data <- merged_data |> mutate(blank_sample_ratio = blank_median/tic_median)
# 
# 
# 
# 
# # Calculating variance and delta rt for each compounds rt
# # I exclude 0's in this calculation because 0 means not detected in this case
# merged_rt_data$rt_variance = apply(merged_rt_data[, 5:ncol(merged_rt_data)],
#                                    1, function(row) {var(row[row != 0])})
# merged_rt_data$avg_rt = apply(merged_rt_data[, 5:ncol(merged_rt_data)],
#                                     1, mean, na.rm = FALSE)
# merged_rt_data <- merged_rt_data |> mutate(avg_delta_rt = avg_rt - rt)
# merged_rt_data <- merged_rt_data |>  select(Name, Formula, rt, mz, avg_delta_rt, rt_variance)
# 
# 
# 
# 
# # Merging TIC and retention time data together
# merged_data <- left_join(merged_data, merged_rt_data, by = c("Name", "Formula", "rt", "mz"))
# 
# 
# 
# 
# # Visualization of tic_median, blank_sample_ratio, rt_variance metrics
# merged_data |>
#   ggplot() +
#   geom_histogram(aes(x = tic_median, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Meidan TIC Value",
#        y = "Density",
#        title = "Distribution of Median TIC values per metabolite identified") +
#   xlim(0, 1000000) +
#   ylim(0, 1.0e-05)
# 
# merged_data |>
#   ggplot() +
#   geom_histogram(aes(x = blank_sample_ratio, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Median Blank Intensity/Median Sample Intensity",
#        y = "Density",
#        title = "Distribution of Ratio of Metabolite to Blank TIC levels")
# 
# merged_data |>
#   ggplot() +
#   geom_histogram(aes(x = rt_variance, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Retention Time Variance",
#        y = "Density",
#        title = "Distribution of Retention Time Variance per Metabolite")
# 
# merged_data |>
#   ggplot() +
#   geom_histogram(aes(x = avg_delta_rt, y = after_stat(density)),
#                  color = 'white', fill = 'steelblue') +
#   labs(x = "Average Delta Retention Time",
#        y = "Density",
#        title = "Distribution of Delta Retention Time Across Compounds")
# 
# merged_data |>
#   ggplot() +
#   geom_histogram(aes(x = avg_delta_rt, y = after_stat(density)),
#                  color = 'white', fill = 'steelblue') +
#   labs(x = "Average Delta Retention Time",
#        y = "Density",
#        title = "Distribution of Delta Retention Time Across Compounds")
# 
# 
# # Filtering the data:
# #1. median_TIC > e4
# #2. rt_variance < 0.05
# #3. blank_sample_ratio < 0.5
# #4. blank_sample_ratio is not NaN or infinite
# #5. abs(avg_delta_rt) < 2
# #6. There are some duplicate compounds - if duplicate choose one with smallest delt_rt
# # The aim here is to not be overly stringent and end up with ~300 compounds
# filtered_data <- merged_data |>
#   filter(tic_median > 10000 & rt_variance < 0.05 & blank_sample_ratio < 0.5
#          & !is.nan(blank_sample_ratio) & abs(avg_delta_rt) < 1.5)
# 
# # Removing duplicated compounds
# filtered_data <- filtered_data |>
#   distinct(Name, .keep_all = TRUE)
# 
# 
# # Visualizing our metrics after filtering
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = tic_median, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Meidan TIC Value",
#        y = "Density",
#        title = "Distribution of Median TIC values per metabolite identified After Filtering") +
#   xlim(0, 1000000) +
#   ylim(0, 1.0e-05)
# 
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = blank_sample_ratio, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Median Blank Intensity/Median Sample Intensity",
#        y = "Density",
#        title = "Distribution of Ratio of Metabolite to Blank TIC levels After Filtering")
# 
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = rt_variance, y = after_stat(density)),
#                  fill = 'steelblue', color = 'white') +
#   labs(x = "Retention Time Variance",
#        y = "Density",
#        title = "Distribution of Retention Time Variance per Metabolite After Filtering")
# 
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = avg_delta_rt, y = after_stat(density)),
#                  color = 'white', fill = 'steelblue') +
#   labs(x = "Average Delta Retention Time",
#        y = "Density",
#        title = "Distribution of Delta Retention Time Across Compounds After Filtering")
# 
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = avg_delta_rt, y = after_stat(density)),
#                  color = 'white', fill = 'steelblue') +
#   labs(x = "Average Delta Retention Time",
#        y = "Density",
#        title = "Distribution of Delta Retention Time Across Compounds After Filtering")



# TIC Normalization --------------------------------------------------------



# # For each sample we want to calculate median TIC value
# # Then we want to normalize all the samples such that median TIC is the same for all of them
# # This will control for differences in concentration of sample etc
# TIC_norm <- merged_data |>
#   select(Name, Formula, !matches("[A-Za-z]")) |>
#   reshape2::melt() |>
#   group_by(variable) |>
#   summarise(tot_IC = sum(value))
# 
# 
# 
# # Initial total intensities per sample relative to the median
# TIC_norm |>
#   ggplot() +
#   geom_point(aes(x = variable, y = tot_IC)) +
#   geom_hline(aes(yintercept = median(TIC_norm$tot_IC)),
#              color = 'red') +
#   labs(x = "Sample",
#        y = "Total IC",
#        title = "Total Ion Current Per Sample")
# 
# 
# TIC_norm <- TIC_norm |>
#   mutate(factor = median(TIC_norm$tot_IC)/tot_IC) |>
#   mutate(new_TIC = tot_IC*factor)
# 
# 
# 
# 
# # After TIC normalization everything gets pushed towards the median
# # As we can see it is not just a line
# TIC_norm |>
#   ggplot() +
#   geom_point(aes(x = variable, y = new_TIC)) +
#   geom_hline(aes(yintercept = median(tot_IC)),
#              color = 'red') +
#   labs(x = "Sample",
#        y = "Total IC",
#        title = "Total Ion Current Per Sample After TIC Normalization")
# 
# 
# 
# 
# # Finally we need to multiply each compound for each sample by its respective factor
# TIC_norm <- TIC_norm |>
#   select(variable, factor)
# 
# filtered_data <- filtered_data |>
#   select(Name, Formula, !matches("[A-Za-z]")) |>
#   reshape2::melt()
# 
# filtered_data <- left_join(filtered_data, TIC_norm, by = 'variable')
# 
# filtered_data <- filtered_data |>
#   mutate(tic_norm_value = value*factor)
# 
# filtered_data <- filtered_data |>
#   rename(sample_id = variable) |>
#   rename(raw_tic = value)
# 
# 
# # Visualization of TIC normalized values
# filtered_data |>
#   ggplot() +
#   geom_histogram(aes(x = tic_norm_value, y = after_stat(density)),
#              fill = 'steelblue', color = 'white') +
#   labs(x = "TIC Normalized Values",
#        y = "Density",
#        title = "Distribution of TIC Normalized Values") +
#   xlim(0, 10000000)
# 
# 
# 
# # Writing the file out so that I don't have to do everything above every single time
# write_rds(filtered_data, file = "TIC Normalized Data.rds")



# Self Normalization -------------------------------------------------

# Excluding samples with processing errors and C7-F, whose numbering is a little messed up and needs to be checked 
TIC_norm <- read_rds("TIC Normalized Data.rds") |> 
  filter(!sample_id %in% c(644, 645, 744, 964,
                           394, 395, 1487:1538,
                           '218_20250410221533',
                           '454_20250310103343')) # Check file sizes for 218 and 454 to make sure you use largest file size

meta_data <- read_csv('Metadata.csv') |> 
  filter(Mouse != "serum blank") |> 
  mutate(sample_id = as.factor(sample_id))


# Merging metadata and TIC_norm data
merged_data <- left_join(TIC_norm, meta_data, by = 'sample_id')


# Inputting values for 0s using Half Minimum (HM) Method to prevent NaN/Inf values in downstream analysis
merged_data <- merged_data |> 
  group_by(Name) |> 
  mutate(half_minimum = min(tic_norm_value[tic_norm_value != 0])/2) |> 
  ungroup() |> 
  mutate(imputted_values = case_when(
    (tic_norm_value == 0) ~ half_minimum,
    (tic_norm_value != 0) ~ tic_norm_value
  ))



merged_data$Timepoint <- factor(merged_data$Timepoint,
                               levels = c("mo00", "mo03", "mo06", "mo09", "mo12"))

# Visualizing missing mo00 values
# merged_data |>
#   group_by(Mouse, Diet, Sex) |>
#   distinct(Timepoint) |>
#   arrange(Timepoint) |>
#   summarise(first_timepoint = first(Timepoint)) |>
#   ungroup() |>
#   ggplot(aes(x = first_timepoint, fill = Diet)) +
#   geom_histogram(stat = "count", position = "dodge", alpha = 0.8) +
#   facet_wrap(Sex ~ .) +
#   labs(x = "First Timepoint", y = "Mouse Count",
#        title = "First Timepoints of Mice in Study")


mice_to_keep <- (merged_data |> 
  group_by(Mouse) |> 
  distinct(Timepoint) |> 
  arrange(Timepoint) |> 
  summarise(first_timepoint = first(Timepoint)) |> 
  filter(first_timepoint == 'mo00') |> 
  ungroup())$Mouse

# Self normalization
merged_data_2 <- merged_data |> 
  filter(Mouse %in% mice_to_keep) |> 
  group_by(Name, Mouse) |> 
  arrange(Timepoint) |> 
  mutate(self_norm_value = imputted_values/first(imputted_values)) |> 
  ungroup()

# Now, only mice data set in the study whose first value is mo00
# merged_data |> 
# merged_data |>
#   group_by(Mouse, Diet, Sex) |>
#   distinct(Timepoint) |>
#   arrange(Timepoint) |>
#   summarise(first_timepoint = first(Timepoint)) |>
#   ungroup() |>
#   ggplot(aes(x = first_timepoint, fill = Diet)) +
#   geom_histogram(stat = "count", position = "dodge", alpha = 0.8) +
#   facet_wrap(Sex ~ .) +
#   labs(x = "First Timepoint", y = "Mouse Count",
#        title = "First Timepoints of Mice in Study")


# Some self normalized values still fall within the range of 10 - 1500
merged_data_2 |> 
  ggplot() +
  geom_histogram(aes(x = self_norm_value, y = after_stat(density)), fill = 'lightgreen', color = 'darkgreen', alpha = 0.8) +
  coord_cartesian(ylim = c(0, 0.00001))
  # stat_function(fun = dnorm, 
  #               args = list(mean = mean(log(merged_data$self_norm_value)), sd = sd(log(merged_data$self_norm_value))), 
  #               aes(color = "Normal distribution"))

# Dealing with outliers using chebyshev
# An alternate approach to try could be chernoff's bounds using exponential distribution
merged_data_2 <- merged_data_2 |>
  mutate(sd = sd(self_norm_value), mean = mean(self_norm_value)) |>
  mutate(k = (self_norm_value - mean)/sd, prob = 1/(k^2)) |>
  mutate(chebyshev_filter = case_when(
    (k >= 1 & prob <= 0.1) ~ (mean + sqrt(10)*sd),
    TRUE ~ self_norm_value
  ))

max(merged_data_2$chebyshev_filter)
merged_data_2 |> 
  ggplot() +
  geom_histogram(aes(x = chebyshev_filter, y = after_stat(density)), fill = 'lightgreen', color = 'darkgreen', alpha = 0.8) +
  coord_cartesian(ylim = c(0, 0.2))

merged_data_2$Timepoint <- (as.numeric(merged_data_2$Timepoint)*3) - 3

# Backward Feature Selection -------------------------------------------------------
# Backward selection is best practice based on ECE364 class because you maintain interacting terms

# Feature selection for control mice using BIC
# Results:  for the majority of metabolites, the Timepoint*Sex combination is the most useful model
feature_selection_cntrl <- function(dat) {
  cntrl_model <- lm(log2(chebyshev_filter) ~ as.numeric(Timepoint)*Sex, dat)
  cntrl_step <- step(cntrl_model, direction = "backward", k = log(length(dat)))
  return(toString(formula(cntrl_step)))
}

merged_data_2 |>
  filter(Diet == 'CON') |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, feature_selection_cntrl)) |>
  unnest_legacy(df) |>
  select(Name, df) |>
  count(df)

# Feature selection for diet using BIC
# Results: The as.numeric(Timepoint) * Diet seems to be the best model
feature_selection_diet <- function(dat) {
  start_model <- lm(log2(chebyshev_filter) ~ as.numeric(Timepoint)*Diet, dat)
  cntrl_step <- step(start_model, direction = "backward", k = log(length(dat)))
  return(toString(formula(cntrl_step)))
}

merged_data_2 |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, feature_selection_diet)) |>
  unnest_legacy(df) |>
  select(Name, df) |>
  count(df)

# Full factorial feature selection for diet and sex using BIC
# Results: the full factorial model does appear to have predictive power
feature_selection_full <- function(dat) {
  start_model <- lm(log2(chebyshev_filter) ~ as.numeric(Timepoint)*Diet + Sex + Sex:Diet + Sex:Diet:as.numeric(Timepoint), dat)
  cntrl_step <- step(start_model, direction = "backward", k = log(length(dat)))
  return(toString(formula(cntrl_step)))
}

merged_data_2 |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, feature_selection_full)) |>
  unnest_legacy(df) |>
  select(Name, df) |>
  count(df) |> 
  arrange(n)


# Age Effect ------------------------------------------------------------------

# Addition of cohort accounts for batch effect
age_fxn <- function(dat){
   lm(log2(chebyshev_filter) ~ Timepoint + as.factor(Cohort), data = dat) |> 
    broom::tidy() #broom:tidy puts output of lm() into a data frame
}


raw_stat <- merged_data_2 |> 
  filter(Diet == 'CON') |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, age_fxn)) |>
  unnest_legacy(df) |>
  filter(term == "Timepoint") |> 
  ungroup()


raw_stat <- raw_stat |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


age_sig_metabolites <- raw_stat |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH) |> 
  filter(Name != "Myristic acid")


write.csv(age_sig_metabolites, file = "Age Significant Metabolites.csv")


top_ten <- age_sig_metabolites |> 
  head(n = 10) |> 
  pull(Name)


# Visualizing significant compounds
filtered_data <- merged_data_2 |>
  mutate(Timepoint = Timepoint + 21) |>
  filter(Name %in% top_ten) |>
  filter(Diet == "CON")

filtered_data |> 
  ggplot() +
  geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              aes(x = Timepoint, y = log2(chebyshev_filter)), alpha = 0.5, size = 0.5) +
  geom_point(data = filtered_data |> filter(Timepoint == 21),
             aes(x = Timepoint, y = log2(chebyshev_filter)), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter)), method = "lm") +
  facet_wrap(~factor(Name, c("Perfluorooctanesulfonic acid","C24:4","C21:5","2-Hydroxyvaleric acid","2-Hydroxy-3-methylbutyric acid",
                             "Taurodeoxycholic acid", "2-hydroxy-b-methylvalerate","palmitoleic acid", "C16:1","Hydroxyisocaproic acid")) ~ ., nrow = 2, labeller = label_wrap_gen(width = 2), scales = "free_y") +
  labs(x = "Age (Months)", y = expression("log"[2]*"(Self Normalized Values)"), title = "Top 10 Metabolites Modulated by Age") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))


colnames(age_sig_metabolites)[colnames(age_sig_metabolites) == "estimate"] <- "\u0394ln(Self-Normalized Value) / \u0394Month"

age_sig_metabolites <- age_sig_metabolites |> 
  rename(`Standard Error` = std.error,
         `Term` = term,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH)

age_sig_metabolites |> 
  gt() |> 
  tab_header(
    title = "Metabolites that Vary Significantly Over Time"
  ) |> 
  fmt_number(
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )
  

merged_data_2  |>
  filter(Diet == "CON", Timepoint == "mo00") |> 
  distinct(Mouse) |> 
  count()

# Volcano Plot 

# Calculating log fold change and significance across each timepoint
raw_stat_volcano <- merged_data_2 |>
  filter(Diet == 'CON') |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, age_fxn)) |>
  unnest_legacy(df) |>
  filter(term == "Timepoint") |>
  rename(log_2_fold_change = estimate) # Here estimate is log(condition 1) - log(condition 2) = log(condition1/condition2) Not super intuitive but this is the log2fold change

# Doing p value adjustment
raw_stat_volcano <- raw_stat_volcano |>
  ungroup() |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))

# Creating volcano plot
raw_stat_volcano <- raw_stat_volcano |>
  mutate(Change = case_when(
    (log_2_fold_change < -0.14 & p_adj_BH < 0.05) ~ "Decrease",
    (log_2_fold_change > 0.14 & p_adj_BH < 0.05) ~ "Increase",
    TRUE ~ "No Significant Change"
  ))

# Significant values that I want to label on the graph
labels <- raw_stat_volcano |>
  filter(p_adj_BH < 0.05 & abs(log_2_fold_change) > 0.14) |>
  arrange(p_adj_BH)

labels <- labels$Name

raw_stat_volcano$labels <- ifelse(raw_stat_volcano$Name %in% labels,
                            raw_stat_volcano$Name, NA)

# Actually plotting the graph
raw_stat_volcano |>
  ggplot(mapping = aes(x = log_2_fold_change, y = -log10(p_adj_BH), col = Change, label = labels)) +
  geom_vline(xintercept = c(-0.14, 0.14), col = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = c(1.30102999566), col = 'gray', linetype = 'dashed') + # 1.3010 .. = -log10(0.05)
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "brown", "grey"),
                     labels = c("Decrease", "Increase", "No Significant Change")) +
  labs(color = 'Log Fold Change',
       x = expression("log"[2]*"(fold change)"),
       y = expression("-log"[10]*"p-value")) +
  ggtitle("Log Fold Change of Metabolites Across Timepoints") +
  geom_text_repel(max.overlaps = Inf, size = 4.4) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 20),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom")

# Making heaatmap

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(circlize)
library(RColorBrewer)

top_25 <- age_sig_metabolites |> 
  arrange(`Adjusted P-Value`) |> 
  head(n = 25) |> 
  pull(Name)
  
# Getting data in proper format and generating heat map
test <- merged_data_2 |>
  filter(Diet == "CON" & Name %in% top_25) |> 
  group_by(Timepoint, Name) |> 
  mutate(median_values = median(chebyshev_filter)) |> # Calculating median self_norm_value for each Timepoint, Sex, Diet, and Pathway combination
  ungroup() |> 
  distinct(Timepoint, Name, median_values) |> 
  mutate(log_fold = log2(median_values))

heatmap_format <- test |> 
  select(Name, Timepoint, log_fold) |> 
  pivot_wider(names_from = Timepoint, values_from = log_fold)

heatmap_format_2 <- heatmap_format[, -1]
rownames(heatmap_format_2) = heatmap_format$Name
heatmap_matrix <- as.matrix(heatmap_format_2)




# Creating Annotations for Heatmap
annotations <- as.data.frame(colnames(heatmap_matrix)) |> 
  mutate(`Age (Months)` = as.numeric(`colnames(heatmap_matrix)`) + 21) |> 
  select(`Age (Months)`)


# Creating and Plotting Heatmap 
age_heatmap <- Heatmap(
  
  heatmap_matrix,
  name = "Log Fold Change",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # Rows
  row_names_gp = gpar(fontsize = 12),
  row_title = "Compound Name",
  row_title_gp = gpar(fontsize = 14),
  row_title_side = "left",
  
  # Columns
  column_labels = annotations$`Age (Months)`,
  column_names_gp = gpar(fontsize = 14),
  column_title = "Age (Months)",
  column_title_gp = gpar(fontsize = 14),
  column_names_rot = 0,
  column_title_side = "bottom",
  
  # Legend
  heatmap_legend_param = list(direction = "horizontal", 
                              at = c(-4, -2, 0, 2, 4)),
  
  # Color
  col = circlize::colorRamp2(c(-4, 0, 4), c("blue","whitesmoke", "red"))
  
)

draw(age_heatmap, heatmap_legend_side = "top")


# Sex Effect -----------------------------------------------

sex_fxn <- function(dat) {
  lm(log2(chebyshev_filter) ~ Timepoint*Sex + as.factor(Cohort), data = dat) |>
    broom::tidy()
}
  

raw_stat_2 <- merged_data_2 |> 
  filter(Diet == 'CON') |>
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, sex_fxn)) |>
  unnest_legacy(df) |>
  filter(term %in% c('SexM','Timepoint:SexM')) |> 
  ungroup()


raw_stat_2 <- raw_stat_2 |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


sex_sig_metabolites <- raw_stat_2 |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH) |> 
  filter(Name != "Myristic acid",
         Name != "PE(22:6/16:0)") #Duplicates


write.csv(sex_sig_metabolites, file = "Sex Significant Metabolites.csv")

top_ten <- sex_sig_metabolites |> 
  head(n = 10) |> 
  pull(Name)


# # Visualizing significant compounds
filtered_data <- merged_data_2 |>
  mutate(Timepoint = Timepoint + 21) |>
  filter(Name %in% top_ten) |>
  filter(Diet == "CON")

sex_colors <- c("F" = "mediumpurple1", "M" = "dodgerblue")

filtered_data |> 
  ggplot() +
  geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), alpha = 0.5, size = 0.5) +
  geom_point(data = filtered_data |> filter(Timepoint == 21),
             aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), method = "lm") +
  facet_wrap(~factor(Name, c("C15:2","15-deoxy-delta-12-14-PGJ2","Perfluorooctanesulfonic acid","PE(16:0/22:6)","Capryloylglycine",
                             "PE(38:6)", "gamma-D-Glutamylglycine","PE(16:0/20:4)", "All-trans-retinoic acid","C14:2")) ~ ., nrow = 2, labeller = label_wrap_gen(width = 2), scales = "free_y") +
  labs(x = "Age (Months)", y = expression("log"[2]*"(Self Normalized Values)"), title = "Top 10 Metabolites Modulated by Sex Over Time") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 12),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 18)) +
  scale_color_manual(values = sex_colors) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))
  
  #coord_cartesian(ylim = c(0,10)) 


colnames(sex_sig_metabolites)[colnames(sex_sig_metabolites) == "estimate"] <- "\u0394Log(Self-Normalized Value) / \u0394Month"

sex_sig_metabolites <- sex_sig_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Term` = term)

sex_sig_metabolites |> 
  gt() |> 
  tab_header(
    title = "Metabolites With Sex By Time Interactions"
  ) |> 
  fmt_number(
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )

merged_data_2  |>
  filter(Diet == "CON", Timepoint == "mo00") |> 
  distinct(Mouse) |> 
  count()

top_25 <- sex_sig_metabolites |> 
  arrange(`Adjusted P-Value`) |> 
  head(n = 25) |> 
  pull(Name)

# Making a heatmap 
# Getting data in proper format and generating heat map
test <- merged_data_2 |>
  filter(Diet == "CON" & Name %in% top_25) |> 
  group_by(Timepoint, Sex, Name) |> 
  mutate(median_values = median(self_norm_value)) |> # Calculating median self_norm_value for each Timepoint, Sex, Diet, and Pathway combination
  ungroup() |> 
  distinct(Timepoint, Sex, Name, median_values) |> 
  mutate(log_fold = log2(median_values))


heatmap_format <- test |> 
  mutate(Time_Sex = paste(Timepoint,Sex)) |> 
  arrange(Sex, Timepoint) |> 
  select(Name, Time_Sex, log_fold) |> 
  pivot_wider(names_from = Time_Sex, values_from = log_fold)

heatmap_format_2 <- heatmap_format[, -1]
rownames(heatmap_format_2) = heatmap_format$Name
heatmap_matrix <- as.matrix(heatmap_format_2)




# Creating Annotations for Heatmap
annotations <- as.data.frame(colnames(heatmap_matrix)) |> 
  separate_wider_delim(`colnames(heatmap_matrix)`, 
                       delim = " ", names = c("Timepoint", "Sex"))

annotations <- annotations |> 
  mutate(`Age (Months)` = as.numeric(Timepoint) + 21)

col_annotation <- HeatmapAnnotation(
  Sex = annotations$Sex,
  col = list(Sex = c("M" = "deepskyblue3", "F" = "mediumorchid3")))



# Creating and Plotting Heatmap 
sex_heatmap <- Heatmap(
  
  heatmap_matrix,
  name = "Log Fold Change",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # Rows
  row_names_gp = gpar(fontsize = 12),
  row_title = "Compound Names",
  row_title_gp = gpar(fontsize = 14),
  row_title_side = "left",
  
  # Columns
  column_labels = annotations$`Age (Months)`,
  column_names_gp = gpar(fontsize = 12),
  column_title = "Age (Months)",
  column_title_gp = gpar(fontsize = 14),
  column_names_rot = 0,
  column_title_side = "bottom",
  
  # Annotations
  top_annotation = col_annotation,
  
)

draw(sex_heatmap)




# Age and Diet Effect ----------------------------------------



merged_data_2$Diet <- relevel(as.factor(merged_data_2$Diet), ref = "CON")

diet_fxn <- function(dat) {
  lm(log2(chebyshev_filter) ~ Timepoint*Diet + as.factor(Cohort), data = dat) |>
    broom::tidy()
}

raw_stat_3 <- merged_data_2 |> 
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, diet_fxn)) |>
  unnest_legacy(df) |>
  filter(grepl("Diet", term)) |> 
  ungroup()


raw_stat_3 <- raw_stat_3 |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


diet_sig_metabolites <- raw_stat_3 |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH) |> 
  filter(Name != "Myristic acid",
         Name != "PE(22:6/16:0)") #Duplicates


write.csv(diet_sig_metabolites, file = "Diet Significant Metabolites.csv")

custom_colors <- c("CON" = "black", "MetR" = "red", "PF" = "orange", "40CR" = "yellow4", 
                    "ZGN201" = "blue", "ZGN1062" = "green4")

# Visualizing Significant Compounds
filtered_data <- merged_data_2 |> 
  filter((Diet %in% c("CON", "PF") & Name == "C16:2") | 
           #(Name == "Citric acid" & Diet %in% c("ZGN201", "CON")) |
           (Name == "Citric acid" & Diet %in% c("ZGN1062", "CON")) |
           (Name == "LysoPE(22:4)"& Diet %in% c("MetR", "CON")) |
           (Name == "Taurodeoxycholic acid" & Diet %in% c("MetR", "CON")) |
           (Name == "LysoPE(18:1(9Z)/0:0)" & Diet %in% c("ZGN201", "CON")) |
           (Name == "N-Glycolylneuraminic acid" & Diet %in% c("ZGN1062", "CON")) |
           (Name == "C24:4" & Diet %in% c("MetR", "CON")) |
           (Name == "C22:4" & Diet %in% c("MetR", "CON")) |
           (Name == "PG(18:2/18:2)" & Diet %in% c("MetR", "CON")) |
           (Name == "C24:5" & Diet %in% c("MetR", "CON"))
         ) |> 
  mutate(Timepoint = Timepoint + 21)

filtered_data |> 
  ggplot() +
  geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet), alpha = 0.5, size = 0.5) +
  geom_point(data = filtered_data |> filter(Timepoint == 21),
              aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet), method = "lm", alpha = 0.5) +
  labs(x = "Age (Months)", y = expression("log"[2]*"(Self Normalized Values)"), title = "Top 10 Metabolites Modulated by Diet Over Time") +
  scale_color_manual(values = custom_colors) +
  facet_wrap(factor(Name, levels = c("C16:2", "Citric acid", "LysoPE(22:4)", "Taurodeoxycholic acid", "LysoPE(18:1(9Z)/0:0)",
                                     "N-Glycolylneuraminic acid", "C24:4", "C22:4", "PG(18:2/18:2)", "C24:5")) ~ ., nrow = 2, labeller = label_wrap_gen(width = 15), scales = "free_y") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 16)) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))
  

merged_data_2 |> 
  filter(Timepoint == "mo00") |> 
  group_by(Diet, Sex) |> 
  summarize(n_distinct(Mouse))
  

colnames(diet_sig_metabolites)[colnames(diet_sig_metabolites) == "estimate"] <- "\u0394Log(Self-Normalized Value) / \u0394Month"

diet_sig_metabolites <- diet_sig_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Term` = term)

diet_sig_metabolites |> 
  filter(Name != "3-[2-(2H-1-3-Benzodioxol-5-yl)ethyl]-5-oxo-2-5-dihydro-1-2-3-oxadiazol-3-ium") |> 
  gt() |> 
  tab_header(
    title = "Metabolites With Time By Diet Interactions"
  ) |> 
  fmt_number(
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )

merged_data_2  |>
  filter(Diet == "CON", Timepoint == "mo00") |> 
  distinct(Mouse) |> 
  count()


# Full Factorial Model --------------------------------------

full_fxn <- function(dat) {
  lm(log2(chebyshev_filter) ~ Timepoint*Diet + Sex + Sex:Timepoint + Sex:Diet + Sex:Diet:Timepoint + Cohort, data = dat) |>
    broom::tidy()
}

raw_stat_4 <- merged_data_2 |> 
  mutate(Diet = factor(Diet) |> relevel(ref = "CON"),
         Sex = factor(Sex) |> relevel(ref = "F")) |> # Ensure Sex reference is also set
  nest(stat_data = -Name) |>
  mutate(df = map(stat_data, full_fxn)) |>
  unnest_legacy(df) |> 
  filter(str_detect(term, "as.numeric(Timepoint):Diet.+:SexM|Diet.+:SexM")) |> 
  ungroup()

raw_stat_4 <- raw_stat_4 |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


full_factorial_metabolites <- raw_stat_4 |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH)


write.csv(full_factorial_metabolites, file = "Full Factorial Metabolites.csv")


# Visualizing significant compounds

custom_colors <- c("CON" = "black", "MetR" = "red", "PF" = "orange", "40CR" = "yellow4", 
                   "ZGN201" = "blue", "ZGN1062" = "green4")

filtered_data <- merged_data_2 |>
  filter((Diet %in% c("CON", "MetR") & Name == "gamma-D-Glutamylglycine") | 
           (Name == "DAG(15:1/22:4)" & Diet %in% c("MetR", "CON")) |
           (Name == "Glutamylserine" & Diet %in% c("MetR", "CON")) |
           (Name == "15-deoxy-delta-12-14-PGJ2" & Diet %in% c("PF", "CON")) |
           (Name == "Phenol sulphate" & Diet %in% c("ZGN201", "CON")) |
           (Name == "Perfluorooctanesulfonic acid" & Diet %in% c("MetR", "CON")) |
           (Name == "PE(22:6/16:0)" & Diet %in% c("ZGN201", "CON")) |
           (Name == "Propenoyl carnitine C3:1" & Diet %in% c("MetR", "CON")) |
           (Name == "C15:2" & Diet %in% c("ZGN1062", "CON")) |
           (Name == "PE(38:6)" & Diet %in% c("ZGN201", "CON")) 
         ) |> 
  mutate(Timepoint = Timepoint + 21)

filtered_data |> 
  ggplot() +
  #geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              #aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet, shape = Sex), alpha = 0.5, size = 0.3) +
  #geom_point(data = filtered_data |> filter(Timepoint == 21),
             #aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet, shape = Sex), alpha = 0.5, size = 0.3) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter), color = Diet, linetype = Sex), method = "lm", size = 1) +
  labs(x = "Age (Months)", y = expression("log"[2]*"(Self Normalized Values)"), title = "Top 10 Metabolites Modulated by Sex and Diet Over Time") +
  scale_color_manual(values = custom_colors) +
  facet_wrap(factor(Name, levels = c("gamma-D-Glutamylglycine", "DAG(15:1/22:4)", "Glutamylserine", "15-deoxy-delta-12-14-PGJ2", "Phenol sulphate", "Perfluorooctanesulfonic acid", "PE(22:6/16:0)", "Propenoyl carnitine C3:1", "C15:2", "PE(38:6)")) ~ ., nrow = 2, labeller = label_wrap_gen(width = 15), scales = "free_y") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 16),
          axis.title.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          strip.text = element_text(size = 12),
          axis.text.x = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 14),
          legend.position = "bottom",
          plot.title = element_text(size = 22)) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))

merged_data_2 |> 
  filter(Timepoint != 0)

# Drawing Table

colnames(full_factorial_metabolites)[colnames(full_factorial_metabolites) == "estimate"] <- "\u0394Log(Self-Normalized Value) / \u0394Month"

full_factorial_metabolites <- full_factorial_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Term` = term)

full_factorial_metabolites |> 
  head(n = 10) |> 
  gt() |> 
  tab_header(
    title = "Metabolites With Sex- and Diet- Specific Changes Over Time"
  ) |> 
  fmt_number(
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )
  


# Outcomes Analysis ------------------------------------------

# For outcome analysis, the goal here is to just use the TIC normalized values to 
# see if absolute values of metabolites correlate with time left to live. There may need 
# to be some normalization between batches. Here I can apply interesting ML techniques,
# decision tree, Bayesian, neural network.


# Looking at batch effect 
# Result: median of each Cohort is systematically lower 
merged_data |> 
  group_by(Cohort, Diet) |> 
  summarise(median = median(imputted_values)) |> 
  ggplot() +
  geom_col(aes(x = Cohort, y = median, fill = Diet), position = "dodge")



# Import of metadata and merge of metadata 
male_frailty <- read_csv("maleFrailty.csv")
female_frailty <- read_csv("femaleFrailty.csv")



# Joining male and female frailty data and removing rows where daysLeft data is not avaliable
frailty <- rbind(male_frailty, female_frailty) |> 
  rename(Timepoint = time, Mouse = ID) |> 
  filter(!is.na(daysLeft)) |> 
  select(Mouse, Timepoint, daysLeft)



# Removing mice that were SACed
wo_SAC <- merged_data |> filter(is.na(SAC))


# Joining frailty and metabolomics data 
outcomes_data <- inner_join(frailty, wo_SAC, by = c("Mouse", "Timepoint")) |> 
  select(Mouse, Sex, Cohort, Diet, Timepoint, daysLeft, Name, sample_id, imputted_values)


# Outcomes analysis
# Including Cohort as a separate variable to hopefully capture systematic difference between Cohorts
outcomes_fxn <- function(dat) {
  lm(log2(imputted_values) ~ daysLeft*Sex + as.factor(Cohort), data = dat) |> 
    broom::tidy()
}


raw_stat_5 <- outcomes_data |> 
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, outcomes_fxn)) |>
  unnest_legacy(df) |> 
  filter(term == "daysLeft") |> 
  ungroup()



raw_stat_5 <- raw_stat_5 |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))



raw_stat_5 <- raw_stat_5 |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH)


write.csv(raw_stat_5, file = "Outcomes.csv")


# Visualizing results
top_ten <- raw_stat_5 |>
  filter(term == 'daysLeft') |>
  arrange(p_adj_BH) |>
  head(n = 9) |>
  pull(Name)


# top_ten <- raw_stat_5 |> 
#   filter(term == 'SexM') |> 
#   arrange(p.value) |> 
#   head(n = 10) |> 
#   pull(Name)

filtered_data <- outcomes_data |> 
  filter(Name %in% top_ten) 

# filtered_data |>
#   filter(Name %in% top_ten) |>
#   ggplot() +
#   geom_jitter(aes(x = daysLeft, y = log2(imputted_values), color = Sex), size = 0.2, alpha = 0.5)  +
#   geom_smooth(aes(x = daysLeft, y = log2(imputted_values), color = Sex), method = "lm") +
#   facet_wrap(Name ~ ., nrow = 2) +
#   labs(title = "Top Ten Metabolites Correlated with Days Left to Live", x = "Days Left", y = "Log2 of TIC Normalized Metabolite Concentration")

filtered_data |> 
  ggplot() +
  geom_jitter(aes(x = daysLeft, y = log2(imputted_values)), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = daysLeft, y = log2(imputted_values)), method = "lm") +
  facet_wrap(factor(Name, levels = c("Perfluorooctanesulfonic acid", "C24:4", "Phosphoglycolic acid", "C24:2", "C20:5", "Neoabietic acid", "C21:3", "Acetylphosphate", "15-deoxy-delta-12-14-PGJ2")) ~ ., nrow = 3, labeller = label_wrap_gen(width = 2), scales = "free_y") +
  labs(x = "Days Left", y = expression("log"[2]*"(TIC Normalized Values)"), title = "Top 10 Metabolites With Absolute Levels Modulated by Sex Overtime") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 17))


# Drawing out table 

outcomes_metabolites <- raw_stat_5 |>
  filter(term == 'daysLeft') |>
  arrange(p_adj_BH)

colnames(outcomes_metabolites)[colnames(outcomes_metabolites) == "estimate"] <- "\u0394Log(Metabolite Abundance)/\u0394Month"

outcomes_metabolites <- outcomes_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Term` = term)

outcomes_metabolites |> 
  head(n = 9) |> 
  gt() |> 
  tab_header(
    title = "Metabolites Whose Abundance Correlates with Days Left to Live"
  ) |> 
  fmt_number(
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )

filtered_data |> 
  group_by(Name, Sex) |> 
  count()

# Deep Learning ---------------------------------------------

# Removing mice that were censored
fustat_filter <- survival |> 
  select(Mouse, fustat)
outcomes_data <- left_join(outcomes_data, fustat_filter, by = "Mouse") |> 
  filter(fustat == 1)

test <- meta_data |> 
  group_by(Cohort, Diet) |> 
  summarize(mice = n_distinct(Mouse)) |> 
  group_by(Cohort) |> 
  mutate(total = sum(mice)) |> 
  ungroup() |> 
  mutate(proportion = mice/total)



# Reformatting Data, and converting character strings to numeric values
long_format_outcomes <- outcomes_data |> 
  pivot_wider(id_cols = c(Mouse, Cohort, Sex, Diet, Timepoint, daysLeft, sample_id, fustat),
              names_from = Name, values_from = imputted_values) |> 
  mutate(Sex = case_when(
    Sex == "M" ~ 0, 
    Sex == "F" ~ 1
  )) |> 
  mutate(Control = case_when(Diet == "CON" ~ 1, TRUE ~ 0),
         MetR = case_when(Diet == "MetR" ~ 1, TRUE ~ 0),
         PF = case_when(Diet == "PF" ~ 1, TRUE ~ 0),
         ZGN1062 = case_when(Diet == "ZGN1062" ~ 1, TRUE ~ 0),
         ZGN201 = case_when(Diet == "ZGN201" ~ 1, TRUE ~ 0
  )) |> 
  mutate(Timepoint = case_when(
    Timepoint == "mo00" ~ 21,
    Timepoint == "mo03" ~ 24,
    Timepoint == "mo06" ~ 27,
    Timepoint == "mo09" ~ 30, 
    Timepoint == "mo12" ~ 33
  )) |> 
  mutate(
    `Cohort 1` = case_when(Cohort == "1" ~ 1, TRUE ~ 0),
    `Cohort 2` = case_when(Cohort == "2" ~ 1, TRUE ~ 0),
    `Cohort 3` = case_when(Cohort == "3" ~ 1, TRUE ~ 0),
    `Cohort 4` = case_when(Cohort == "4" ~ 1, TRUE ~ 0),
    `Cohort 5` = case_when(Cohort == "5" ~ 1, TRUE ~ 0),
    `Cohort 6` = case_when(Cohort == "6" ~ 1, TRUE ~ 0),
    `Cohort 7` = case_when(Cohort == "7" ~ 1, TRUE ~ 0),
    `Cohort 8` = case_when(Cohort == "8" ~ 1, TRUE ~ 0),
    `Cohort 9` = case_when(Cohort == "9" ~ 1, TRUE ~ 0),
    `Cohort C1-F` = case_when(Cohort == "1-F" ~ 1, TRUE ~ 0),
    `Cohort C2-F` = case_when(Cohort == "C2-F" ~ 1, TRUE ~ 0),
    `Cohort C3-F` = case_when(Cohort == "C3-F" ~ 1, TRUE ~ 0),
    `Cohort C4-F` = case_when(Cohort == "C4-F" ~ 1, TRUE ~ 0),
    `Cohort C5-F` = case_when(Cohort == "C5-F" ~ 1, TRUE ~ 0),
    `Cohort C6-F` = case_when(Cohort == "C6-F" ~ 1, TRUE ~ 0)
  )




# Looking at distribution of outcomes for binning
long_format_outcomes |> 
  ggplot() +
  geom_histogram(aes(x = daysLeft, y = after_stat(density)))


# Creating bins for daysLeft
long_format_outcomes <- long_format_outcomes |> 
  mutate(bins = ntile(daysLeft, n = 10)) |> 
  relocate(bins, .after = daysLeft)


write.csv(long_format_outcomes, file = "ML_data_2.csv")

long_format_outcomes |> 
  summarize(min(daysLeft))

724*0.01


# See Colab notebook for the rest of the ML stuff 


# Survival Analysis  --------------------------------------------------

# Reading in data
f_survival <- read.csv("SurvivalF.csv") |> 
  rename(Mouse = ID)
m_survival <- read.csv("SurvivalM.csv") |> 
  rename(Mouse = ID)
survival <- rbind(f_survival, m_survival)


# Only mouse, sex, and cohort metadata is necessary for baseline measurement comparison
small_meta_data <- meta_data |> 
  select(Mouse, Sex, Cohort) |> 
  distinct(Mouse, Sex, Cohort)


# Joining survival and metadata 
survival <- left_join(survival, small_meta_data, by = "Mouse") |> 
  filter(!is.na(Sex)) |> 
  mutate(Cohort = case_when(
    (Cohort == "1-F") ~ "C1-F",
    TRUE ~ Cohort
  ))


# Looking at diets relation to hazard 
# Results: ZGN1062, PF, and MetR appear to increase survival
surv <- Surv(survival$futime, survival$fustat)
fit <- survfit(surv ~ Tx, data = survival)  
ggsurvplot(fit, data = survival, pval = T, 
          risk.table = F, break.time.by = 100,
          risk.table.y.test.col = T, #ggtheme = survTheme,
          risk.table.fontsize = 3,
          font.tickslab = c(16),
          font.x = 16,
          font.y = 16,
          tables.theme = theme_survminer(font.tickslab = c(8)),
          xlab = "Age (days)", xlim = c(500, 1250),
          title = "Survival Curves by Diet",
          font.main = 20,
          font.legend = 14,
          legend.labs = c("Control", "MetR", "PF", "ZGN1062", "ZGN201"),
          palette = c("#000000","red", "orange", "darkgreen","blue"))


# Looking at Sex relation to hazard
# Results: Females live shorter time than males
surv <- Surv(survival$futime, survival$fustat)
fit <- survfit(surv ~ Sex, data = survival)
ggsurvplot(fit, data = survival, pval = T, 
           risk.table = F, break.time.by = 100,
           risk.table.y.test.col = T, #ggtheme = survTheme,
           font.tickslab = c(16),
           font.x = 16,
           font.y = 16,
           font.main = 20,
           font.legend = 16,
           xlab = "Age (days)", xlim = c(500, 1250),
           title = "Survival Curves by Sex",
           legend.labs = c("Female", "Male"),
           palette = c("mediumpurple1", "dodgerblue"))


# Looking at Cohort relation to hazard
# Results: There appears to be variation but not insane systematic variation
# The major confounding variable in terms of survival appears to be sex
surv <- Surv(survival$futime, survival$fustat)
fit <- survfit(surv ~ Cohort, data = survival)
ggsurvplot(fit, data = survival, pval = T, break.time.by = 100,
           risk.table.y.test.col = T, risk.table = T,
           #ggtheme = survTheme,
           xlab = "Age (days)", xlim = c(500, 1250),
           legend.labs = c("Cohort 1", "Cohort 2", "Cohort 3", "Cohort 4",
                           "Cohort 5", "Cohort 6", "Cohort 7", "Cohort 8", "Cohort 9", "Cohort 1-F",
                           "Cohort 2-F", "Cohort 3-F", "Cohort 4-F", "Cohort 5-F", "Cohort 6-F", "Cohort 7-F"),
           palette = c("deepskyblue", "dodgerblue", "darkturquoise","powderblue", "blue", "blue4", "skyblue",
                       "seagreen", "slateblue", "orange", "orange4", "brown", "coral", "tomato", "red", "gold"))




# Analysis of Males
baseline_measures_m <- merged_data |> 
  filter(Timepoint == "mo00", Sex == "M")


joined_cox_m <- left_join(baseline_measures_m, m_survival, by = "Mouse") |> 
  group_by(Name) |> 
  mutate(z_scaled = scale(imputted_values, center = T)) |> 
  filter(!is.na(futime))

# Function that extracts hazard ratio and p value for each metabolite
# Model extracts hazard ratio of metabolites controlling for confounding variables Cohort 
# Here I used z_scaled versions of metabolite levels such that hazard ratio number varied more from 1
cox_function <- function(dat) {
  model <- coxph(Surv(futime, fustat) ~ z_scaled + as.factor(Diet) + as.factor(Cohort), data=dat, x=TRUE)
  pvalue <- c(summary(model)$coefficients[1, 5])
  m_h_ratio <- c(exp(model$coefficients[1]))
  return(data.frame(m_h_ratio, pvalue))
}

# Applying function to each metabolite in the data set
raw_stat_6 <- joined_cox_m |> 
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, cox_function)) |> 
  unnest(df) |> 
  ungroup()


raw_stat_6 <- raw_stat_6 |>
  select(Name, m_h_ratio, pvalue) |> 
  mutate(p_adj_BH = p.adjust(pvalue, method = "BH"))

raw_stat_6 |> 
  filter(p_adj_BH <= 0.05) 

write.csv(raw_stat_6, file = "Survival_metabolites_male.csv")


raw_stat_6 <- raw_stat_6 |> 
  mutate(label = case_when(
    p_adj_BH <= 0.05 ~ Name,
    abs(1 - m_h_ratio) > 0.3 ~ Name,
    TRUE ~ NA
  ))
  

# Visualizing results
# Results: Hazard Ratios > 1 suggests that the hazard (in this case death) is more probable
# if the compound is elevated at baseline measurment 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5388384/
raw_stat_6 |> 
  ggplot(aes(m_h_ratio, -log10(p_adj_BH), label = label)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 8, min.segment.length = unit(0, 'lines')) +
  labs(x = "Hazard Ratio", y = expression("-Log"[10] * "(p-value)"),
       title = "Hazard Ratio Of Baseline Metabolite Levels in Males") +
  theme_minimal() +
       #caption = "Figure 1b. Results from Cox proportional hazard model on male mice, where the impact of baseline metabolite levels (standardized) on 
#survival was analyzed. Compounds with a p-value less than or equal to 0.1 are displayed. A hazard ratio > 1 indicates that
#the event of interest, in this case death, is more likely to occur when baseline levels of that metabolite are 
#elevated. A hazard ratio < 1 suggests a smaller risk of death if baseline levels of that metabolite are elevated.") 
  theme(plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")



# Visualizing survival curves of top 10 significant metabolites 
# PE(38:6), PG(18:1/22:6), PE(36:4), PE(38:4), PE(22:6/16:0), Linoleic acid, PE(36:2), 
# 2-hydroxy-b-methylvalerate, 2-hydroxy-b-methylvalerate, Hydroxyisocaproic acid
# Example
compound_specific_filter <- joined_cox_m |> 
  filter(Name == "PE(38:6)") |> 
  filter(!is.na(futime)) |> 
  mutate(quantile = as.factor(ntile(z_scaled, 4)))


fit <- survfit(Surv(futime, fustat) ~ quantile, data = compound_specific_filter)
ggsurvplot(fit, data = as.data.frame(compound_specific_filter), risk.table = F, pval = F,
           risk.table.fontsize = 3,
           font.tickslab = c(16),
           font.x = 16,
           font.y = 16,
           font.main = 20,
           font.legend = 14,
          legend.labs = c("1st Quantile", "2nd Quantile", "3rd Quantile", "4th Quantile"),
          xlab = "Age (days)", title = "Male Survival Curves\nStratified Along PE(38:6) Baseline Levels")



# Survival Analysis on Female Mice
baseline_measures_f <- merged_data |> 
  filter(Timepoint == "mo00", Sex == "F")


joined_cox_f <- left_join(baseline_measures_f, f_survival, by = "Mouse") |> 
  group_by(Name) |> 
  mutate(z_scaled = scale(imputted_values, center = T)) |> 
  filter(!is.na(futime))


# Cox_function modified for females
cox_function <- function(dat) {
  model <- coxph(Surv(futime, fustat) ~ z_scaled + as.factor(Diet) + as.factor(Cohort), data=dat, x=TRUE)
  pvalue <- c(summary(model)$coefficients[1, 5])
  f_h_ratio <- c(exp(model$coefficients[1]))
  return(data.frame(f_h_ratio, pvalue))
}


# Applying function to each metabolite in the data set for females
raw_stat_7 <- joined_cox_f |> 
  nest(stat_data = -Name) |>
  mutate(df = purrr::map(stat_data, cox_function)) |> 
  unnest(df) |> 
  ungroup()


raw_stat_7 <- raw_stat_7 |>
  select(Name, f_h_ratio, pvalue) |> 
  filter(Name != "Myristic acid") |> # Two myristic acid's
  mutate(p_adj_BH = p.adjust(pvalue, method = "BH"))


write.csv(raw_stat_7, file = "Survival_metabolites_female.csv")


# Making compound labels
raw_stat_7 <- raw_stat_7 |> 
  mutate(label = case_when(
    p_adj_BH <= 0.54 ~ Name,
    TRUE ~ NA
  ))


# Visualizing compounds associated with hazard in females
raw_stat_7 |> 
  ggplot(aes(f_h_ratio, -log10(p_adj_BH), label = label)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 14, min.segment.length = unit(0, 'lines')) +
  labs(x = "Hazard Ratio", y = expression("-Log"[10] * "(p-value)"),
       title = "Hazard Ratio Of Baseline Metabolite Levels in Females") +
  theme(plot.caption = element_text(hjust = 0)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey")


# Looking at survival curves of top 10 significant metabolites 
# C26:4, C14:0(Myristic acid), C22:4, PE(18:1/22:6), Indole-3-carboxylic acid,  
# C24:2, Carnosine, D-Ribose 5-phosphate, C24:4, C24:1
# Example
compound_specific_filter <- joined_cox_f |> 
  filter(Name == "C26:4") |> 
  filter(!is.na(futime)) |> 
  mutate(quantile = as.factor(ntile(z_scaled, 4)))


fit <- survfit(Surv(futime, fustat) ~ quantile, data = compound_specific_filter)
ggsurvplot(fit, data = as.data.frame(compound_specific_filter), risk.table = F, pval = F,
           risk.table.fontsize = 3,
           font.tickslab = c(16),
           font.x = 16,
           font.y = 16,
           font.main = 20,
           font.legend = 14,
           legend.labs = c("1st Quantile", "2nd Quantile", "3rd Quantile", "4th Quantile"),
           xlab = "Age (days)", title = "Female Survival Curves\nStratified Along C26:4 Baseline Levels")




# Comparing hazard ratio relative to beta value (slope) from linear models
slope_values <- raw_stat_2 |> 
  select(Name, term, estimate) |> 
  pivot_wider(names_from = term, values_from = estimate) |> 
  rename(f_slope = Timepoint) |> 
  mutate(m_slope = f_slope + `Timepoint:SexM`) |> 
  select(Name, f_slope, m_slope)


female_hazard <- raw_stat_7 |> 
  select(Name, f_h_ratio)


male_hazard <- raw_stat_6 |> 
  select(Name, m_h_ratio)


hazard_beta_merge <- left_join(left_join(slope_values, female_hazard, by = "Name"), 
                               male_hazard, by = "Name")


# Visualizing results for males
hazard_beta_merge |> 
  ggplot(aes(x = m_slope, y = m_h_ratio, label = Name)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 10, min.segment.length = unit(0, 'lines')) +
  geom_vline(xintercept = 0, alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, alpha = 0.4, linetype = "dashed") +
  labs(x = "Slope", y = "Hazard Ratio", title = "Hazard Ratio vs. Slope of Linear Regression Model\nIn Male Mice") +
  theme_minimal()+
  stat_cor()+
  theme(axis.title.y = element_text(size = 16),
                axis.title.x = element_text(size = 16),
                axis.text.y = element_text(size = 16),
                axis.text.x = element_text(size = 16),
                plot.title = element_text(size = 20))


# Visualizing results for females
hazard_beta_merge |> 
  ggplot(aes(x = f_slope, y = f_h_ratio, label = Name)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 8, min.segment.length = unit(0, 'lines')) +
  geom_vline(xintercept = 0, alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, alpha = 0.4, linetype = "dashed") +
  labs(x = "Slope", y = "Hazard Ratio", title = "Hazard Ratio vs. Slope of Linear Regression Model\nIn Female Mice") +
  theme_minimal()+
  stat_cor() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20))


# Pathway Analysis I ------------------------------------------------------

# Installing necessary packages
library(stringr)
library(dplyr)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("KEGGREST")
library(KEGGREST)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(circlize)
library(RColorBrewer)


# Loading compound database
c_database <- read.csv("Knowns Compound Database Mapping.csv")



# Finding compounds that actually have KEGG ID's and extracting compound pathways from those
pathway_merge <- left_join(merged_data_2, c_database, by = "Name") |> 
  filter(!is.na(KEGG) & KEGG != "")

pathways <- pathway_merge |> 
  distinct(KEGG) |> 
  rowwise() |> 
  mutate(pathways = list(keggGet(KEGG)[[1]]$PATHWAY)) |> 
  ungroup()

pathway_merge <- left_join(pathway_merge, pathways, by = "KEGG") 



# Un-nesting pathways
pathway_merge <- pathway_merge |> 
  tidyr::unnest(pathways)

# Counting pathways and filtering out those without many compounds within them
pathway_merge_2 <- pathway_merge |> 
  add_count(pathways, name = "p_count")

# Pathways with less than four compounds within it are excluded from analysis
pathway_merge_2 |> 
  filter(p_count <= 4132) |> 
  group_by(pathways) |> 
  summarise(n_distinct(Name)) |> 
  print(n = 143)
  
pathway_merge_2 <- pathway_merge_2 |> 
  filter(p_count > 4132) |> 
  filter(Diet != "40CR") |> 
  filter(!(pathways %in% c("Carbon fixation by Calvin cycle", "Biosynthesis of plant hormones","Biosynthesis of plant secondary metabolites","Biosynthesis of alkaloids derived from ornithine, lysine and nicotinic acid",
          "Biosynthesis of alkaloids derived from shikimate pathway", "Biosynthesis of various antibiotics","Biosynthesis of phenylpropanoids","Biosynthesis of various other secondary metabolites",
          "Metabolic pathways", "Microbial metabolism in diverse environments","Diterpenoid biosynthesis",
          "Leishmaniasis", "Monobactam biosynthesis", "Other carbon fixation pathways", "Sulfoquinovose metabolism",
          "Biosynthesis of alkaloids derived from terpenoid and polyketide", "Carbon metabolism",
          "Vancomycin resistance", "Two-component system", "Biosynthesis of various siderophores",
          "Galactose metabolism")))

pathway

# Getting data in proper format and generating heat map
raw_stat_pathways <- pathway_merge_2 |>
  nest(data = -c(Sex, Diet, pathways)) |>
  mutate(df = purrr::map(data, age_fxn)) |> 
  unnest_legacy(df) |> 
  filter(term == "Timepoint")
  
heatmap_format <- raw_stat_pathways |> 
  mutate(Sex_Diet = paste(Sex,Diet)) |> 
  arrange(Diet, Sex) |> 
  select(pathways, Sex_Diet, estimate) |> 
  pivot_wider(names_from = Sex_Diet, values_from = estimate)

heatmap_format_2 <- heatmap_format[, -1]
rownames(heatmap_format_2) = heatmap_format$pathways
heatmap_matrix <- as.matrix(heatmap_format_2)




# Creating Annotations for Heatmap
annotations <- as.data.frame(colnames(heatmap_matrix)) |> 
  separate_wider_delim(`colnames(heatmap_matrix)`, 
                       delim = " ", names = c("Sex", "Diet"))

col_annotation <- HeatmapAnnotation(
  Diet = annotations$Diet,
  Sex = annotations$Sex,
  col = list(Sex = c("M" = "deepskyblue3", "F" = "mediumorchid3"),
             Diet = c("CON" = "black", "MetR" = "darkred",  "ZGN1062" = "forestgreen", "ZGN201" = "blue3", "PF" = "orange2")),
  annotation_legend_param = list(
    Diet = list(direction = "horizontal", nrow = 1),
    Sex = list(direction = "horizontal", nrow = 1)
  )
  )




# Creating and Plotting Heatmap 
pathway_heatmap <- Heatmap(
  
  heatmap_matrix,
  name = "\u0394Log Fold Change/Time",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # Rows
  row_names_gp = gpar(fontsize = 12),
  row_title = "Metabolic Pathways",
  row_title_gp = gpar(fontsize = 12),
  row_title_side = "left",
  
  # Columns
  show_column_names = F,
  column_title = "\u0394Log Fold Change/Time",
  column_title_gp = gpar(fontsize = 12),
  column_names_rot = 0,
  column_title_side = "bottom",
  
  # Annotations
  top_annotation = col_annotation,
  column_split = c(annotations$Diet),
  
  # Legend
  heatmap_legend_param = list(direction = "horizontal"),
  
  # Color
  col = circlize::colorRamp2(c(-0.1, 0, 0.1), c("blue","whitesmoke", "red"))
  
  )
  

draw(pathway_heatmap, heatmap_legend_side = "top", annotation_legend_side = "top", padding = unit(c(2, 2, 2, 100), "mm"))

# Summarizing metabolites in each pathway 
mets_in_paths <- pathway_merge_2 |> 
  select(pathways, Name) |> 
  group_by(pathways) |> 
  summarize(Metabolites = paste(unique(Name), collapse = ", ")) |> 
  ungroup() |> 
  rename(Pathways = pathways)

write.csv(mets_in_paths, file = "Metabolites In Each Pathway.csv")

mets_in_paths |> 
  gt() |> 
  tab_header(
    title = "Metabolites in Each Pathway"
  )


# Pathway Analysis II --------------------------------------------------------

# Glyoxylate and dicarboxylate metabolism in males compared in CON vs. ZGN 1062
# Here we are attempting to answer the question - why do the ZGN 1062 mice live so much longer
# Specifically we are looking at mo06 mice for both groups

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pathview")
library(pathview)

test <- pathway_merge |> 
  filter(Diet == "ZGN1062" & Sex == "F" & pathways == "Histidine metabolism" ) |> 
  nest(data = -Name) |> 
  mutate(df = purrr::map(data, age_fxn)) |> 
  unnest_legacy(df) |> 
  filter(term == "Timepoint") |> 
  rename(ZGN1062 = estimate) |> 
  select(Name, ZGN1062)
  
test2 <- pathway_merge |> 
  filter(Diet == "CON" & Sex == "F" & pathways == "Histidine metabolism" ) |> 
  nest(data = -Name) |> 
  mutate(df = purrr::map(data, age_fxn)) |> 
  unnest(df) |> 
  filter(term == "Timepoint") |> 
  rename(Control = estimate) |> 
  select(Name, Control)

big_test <- inner_join(test, test2, by = "Name")
big_test <- left_join(big_test, c_database, by = "Name") 
big_test <- big_test |> 
  select(KEGG, Control, ZGN1062)
  

matrix_format <- big_test[, -1]
rownames(matrix_format) = big_test$KEGG
matrix <- as.matrix(matrix_format)
colnames(matrix) <- c("Control", "ZGN1062")



pv.out <- pathview(
  cpd.data = matrix,
  pathway.id = "mmu00340",
  species = "mmu",        
  cpd.idtype = "kegg",    
  out.suffix = "metabolomics_plot",
  multi.state = T,
  same.layer = F,
  limit = 0.1
)



# Physiology Section ----

female_frailty <- read_csv("femaleFrailty.csv")
male_frailty <- read_csv("maleFrailty.csv")
frailty <- rbind(male_frailty, female_frailty) |> 
  mutate(Diet = case_when(
    Tx == "Con" ~ "Control",
    Tx == "Zgn201" ~ "ZGN201",
    TRUE ~ Tx
  )) |> 
  filter(Diet != "40CR")

frailty_summary <- frailty |> 
  filter(Diet != "40CR") |> 
  group_by(Diet, time) |> 
  summarise(median_frailty = median(TotalScore, na.rm = TRUE),
            median_bweight = median(BodyWeightRaw, na.rm = TRUE))

colors <- c("Control" = "black", "MetR" = "red", "PF" = "orange", "ZGN1062" = "forestgreen", 
            "ZGN201" = "dodgerblue")
ggplot() +
  geom_jitter(data = frailty, mapping = aes(x = as.numeric(str_sub(time, 3, 4)) + 21, y = TotalScore, color = Diet), alpha = 0.15, size = 1) +
  geom_line(data = frailty_summary, mapping = aes(x = as.numeric(str_sub(time, 3, 4)) + 21, y = median_frailty, color = Diet), size = 1) +
  labs(x = "Age (Months)", y = "Median Frailty Score", title = "Frailty Score Trajectories By Diet") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  scale_colour_manual(values = colors)


ggplot() +
  geom_jitter(data = frailty, mapping = aes(x = as.numeric(str_sub(time, 3, 4)) + 21, y = BodyWeightRaw, color = Diet), alpha = 0.15, size = 1) +
  geom_line(data = frailty_summary, mapping = aes(x = as.numeric(str_sub(time, 3, 4)) + 21, y = median_bweight, color = Diet), size = 1) +
  labs(x = "Age (Months)", y = "Median Body Weight", title = "Raw Body Weight Trajectories By Diet") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        strip.text = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        legend.text = element_text(size = 16),
        legend.title = element_text(size = 16),
        plot.title = element_text(size = 18)) +
  coord_cartesian(xlim = c(20, 35)) +
  scale_colour_manual(values = colors)

frailty |> 
  filter(time == "mo00") |> 
  group_by(Sex) |> 
  count()


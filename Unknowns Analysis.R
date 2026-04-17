# Unknowns Compound Analysis

library(tidyverse)
library(ggrepel)
library(survminer)
library(gt)


# Loading Data ------

#list.files(path = default working directory, pattern = regular expression to match, ...)
aging_files <- list.files(pattern = "*intens.csv",
                         full.names = TRUE,
                         ignore.case = TRUE)


#lapply applies read_csv over all elements of the list
aging_data <- lapply(aging_files, read_csv)


# Merging Data
merged_data <- aging_data[[1]]
for (i in 2:length(aging_data)) {
  merged_data <- left_join(merged_data, aging_data[[i]] |> 
                             select(-c(id, log10_inten, class, ppm_error, score_known_rt, 
                                       score_MS2, score_MS2, level, accession, compound, isFrag)), 
                           by = c("peak_id", "formula", "rt", "mz"))
}


# Peak Picking ------

# apply(data, 1 indicates that we are iterating over rows, function)
# Calculation of tic median per metabolite (per row)
merged_data$tic_median = apply(merged_data[, 5:ncol(merged_data)] |> select(!matches("[A-Za-z]")),
                               1, median, na.rm=FALSE)




# Calculation of blank to sample ratio per metabolite (per row)
merged_data$blank_median = apply(merged_data |> select(matches("[Bb]...k.*")),
                                     1, median, na.rm = FALSE)
merged_data <- merged_data |> mutate(blank_sample_ratio = blank_median/tic_median)



# Visualization of tic_median, blank_sample_ratio
merged_data |>
  ggplot() +
  geom_histogram(aes(x = tic_median, y = after_stat(density)),
                 fill = 'steelblue', color = 'white') +
  labs(x = "Meidan TIC Value",
       y = "Density",
       title = "Distribution of Median TIC values per metabolite identified") +
  xlim(0, 1000000) +
  ylim(0, 1.0e-05)

merged_data |>
  ggplot() +
  geom_histogram(aes(x = blank_sample_ratio, y = after_stat(density)),
                 fill = 'steelblue', color = 'white') +
  labs(x = "Median Blank Intensity/Median Sample Intensity",
       y = "Density",
       title = "Distribution of Ratio of Metabolite to Blank TIC levels")


# Filtering the data:
#1. median_TIC > e4
#2. blank_sample_ratio < 0.5
#3. blank_sample_ratio is not NaN or infinite
#4. There are some duplicate compounds - if duplicate choose one with smallest delt_rt
# The aim here is to not be overly stringent and end up with ~300 compounds
filtered_data <- merged_data |>
  filter(tic_median > 10000 & blank_sample_ratio < 0.5 & !is.nan(blank_sample_ratio))


# Removing duplicated compounds
filtered_data <- filtered_data |>
  distinct(id, .keep_all = TRUE)


# Visualizing our metrics after filtering
filtered_data |>
  ggplot() +
  geom_histogram(aes(x = tic_median, y = after_stat(density)),
                 fill = 'steelblue', color = 'white') +
  labs(x = "Meidan TIC Value",
       y = "Density",
       title = "Distribution of Median TIC values per metabolite identified After Filtering") +
  xlim(0, 1000000) +
  ylim(0, 1.0e-05)

filtered_data |>
  ggplot() +
  geom_histogram(aes(x = blank_sample_ratio, y = after_stat(density)),
                 fill = 'steelblue', color = 'white') +
  labs(x = "Median Blank Intensity/Median Sample Intensity",
       y = "Density",
       title = "Distribution of Ratio of Metabolite to Blank TIC levels After Filtering")


# TIC Normalization ------ 

# For each sample we want to calculate median TIC value
# Then we want to normalize all the samples such that median TIC is the same for all of them
# This will control for differences in concentration of sample etc
TIC_norm <- filtered_data |>
  select(peak_id, !matches("[A-Za-z]")) |>
  reshape2::melt(id.vars = "peak_id", variable.name = "Sample", value.name = "TIC") |> 
  group_by(Sample) |>
  summarise(tot_IC = sum(TIC))


# Initial total intensities per sample relative to the median
# Many high outliers
TIC_norm |>
  ggplot() +
  geom_point(aes(x = Sample, y = tot_IC)) +
  geom_hline(aes(yintercept = median(TIC_norm$tot_IC)),
             color = 'red') +
  labs(x = "Sample",
       y = "Total IC",
       title = "Total Ion Current Per Sample")


TIC_norm <- TIC_norm |>
  mutate(factor = median(TIC_norm$tot_IC)/tot_IC) |>
  mutate(new_TIC = tot_IC*factor)




# After TIC normalization everything gets pushed towards the median
# As we can see it is not just a line
TIC_norm |>
  ggplot() +
  geom_point(aes(x = Sample, y = new_TIC)) +
  geom_hline(aes(yintercept = median(tot_IC)),
             color = 'red') +
  labs(x = "Sample",
       y = "Total IC",
       title = "Total Ion Current Per Sample After TIC Normalization")




# Finally we need to multiply each compound for each sample by its respective factor
TIC_norm <- TIC_norm |>
  select(Sample, factor)

filtered_data <- filtered_data |>
  select(peak_id, !matches("[A-Za-z]")) |>
  reshape2::melt(id.vars = "peak_id", variable.name = "Sample", value.name = "TIC")

filtered_data <- left_join(filtered_data, TIC_norm, by = 'Sample')

filtered_data <- filtered_data |>
  mutate(tic_norm_value = TIC*factor)


# Visualization of TIC normalized values
filtered_data |>
  ggplot() +
  geom_histogram(aes(x = tic_norm_value, y = after_stat(count)),
             fill = 'steelblue', color = 'white', bins = 10) +
  labs(x = "TIC Normalized Values",
       y = "Density",
       title = "Distribution of TIC Normalized Values") +
  coord_cartesian(ylim = c(0, 10000))
  


# Writing the file out so that I don't have to do everything above every single time
write_rds(filtered_data, file = "TIC Normalized Data.rds")


# Self Normalization -----

# Excluding samples with processing errors and C7-F, whose numbering is a little messed up and needs to be checked 
TIC_norm <- read_rds("TIC Normalized Data.rds") |> 
  filter(!Sample %in% c(644, 645, 744, 964,
                           394, 395, 1487:1538,
                           '218_20250410221533',
                           '454_20250310103343')) # Check file sizes for 218 and 454 to make sure you use largest file size

meta_data <- read_csv('Metadata.csv') |> 
  filter(Mouse != "serum blank") |> 
  mutate(Sample = as.factor(sample_id))


# Merging metadata and TIC_norm data
merged_data <- left_join(TIC_norm, meta_data, by = 'Sample')

merged_data |> 
  distinct(peak_id) |> 
  count()

# Inputting values for 0s using Half Minimum (HM) Method to prevent NaN/Inf values in downstream analysis
merged_data <- merged_data |> 
  group_by(peak_id) |> 
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
  group_by(peak_id, Mouse) |> 
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
  geom_histogram(aes(x = self_norm_value, y = after_stat(count)), fill = 'lightgreen', color = 'darkgreen', alpha = 0.8) +
  coord_cartesian(ylim = c(0, 100))
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

merged_data_2 |> 
  ggplot() +
  geom_histogram(aes(x = chebyshev_filter, y = after_stat(count)), fill = 'lightgreen', color = 'darkgreen', alpha = 0.8) +
  coord_cartesian(ylim = c(0, 800)) # Not the best result but like still okay 

 
merged_data_2 <- merged_data_2 |> 
  mutate(Timepoint = as.numeric(str_sub(Timepoint, 3, 4)))

# Age Effect ----- 


library(ggrepel)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)

library(circlize)
library(RColorBrewer)

# Addition of cohort accounts for batch effect
age_fxn <- function(dat){
  lm(log2(chebyshev_filter) ~ Timepoint + as.factor(Cohort), data = dat) |> 
    broom::tidy() #broom:tidy puts output of lm() into a data frame
}


raw_stat <- merged_data_2 |> 
  filter(Diet == 'CON') |>
  nest(stat_data = -peak_id) |>
  mutate(df = purrr::map(stat_data, age_fxn)) |>
  unnest_legacy(df) |>
  filter(term == "Timepoint") |> 
  ungroup()


raw_stat <- raw_stat |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


unknowns_age_sig_metabolites <- raw_stat |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH)


write.csv(unknowns_age_sig_metabolites, file = "Unknowns Age Significant Metabolites.csv")


top_ten <- unknowns_age_sig_metabolites |> 
  head(n = 10) |> 
  pull(peak_id)

filtered_data <- merged_data_2 |> 
  mutate(Timepoint = Timepoint + 21) |> 
  filter(Diet == "CON") |> 
  filter(peak_id %in% top_ten)

## Visualizing significant compounds
filtered_data |> 
  ggplot() +
  geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              aes(x = Timepoint, y = log2(chebyshev_filter)), alpha = 0.5, size = 0.5) +
  geom_point(data = filtered_data |> filter(Timepoint == 21),
             aes(x = Timepoint, y = log2(chebyshev_filter)), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter)), method = "lm") +
  facet_wrap(factor(peak_id, levels = c("2935", "3877", "399", "3887", "3878", "3883", "2706", "955", "5193", "3892")) ~ ., nrow = 2, scales = "free_y") +
  labs(x = "Age (Months)", y = expression("log"[2]*" (Self Normalized Values)"), title = "Top 10 Most Significant Unknown Metabolites Associated with Age") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 13),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom",
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))


# Making table prettier 

colnames(unknowns_age_sig_metabolites)[colnames(unknowns_age_sig_metabolites) == "estimate"] <- "\u0394Log(Self-Normalized Value) / \u0394Month"

unknowns_age_sig_metabolites <- unknowns_age_sig_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Peak ID` = peak_id,
         `Term` = term)

unknowns_age_sig_metabolites |> 
  head(10) |> 
  gt() |> 
  tab_header(
    title = "Unknown Metabolites that Vary Significantly Over Time"
  ) |> 
  fmt_number(
    columns = -`Peak ID`,
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )


# Volcano Plot 

# Calculating log fold change and significance across each timepoint
raw_stat_volcano <- merged_data_2 |>
  filter(Diet == 'CON') |>
  nest(stat_data = -peak_id) |>
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
    (log_2_fold_change < -0.2 & p_adj_BH < 0.05) ~ "Decrease",
    (log_2_fold_change > 0.2 & p_adj_BH < 0.05) ~ "Increase",
    TRUE ~ "No Significant Change"
  ))

# Significant values that I want to label on the graph
labels <- raw_stat_volcano |>
  filter(p_adj_BH < 0.05 & abs(log_2_fold_change) > 0.2) |>
  arrange(p_adj_BH)

labels <- labels$peak_id

raw_stat_volcano$labels <- ifelse(raw_stat_volcano$peak_id %in% labels,
                                  raw_stat_volcano$peak_id, NA)

# Actually plotting the graph
raw_stat_volcano |>
  ggplot(mapping = aes(x = log_2_fold_change, y = -log10(p_adj_BH), col = Change, label = labels)) +
  geom_vline(xintercept = c(-0.2, 0.2), col = 'gray', linetype = 'dashed') +
  geom_hline(yintercept = c(1.30102999566), col = 'gray', linetype = 'dashed') + # 1.3010 .. = -log10(0.05)
  geom_point(size = 1) +
  scale_color_manual(values = c("#00AFBB", "brown", "grey"),
                     labels = c("Decrease", "Increase", "No Significant Change")) +
  labs(color = 'Log Fold Change',
       x = expression("log"[2]*"FC"),
       y = expression("-log"[10]*"p-value")) +
  ggtitle(expression("Log"[2]*" Fold Change of Unknown Metabolites Correlated with Age")) +
  geom_text_repel(max.overlaps = 10, size = 6, min.segment.length = unit(0, 'lines')) +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 19),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        legend.position = "bottom")


# Creating heat map of labels 

# Getting data in proper format and generating heat map
heatmap_1 <- merged_data_2 |>
  filter(Diet == "CON" & peak_id %in% labels) |> 
  group_by(Timepoint, peak_id) |> 
  mutate(median_values = median(chebyshev_filter)) |> # Calculating median self_norm_value for each Timepoint and compound combination
  ungroup() |> 
  distinct(peak_id, Timepoint, median_values) |> 
  mutate(log_fold = log2(median_values)) |> 
  select(peak_id, Timepoint, log_fold)

heatmap_format <- heatmap_1 |> 
  pivot_wider(names_from = Timepoint, values_from = log_fold)

heatmap_format_2 <- heatmap_format[, -1]
rownames(heatmap_format_2) = heatmap_format$peak_id
heatmap_matrix <- as.matrix(heatmap_format_2)




# Creating Annotations for Heatmap
annotations <- as.data.frame(colnames(heatmap_matrix)) |> 
  mutate(`Age (Months)` = as.numeric(`colnames(heatmap_matrix)`) + 21)


# Creating and Plotting Heatmap 
unknowns_age_heatmap <- Heatmap(
  
  heatmap_matrix,
  name = "Log Fold Change",
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  
  # Rows
  row_names_gp = gpar(fontsize = 14),
  row_title = "Compound IDs",
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
  
)

draw(unknowns_age_heatmap, heatmap_legend_side = "top")


# Sex Effect -----

sex_fxn <- function(dat) {
  lm(log2(chebyshev_filter) ~ Timepoint*Sex + as.factor(Cohort), data = dat) |>
    broom::tidy()
}
  

raw_stat_2 <- merged_data_2 |> 
  filter(Diet == 'CON') |>
  nest(stat_data = -peak_id) |>
  mutate(df = purrr::map(stat_data, sex_fxn)) |>
  unnest_legacy(df) |>
  filter(term %in% c('SexM','Timepoint:SexM')) |> 
  ungroup()


raw_stat_2 <- raw_stat_2 |>
  mutate(p_adj_BH = p.adjust(p.value, method = "BH"))


sex_sig_metabolites <- raw_stat_2 |>
  filter(p_adj_BH < 0.05) |>
  arrange(p_adj_BH)


write.csv(sex_sig_metabolites, file = "Unknown Sex Significant Metabolites.csv")

top_ten <- sex_sig_metabolites |> 
  head(n = 10) |> 
  pull(peak_id)


filtered_data <- merged_data_2 |> 
  mutate(Timepoint = Timepoint + 21) |> 
  filter(Diet == "CON") |> 
  filter(peak_id %in% top_ten)

# # Visualizing significant compounds
filtered_data |> 
  ggplot() +
  geom_jitter(data = filtered_data |> filter(Timepoint != 21),
              aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), alpha = 0.5, size = 0.5) +
  geom_point(data = filtered_data |> filter(Timepoint == 21),
             aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), alpha = 0.5, size = 0.5) +
  geom_smooth(aes(x = Timepoint, y = log2(chebyshev_filter), color = Sex), method = "lm") +
  facet_wrap(factor(peak_id, levels = c("528", "1985", "536", "2891", "5468", "766", "526", "1875", "2874", "2735")) ~ ., nrow = 2, scales = "free_y") +
  labs(x = "Age (Months)", y = expression("log"[2]*" (Self Normalized Values)"), title = "Top 10 Most Significant Unknown Metabolites That Vary With Sex") +
  theme_minimal() +
  theme(axis.title.y = element_text(size = 14),
        axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        plot.title = element_text(size = 20)) +
  scale_x_continuous(breaks = seq(21, 34, by = 3))

# Visualizing the table nicely 
colnames(sex_sig_metabolites)[colnames(sex_sig_metabolites) == "estimate"] <- "\u0394ln(Self-Normalized Value) / \u0394Month"

sex_sig_metabolites <- sex_sig_metabolites |> 
  rename(`Standard Error` = std.error,
         `Statistic` = statistic,
         `P-value` = p.value,
         `Adjusted P-Value` = p_adj_BH,
         `Peak ID` = peak_id)

sex_sig_metabolites |> 
  head(10) |> 
  gt() |> 
  tab_header(
    title = "Unknown Metabolites That Vary Significantly Between Sexes"
  ) |> 
  fmt_number(
    columns = -`Peak ID`,
    decimals = 4
  ) |> 
  fmt_scientific(
    columns = c(`Adjusted P-Value`, `P-value`),
    decimals = 2
  )

top_25 <- sex_sig_metabolites |> 
  head(25) |> 
  pull(`Peak ID`)


# Creating heat map of labels 

# Getting data in proper format and generating heat map
test <- merged_data_2 |>
  filter(Diet == "CON" & peak_id %in% top_25) |> 
  group_by(Timepoint, Sex, peak_id) |> 
  mutate(median_values = median(chebyshev_filter)) |> # Calculating median self_norm_value for each Timepoint, Sex, Diet, and Pathway combination
  ungroup() |> 
  distinct(Timepoint, Sex, peak_id, median_values) |> 
  mutate(log_fold = log2(median_values))


heatmap_format <- test |> 
  mutate(Time_Sex = paste(Timepoint,Sex)) |> 
  arrange(Sex, Timepoint) |> 
  select(peak_id, Time_Sex, log_fold) |> 
  pivot_wider(names_from = Time_Sex, values_from = log_fold)

heatmap_format_2 <- heatmap_format[, -1]
rownames(heatmap_format_2) = heatmap_format$peak_id
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
  row_names_gp = gpar(fontsize = 13),
  row_title = "Compound Names",
  row_title_gp = gpar(fontsize = 14),
  row_title_side = "left",
  
  # Columns
  column_labels = annotations$`Age (Months)`,
  column_names_gp = gpar(fontsize = 14),
  column_title = "Age (Months)",
  column_title_gp = gpar(fontsize = 14),
  column_names_rot = 0,
  column_title_side = "bottom",
  
  # Annotations
  top_annotation = col_annotation,
  
)

draw(sex_heatmap)


# Hazard Analysis on Unknown Compounds -----------------------

library(survival)

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

# Analysis of Males
baseline_measures_m <- merged_data |> 
  filter(Timepoint == "mo00", Sex == "M")


joined_cox_m <- left_join(baseline_measures_m, m_survival, by = "Mouse") |> 
  group_by(peak_id) |> 
  mutate(z_scaled = scale(imputted_values, center = T)) |> 
  filter(!is.na(futime))


# Function that extracts hazard ratio and p value for each metabolite
# Model extracts hazard ratio of metabolites controlling for confounding variables Cohort 
# Here I used z_scaled versions of metabolite levels such that hazard ratio number varied more from 1
cox_function <- function(dat) {
  model <- coxph(Surv(futime, fustat) ~ z_scaled + as.factor(Cohort) + as.factor(Diet), data=dat, x=TRUE)
  pvalue <- c(summary(model)$coefficients[1, 5])
  m_h_ratio <- c(exp(model$coefficients[1]))
  return(data.frame(m_h_ratio, pvalue))
}

# Applying function to each metabolite in the data set
raw_stat_6 <- joined_cox_m |> 
  nest(stat_data = -peak_id) |>
  mutate(df = purrr::map(stat_data, cox_function)) |> 
  unnest(df) |> 
  ungroup()


raw_stat_6 <- raw_stat_6 |>
  select(peak_id, m_h_ratio, pvalue) |> 
  mutate(p_adj_BH = p.adjust(pvalue, method = "BH"))


write.csv(raw_stat_6, file = "U_Survival_metabolites_male.csv")


raw_stat_6 <- raw_stat_6 |> 
  mutate(label = case_when(
    p_adj_BH <= 0.05 ~ peak_id,
    TRUE ~ NA
  ))


# Visualizing results
# Results: Hazard Ratios > 1 suggests that the hazard (in this case death) is more probable
# if the compound is elevated at baseline measurment 
# https://pmc.ncbi.nlm.nih.gov/articles/PMC5388384/
raw_stat_6 |> 
  ggplot(aes(m_h_ratio, -log10(p_adj_BH), label = label)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 15, min.segment.length = unit(0, 'lines')) +
  labs(x = "Hazard Ratio", y = expression("-Log"[10] * "(p-value)"),
       title = "Hazard Ratio Of Unknown Metabolites in Males") +
  theme_minimal() +
  #caption = "Figure 1b. Results from Cox proportional hazard model on male mice, where the impact of baseline metabolite levels (standardized) on 
  #survival was analyzed. Compounds with a p-value less than or equal to 0.1 are displayed. A hazard ratio > 1 indicates that
  #the event of interest, in this case death, is more likely to occur when baseline levels of that metabolite are 
  #elevated. A hazard ratio < 1 suggests a smaller risk of death if baseline levels of that metabolite are elevated.") 
  theme(plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 18),
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
  filter(peak_id == "66") |> 
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
           xlab = "Age (days)", title = "Male Survival Curves\nStratified Along 66 Baseline Levels")



# Survival Analysis on Female Mice
baseline_measures_f <- merged_data |> 
  filter(Timepoint == "mo00", Sex == "F")


joined_cox_f <- left_join(baseline_measures_f, f_survival, by = "Mouse") |> 
  group_by(peak_id) |> 
  mutate(z_scaled = scale(imputted_values, center = T)) |> 
  filter(!is.na(futime))


# Cox_function modified for females
cox_function <- function(dat) {
  model <- coxph(Surv(futime, fustat) ~ z_scaled + as.factor(Cohort) + as.factor(Diet), data=dat, x=TRUE)
  pvalue <- c(summary(model)$coefficients[1, 5])
  f_h_ratio <- c(exp(model$coefficients[1]))
  return(data.frame(f_h_ratio, pvalue))
}


# Applying function to each metabolite in the data set for females
raw_stat_7 <- joined_cox_f |> 
  nest(stat_data = -peak_id) |>
  mutate(df = purrr::map(stat_data, cox_function)) |> 
  unnest(df) |> 
  ungroup()


raw_stat_7 <- raw_stat_7 |>
  select(peak_id, f_h_ratio, pvalue) |> 
  mutate(p_adj_BH = p.adjust(pvalue, method = "BH"))


write.csv(raw_stat_7, file = "U_Survival_metabolites_female.csv")


# Making compound labels
raw_stat_7 <- raw_stat_7 |> 
  mutate(label = case_when(
    p_adj_BH <= 0.6 ~ peak_id,
    f_h_ratio < 0.4 ~ peak_id,
    TRUE ~ NA
  ))


# Visualizing compounds associated with hazard in females
raw_stat_7 |> 
  ggplot(aes(f_h_ratio, -log10(p_adj_BH), label = label)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 10, min.segment.length = unit(0, 'lines')) +
  labs(x = "Hazard Ratio", y = expression("-Log"[10] * "(p-value)"),
       title = "Hazard Ratio Of Unknown Metabolites in Females") +
  theme(plot.caption = element_text(hjust = 0)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  theme_minimal() +
  theme(plot.caption = element_text(hjust = 0),
        axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20))


# Looking at survival curves of top 10 significant metabolites 
# C26:4, C14:0(Myristic acid), C22:4, PE(18:1/22:6), Indole-3-carboxylic acid,  
# C24:2, Carnosine, D-Ribose 5-phosphate, C24:4, C24:1
# Example
compound_specific_filter <- joined_cox_f |> 
  filter(peak_id == "4680") |> 
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
           xlab = "Age (days)", title = "Female Survival Curves\nStratified Along 4680 Baseline Levels")




# Comparing hazard ratio relative to beta value (slope) from linear models
slope_values <- raw_stat_2 |> 
  select(peak_id, term, estimate) |> 
  pivot_wider(names_from = term, values_from = estimate) |> 
  rename(f_slope = Timepoint) |> 
  mutate(m_slope = f_slope + `Timepoint:SexM`) |> 
  select(peak_id, f_slope, m_slope)


female_hazard <- raw_stat_7 |> 
  select(peak_id, f_h_ratio)


male_hazard <- raw_stat_6 |> 
  select(peak_id, m_h_ratio)


hazard_beta_merge <- left_join(left_join(slope_values, female_hazard, by = "peak_id"), 
                               male_hazard, by = "peak_id")


# Visualizing results for males
hazard_beta_merge |> 
  ggplot(aes(x = m_slope, y = m_h_ratio, label = peak_id)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 15, min.segment.length = unit(0, 'lines')) +
  geom_vline(xintercept = 0, alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, alpha = 0.4, linetype = "dashed") +
  labs(x = "Slope", y = "Hazard Ratio", title = "Hazard Ratio vs. Slope of Linear Regression Model\nIn Male Mice Unknown Compounds") +
  theme_minimal()+
  stat_cor() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20))


# Visualizing results for females
hazard_beta_merge |> 
  ggplot(aes(x = f_slope, y = f_h_ratio, label = peak_id)) +
  geom_point(size = 1, alpha = 0.7, color = "darkblue") +
  geom_text_repel(size = 5, max.overlaps = 8, min.segment.length = unit(0, 'lines')) +
  geom_vline(xintercept = 0, alpha = 0.4, linetype = "dashed") +
  geom_hline(yintercept = 1, alpha = 0.4, linetype = "dashed") +
  labs(x = "Slope", y = "Hazard Ratio", title = "Hazard Ratio vs. Slope of Linear Regression Model\nIn Female Mice Unknown Compounds") +
  theme_minimal()+
  stat_cor() +
  theme(axis.title.y = element_text(size = 16),
        axis.title.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        plot.title = element_text(size = 20))

merged_data_2 |> 
  filter(Diet == "CON") |> 
  distinct(Mouse, Sex) |> 
  count(Sex)


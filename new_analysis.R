getwd()
setwd("C:/Users/fishm/OneDrive/Documents/MSc EEC/Thesis/Results/all_analysis/Tyler_Christian_MSc_EEC_Thesis")

library(tidyverse)
library(factoextra)
library(cluster)
library(lme4)
library(patchwork)
library(MCMCglmm)
library(stats)
library(parallel)

# data --------------------------------------------------------------------

best_calls <- read.csv("avisoft_data.csv") %>%
  separate(audio_marker,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  mutate(
    audio_marker = as.factor(audio_marker),
    before_after = as.factor(before_after),
    site = as.factor(site),
    call = as.factor(call),
    segment = as.factor(segment)
  )

koe_data <- read.csv("koe_song_details.csv", na.strings = "") %>%
  mutate(call_type = sequence) %>%
  filter(quality != "INVALID") %>%
  group_by(filename) %>%
  separate(filename,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  separate_longer_delim(cols = call_type, delim = "-") %>%
  mutate(across(call_type, ~ str_remove_all(., '"')),
         across(call_type, ~ str_remove_all(., '_'))) %>%
  group_by(audio_marker, before_after, site, duration) %>%
  count(call_type) %>%
  pivot_wider(names_from = call_type, values_from = n) %>%
  replace(is.na(.), 0) %>%
  select(!c("NA")) %>%
  mutate(
    audio_marker = as.factor(audio_marker),
    before_after = as.factor(before_after),
    site = as.factor(site),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
    site = gsub(x = site, pattern = "nan", replacement = "mainland")
  ) %>%
  mutate(
    CC_per_time = CC / duration,
    DW_per_time = DW / duration,
    HA_per_time = HA / duration,
    CA_per_time = CA / duration,
    NC_per_time = NC / duration,
    TQ_per_time = TQ / duration
  )

behaviour_data <- read.csv("behavioural_data.csv") %>%
  replace(is.na(.), 0) %>%
  mutate(
    date = as.factor(date),
    site = as.factor(site),
    taxidermy = as.factor(taxidermy),
    audio_marker = as.factor(audio_marker),
    time = as.factor(time),
    taxidermy_present = as.factor(taxidermy_present),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
    site = gsub(x = site, pattern = "nan", replacement = "mainland")
  ) %>%
  filter(usable_stimuli == "x") %>%
  select(
    !c(
      trial,
      observation,
      number_male,
      number_female,
      audio_start:angle_video_finish,
      wild_pigeon:wild_other_passerine,
      usable_stimuli,
      notes
    )
  )

combined_data <- merge(behaviour_data, koe_data)

calls_long <- read.csv("koe_song_details.csv", na.strings = "") %>%
  mutate(call_type = sequence) %>%
  filter(quality != "INVALID") %>%
  group_by(filename) %>%
  separate(filename,
           into = c("audio_marker", "before_after"),
           sep = "_(?=a|b)") %>%
  separate_longer_delim(cols = call_type, delim = "-") %>%
  mutate(
    across(call_type, ~ str_remove_all(., '"')),
    across(call_type, ~ str_remove_all(., '_')),
    site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
    site = gsub(x = site, pattern = "nan", replacement = "mainland")
  ) %>%
  filter(site == "lundy") %>% 
  select(audio_marker, before_after, site, call_type, duration) %>%
  group_by(audio_marker, before_after, site) %>%
  count(call_type) %>%
  ungroup()

calls_long_stimuli <- merge(calls_long, behaviour_data) %>% 
  mutate(wild_starling = ifelse(wild_starling == 0,
                                NA,
                                "wild_starling"),
         wild_corvid = ifelse(wild_corvid == 0,
                              NA,
                              "wild_corvid"),
         wild_gull = ifelse(wild_gull == 0,
                            NA,
                            "wild_gull"),
         wild_raptor = ifelse(wild_raptor == 0,
                              NA,
                              "wild_raptor"),
         human_disturbance = ifelse(human_disturbance == 0,
                                    NA,
                                    "human_disturbance")) %>% 
  pivot_longer(cols = c(wild_starling:human_disturbance), names_to = "remove_me", values_to = "stimuli") %>% 
  select(c(site, before_after, call_type, stimuli, n)) %>% 
  na.omit() %>% 
  mutate(site = as.factor(site),
         before_after = as.factor(before_after),
         call_type = as.factor(call_type),
         stimuli = as.factor(stimuli)) %>% 
  filter(stimuli != "closed")

calls_long_taxidermy <- merge(calls_long, behaviour_data) %>% 
  filter(taxidermy_present != "closed") %>% 
  select(c(site, before_after, call_type, n, audio_marker)) %>% 
  na.omit() %>% 
  mutate(site = as.factor(site),
         before_after = as.factor(before_after),
         call_type = as.factor(call_type))

calls_long_mean_no_stimuli <- calls_long_stimuli %>% 
  group_by(before_after, call_type) %>% 
  summarise("mean" = mean(n)) %>% 
  mutate(call_type = replace_na(call_type, "none"),
         before_after = gsub(pattern = "a", replacement = "After", x = before_after),
         before_after = gsub(pattern = "b", replacement = "Before", x = before_after)) %>% 
  replace(is.na(.), 0) %>% 
  ungroup() %>% 
  add_row(before_after = "Before", call_type = "HA", mean = 0) %>% 
  add_row(before_after = "Before", call_type = "TQ", mean = 0)

combined_data_long <- combined_data %>% 
  group_by(audio_marker, before_after, site) %>% 
  mutate(wild_starling = ifelse(wild_starling == 0,
                                NA,
                                "starling"),
         wild_corvid = ifelse(wild_corvid == 0,
                              NA,
                              "corvid"),
         wild_gull = ifelse(wild_gull == 0,
                            NA,
                            "gull"),
         wild_raptor = ifelse(wild_raptor == 0,
                              NA,
                              "raptor"),
         human_disturbance = ifelse(human_disturbance == 0,
                                    NA,
                                    "human")) %>% 
  pivot_longer(cols = c(wild_starling:human_disturbance, taxidermy_present), 
               names_to = "stimuli",
               values_to = "stimuli_type") %>%
  na.omit() %>% 
  mutate(taxidermy_wild = ifelse(stimuli_type == "sparrow_taxidermy",
                                 "taxidermy",
                                 ifelse(stimuli_type == "woodpecker_taxidermy",
                                        "taxidermy",
                                        ifelse(stimuli_type == "starling_taxidermy",
                                               "taxidermy",
                                               ifelse(stimuli_type == "sparrowhawk_taxidermy",
                                                      "taxidermy",
                                                      ifelse(stimuli_type == "empty_taxidermy",
                                                             "taxidermy",
                                                             "wild")))))) %>% 
  ungroup()

# Dendrogram All Columns --------------------------------------------------------------

# Aggregate the data, filtering for lundy and removing HQ (only 1 sample but I don't want to delete data just incase)
aggregated_data <- best_calls %>% 
  filter(site == "lundy")

# View the aggregated data
print(aggregated_data)

# Standardize the data
data_scaled <- scale(aggregated_data %>%
                       select(c(7,9:14,16:21,23:27)))

# Calculate the distance matrix
dist_matrix <- dist(data_scaled, 
                        method = "euclidean")

# Perform hierarchical clustering
hc <- hclust(dist_matrix, 
                 method = "ward.D2")

# cophenetic correlation coefficient: 
coph_corr <- cor(cophenetic(hc), dist(aggregated_data))
coph_corr # acceptable

### The cophenetic correlation coefficient measures how faithfully a dendrogram preserves the pairwise distances between the original data points. A higher value indicates a better representation of the data.
# 0-0.5 bad
# 0.5-0.75 acceptable
# 0.75-1.0 very good

# Dendrogram Iterations --------------------------------------------------

### May take a while to run

# Define the column groups
column_groups <- list(
  c(7, 14, 21),
  c(9, 16, 23),
  c(10, 17, 24),
  c(11, 18, 25),
  c(12, 19, 26),
  c(13, 20, 27)
)

# Optional column
optional_column <- 6

# Function to calculate the cophenetic correlation coefficient for a given set of columns
calculate_cophenetic_corr <- function(data, columns) {
  data_scaled <- scale(data %>% select(all_of(columns)))
  dist_matrix <- dist(data_scaled, method = "euclidean")
  hc <- hclust(dist_matrix, method = "ward.D2")
  cophenetic_corr <- cor(cophenetic(hc), dist(dist_matrix))
  return(cophenetic_corr)
}

# Generate all combinations of columns from the groups
generate_combinations <- function(groups, optional) {
  all_combinations <- expand.grid(lapply(groups, function(group) group))
  colnames(all_combinations) <- paste0("V", 1:ncol(all_combinations))
  combinations_with_optional <- rbind(
    cbind(Optional = optional, all_combinations),
    cbind(Optional = NA, all_combinations)
  )
  return(combinations_with_optional)
}

# Find the best combination of columns
best_combination <- NULL
best_cophenetic_corr <- -Inf

# Generate all possible combinations
column_combinations <- generate_combinations(column_groups, optional_column)
total_combinations <- nrow(column_combinations)
current_combination <- 0

# Function to process each combination
process_combination <- function(i) {
  combination <- na.omit(as.numeric(column_combinations[i, ]))
  cophenetic_corr <- calculate_cophenetic_corr(aggregated_data, combination)
  return(list(combination = combination, cophenetic_corr = cophenetic_corr))
}

# Counter function to print progress
counter_function <- function(i) {
  current_combination <<- current_combination + 1
  cat("Processing combination", current_combination, "of", total_combinations, "\n")
  return(process_combination(i))
}

# Run the combinations in parallel using parLapply
num_cores <- detectCores() - 1
cl <- makeCluster(num_cores)

# Load necessary libraries and export objects to each node in the cluster
clusterEvalQ(cl, {
  library(dplyr)
  library(stats)
})
clusterExport(cl, list("column_combinations", "calculate_cophenetic_corr", "aggregated_data", "process_combination", "current_combination", "total_combinations"))

# Run the combinations in parallel
results <- parLapply(cl, 1:total_combinations, function(i) {
  current_combination <<- current_combination + 1
  cat("Processing combination", current_combination, "of", total_combinations, "\n")
  return(process_combination(i))
})

stopCluster(cl)

# Find the best result
for (result in results) {
  if (result$cophenetic_corr > best_cophenetic_corr) {
    best_cophenetic_corr <- result$cophenetic_corr
    best_combination <- result$combination
  }
}

# Print the best combination of columns and the corresponding cophenetic correlation coefficient
cat("Best combination of columns:", best_combination, "\n")
cat("Best cophenetic correlation coefficient:", best_cophenetic_corr, "\n")

### Output:
# Best combination of columns: 14 23 17 18 26 27 
# Best cophenetic correlation coefficient: 0.8334661 

# Optimised Dendrogram (Final) ----------------------------------------------------

# Aggregate the data, filtering for lundy and removing HQ (only 1 sample but I don't want to delete data just incase)
aggregated_data_optimised <- best_calls %>% 
  filter(site == "lundy")

# View the aggregated data
print(aggregated_data_optimised)

# Standardize the data
data_scaled_optimised <- scale(aggregated_data_optimised %>%
                       select(c(14,23,17,18,26,27)))

# Calculate the distance matrix
dist_matrix_optimised <- dist(data_scaled_optimised, 
                    method = "euclidean")

# Perform hierarchical clustering
hc_optimised <- hclust(dist_matrix_optimised, 
             method = "ward.D2")

# Set colour groups per call
aggregated_data_optimised <- aggregated_data_optimised %>% 
  mutate(color_call = ifelse(call == "CA",
                             "#88ccee",
                             ifelse(call == "CC",
                                    "#661100",
                                    ifelse(call == "DW",
                                           "#117733",
                                           ifelse(call == "HA",
                                                  "#999933",
                                                  ifelse(call == "NC",
                                                         "#e31a1c",
                                                         "#aa4499" ### "TQ"
                                                  )
                                           )
                                    )
                             )
  )
  ) %>% 
  mutate(color_call = as.factor(color_call))

# Enhanced dendrogram visualization with color-coded nodes by call 
# (only the colours at the bottom have any meaning, see attached ppt for interpretation)
dendrogram_optimised <- 
  fviz_dend(
    hc_optimised,
    rect = FALSE,
    cex = 0.5,
    lwd = 0.5,
    k_colors = "black",
    labels_track_height = 0.8,
    show_labels = TRUE,
    label_cols = aggregated_data_optimised$color_call
  ) +
  theme(plot.title = element_blank()) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 44,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#88ccee",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 22,
    y = -1.55,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 44,
    xmax = 95,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#661100",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 69.5,
    y = -1.55,
    size = 2.5,
    label = "CC"
  ) +
  annotate(
    "rect",
    xmin = 95,
    xmax = 103,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 99,
    y = -1.55,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 103,
    xmax = 118,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#88ccee",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 110.5,
    y = -1.55,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 118,
    xmax = 125,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#117733",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 121.5,
    y = -1.55,
    size = 2.5,
    label = "DW"
  ) +
  annotate(
    "rect",
    xmin = 125,
    xmax = 137,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#999933",
    color = "black",
    lwd = 0.3
  )  +
  annotate(
    "text",
    x = 131,
    y = -1.55,
    size = 2.5,
    label = "HA"
  ) +
  annotate(
    "rect",
    xmin = 137,
    xmax = 142,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#aa4499",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 139.5,
    y = -1.55,
    size = 2.5,
    label = "TQ"
  ) +
  annotate(
    "rect",
    xmin = 142,
    xmax = 146,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 144,
    y = -1.55,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 146,
    xmax = 159,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#aa4499",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 152.5,
    y = -1.55,
    size = 2.5,
    label = "TQ"
  ) +
  annotate(
    "rect",
    xmin = 159,
    xmax = 167,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#e31a1c",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 163,
    y = -1.55,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 50,
    ymin = -3.2,
    ymax = -2.3,
    alpha = .5,
    fill = "white",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 25,
    y = -2.75,
    size = 2.5,
    label = (bold("1:") ~ " CA")
  ) +
  annotate(
    "rect",
    xmin = 50,
    xmax = 109,
    ymin = -3.2,
    ymax = -2.3,
    alpha = .5,
    fill = "white",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 79.5,
    y = -2.75,
    size = 2.5,
    label = (bold("2:") ~ " CC")
  ) +
  annotate(
    "rect",
    xmin = 109,
    xmax = 167,
    ymin = -3.2,
    ymax = -2.3,
    alpha = .5,
    fill = "white",
    color = "black",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 138,
    y = -2.75,
    size = 2.5,
    label = (bold("3:") ~ " Other")
  ) +
  annotate(
    "rect",
    xmin = 158,
    xmax = 169,
    ymin = 17.5,
    ymax = 24.5,
    color = "black",
    fill = "white",
    lwd = 0.3
  ) +
  annotate(
    "text",
    x = 163,
    y = 24,
    size = 3.7,
    label = 'bold(" Key:")',
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 163,
    y = 23,
    size = 3.7,
    label = "CC",
    color = "#661100"
  ) +
  annotate(
    "text",
    x = 163,
    y = 22,
    size = 3.7,
    label = "NC",
    color = "#e31a1c"
  ) +
  annotate(
    "text",
    x = 163,
    y = 21,
    size = 3.7,
    label = "CA",
    color = "#88ccee"
  ) +
  annotate(
    "text",
    x = 163,
    y = 20,
    size = 3.7,
    label = "TQ",
    color = "#aa4499"
  ) +
  annotate(
    "text",
    x = 163,
    y = 19,
    size = 3.7,
    label = "HA",
    color = "#999933"
  ) +
  annotate(
    "text",
    x = 163,
    y = 18,
    size = 3.7,
    label = "DW",
    color = "#117733"
  ) +
  annotate(
    "text",
    x = -5.1,
    y = -1.55,
    size = 2.5,
    label = "Minor"
  ) +
  annotate(
    "text",
    x = -5.1,
    y = -2.75,
    size = 2.5,
    label = "Major"
  )

dendrogram_optimised

ggsave("dendrogram_optimised.jpg", width = 6.8, height = 6.5, units = "in")

# LM ----------------------------------------------------------------------

plot(combined_data$CC_per_time)
plot(combined_data$NC_per_time)
plot(combined_data$CA_per_time) 
plot(combined_data$HA_per_time)
plot(combined_data$TQ_per_time)
plot(combined_data$DW_per_time)
# not normally distributed, therefore use GLM

### use these ###
# explain taxidermy and stimuli may impact behavior and are thus random variables
# site may have some impacts and is thus random (link to potential island effect but not too much detail)

CC_model <- lmer(CC_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)
NC_model <- lmer(NC_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)
TQ_model <- lmer(TQ_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)
HA_model <- lmer(HA_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)
CA_model <- lmer(CA_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)
DW_model <- lmer(DW_per_time ~ before_after + (1|taxidermy_wild) * (1 | stimuli_type) + (1|site), data = combined_data_long)

summary(CC_model)
summary(NC_model)
summary(TQ_model)
summary(HA_model)
summary(CA_model)
summary(DW_model)

plot(CC_model)
plot(NC_model)
plot(TQ_model)
plot(HA_model)
plot(CA_model)
plot(DW_model)

# Before/After plot -------------------------------------------------------

before_after_plot <- calls_long_mean_no_stimuli %>% 
  ggplot(aes(
    x = before_after,
    y = mean,
    group = call_type,
    color = call_type
  )) +
  geom_point() + 
  geom_line() +
  ylab("Mean Call Frequency") +
  labs(color = "Call Type") +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  scale_color_manual(values = c("#88ccee","#661100","#117733","#999933","#e31a1c","#aa4499" ))

before_after_plot

ggsave("before_after_plot.jpg", height = 7.5, width = 6.8, units = "in")

# Taxidermy plot (Methods) ------------------------------------------------

change_data <- combined_data_long %>% 
  filter(stimuli_type == c("starling", "starling_taxidermy")) %>% 
  group_by(before_after, audio_marker) %>% 
  pivot_longer(cols = CC:TQ, names_to = "call_type", values_to = "call_frequency") %>% 
  pivot_wider(names_from = before_after, values_from = call_frequency) %>% 
  replace(is.na(.), 0) %>% 
  mutate(difference_before_after = a-b) %>% 
  ungroup() %>% 
  group_by(site, taxidermy_wild, call_type) %>% 
  summarise(mean_difference_before_after = mean(difference_before_after)) %>% 
  mutate(call_type = gsub(x = call_type, pattern = "_per_time", replacement = "")) %>% 
  ungroup()

change_plot <- change_data %>% 
  group_by(taxidermy_wild, call_type) %>% 
  summarise(mean_difference_before_after = mean(mean_difference_before_after)) %>% 
  mutate(taxidermy_wild = gsub(x = taxidermy_wild, pattern = "t", replacement = "T"),
         taxidermy_wild = gsub(x = taxidermy_wild, pattern = "w", replacement = "W")) %>% 
  ggplot(aes(x = call_type,
             y = mean_difference_before_after,
             group = taxidermy_wild,
             fill = taxidermy_wild)) +
  geom_col(position = "dodge",
           color = "black",
           lwd = 0.4) +
  theme_bw() +
  xlab("Call Type") +
  ylab("Change in Call Frequency After Stimuli Reveal") +
  labs(fill = "Stimuli Type") +
  geom_hline(yintercept = 0, lwd = 0.4) +
  scale_fill_manual(values=c("#661100","#117733"))

change_plot

ggsave("change_plot.jpg", height = 6.8, width = 6.8, units = "in")

change_model <- lm(mean_difference_before_after ~ call_type * taxidermy_wild, data = change_data)

summary(change_model)

# explain that whilst there are no significant differences, taxidermy has a visual impact on calls produced before/after with 50% being the opposite effect to wild stimuli. 
# therefore taxidermy is included as a random effect but the reality of its responses should be tested before being used in similar studies.
# stronger evidence is required to validate this

# Island Effect PCA (SI) -------------------------------------------------------

### There seems to be no island effect (PCA clustered through calls, not location). 
### Maybe don't mention island effect at all???

call_pca_data <- best_calls %>% 
  #filter(call == "CC") %>% 
  mutate(site = gsub(x = site, pattern = "ascot", replacement = "mainland"),
         site = gsub(x = site, pattern = "nan", replacement = "mainland"))

call_pca <- princomp(call_pca_data[, c(7,9:14,16:21,23:27)], # remove duration and peak amplitude
                          cor = TRUE)

eigenvalues <- call_pca$sdev^2
eigenvalues
# 2 components

### The Scree Plot Criterion:

plot(call_pca, type="lines", ylim=c(0,10))
# 2 or 3 components

### The Relative Percent Variance Criterion:

summary(call_pca)
# 2 components represent 55.9% of variance

### Interpreting the component structure:
loadings(call_pca)

call_pca_plot <- fviz_pca_biplot(call_pca, 
                           repel = TRUE,
                           geom.ind = "point",
                           geom.var = "none",
                           ellipse.level = 0.95,
                           col.var = "black", 
                           labelsize = 4,
                           habillage = call_pca_data$call,
                           addEllipses = TRUE,
                           alpha.ind = 0.5,
                           title = "Calls") # change group colours


call_pca_plot

island_pca_plot <- fviz_pca_biplot(call_pca, 
                                 repel = TRUE,
                                 geom.ind = "point",
                                 geom.var = "none",
                                 ellipse.level = 0.95,
                                 col.var = "black", 
                                 labelsize = 4,
                                 habillage = call_pca_data$site,
                                 addEllipses = TRUE,
                                 alpha.ind = 0.5,
                                 title = "Site") # change group colours


island_pca_plot

call_pca_plot + island_pca_plot


####

# Extract PCA results
pca_results <- as.data.frame(call_pca$scores)
pca_results <- cbind(call_pca_data, pca_results)

# Create the base PCA biplot
base_pca_plot <- ggplot(pca_results, aes(x = Comp.1, y = Comp.2)) +
  geom_point(aes(color = call)) +
  theme_minimal() +
  labs(title = "PCA Biplot", x = "PC1", y = "PC2")

# Add ellipses for site
site_ellipses <- base_pca_plot +
  stat_ellipse(aes(x = Comp.1, y = Comp.2, color = site, linetype = "solid"), level = 0.95)

# Add ellipses for call
final_pca_plot <- site_ellipses +
  stat_ellipse(aes(x = Comp.1, y = Comp.2, color = call, linetype = "dashed"), level = 0.95) +
  scale_color_manual(values = c("black", "grey","blue", "red", "green", "purple", "orange", "yellow")) # Adjust colors as needed

# Display the plot
print(final_pca_plot)


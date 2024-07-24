getwd()
setwd("C:/Users/fishm/OneDrive/Documents/MSc EEC/Thesis/Results/all_analysis/Tyler_Christian_MSc_EEC_Thesis")

library(tidyverse)
library(factoextra)
library(cluster)
library(lme4)
library(patchwork)
library(MCMCglmm)

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
    site = as.factor(site)
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
    taxidermy_present = as.factor(taxidermy_present)
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

# Dendrogram --------------------------------------------------------------

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

# Set colour groups per call
aggregated_data <- aggregated_data %>% 
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
dendrogram <- 
  fviz_dend(
  hc,
  rect = FALSE,
  cex = 0.5,
  lwd = 0.5,
  k_colors = "black",
  labels_track_height = 0.8,
  show_labels = TRUE,
  label_cols = aggregated_data$color_call
) +
  theme(plot.title = element_blank()) +
  annotate(
    "rect",
    xmin = 0,
    xmax = 45,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#88ccee"
  ) +
  annotate(
    "text",
    x = 22.5,
    y = -1.55,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 45,
    xmax = 56,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#661100"
  ) +
  annotate(
    "text",
    x = 50.5,
    y = -1.55,
    size = 2.5,
    label = "CC"
  ) +
  annotate(
    "rect",
    xmin = 56,
    xmax = 60,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#88ccee"
  ) +
  annotate(
    "text",
    x = 58,
    y = -1.55,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 60,
    xmax = 96,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#661100"
  ) +
  annotate(
    "text",
    x = 78,
    y = -1.55,
    size = 2.5,
    label = "CC"
  ) +
  annotate(
    "rect",
    xmin = 96,
    xmax = 101,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#e31a1c"
  ) +
  annotate(
    "text",
    x = 98.5,
    y = -1.55,
    size = 2.5,
    label = "NC"
  ) +
  annotate(
    "rect",
    xmin = 101,
    xmax = 117,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#88ccee"
  ) +
  annotate(
    "text",
    x = 109,
    y = -1.55,
    size = 2.5,
    label = "CA"
  ) +
  annotate(
    "rect",
    xmin = 117,
    xmax = 129,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#117733"
  ) +
  annotate(
    "text",
    x = 123,
    y = -1.55,
    size = 2.5,
    label = "DW"
  ) +
  annotate(
    "rect",
    xmin = 129,
    xmax = 137,
    ymin = -2,
    ymax = -1.1,
    alpha = .3,
    fill = "#999933"
  )  +
  annotate(
    "text",
    x = 133,
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
    fill = "#aa4499"
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
    fill = "#e31a1c"
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
    fill = "#aa4499"
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
    fill = "#e31a1c"
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
    xmin = 158,
    xmax = 169,
    ymin = 26.5,
    ymax = 33.5,
    color = "black",
    fill = "white"
  ) +
  annotate(
    "text",
    x = 163,
    y = 33,
    size = 3.7,
    label = 'bold(" Key:")',
    parse = TRUE
  ) +
  annotate(
    "text",
    x = 163,
    y = 32,
    size = 3.7,
    label = "CC",
    color = "#661100"
  ) +
  annotate(
    "text",
    x = 163,
    y = 31,
    size = 3.7,
    label = "NC",
    color = "#e31a1c"
  ) +
  annotate(
    "text",
    x = 163,
    y = 30,
    size = 3.7,
    label = "CA",
    color = "#88ccee"
  ) +
  annotate(
    "text",
    x = 163,
    y = 29,
    size = 3.7,
    label = "TQ",
    color = "#aa4499"
  ) +
  annotate(
    "text",
    x = 163,
    y = 28,
    size = 3.7,
    label = "HA",
    color = "#999933"
  ) +
  annotate(
    "text",
    x = 163,
    y = 27,
    size = 3.7,
    label = "DW",
    color = "#117733"
  ) 

dendrogram

ggsave("dendrogram.jpg", width = 6.8, height = 7.5, units = "in")

coph_corr <- cor(cophenetic(hc), dist(aggregated_data))
coph_corr # acceptable

# The cophenetic correlation coefficient measures how faithfully a dendrogram preserves the pairwise distances between the original data points. A higher value indicates a better representation of the data.

# 0-0.5 bad
# 0.5-0.75 acceptable
# 0.75-1.0 very good




# LM ----------------------------------------------------------------------

plot(combined_data$CC_per_time)
plot(combined_data$NC_per_time)
plot(combined_data$CA_per_time) 
plot(combined_data$HA_per_time)
plot(combined_data$TQ_per_time)
plot(combined_data$DW_per_time)
# not normally distributed, therfore use GLM

combined_data_no_tax <- combined_data %>% 
  filter(open == 0)

summary(glm(combined_data$CC_per_time ~ combined_data$before_after))
summary(glm(combined_data$NC_per_time ~ combined_data$before_after))
summary(glm(combined_data$CA_per_time ~ combined_data$before_after))
summary(glm(combined_data$HA_per_time ~ combined_data$before_after))
summary(glm(combined_data$TQ_per_time ~ combined_data$before_after))
summary(glm(combined_data$DW_per_time ~ combined_data$before_after))

combined_data_long <- combined_data %>% 
  group_by(audio_marker, before_after) %>% 
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
  pivot_longer(cols = wild_starling:human_disturbance, 
               names_to = "stimuli",
               values_to = "stimuli_type") %>% 
  na.omit()

combined_data_numeric_stimuli <- combined_data %>% 
  mutate(taxidermy_present = ifelse(taxidermy_present == "starling_taxidermy",
                                    y = 0.4,
                                    n = ifelse(
                                      taxidermy_present == "woodpecker_taxidermy",
                                      y = 0.15,
                                      n = ifelse(
                                        taxidermy_present == "sparrow_taxidermy",
                                        y = 0.15,
                                        n = ifelse(
                                          taxidermy_present == "sparrowhawk_taxidermy",
                                          y = 0.5,
                                          n = ifelse(
                                            taxidermy_present == "empty_taxidermy",
                                            y = 0.1,
                                            n = 0
                                    ))))),
         wild_starling = ifelse(wild_starling == 0,
                                0,
                                0.8),
         wild_corvid = ifelse(wild_corvid == 0,
                              0,
                              0.9),
         wild_gull = ifelse(wild_gull == 0,
                            0,
                            0.9),
         wild_raptor = ifelse(wild_raptor == 0,
                              0,
                              1),
         total_threat = wild_starling + wild_corvid + wild_gull + wild_raptor + taxidermy_present)

CC_model1 <- lmer(CC_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type), data = combined_data_long)
CC_model2 <- lmer(CC_per_time ~ before_after + (1|total_threat), data = combined_data_numeric_stimuli)
CC_model3 <- lmer(CC_per_time ~ before_after + (1|taxidermy_present), data = combined_data)
CC_model4 <- lmer(CC_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
CC_model5 <- lmer(CC_per_time ~ before_after + (1|total_threat)+(1|taxidermy_present), data = combined_data_numeric_stimuli)

CC_model <- lmer(CC_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
NC_model <- lmer(NC_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
TQ_model <- lmer(TQ_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
HA_model <- lmer(HA_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
CA_model <- lmer(CA_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)
DW_model <- lmer(DW_per_time ~ before_after + (1|total_threat)*(1|taxidermy_present), data = combined_data_numeric_stimuli)


CC_model <- lmer(CC_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
NC_model <- lmer(NC_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
TQ_model <- lmer(TQ_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
HA_model <- lmer(HA_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
CA_model <- lmer(CA_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
DW_model <- lmer(DW_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)

combined_data_numeric_stimuli <- combined_data_numeric_stimuli %>% 
  filter(total_threat > 0.5)

CC_model <- glm(CC_per_time ~ before_after, data = combined_data_numeric_stimuli)
NC_model <- glm(NC_per_time ~ before_after, data = combined_data_numeric_stimuli)
TQ_model <- glm(TQ_per_time ~ before_after, data = combined_data_numeric_stimuli)
HA_model <- glm(HA_per_time ~ before_after, data = combined_data_numeric_stimuli)
CA_model <- glm(CA_per_time ~ before_after, data = combined_data_numeric_stimuli)
DW_model <- glm(DW_per_time ~ before_after, data = combined_data_numeric_stimuli)


summary(CC_model)
summary(NC_model)
summary(TQ_model)
summary(HA_model)
summary(CA_model)
summary(DW_model)

lrtest(CC_model2, CC_model4)

### use these ###
# explain taxidermy and stimuli may impact behavior and are thus random variables
# site may have some impacts and is thus random (link to potential island effect but not too much detail)

CC_model <- lmer(CC_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
NC_model <- lmer(NC_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
TQ_model <- lmer(TQ_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
HA_model <- lmer(HA_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
CA_model <- lmer(CA_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)
DW_model <- lmer(DW_per_time ~ before_after + (1|taxidermy_present) * (1 | stimuli_type) + (1|site), data = combined_data_long)

summary(CC_model)
summary(NC_model)
summary(TQ_model)
summary(HA_model)
summary(CA_model)
summary(DW_model)

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
  theme(axis.title.x = element_blank())

before_after_plot

# takes mean call frequency for all wild stimuli on lundy before and after appearance

### repeat for each type of wild stimuli???

# Calls over time ---------------------------------------------------------

time_data_tax <- calls_long_taxidermy %>% 
  separate(audio_marker,
           into = c("trial", "observation"),
           sep = "_") %>%
  mutate(trial = as.numeric(trial)) %>% 
  group_by(trial, call_type) 

time_plot <- time_data_tax %>% 
  ggplot(aes(x = trial,
             y = n,
             group = call_type,
             color = call_type)) +
  geom_smooth(se = FALSE, method = "lm") +
  geom_point()

time_plot # do calls decrease with time? 

summary(lm(time_data_tax$n[time_data_tax$call_type=="CC"] ~ time_data_tax$trial[time_data_tax$call_type=="CC"]))
summary(lm(time_data_tax$n[time_data_tax$call_type=="TQ"] ~ time_data_tax$trial[time_data_tax$call_type=="TQ"]))
summary(lm(time_data_tax$n[time_data_tax$call_type=="HA"] ~ time_data_tax$trial[time_data_tax$call_type=="HA"]))
summary(lm(time_data_tax$n[time_data_tax$call_type=="CA"] ~ time_data_tax$trial[time_data_tax$call_type=="CA"]))
summary(lm(time_data_tax$n[time_data_tax$call_type=="DW"] ~ time_data_tax$trial[time_data_tax$call_type=="DW"]))

# CC has a significant decrease over time (all others insignificant)
# y = -0.796x + 21.080 
# R = 0.348 (poor fit)


# Taxidermy plot (Methods) ------------------------------------------------

change_data <- read.csv("change_in_calls.csv") %>% 
  pivot_longer(cols = Taxidermy:Wild,
               values_to = "change",
               names_to = "Specimen") %>% 
  filter(Call != "Total") %>% 
  mutate(Call = as.factor(Call),
         Specimen = as.factor(Specimen)) 

change_plot <- change_data %>% 
  ggplot(aes(x = fct_inorder(Call),
             y = change,
             group = Specimen,
             fill = Specimen)) +
  geom_col(position = "dodge") +
  theme_bw() +
  xlab("Call Type") +
  ylab("Change in Frequency After Stimuli Reveal")

change_plot

summary(lm(change_data$change~ change_data$Specimen))
t.test(change_data$change~ change_data$Specimen)

# Island Effect PCA (Methods) -------------------------------------------------------

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



## Competition constrains parasite adaptation to thermal heterogeneity
## Samuel TE Greenrod, Daniel Cazares, Weronika Ślesak, Tobias E Hector, R. Craig MacLean, Kayla C King

## This R script was used to generate all figures and run all analyses for the manuscript labelled above.




library("car")
library("ggplot2")
library("RColorBrewer")
library("patchwork")
library("zoo") 
library("broom")
library("ggbeeswarm")
library("emmeans")
library("ggsignif")
library("FSA")
library("pgirmess")
library("DescTools")
library("ggpmisc")
library("plyr")
library("devtools")
library("ggradar")
library("scales")
library("tidyverse") 
library("patchwork")
library("forcats")
library("cowplot")
library("performance")
library("lme4")
library("read.gb")
library("ape")
library("gggenes")
library("vegan")
library("ggtree")
library("dplyr")
library("paco")
library("ggh4x")
library("ggpattern")

brewer.pal(12, "Paired")


# Figure 1

# Figure 1A

pfu <- read_csv("Data/PFU_values_each_passage.csv")
pfu$phage_measured <- factor(pfu$phage_measured, levels = c("phage14_1", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))
pfu$phage_treatment <- factor(pfu$phage_treatment, levels = c("Single_phage", "Competition"), labels = c("Monoculture", "Co-culture"))

pfu <- pfu %>%
  filter(phage_treatment == "Monoculture") %>%
  filter(passage %% 1 == 0)

x_min <- floor(min(pfu$passage, na.rm = TRUE))
x_max <- ceiling(max(pfu$passage, na.rm = TRUE))
rects <- tibble(
  xmin = seq(x_min, x_max - 1, by = 1),
  xmax = xmin + 1,
  fill = rep(c("tomato", "lightblue"), length.out = length(xmin))
)

ggplot(pfu, aes(passage, pfu_ml, group = interaction(phage_treatment, replicate), 
                color = phage_treatment)) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = 1e7, ymax = 1e11, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_line(linewidth = 0.2, linetype = "dashed") +
  geom_line(aes(group = phage_treatment), 
            stat = "summary", fun = mean, linewidth = 1.2, alpha = 0.5) +
  facet_wrap(~phage_measured) +  
  scale_y_continuous(trans = 'log10', limits = c(1e7, 1e11)) +
  annotation_logticks(sides = "l") +
  scale_color_manual(values = c("Monoculture" = "#6A3D9A", "Co-culture" = "grey10")) +
  labs(x = "Passage", y = "Phage density (PFU/ml)", colour = "Temperature treatment") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )

ggsave("Figures/Fig_1a.png", width = 8, height = 3, dpi = 300)



# Figure 1B

phage_fitness_pfu <- read.csv("Data/Phage_population_growth_rates.csv")

phage_fitness_pfu_update <- phage_fitness_pfu %>% 
  mutate(competition = ifelse(grepl("anc", line), "anc",
                              ifelse(grepl("C", line), "Competition", "No_competition"))) %>%  # Label lines which evolved with competitors
  mutate(evolution_treatment = case_when(   # Label lines based on their temperature evolution treatment
    grepl("anc", line) ~ "anc",
    grepl("37", line) ~ "37_static",
    grepl("42", line) ~ "42_static",
    grepl("F", line) ~ "Fluctuating",
    TRUE ~ NA_character_
  )) %>%
  filter(!(phage == "phage14" & time == "2"))

phage_fitness_pfu_update$temp <-  factor(phage_fitness_pfu_update$temp, levels = c("37","42"), labels = c("37\u00B0C", "42\u00B0C"))
phage_fitness_pfu_update$line <- factor(phage_fitness_pfu_update$line, levels = c("anc","37_evo","37_C_evo", "42_evo", "42_C_evo", "F_evo", "F_C_evo"), labels = c("Ancestor", "37\u00B0C static", "37\u00B0C static (comp)", "42\u00B0C static", "42\u00B0C static (comp)","Fluctuating", "Fluctuating (comp)"))
phage_fitness_pfu_update$phage <- factor(phage_fitness_pfu_update$phage, levels = c("phage14", "LUZ19"), labels = c("ϕ14-1", "ϕLUZ19"))

phage_fitness_pfu_update_single <- phage_fitness_pfu_update %>%
  filter(line %in% c("Ancestor", "37\u00B0C static", "42\u00B0C static", "Fluctuating"))

phage_fitness_pfu_update_single <- phage_fitness_pfu_update_single %>%
  select(-c(dilution, pfu)) %>%
  group_by(across(-c(batch, pfu_ml))) %>%  # Group by all columns except `value_column`
  dplyr::summarize(pfu_ml = mean(pfu_ml, na.rm = TRUE)) %>%  # Average duplicates in `value_column`
  ungroup()

phage_fitness_pfu_update_boxplots = phage_fitness_pfu_update_single %>%
  filter((phage == "ϕ14-1" & time == "4") |
           (phage == "ϕLUZ19" & time == "2"))   # Remove T2 values from 14-1 dataset as only carried out for a single rep


phage_fitness_pfu_update_boxplots_single = phage_fitness_pfu_update_boxplots %>%
  mutate(phage_temp = paste(phage, temperature, sep = "_"))

phage_fitness_pfu_update_boxplots_single$phage_temp = factor(phage_fitness_pfu_update_boxplots_single$phage_temp, levels = c("ϕ14-1_37", "ϕ14-1_42", "ϕLUZ19_37", "ϕLUZ19_42"), labels = c("ϕ14-1 (37°C)", "ϕ14-1 (42°C)", "ϕLUZ19 (37°C)", "ϕLUZ19 (42°C)"))



phage_fitness_pfu_update_boxplots_single_normalised <- phage_fitness_pfu_update_boxplots_single %>%
  dplyr::group_by(phage_temp) %>%
  dplyr::mutate(
    ancestor_mean = mean(pfu_ml[line == "Ancestor"], na.rm = TRUE),
    PFU_relative = pfu_ml / ancestor_mean
  ) %>%
  ungroup()

ancestor_se <- phage_fitness_pfu_update_boxplots_single_normalised %>%
  dplyr::filter(line == "Ancestor") %>%
  dplyr::group_by(phage_temp) %>%
  dplyr::summarise(se = sd(PFU_relative, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  dplyr::mutate(
    ymin = 1 - se,
    ymax = 1 + se,
    xmin = -Inf,
    xmax = Inf
  )

phage_fitness_pfu_update_boxplots_single_normalised = phage_fitness_pfu_update_boxplots_single_normalised %>%
  dplyr::filter(line != "Ancestor")

cbPalette2 <- c("#377EB8", "#6A3D9A", "#E41A1C")

phage_fitness_pfu_update_boxplots_single_normalised$line = factor(phage_fitness_pfu_update_boxplots_single_normalised$line, levels = c("37°C static", "Fluctuating", "42°C static"), labels = c("37°C (static)", "Fluctuating", "42°C (static)"))

label_map <- c(
  "ϕ14-1 (37°C)" = "ϕ14-1",
  "ϕ14-1 (42°C)" = "ϕ14-1",
  "ϕLUZ19 (37°C)" = "ϕLUZ19",
  "ϕLUZ19 (42°C)" = "ϕLUZ19"
)

my_labeller <- as_labeller(label_map)


ggplot(phage_fitness_pfu_update_boxplots_single_normalised, aes(line, PFU_relative, colour = as.factor(line))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap2(~ phage_temp, nrow = 1, strip.position = "top", strip = strip_themed(
    background_x = list(
      A = element_rect(fill = "lightblue"),
      B = element_rect(fill = "tomato"),
      C = element_rect(fill = "lightblue"),
      D = element_rect(fill = "tomato")
    ),
    
    text_x = list(
      A = element_text(color = "black", face = "bold"),
      B = element_text(color = "black", face = "bold"),
      C = element_text(color = "black", face = "bold"),
      D = element_text(color = "black", face = "bold")
    )
  ),
  labeller = my_labeller 
  )+
  ylab("Phage growth rel. to ancestor") +
  xlab("Evolution treatment")+
  labs(colour="Evolution treatment")+
  scale_y_continuous(trans='log10', limits = c(0.01,10^6)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "black", alpha = 0.2, inherit.aes = FALSE)+
  scale_y_continuous(trans='log10', limits = c(0.01,10^6)) +
  scale_colour_manual(values=cbPalette2) +
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Figures/Fig_1b.png", width = 8, height = 4, dpi = 300)


luz19_stats = phage_fitness_pfu_update_boxplots %>%
  filter(Phage == "\U03D5LUZ19")

p14_stats = phage_fitness_pfu_update_boxplots %>%
  filter(Phage == "\U03D5 14-1")

mod = lmer(log(pfu_ml) ~ line * temp + (1|rep), data = luz19_stats)
joint_tests(mod)
summary(mod)
plot(mod)
emmeans(mod, pairwise ~ line | temp )


mod = lmer(log(pfu_ml) ~ line * temp + (1|rep), p14_stats)
joint_tests(mod)

summary(mod)
anova(mod)
plot(mod)
emmeans(mod, pairwise ~ line | temp)



# Figure 2

# Figure 2A

pfu <- read_csv("Data/PFU_values_each_passage.csv")
pfu$phage_measured <- factor(pfu$phage_measured, levels = c("phage14_1", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))
pfu$phage_treatment <- factor(pfu$phage_treatment, levels = c("Single_phage", "Competition"), labels = c("Fluctuating (Monoculture)", "Fluctuating (Co-culture)"))

pfu <- pfu %>%
  filter(passage %% 1 == 0)

x_min <- floor(min(pfu$passage, na.rm = TRUE))
x_max <- ceiling(max(pfu$passage, na.rm = TRUE))
rects <- tibble(
  xmin = seq(x_min, x_max - 1, by = 1),
  xmax = xmin + 1,
  fill = rep(c("tomato", "lightblue"), length.out = length(xmin))
)

ggplot(pfu, aes(passage, pfu_ml, group = interaction(phage_treatment, replicate), 
                color = phage_treatment)) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = 1e7, ymax = 1e11, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_line(linewidth = 0.2, linetype = "dashed") +
  geom_line(aes(group = phage_treatment), 
            stat = "summary", fun = mean, linewidth = 1.2, alpha = 0.5) +
  facet_wrap(~phage_measured) +  
  scale_y_continuous(trans = 'log10', limits = c(1e7, 1e11)) +
  annotation_logticks(sides = "l") +
  scale_color_manual(values = c("Fluctuating (Monoculture)" = "#6A3D9A", "Fluctuating (Co-culture)" = "#CAB2D6")) +
  labs(x = "Passage", y = "Phage density (PFU/ml)", colour = "Evolution treatment") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )

ggsave("Figures/Fig_2a.png", width = 8, height = 3, dpi = 300)



# Figure 2B

phage_fitness_pfu <- read_csv("Data/Phage_population_growth_rates.csv")

phage_fitness_pfu_update <- phage_fitness_pfu %>% 
  dplyr::mutate(competition = ifelse(grepl("anc", line), "anc",
                                     ifelse(grepl("C", line), "Competition", "No_competition"))) %>%  # Label lines which evolved with competitors
  dplyr::mutate(evolution_treatment = case_when(   # Label lines based on their temperature evolution treatment
    grepl("anc", line) ~ "anc",
    grepl("37", line) ~ "37_static",
    grepl("42", line) ~ "42_static",
    grepl("F", line) ~ "Fluctuating",
    TRUE ~ NA_character_
  )) %>%
  dplyr::filter(!(phage == "phage14" & time == "2"))

phage_fitness_pfu_update$temperature <-  factor(phage_fitness_pfu_update$temperature, levels = c("37","42"), labels = c("37\u00B0C", "42\u00B0C"))
phage_fitness_pfu_update$line <- factor(phage_fitness_pfu_update$line, levels = c("anc","37_evo","37_C_evo", "42_evo", "42_C_evo", "F_evo", "F_C_evo"), labels = c("Ancestor", "37\u00B0C (mono.)", "37\u00B0C (comp.)", "42\u00B0C (mono.)", "42\u00B0C (comp.)","Fluctuating (mono.)", "Fluctuating (comp.)"))
phage_fitness_pfu_update$phage <- factor(phage_fitness_pfu_update$phage, levels = c("phage14", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))

phage_fitness_pfu_update_comp <- phage_fitness_pfu_update %>%
  dplyr::select(-c(dilution, pfu)) %>%
  dplyr::group_by(across(-c(batch, pfu_ml))) %>%  # Group by all columns except `value_column`
  dplyr::summarize(pfu_ml = mean(pfu_ml, na.rm = TRUE)) %>%  # Average duplicates in `value_column`
  ungroup()

phage_fitness_pfu_update_boxplots = phage_fitness_pfu_update_comp %>%
  dplyr::filter((phage == "ϕ 14-1" & time == "4") |
                  (phage == "ϕLUZ19" & time == "2"))   # Remove T2 values from 14-1 dataset as only carried out for a single rep


phage_fitness_pfu_update_boxplots_normalised <- phage_fitness_pfu_update_boxplots %>%
  dplyr::group_by(phage, temperature) %>%
  dplyr::mutate(
    ancestor_mean = mean(pfu_ml[line == "Ancestor"], na.rm = TRUE),
    PFU_relative = pfu_ml / ancestor_mean
  ) %>%
  ungroup()

phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::mutate(phage_temp = paste(phage, temperature, sep = "_"))

phage_fitness_pfu_update_boxplots_normalised$phage_temp = factor(phage_fitness_pfu_update_boxplots_normalised$phage_temp, levels = c("ϕ 14-1_37°C", "ϕ 14-1_42°C", "ϕLUZ19_37°C", "ϕLUZ19_42°C"), labels = c("ϕ14-1 (37°C)", "ϕ14-1 (42°C)", "ϕLUZ19 (37°C)", "ϕLUZ19 (42°C)"))


ancestor_se <- phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::filter(line == "Ancestor") %>%
  dplyr::group_by(phage_temp) %>%
  dplyr::summarise(se = sd(PFU_relative, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  dplyr::mutate(
    ymin = 1 - se,
    ymax = 1 + se,
    xmin = -Inf,
    xmax = Inf
  )


phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::filter(line != "Ancestor")


phage_fitness_pfu_update_boxplots_normalised$Legend_group1 <- dplyr::case_when(
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (mono.)" ~ "37°C\nStatic\nMono.",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (mono.)" ~ "42°C\nStatic\nMono.",
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (comp.)" ~ "37°C\nStatic\nCo.",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (comp.)" ~ "42°C\nStatic\nCo.",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (mono.)" ~ "Fluc.",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (comp.)" ~ "Fluc.\nCo.",
)
phage_fitness_pfu_update_boxplots_normalised$Legend_group2 <- dplyr::case_when(
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (mono.)" ~ "37°C Static (Monoculture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (mono.)" ~ "42°C Static (Monoculture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (comp.)" ~ "37°C Static (Co-culture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (comp.)" ~ "42°C Static (Co-culture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (mono.)" ~ "Fluctuating (Monoculture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (comp.)" ~ "Fluctuating (Co-culture)",
)


label_map <- c(
  "ϕ14-1 (37°C)" = "ϕ14-1",
  "ϕ14-1 (42°C)" = "ϕ14-1",
  "ϕLUZ19 (37°C)" = "ϕLUZ19",
  "ϕLUZ19 (42°C)" = "ϕLUZ19"
)

my_labeller <- as_labeller(label_map)

phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  filter(
    !(
      (temperature == "37°C" & evolution_treatment == "42_static" |
         (temperature == "42°C" & evolution_treatment == "37_static"
         ))))


phage_fitness_pfu_update_boxplots_normalised$Legend_group1 = factor(phage_fitness_pfu_update_boxplots_normalised$Legend_group1, levels = c("37°C\nStatic\nMono.", "37°C\nStatic\nCo.", "42°C\nStatic\nMono.", "42°C\nStatic\nCo.", "Fluc.", "Fluc.\nCo."))

ggplot(phage_fitness_pfu_update_boxplots_normalised, aes(Legend_group1, PFU_relative, color = as.factor(Legend_group2))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap2(~ phage_temp, nrow = 1, scales = "free_x", strip.position = "top", strip = strip_themed(
    background_x = list(
      A = element_rect(fill = "lightblue"),
      B = element_rect(fill = "tomato"),
      C = element_rect(fill = "lightblue"),
      D = element_rect(fill = "tomato")
    ),
    
    text_x = list(
      A = element_text(color = "black", face = "bold"),
      B = element_text(color = "black", face = "bold"),
      C = element_text(color = "black", face = "bold"),
      D = element_text(color = "black", face = "bold")
    )
  ),
  labeller = my_labeller 
  ) +
  scale_color_manual(values = c("37°C Static (Monoculture)" = "#1F78B4", 
  "37°C Static (Co-culture)" = "#A6CEE3", 
  "42°C Static (Monoculture)" = "#E31A1C", 
  "42°C Static (Co-culture)" = "#FB9A99",
  "Fluctuating (Monoculture)" = "#6A3D9A", 
  "Fluctuating (Co-culture)" = "#CAB2D6")) +

  ylab("Phage growth rel. to ancestor") +
  xlab("Evolution treatment")+
  labs(colour="Evolution treatment")+
  scale_y_continuous(trans='log10', limits = c(0.1,5*10^5)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "black")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "black", alpha = 0.2, inherit.aes = FALSE)+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )


ggsave("Figures/Fig_2b.png", width = 8, height = 4, dpi = 300)


luz19_stats = phage_fitness_pfu_update_boxplots %>%
  filter(Phage == "\U03D5LUZ19")

p14_stats = phage_fitness_pfu_update_boxplots %>%
  filter(Phage == "\U03D5 14-1")

mod = lmer(log(pfu_ml) ~ line * temp + (1|rep), luz19_stats)
joint_tests(mod)

summary(mod)
plot(mod)
emmeans(mod, pairwise ~ line * temp)

mod = lmer(log(pfu_ml) ~ line * temp + (1|rep), p14_stats)
joint_tests(mod)

summary(mod)
plot(mod)
emmeans(mod, pairwise ~ line * temp)





## Figure 3

# Figure 3A

# Euclidean distance trees

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

# identify ancestral mutations (any mutation present in any "anc" phage)
ancestor_mutations <- breseq_14 %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# drop any rows carrying those ancestral mutations
breseq_14_filtered <- breseq_14 %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))

# cast to factors if needed (not strictly required for distance)
breseq_14_filtered <- breseq_14_filtered %>%
  mutate(
    position = as.factor(position),
    mutation = as.factor(mutation),
    phage = as.factor(phage)
  )

breseq_14_filtered <- breseq_14_filtered %>%
  filter(frequency > 0.1)

# reshape wide: mutation identifier
df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  pivot_wider(
    names_from = mutation_id,
    values_from = frequency,
    values_fill = 0
  )

# ensure ancestor row exists with zeros (ancestor has no SNPs after filtering)
ancestor_name <- "14-1-anc"
if (!ancestor_name %in% df_wide$phage) {
  mut_cols <- setdiff(colnames(df_wide), "phage")
  zero_row <- tibble(phage = ancestor_name)
  zero_row[, mut_cols] <- 0
  df_wide <- bind_rows(df_wide, zero_row)
}

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

remove_labels <- grepl("C", rownames(distance_matrix)) |grepl("C", colnames(distance_matrix))
phage14_distance_matrix_nc <- distance_matrix[!remove_labels, !remove_labels]

tree_nc_14 <- nj(phage14_distance_matrix_nc)
tree_rooted_nc_14 <- root(tree_nc_14, outgroup = "14-1-anc", resolve.root = TRUE)

# Plot the tree

metadata <- data.frame(
  label = tree_rooted_nc_14$tip.label)  # Assign groups

metadata <- metadata %>%
  mutate(group = sub("-[^-]+$", "", label))

metadata$group <- factor(metadata$group, levels = c("14-1", "14-1-37", "14-1-F", "14-1-42"), labels = c("Ancestor", "37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))

cbPalette4 <- c("black", "#377EB8", "#6A3D9A", "#E41A1C")

phage14_eucl_tree = ggtree(tree_rooted_nc_14, layout = "rectangular") %>%
  ggtree::`%<+%`(metadata) +  # Attach metadata
  geom_tiplab(aes(label = group, color = group), size = 3, show.legend=FALSE, offset = 0.001) +  # Add tip labels colored by group
  geom_tippoint(aes(color = group), size = 0) +
  scale_colour_manual(values=cbPalette4) +
  theme_tree2() +
  xlim(0, 2.5)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  labs(color = "Evolved line")+
  xlab("Euclidean genetic distance")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=10,face="bold"))+
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6)))


## LUZ19

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

# identify ancestral mutations (any mutation present in any "anc" phage)
ancestor_mutations <- breseq_luz19 %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# drop any rows carrying those ancestral mutations
breseq_luz19_filtered <- breseq_luz19 %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))

# cast to factors if needed (not strictly required for distance)
breseq_luz19_filtered <- breseq_luz19_filtered %>%
  mutate(
    position = as.factor(position),
    mutation = as.factor(mutation),
    phage = as.factor(phage)
  )

breseq_luz19_filtered <- breseq_luz19_filtered %>%
  filter(frequency > 0.1)

# reshape wide: mutation identifier
df_wide <- breseq_luz19_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  pivot_wider(
    names_from = mutation_id,
    values_from = frequency,
    values_fill = 0
  )

# ensure ancestor row exists with zeros (ancestor has no SNPs after filtering)
ancestor_name <- "LUZ19-anc"
if (!ancestor_name %in% df_wide$phage) {
  mut_cols <- setdiff(colnames(df_wide), "phage")
  zero_row <- tibble(phage = ancestor_name)
  zero_row[, mut_cols] <- 0
  df_wide <- bind_rows(df_wide, zero_row)
}

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

remove_labels <- grepl("C", rownames(distance_matrix)) |grepl("C", colnames(distance_matrix))
luz19_distance_matrix_nc <- distance_matrix[!remove_labels, !remove_labels]

tree_nc_L <- nj(luz19_distance_matrix_nc)
tree_rooted_nc_L <- root(tree_nc_L, outgroup = "LUZ19-anc", resolve.root = TRUE)

# Plot the tree

metadata <- data.frame(
  label = tree_rooted_nc_L$tip.label)  # Assign groups

metadata <- metadata %>%
  mutate(group = sub("-[^-]+$", "", label))

metadata$group <- factor(metadata$group, levels = c("LUZ19", "LUZ19-37", "LUZ19-F","LUZ19-42"), labels = c("Ancestor", "37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))


# Convert to a ggtree object and join metadata
cbPalette4 <- c("black", "#377EB8", "#6A3D9A", "#E41A1C")

LUZ19_eucl_tree = ggtree(tree_rooted_nc_L, layout = "rectangular") %>%
  ggtree::`%<+%`(metadata) +  # Attach metadata
  geom_tiplab(aes(label = group, color = group), size = 3, show.legend=FALSE, offset = 0.001) +  # Add tip labels colored by group
  geom_tippoint(aes(color = group), size = 0) +
  scale_colour_manual(values=cbPalette4) +
  theme_tree2() +
  xlim(0, 2.5)+
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())+
  labs(color = "Evolved line")+
  xlab("Euclidean genetic distance")+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=10,face="bold"))+
  guides(color = guide_legend(override.aes = list(shape = 15, size = 6)))

phage14_eucl_tree + LUZ19_eucl_tree + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figures/Fig_3a.png", width = 8, height = 4, dpi = 300)




# Figure 3B

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.1)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

breseq_14_filtered = breseq_14_filtered %>%
  filter(!grepl("C", phage))

# Anc comparison

df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
phage14_distance_matrix_nc_anc <- as.matrix(distance_matrix)
rownames(phage14_distance_matrix_nc_anc) <- df_wide$phage
colnames(phage14_distance_matrix_nc_anc) <- df_wide$phage 

distance_long <- as.data.frame(phage14_distance_matrix_nc_anc) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  dplyr::mutate(evol_group_first = sub("-[^-]+$", "", phage1),
                evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  dplyr::mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  dplyr::arrange(pair, other) %>%
  dplyr::distinct(pair, other, .keep_all = TRUE) %>%
  dplyr::select(-pair, -other) %>%
  dplyr::filter(phage2 == "14-1-anc-2" | phage1 == "14-1-anc-2") %>%
  dplyr::filter(!(phage2 == "14-1-anc-2" & phage1 == "14-1-anc-2"))

distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "14-1-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "14-1-37"           ~ "Static\n37\u00B0C",
    evol_group  == "14-1-42"           ~ "Static\n42\u00B0C",
    evol_group == "14-1-37-C"          ~ "Static\n37\u00B0C\n(Co-culture)",
    evol_group  == "14-1-42-C"         ~ "Static\n42\u00B0C\n(Co-culture)",
    evol_group == "14-1-F"             ~ "Fluctuating",
    evol_group == "14-1-F-C"           ~ "Fluctuating\n(Co-culture)",
    TRUE                                       ~ NA_character_
  ))


distance_long$evol_type = factor(distance_long$evol_type, levels = c("Static\n37\u00B0C", "Fluctuating", "Static\n42\u00B0C"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))

cbPalette3 <- c("#1f78b4", "#6A3D9A", "#e31a1c")

phage14_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  labs(colour="Evolution treatment")+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"))

mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)
summary(mod)

emmeans(mod, pairwise ~ evol_type)

distance_long$dev <- with(distance_long, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)



# LUZ19 

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.1)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

breseq_luz19_filtered = breseq_luz19_filtered %>%
  filter(!grepl("C", phage))

# Anc comparison

df_wide <- breseq_luz19_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

distance_long <- as.data.frame(distance_matrix) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  mutate(evol_group_first = sub("-[^-]+$", "", phage1),
         evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  arrange(pair, other) %>%
  distinct(pair, other, .keep_all = TRUE) %>%
  select(-pair, -other) %>%
  filter(phage2 == "LUZ19-anc-3" | phage1 == "LUZ19-anc-3") %>%
  filter(!(phage2 == "LUZ19-anc-3" & phage1 == "LUZ19-anc-3"))


distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "LUZ19-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "LUZ19-37"           ~ "Static\n37\u00B0C",
    evol_group  == "LUZ19-42"           ~ "Static\n42\u00B0C",
    evol_group == "LUZ19-37-C"          ~ "Static\n37\u00B0C\n(Co-culture)",
    evol_group  == "LUZ19-42-C"         ~ "Static\n42\u00B0C\n(Co-culture)",
    evol_group == "LUZ19-F"             ~ "Fluctuating",
    evol_group == "LUZ19-F-C"           ~ "Fluctuating\n(Co-culture)",
    TRUE                                       ~ NA_character_
  ))

distance_long$evol_type = factor(distance_long$evol_type, levels = c("Static\n37\u00B0C", "Fluctuating", "Static\n42\u00B0C"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))

luz19_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  labs(colour="Evolution treatment")+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"), axis.title.y = element_blank())

phage14_dist_from_anc + luz19_dist_from_anc + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figures/Fig_3b.png", width = 8, height = 3, dpi = 300)


mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)

emmeans(mod, pairwise ~ evol_type)


leveneTest(distance ~ evol_type, data = distance_long)

distance_long$dev <- with(distance_long, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)




# Figure 3C

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_14_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_14_filtered <- breseq_14_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation")) %>%
  filter(!grepl("C", phage))



df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  mutate(present = ifelse(frequency > 0, 1, 0)) %>%  # presence/absence
  select(-frequency) %>%
  spread(key = mutation_id, value = present, fill = 0)

df_wide = df_wide %>%
  mutate(phage_type = sub("^(([^-]*-){2}[^-]*).*", "\\1", phage)) %>%
  group_by(phage_type) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


# Look at unique v shared without accounting for which groups are sharing
mutation_summary <- df_wide %>%
  pivot_longer(
    cols = -phage_type, 
    names_to = "mutation_id", 
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%   # only keep present mutations
  group_by(mutation_id) %>%
  summarise(
    phage_types = paste(unique(phage_type), collapse = ", "),
    n_phage_types = n_distinct(phage_type),
    .groups = "drop"
  ) %>%
  mutate(category = ifelse(n_phage_types == 1, "unique", "shared"))

p14_mutation_counts <- mutation_summary %>%
  separate_rows(phage_types, sep = ", ") %>%
  group_by(phage_types, category) %>%
  summarise(n_mutations = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = category,
    values_from = n_mutations,
    values_fill = 0
  ) %>%
  rename(Evolved_line = phage_types)

p14_mutation_counts

p14_mutation_counts_long <- p14_mutation_counts %>%
  pivot_longer(cols = c(unique, shared), names_to = "Mutation", values_to = "count")

p14_tbl <- xtabs(count ~ Evolved_line + Mutation, data = p14_mutation_counts_long)
chisq.test(p14_tbl)
fisher.test(p14_tbl)


p14_mutation_counts_long$Evolved_line = factor(p14_mutation_counts_long$Evolved_line, levels = c("14-1-37", "14-1-F", "14-1-42"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))
p14_mutation_counts_long$Mutation = factor(p14_mutation_counts_long$Mutation, levels = c("unique", "shared"), labels = c("Unique", "Shared"))

cbPalette4 <- c("#377EB8", "#6A3D9A", "#E41A1C")

phage14_mutation_counts = ggplot(p14_mutation_counts_long, aes(x = Evolved_line, y = count, fill = Evolved_line, pattern = Mutation)) +
  labs(x = "Evolution treatment", y = "Number of mutations", fill = "Evolution treatment") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,25)+
  labs(colour="Evolved line")+
  geom_bar_pattern(
    stat = "identity",
    color = "black",             # bar outline
    pattern_type = "stripe",      # explicitly set pattern type
    pattern_fill = "black",      # color of hatch
    pattern_angle = 45,          # angle of lines
    pattern_density = 0.1,       # line density
    pattern_spacing = 0.05       # space between lines
  ) +
  scale_pattern_manual(values = c("Unique" = "none", "Shared" = "stripe")) +
  scale_fill_manual(values=cbPalette3) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"))



# LUZ19

breseq_L <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_L_filtered <- breseq_L %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_L_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_L_filtered <- breseq_L_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation")) %>%
  filter(!grepl("C", phage))



df_wide <- breseq_L_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  mutate(present = ifelse(frequency > 0, 1, 0)) %>%  # presence/absence
  select(-frequency) %>%
  spread(key = mutation_id, value = present, fill = 0)

df_wide = df_wide %>%
  mutate(phage_type = sub("^(([^-]*-){1}[^-]*).*", "\\1", phage)) %>%
  group_by(phage_type) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


# Assess unique v shared without breaking down by who shares with who

mutation_summary <- df_wide %>%
  pivot_longer(
    cols = -phage_type, 
    names_to = "mutation_id", 
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%   # only keep present mutations
  group_by(mutation_id) %>%
  summarise(
    phage_types = paste(unique(phage_type), collapse = ", "),
    n_phage_types = n_distinct(phage_type),
    .groups = "drop"
  ) %>%
  mutate(category = ifelse(n_phage_types == 1, "unique", "shared"))

L_mutation_counts <- mutation_summary %>%
  separate_rows(phage_types, sep = ", ") %>%
  group_by(phage_types, category) %>%
  summarise(n_mutations = n(), .groups = "drop") %>%
  pivot_wider(
    names_from = category,
    values_from = n_mutations,
    values_fill = 0
  ) %>%
  rename(Evolved_line = phage_types)

L_mutation_counts

L_mutation_counts_long <- L_mutation_counts %>%
  pivot_longer(cols = c(unique, shared), names_to = "Mutation", values_to = "count")

L_tbl <- xtabs(count ~ Evolved_line + Mutation, data = L_mutation_counts_long)
chisq.test(L_tbl)
fisher.test(L_tbl)



L_mutation_counts_long$Evolved_line = factor(L_mutation_counts_long$Evolved_line, levels = c("LUZ19-37", "LUZ19-F", "LUZ19-42"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))
L_mutation_counts_long$Mutation = factor(L_mutation_counts_long$Mutation, levels = c("unique", "shared"), labels = c("Unique", "Shared"))

luz19_mutation_count = ggplot(L_mutation_counts_long, aes(x = Evolved_line, y = count, fill = Evolved_line, pattern = Mutation)) +
  labs(x = "Evolution treatment", y = "Number of mutations", fill = "Evolution treatment") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,25)+
  geom_bar_pattern(
    stat = "identity",
    color = "black",             # bar outline
    pattern_type = "stripe",      # explicitly set pattern type
    pattern_fill = "black",      # color of hatch
    pattern_angle = 45,          # angle of lines
    pattern_density = 0.1,       # line density
    pattern_spacing = 0.05       # space between lines
  ) +
  scale_pattern_manual(values = c("Unique" = "none", "Shared" = "stripe")) +
  scale_fill_manual(values=cbPalette3) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"), axis.title.y = element_blank())


phage14_mutation_counts + luz19_mutation_count + plot_layout(guides = "collect") & theme(legend.position = "bottom") &
  guides(
    pattern = guide_legend(override.aes = list(fill = "white", color = "black"), order = 2),
    fill = guide_legend(override.aes = list(pattern = "none"), order = 1),  # removes pattern from fill legend
  )

ggsave("Figures/Fig_3C.png", width = 8, height = 3, dpi = 300)


# Combined stats

all_snps = p14_mutation_counts_long %>%
  left_join(L_mutation_counts_long, by = c("Evolved_line", "Mutation")) %>%
  mutate(all_snps_count = rowSums(across(c(count.x, count.y)), na.rm = TRUE))

all_tbl <- xtabs(all_snps_count ~ Evolved_line + Mutation, data = all_snps)
chisq.test(all_tbl)
fisher.test(all_tbl)




# Figure 4

# Figure 4A

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.1)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

breseq_14_filtered = breseq_14_filtered %>%
  filter(!grepl("37", phage)) %>%
  filter(!grepl("42", phage))
  

# Anc comparison

df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
phage14_distance_matrix_nc_anc <- as.matrix(distance_matrix)
rownames(phage14_distance_matrix_nc_anc) <- df_wide$phage
colnames(phage14_distance_matrix_nc_anc) <- df_wide$phage 

distance_long <- as.data.frame(phage14_distance_matrix_nc_anc) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  dplyr::mutate(evol_group_first = sub("-[^-]+$", "", phage1),
                evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  dplyr::mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  dplyr::arrange(pair, other) %>%
  dplyr::distinct(pair, other, .keep_all = TRUE) %>%
  dplyr::select(-pair, -other) %>%
  dplyr::filter(phage2 == "14-1-anc-2" | phage1 == "14-1-anc-2") %>%
  dplyr::filter(!(phage2 == "14-1-anc-2" & phage1 == "14-1-anc-2"))

distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "14-1-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "14-1-F"             ~ "Fluctuating\nMono.",
    evol_group == "14-1-F-C"           ~ "Fluctuating\nCo.",
    TRUE                                       ~ NA_character_
  ))



distance_long$evol_type = factor(distance_long$evol_type, levels = c("Fluctuating\nMono.", "Fluctuating\nCo."))

cbPalette2 <- c("#6A3D9A", "#CAB2D6")

phage14_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  labs(colour="Evolution treatment")+
  scale_colour_manual(values=cbPalette2) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"))

distance_long2 = distance_long %>%
  filter(phage1 != "14-1-F-C-2")


mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)
summary(mod)

emmeans(mod, pairwise ~ evol_type)

distance_long2$dev <- with(distance_long2, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long2)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)



# LUZ19 

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.1)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

breseq_luz19_filtered = breseq_luz19_filtered %>%
  filter(!grepl("37", phage)) %>%
  filter(!grepl("42", phage))

# Anc comparison

df_wide <- breseq_luz19_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

distance_long <- as.data.frame(distance_matrix) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  mutate(evol_group_first = sub("-[^-]+$", "", phage1),
         evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  arrange(pair, other) %>%
  distinct(pair, other, .keep_all = TRUE) %>%
  select(-pair, -other) %>%
  filter(phage2 == "LUZ19-anc-3" | phage1 == "LUZ19-anc-3") %>%
  filter(!(phage2 == "LUZ19-anc-3" & phage1 == "LUZ19-anc-3"))


distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "LUZ19-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "LUZ19-F"             ~ "Fluctuating\nMono.",
    evol_group == "LUZ19-F-C"           ~ "Fluctuating\nCo.",
    TRUE                                       ~ NA_character_
  ))

distance_long$evol_type = factor(distance_long$evol_type, levels = c("Fluctuating\nMono.", "Fluctuating\nCo."))


luz19_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  scale_colour_manual(values=cbPalette2) +
  labs(colour="Evolution treatment")+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"), axis.title.y = element_blank())


phage14_dist_from_anc + luz19_dist_from_anc + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figures/Fig_4A.png", width = 4, height = 3, dpi = 300)


mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)

emmeans(mod, pairwise ~ evol_type)


leveneTest(distance ~ evol_type, data = distance_long)

distance_long$dev <- with(distance_long, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)




# Figure 4B

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.1)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

breseq_14_filtered = breseq_14_filtered %>%
  filter(grepl("anc", phage) | grepl("C", phage))


# Anc comparison

df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
phage14_distance_matrix_nc_anc <- as.matrix(distance_matrix)
rownames(phage14_distance_matrix_nc_anc) <- df_wide$phage
colnames(phage14_distance_matrix_nc_anc) <- df_wide$phage 

distance_long <- as.data.frame(phage14_distance_matrix_nc_anc) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  dplyr::mutate(evol_group_first = sub("-[^-]+$", "", phage1),
                evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  dplyr::mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  dplyr::arrange(pair, other) %>%
  dplyr::distinct(pair, other, .keep_all = TRUE) %>%
  dplyr::select(-pair, -other) %>%
  dplyr::filter(phage2 == "14-1-anc-2" | phage1 == "14-1-anc-2") %>%
  dplyr::filter(!(phage2 == "14-1-anc-2" & phage1 == "14-1-anc-2"))

distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "14-1-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "14-1-37-C"          ~ "37\u00B0C (static)\nCo.",
    evol_group  == "14-1-42-C"         ~ "42\u00B0C (static)\nCo.",
    evol_group == "14-1-F-C"           ~ "Fluctuating\nCo.",
    TRUE                                       ~ NA_character_
  ))


distance_long$evol_type = factor(distance_long$evol_type, levels = c("37\u00B0C (static)\nCo.", "Fluctuating\nCo.", "42\u00B0C (static)\nCo."))

cbPalette3 <- c("#A6CEE3", "#CAB2D6", "#FB9A99")

  
phage14_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  labs(colour="Evolution treatment")+
  scale_colour_manual(values=cbPalette3) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"))

mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)
summary(mod)

emmeans(mod, pairwise ~ evol_type)

leveneTest(distance ~ evol_type, data = distance_long)

distance_long$dev <- with(distance_long, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)



# LUZ19 

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.1)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

breseq_luz19_filtered = breseq_luz19_filtered %>%
  filter(grepl("anc", phage) | grepl("C", phage))


# Anc comparison

df_wide <- breseq_luz19_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

distance_long <- as.data.frame(distance_matrix) %>%
  rownames_to_column(var = "phage1") %>%
  pivot_longer(-phage1, names_to = "phage2", values_to = "distance") %>%
  mutate(evol_group_first = sub("-[^-]+$", "", phage1),
         evol_group_second = sub("-[^-]+$", "", phage2)) %>%
  mutate(
    pair = pmin(phage1, phage2),
    other = pmax(phage1, phage2)
  ) %>%
  arrange(pair, other) %>%
  distinct(pair, other, .keep_all = TRUE) %>%
  select(-pair, -other) %>%
  filter(phage2 == "LUZ19-anc-3" | phage1 == "LUZ19-anc-3") %>%
  filter(!(phage2 == "LUZ19-anc-3" & phage1 == "LUZ19-anc-3"))


distance_long <- distance_long %>%
  mutate(
    phage1_new = ifelse(str_detect(phage1, "anc"), phage2, phage1),
    phage2_new = ifelse(str_detect(phage1, "anc"), phage1, phage2)
  )%>%
  select(phage1 = phage1_new, phage2 = phage2_new, distance) %>%
  mutate(evol_group = sub("-[^-]+$", "", phage1)) %>%
  filter(evol_group != "LUZ19-anc") %>%
  mutate(evol_type = case_when(
    evol_group == "LUZ19-37-C"          ~ "37\u00B0C (static)\nCo.",
    evol_group  == "LUZ19-42-C"         ~ "42\u00B0C (static)\nCo.",
    evol_group == "LUZ19-F-C"           ~ "Fluctuating\nCo.",
    TRUE                                       ~ NA_character_
  ))


distance_long$evol_type = factor(distance_long$evol_type, levels = c("37\u00B0C (static)\nCo.", "Fluctuating\nCo.", "42\u00B0C (static)\nCo."))

luz19_dist_from_anc = ggplot(distance_long, aes(evol_type, distance, color = evol_type)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_colour_manual(values=cbPalette3) +
  
  ylab("Genetic distance\n from ancestor") +
  xlab("Evolution treatment")+
  ylim(0,2.7)+
  scale_colour_manual(values=cbPalette3) +
  labs(colour="Evolution treatment")+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"), axis.title.y = element_blank())


phage14_dist_from_anc + luz19_dist_from_anc + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figures/Fig_4B.png", width = 6, height = 3, dpi = 300)


mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)

emmeans(mod, pairwise ~ evol_type)


leveneTest(distance ~ evol_type, data = distance_long)

distance_long$dev <- with(distance_long, abs(distance - ave(distance, evol_type, FUN = median)))

# Step 2: run ANOVA on these deviations
levene_aov <- aov(dev ~ evol_type, data = distance_long)

# Step 3: extract test result
summary(levene_aov)
emmeans(levene_aov, pairwise ~ evol_type)





# Supplementary figures

# Figure S1

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_14_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_14_filtered <- breseq_14_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation")) %>%
  filter(!grepl("C", phage))


df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  mutate(present = ifelse(frequency > 0, 1, 0)) %>%  # presence/absence
  select(-frequency) %>%
  spread(key = mutation_id, value = present, fill = 0)

df_wide = df_wide %>%
  mutate(phage_type = sub("^(([^-]*-){2}[^-]*).*", "\\1", phage)) %>%
  group_by(phage_type) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


mutation_summary <- df_wide %>%
  pivot_longer(
    cols = -phage_type, 
    names_to = "mutation_id", 
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  group_by(mutation_id) %>%
  summarise(
    phage_types = paste(unique(phage_type), collapse = ", "),
    n_phage_types = n_distinct(phage_type),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(category = {
    # split phage_types string into a vector
    types <- unlist(strsplit(phage_types, ",\\s*"))
    
    # map to labels
    mapped <- c()
    if ("14-1-37" %in% types) mapped <- c(mapped, "37\u00B0C")
    if ("14-1-42" %in% types) mapped <- c(mapped, "42\u00B0C")
    if ("14-1-F"  %in% types) mapped <- c(mapped, "Fluc.")
    
    # decide output
    if (length(mapped) == 1) {
      paste0(mapped, " alone")
    } else if (length(mapped) == 3) {
      "Shared by all"
    } else {
      paste(mapped, collapse = " + ")
    }
  }) %>%
  ungroup()


p14_mutation_counts <- mutation_summary %>%
  group_by(category) %>%
  summarise(n_mutations = n(), .groups = "drop") 

p14_mutation_counts <- rbind(p14_mutation_counts, data.frame(category = "37\u00B0C + Fluc.", n_mutations = 0))

p14_mutation_counts$category = factor(p14_mutation_counts$category, levels = c("37\u00B0C alone", "42\u00B0C alone", "37\u00B0C + 42\u00B0C", "Fluc. alone", "37\u00B0C + Fluc.", "42\u00B0C + Fluc.", "Shared by all"), labels = c("37\u00B0C\nalone", "42\u00B0C\nalone", "37\u00B0C +\n42\u00B0C", "Fluc.\nalone", "37\u00B0C +\nFluc.", "42\u00B0C +\nFluc.", "Shared\nby all"))

phage14_mutation_breakdown = ggplot(p14_mutation_counts, aes(x = category, y = n_mutations)) +
  labs(x = "Evolution treatment", y = "Number of mutations") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,15)+
  labs(colour="Evolved line")+
  geom_col(fill = "lightgrey", col = "black")+
  #scale_fill_manual(values=cbPalette3) +
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"))


# LUZ19


breseq_L <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_L_filtered <- breseq_L %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_L_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_L_filtered <- breseq_L_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation")) %>%
  filter(!grepl("C", phage))


df_wide <- breseq_L_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%
  select(phage, mutation_id, frequency) %>%
  mutate(present = ifelse(frequency > 0, 1, 0)) %>%  # presence/absence
  select(-frequency) %>%
  spread(key = mutation_id, value = present, fill = 0)

df_wide = df_wide %>%
  mutate(phage_type = sub("^(([^-]*-){1}[^-]*).*", "\\1", phage)) %>%
  group_by(phage_type) %>%
  summarise(across(where(is.numeric), sum), .groups = "drop")


mutation_summary <- df_wide %>%
  pivot_longer(
    cols = -phage_type, 
    names_to = "mutation_id", 
    values_to = "presence"
  ) %>%
  filter(presence > 0) %>%
  group_by(mutation_id) %>%
  summarise(
    phage_types = paste(unique(phage_type), collapse = ", "),
    n_phage_types = n_distinct(phage_type),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(category = {
    # split phage_types string into a vector
    types <- unlist(strsplit(phage_types, ",\\s*"))
    
    # map to labels
    mapped <- c()
    if ("LUZ19-37" %in% types) mapped <- c(mapped, "37\u00B0C")
    if ("LUZ19-42" %in% types) mapped <- c(mapped, "42\u00B0C")
    if ("LUZ19-F"  %in% types) mapped <- c(mapped, "Fluc.")
    
    # decide output
    if (length(mapped) == 1) {
      paste0(mapped, " alone")
    } else if (length(mapped) == 3) {
      "Shared by all"
    } else {
      paste(mapped, collapse = " + ")
    }
  }) %>%
  ungroup()


LUZ19_mutation_counts <- mutation_summary %>%
  group_by(category) %>%
  summarise(n_mutations = n(), .groups = "drop") 

LUZ19_mutation_counts <- rbind(LUZ19_mutation_counts, data.frame(category = "37\u00B0C + Fluc.", n_mutations = 0))

LUZ19_mutation_counts$category = factor(LUZ19_mutation_counts$category, levels = c("37\u00B0C alone", "42\u00B0C alone", "37\u00B0C + 42\u00B0C", "Fluc. alone", "37\u00B0C + Fluc.", "42\u00B0C + Fluc.", "Shared by all"), labels = c("37\u00B0C\nalone", "42\u00B0C\nalone", "37\u00B0C +\n42\u00B0C", "Fluc.\nalone", "37\u00B0C +\nFluc.", "42\u00B0C +\nFluc.", "Shared\nby all"))

LUZ19_mutation_breakdown = ggplot(LUZ19_mutation_counts, aes(x = category, y = n_mutations)) +
  labs(x = "Evolution treatment", y = "Number of mutations") +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  ylim(0,15)+
  geom_col(fill = "lightgrey", col = "black")+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"),
        axis.title.y = element_blank())


phage14_mutation_breakdown + LUZ19_mutation_breakdown + plot_layout(guides = "collect")

ggsave("Figures/Fig_S1.png", width = 8, height = 3, dpi = 300)



# Figure S2

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

# identify ancestral mutations (any mutation present in any "anc" phage)
ancestor_mutations <- breseq_14_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# drop any rows carrying those ancestral mutations
breseq_14_filtered <- breseq_14_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))


df <- read_tsv("Data/14_1_prokka.gff", skip = 2,col_names = FALSE)

df_subset <- df %>%
  dplyr::select(4,5,9) %>%
  dplyr::rename(
    start = X4,
    end = X5,
    annotation = X9
  ) %>%
  mutate(annotation = str_extract(annotation, "(?<=product=).*")) %>%
  drop_na()


# Prepare phage lines and allele frequencies
phage_data <- breseq_14_filtered %>%
  dplyr::select(phage, position, type, mutation, frequency, evol_group,cleaned_function) %>%
  mutate(
    mutation_length = case_when(
      str_detect(mutation, "\\+") ~ nchar(str_remove(mutation, "\\+")),  # Remove "+" and count remaining chars
      str_detect(mutation, "bp") ~ as.numeric(str_extract(mutation, "\\d+")),  # Extract numeric part from "delta"
      TRUE ~ NA_real_  # Assign NA if neither condition is met
    )
  )

phage_data_F <- phage_data %>%
  filter(evol_group == "14-1-F") %>%
  mutate(phage = factor(phage, levels = c("14-1-F-1","14-1-F-2","14-1-F-3","14-1-F-4","14-1-F-5","14-1-F-6"))) %>%
  complete(phage, fill = list(value = NA))

phage_data_F$phage = factor(phage_data_F$phage, levels = c("14-1-F-1","14-1-F-2","14-1-F-3","14-1-F-4","14-1-F-5","14-1-F-6"))


phage_data_42 <- phage_data %>%
  filter(evol_group == "14-1-F-C")

phage_data_42$phage <- factor(phage_data_42$phage, levels = c("14-1-F-C-1","14-1-F-C-2","14-1-F-C-3","14-1-F-C-4","14-1-F-C-5","14-1-F-C-6"))


snp_F <-  ggplot() +
  geom_segment(data = phage_data_F, 
               aes(x = 0, xend = 66298, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Insertion points
  geom_point(data = phage_data_F %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F %>% filter(type == "SNP"),
             aes(x = position, y = phage, size = frequency, shape = "SNP"), 
             fill = "darkorange", color = "black", alpha = 0.6) +
  
  # Customize the legend for different aesthetics
  scale_size_continuous(limits = c(0.2, 1), name = "Frequency") +
  scale_shape_manual(values = c("Insertion" = 25, "SNP" = 21)) +  # Different shapes for INS and SNP
  scale_fill_manual(values = c("Deletion" = "blue")) +  # Different colors for deletions
  
  labs(x = "Genome Position", y = "Fluctuating\n(Monoculture)", size = "Allele Freq") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(legend.position = "none")+ 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

genes <- ggplot() +
  geom_gene_arrow(data = df_subset,
                  aes(xmin = start, xmax = end, y = as.factor(-1)),
                  arrowhead_height = unit(4, "mm"),
                  arrowhead_width = unit(2, "mm"))+
  theme_minimal()+
  labs(x = "Genome Position") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))



snp_F / genes  + plot_layout(guides = "collect", heights = c(1, 0.3))

ggsave("Figures/FigS2A.png", width = 4, height = 3, dpi = 300)




#LUZ19

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_luz19_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_luz19_filtered <- breseq_luz19_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))



df <- read_tsv("Data/luz19_prokka.gff", skip = 2,col_names = FALSE)

df_subset <- df %>%
  dplyr::select(4,5,9) %>%
  dplyr::rename(
    start = X4,
    end = X5,
    annotation = X9
  ) %>%
  mutate(annotation = str_extract(annotation, "(?<=product=).*")) %>%
  drop_na()


# Prepare phage lines and allele frequencies
phage_data <- breseq_luz19_filtered %>%
  dplyr::select(phage, position, type, mutation, frequency, evol_group,cleaned_function) %>%
  mutate(
    mutation_length = case_when(
      str_detect(mutation, "\\+") ~ nchar(str_remove(mutation, "\\+")),  # Remove "+" and count remaining chars
      str_detect(mutation, "bp") ~ as.numeric(str_extract(mutation, "\\d+")),  # Extract numeric part from "delta"
      TRUE ~ NA_real_  # Assign NA if neither condition is met
    )
  )

phage_data_F <- phage_data %>%
  filter(evol_group == "LUZ19-F") %>%
  mutate(phage = factor(phage, levels = c("LUZ19-F-1","LUZ19-F-2","LUZ19-F-3","LUZ19-F-4","LUZ19-F-5","LUZ19-F-6"))) %>%
  complete(phage, fill = list(value = NA))

phage_data_F$phage = factor(phage_data_F$phage, levels = c("LUZ19-F-1","LUZ19-F-2","LUZ19-F-3","LUZ19-F-4","LUZ19-F-5","LUZ19-F-6"))

snp_F <-  ggplot() +
  geom_segment(data = phage_data_F, 
               aes(x = 0, xend = 43187, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.6) +
  
  # Insertion points
  geom_point(data = phage_data_F %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F %>% filter(type == "SNP"),
             aes(x = position, y = phage, size = frequency, shape = "SNP"), 
             fill = "darkorange", color = "black", alpha = 0.6) +
  
  # Customize the legend for different aesthetics
  scale_size_continuous(limits = c(0.2, 1), name = "Frequency") +
  scale_shape_manual(values = c("Insertion" = 25, "SNP" = 21)) +  # Different shapes for INS and SNP
  scale_fill_manual(values = c("Deletion" = "blue")) +  # Different colors for deletions
  
  labs(x = "Genome Position", y = "Fluctuating\n(Monoculture)", size = "Allele Freq") +
  
  # Adjust the legends
  guides(
    fill = guide_legend(title = "Deletion Type"),
    shape = guide_legend(title = "Mutation Type"),
    size = guide_legend(title = "Frequency")
  ) +
  
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(legend.title = element_text(colour="black", size=11, face="bold"), 
        legend.text = element_text(size=10),
        legend.position = "none") +
  theme(axis.title = element_text(size=11, face="bold"))+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


genes <- ggplot() +
  geom_gene_arrow(data = df_subset,
                  aes(xmin = start, xmax = end, y = as.factor(-1)),
                  arrowhead_height = unit(4, "mm"),
                  arrowhead_width = unit(2, "mm"))+
  theme_minimal()+
  labs(x = "Genome Position") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))



snp_F / genes + plot_layout(guides = "collect", heights = c(1, 0.3))

ggsave("Figures/FigS2B.png", width = 4, height = 3, dpi = 300)




# Figure S3

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

# identify ancestral mutations (any mutation present in any "anc" phage)
ancestor_mutations <- breseq_14_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# drop any rows carrying those ancestral mutations
breseq_14_filtered <- breseq_14_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))


df <- read_tsv("Data/14_1_prokka.gff", skip = 2,col_names = FALSE)

df_subset <- df %>%
  dplyr::select(4,5,9) %>%
  dplyr::rename(
    start = X4,
    end = X5,
    annotation = X9
  ) %>%
  mutate(annotation = str_extract(annotation, "(?<=product=).*")) %>%
  drop_na()


# Prepare phage lines and allele frequencies
phage_data <- breseq_14_filtered %>%
  dplyr::select(phage, position, type, mutation, frequency, evol_group,cleaned_function) %>%
  mutate(
    mutation_length = case_when(
      str_detect(mutation, "\\+") ~ nchar(str_remove(mutation, "\\+")),  # Remove "+" and count remaining chars
      str_detect(mutation, "bp") ~ as.numeric(str_extract(mutation, "\\d+")),  # Extract numeric part from "delta"
      TRUE ~ NA_real_  # Assign NA if neither condition is met
    )
  )


phage_data_F <- phage_data %>%
  filter(evol_group == "14-1-F-C") %>%
  mutate(phage = factor(phage, levels = c("14-1-F-C-1","14-1-F-C-2","14-1-F-C-3","14-1-F-C-4","14-1-F-C-5","14-1-F-C-6"))) %>%
  complete(phage, fill = list(value = NA))

phage_data_F$phage <- factor(phage_data_F$phage, levels = c("14-1-F-C-1","14-1-F-C-2","14-1-F-C-3","14-1-F-C-4","14-1-F-C-5","14-1-F-C-6"))

snp_F <-  ggplot() +
  geom_segment(data = phage_data_F, 
               aes(x = 0, xend = 66298, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Insertion points
  geom_point(data = phage_data_F %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F %>% filter(type == "SNP"),
             aes(x = position, y = phage, size = frequency, shape = "SNP"), 
             fill = "darkorange", color = "black", alpha = 0.6) +
  
  # Customize the legend for different aesthetics
  scale_size_continuous(limits = c(0.2, 1), name = "Frequency") +
  scale_shape_manual(values = c("Insertion" = 25, "SNP" = 21)) +  # Different shapes for INS and SNP
  scale_fill_manual(values = c("Deletion" = "blue")) +  # Different colors for deletions
  
  labs(x = "Genome Position", y = "Fluctuating\n(Co-culture)", size = "Allele Freq") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(legend.position = "none")+ 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

genes <- ggplot() +
  geom_gene_arrow(data = df_subset,
                  aes(xmin = start, xmax = end, y = as.factor(-1)),
                  arrowhead_height = unit(4, "mm"),
                  arrowhead_width = unit(2, "mm"))+
  theme_minimal()+
  labs(x = "Genome Position") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))



snp_F / genes  + plot_layout(guides = "collect", heights = c(1, 0.3))

ggsave("Figures/FigS3A.png", width = 4, height = 3, dpi = 300)




#LUZ19

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency >= 0.2) %>%
  mutate(rep = substr(phage, nchar(phage), nchar(phage)),
         evol_group = sub("-[^-]+$", "", phage),
         competition_status = ifelse(grepl("C", phage), "competition", "no_competition"),
         cleaned_function = case_when(
           grepl("hypothetical", description, ignore.case = TRUE) ~ "hypothetical protein",
           grepl("IGINKLAI", description, ignore.case = TRUE) ~ "Multiple genes",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ description  # Keep original if no match
         ),
         cleaned_function_shortened = case_when(
           grepl("DNA", cleaned_function, ignore.case = TRUE) ~ "DNA replication",
           grepl("tail", cleaned_function, ignore.case = TRUE) ~ "tail-related protein",
           #grepl("IGINKLAI", description, ignore.case = TRUE) ~ gene,
           TRUE ~ cleaned_function  # Keep original if no match
         )
  )

ancestor_mutations <- breseq_luz19_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_luz19_filtered <- breseq_luz19_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))



df <- read_tsv("Data/luz19_prokka.gff", skip = 2,col_names = FALSE)

df_subset <- df %>%
  dplyr::select(4,5,9) %>%
  dplyr::rename(
    start = X4,
    end = X5,
    annotation = X9
  ) %>%
  mutate(annotation = str_extract(annotation, "(?<=product=).*")) %>%
  drop_na()


# Prepare phage lines and allele frequencies
phage_data <- breseq_luz19_filtered %>%
  dplyr::select(phage, position, type, mutation, frequency, evol_group,cleaned_function) %>%
  mutate(
    mutation_length = case_when(
      str_detect(mutation, "\\+") ~ nchar(str_remove(mutation, "\\+")),  # Remove "+" and count remaining chars
      str_detect(mutation, "bp") ~ as.numeric(str_extract(mutation, "\\d+")),  # Extract numeric part from "delta"
      TRUE ~ NA_real_  # Assign NA if neither condition is met
    )
  )

phage_data_F <- phage_data %>%
  filter(evol_group == "LUZ19-F-C") %>%
  mutate(phage = factor(phage, levels = c("LUZ19-F-C-1","LUZ19-F-C-2","LUZ19-F-C-3","LUZ19-F-C-4","LUZ19-F-C-5","LUZ19-F-C-6"))) %>%
  complete(phage, fill = list(value = NA))

phage_data_F$phage = factor(phage_data_F$phage, levels = c("LUZ19-F-C-1","LUZ19-F-C-2","LUZ19-F-C-3","LUZ19-F-C-4","LUZ19-F-C-5","LUZ19-F-C-6"))


snp_F <-  ggplot() +
  geom_segment(data = phage_data_F, 
               aes(x = 0, xend = 43187, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.6) +
  
  # Insertion points
  geom_point(data = phage_data_F %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F %>% filter(type == "SNP"),
             aes(x = position, y = phage, size = frequency, shape = "SNP"), 
             fill = "darkorange", color = "black", alpha = 0.6) +
  
  # Customize the legend for different aesthetics
  scale_size_continuous(limits = c(0.2, 1), name = "Frequency") +
  scale_shape_manual(values = c("Insertion" = 25, "SNP" = 21)) +  # Different shapes for INS and SNP
  scale_fill_manual(values = c("Deletion" = "blue")) +  # Different colors for deletions
  
  labs(x = "Genome Position", y = "Fluctuating\n(Co-culture)", size = "Allele Freq") +
  
  # Adjust the legends
  guides(
    fill = guide_legend(title = "Deletion Type"),
    shape = guide_legend(title = "Mutation Type"),
    size = guide_legend(title = "Frequency")
  ) +
  
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  theme(legend.title = element_text(colour="black", size=11, face="bold"), 
        legend.text = element_text(size=10),
        legend.position = "none") +
  theme(axis.title = element_text(size=11, face="bold"))+
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


genes <- ggplot() +
  geom_gene_arrow(data = df_subset,
                  aes(xmin = start, xmax = end, y = as.factor(-1)),
                  arrowhead_height = unit(4, "mm"),
                  arrowhead_width = unit(2, "mm"))+
  theme_minimal()+
  labs(x = "Genome Position") +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank())+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10))+
  theme(axis.title=element_text(size=11,face="bold"))



snp_F / genes + plot_layout(guides = "collect", heights = c(1, 0.3))

ggsave("Figures/FigS3B.png", width = 4, height = 3, dpi = 300)




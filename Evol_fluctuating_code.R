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


# Population dynamics Figure 1

pfu <- read_csv("Data/PFU_values_each_passage.csv")
pfu$Temperature_treatment <-  factor(pfu$Temperature_treatment, levels = c("37C_static","42C_static","Fluctuating"), labels = c("37\u00B0C", "42\u00B0C", "Fluctuating"))
pfu$Phage_measured <- factor(pfu$Phage_measured, levels = c("phage14_1", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))
pfu$Phage_treatment <- factor(pfu$Phage_treatment, levels = c("Single_phage", "Competition"), labels = c("Monoculture", "Co-culture"))

pfu <- pfu %>%
  filter(Temperature_treatment == "Fluctuating") %>%
  filter(Phage_treatment == "Monoculture") %>%
  filter(Passage %% 1 == 0)

x_min <- floor(min(pfu$Passage, na.rm = TRUE))
x_max <- ceiling(max(pfu$Passage, na.rm = TRUE))
rects <- tibble(
  xmin = seq(x_min, x_max - 1, by = 1),
  xmax = xmin + 1,
  fill = rep(c("tomato", "lightblue"), length.out = length(xmin))
)

ggplot(pfu, aes(Passage, PFU_ML, group = interaction(Phage_treatment, Rep), 
                color = Phage_treatment)) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = 1e7, ymax = 1e11, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_line(linewidth = 0.2, linetype = "dashed") +
  geom_line(aes(group = Phage_treatment), 
            stat = "summary", fun = mean, linewidth = 1.2, alpha = 0.5) +
  facet_wrap(~Phage_measured) +  
  scale_y_continuous(trans = 'log10', limits = c(1e7, 1e11)) +
  annotation_logticks(sides = "l") +
  scale_color_manual(values = c("Monoculture" = "#FF7F00", "Co-culture" = "grey10")) +
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

ggsave("Figures/Fig1a.png", width = 8, height = 3, dpi = 300)


## Fitness assay

phage_fitness_pfu <- read.csv("Data/All_phage_PFU.csv")

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
  filter(!(Phage == "phage14" & Time == "2"))

phage_fitness_pfu_update$temp <-  factor(phage_fitness_pfu_update$temp, levels = c("37","42"), labels = c("37\u00B0C", "42\u00B0C"))
phage_fitness_pfu_update$line <- factor(phage_fitness_pfu_update$line, levels = c("anc","37_evo","37_C_evo", "42_evo", "42_C_evo", "F_evo", "F_C_evo"), labels = c("Ancestor", "37\u00B0C static", "37\u00B0C static (comp)", "42\u00B0C static", "42\u00B0C static (comp)","Fluctuating", "Fluctuating (comp)"))
phage_fitness_pfu_update$Phage <- factor(phage_fitness_pfu_update$Phage, levels = c("phage14", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))

phage_fitness_pfu_update_single <- phage_fitness_pfu_update %>%
  filter(line %in% c("Ancestor", "37\u00B0C static", "42\u00B0C static", "Fluctuating"))

phage_fitness_pfu_update_single <- phage_fitness_pfu_update_single %>%
  select(-c(dil, pfu)) %>%
  group_by(across(-c(Batch, pfu_ml))) %>%  # Group by all columns except `value_column`
  dplyr::summarize(pfu_ml = mean(pfu_ml, na.rm = TRUE)) %>%  # Average duplicates in `value_column`
  ungroup()

phage_fitness_pfu_update_boxplots = phage_fitness_pfu_update_single %>%
  filter((Phage == "\U03D5 14-1" & Time == "4") |
           (Phage == "\U03D5LUZ19" & Time == "2"))   # Remove T2 values from 14-1 dataset as only carried out for a single rep


phage_fitness_pfu_update_single$Time = as.numeric(phage_fitness_pfu_update_single$Time)

phage_fitness_pfu_update_single = phage_fitness_pfu_update_single %>%
  mutate(Phage_temp = paste(Phage, temp, sep = "_"))

phage_fitness_pfu_update_single$Phage_temp = factor(phage_fitness_pfu_update_single$Phage_temp, levels = c("ϕ 14-1_37°C", "ϕ 14-1_42°C", "ϕLUZ19_37°C", "ϕLUZ19_42°C"), labels = c("ϕ14-1 (37°C)", "ϕ14-1 (42°C)", "ϕLUZ19 (37°C)", "ϕLUZ19 (42°C)"))

phage_fitness_pfu_update_boxplots_single = phage_fitness_pfu_update_boxplots %>%
  mutate(Phage_temp = paste(Phage, temp, sep = "_"))

unique(phage_fitness_pfu_update_boxplots_single$Phage_temp)

phage_fitness_pfu_update_boxplots_single$Phage_temp = factor(phage_fitness_pfu_update_boxplots_single$Phage_temp, levels = c("ϕ 14-1_37°C", "ϕ 14-1_42°C", "ϕLUZ19_37°C", "ϕLUZ19_42°C"), labels = c("ϕ14-1 (37°C)", "ϕ14-1 (42°C)", "ϕLUZ19 (37°C)", "ϕLUZ19 (42°C)"))


phage_fitness_pfu_update_boxplots_single_normalised <- phage_fitness_pfu_update_boxplots_single %>%
  dplyr::group_by(Phage, temp) %>%
  dplyr::mutate(
    ancestor_mean = mean(pfu_ml[line == "Ancestor"], na.rm = TRUE),
    PFU_relative = pfu_ml / ancestor_mean
  ) %>%
  ungroup()

ancestor_se <- phage_fitness_pfu_update_boxplots_single_normalised %>%
  dplyr::filter(line == "Ancestor") %>%
  dplyr::group_by(Phage_temp) %>%
  dplyr::summarise(se = sd(PFU_relative, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  dplyr::mutate(
    ymin = 1 - se,
    ymax = 1 + se,
    xmin = -Inf,
    xmax = Inf
  )

phage_fitness_pfu_update_boxplots_single_normalised = phage_fitness_pfu_update_boxplots_single_normalised %>%
  dplyr::filter(line != "Ancestor")

cbPalette2 <- c("#377EB8", "#FF7F00", "#E41A1C")

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
  facet_wrap2(~ Phage_temp, nrow = 1, strip.position = "top", strip = strip_themed(
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
  xlab("Evolved line")+
  labs(colour="Evolved line")+
  scale_y_continuous(trans='log10', limits = c(0.01,10^6)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE)+
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

ggsave("Figures/Fig1b.png", width = 8, height = 4, dpi = 300)


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



# Figure 2 - monoculture v co-culture
pfu <- read_csv("Data/PFU_values_each_passage.csv")
pfu$Temperature_treatment <-  factor(pfu$Temperature_treatment, levels = c("37C_static","42C_static","Fluctuating"), labels = c("37\u00B0C", "42\u00B0C", "Fluctuating"))
pfu$Phage_measured <- factor(pfu$Phage_measured, levels = c("phage14_1", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))
pfu$Phage_treatment <- factor(pfu$Phage_treatment, levels = c("Single_phage", "Competition"), labels = c("Monoculture", "Co-culture"))

pfu <- pfu %>%
  filter(Temperature_treatment == "Fluctuating") %>%
  filter(Passage %% 1 == 0)

x_min <- floor(min(pfu$Passage, na.rm = TRUE))
x_max <- ceiling(max(pfu$Passage, na.rm = TRUE))
rects <- tibble(
  xmin = seq(x_min, x_max - 1, by = 1),
  xmax = xmin + 1,
  fill = rep(c("tomato", "lightblue"), length.out = length(xmin))
)

ggplot(pfu, aes(Passage, PFU_ML, group = interaction(Phage_treatment, Rep), 
                color = Phage_treatment)) +
  geom_rect(
    data = rects,
    aes(xmin = xmin, xmax = xmax, ymin = 1e7, ymax = 1e11, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.1
  ) +
  geom_line(linewidth = 0.2, linetype = "dashed") +
  geom_line(aes(group = Phage_treatment), 
            stat = "summary", fun = mean, linewidth = 1.2, alpha = 0.5) +
  facet_wrap(~Phage_measured) +  
  scale_y_continuous(trans = 'log10', limits = c(1e7, 1e11)) +
  annotation_logticks(sides = "l") +
  scale_color_manual(values = c("Monoculture" = "#FF7F00", "Co-culture" = "grey10")) +
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

ggsave("Figures/Fig2a.png", width = 8, height = 3, dpi = 300)





## Phage fitness assay ##

phage_fitness_pfu <- read_csv("Data/All_phage_PFU.csv")

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
  dplyr::filter(!(Phage == "phage14" & Time == "2"))

phage_fitness_pfu_update$temp <-  factor(phage_fitness_pfu_update$temp, levels = c("37","42"), labels = c("37\u00B0C", "42\u00B0C"))
phage_fitness_pfu_update$line <- factor(phage_fitness_pfu_update$line, levels = c("anc","37_evo","37_C_evo", "42_evo", "42_C_evo", "F_evo", "F_C_evo"), labels = c("Ancestor", "37\u00B0C (mono.)", "37\u00B0C (comp.)", "42\u00B0C (mono.)", "42\u00B0C (comp.)","Fluctuating (mono.)", "Fluctuating (comp.)"))
phage_fitness_pfu_update$Phage <- factor(phage_fitness_pfu_update$Phage, levels = c("phage14", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))

phage_fitness_pfu_update_comp <- phage_fitness_pfu_update %>%
  dplyr::select(-c(dil, pfu)) %>%
  dplyr::group_by(across(-c(Batch, pfu_ml))) %>%  # Group by all columns except `value_column`
  dplyr::summarize(pfu_ml = mean(pfu_ml, na.rm = TRUE)) %>%  # Average duplicates in `value_column`
  ungroup()

phage_fitness_pfu_update_boxplots = phage_fitness_pfu_update_comp %>%
  dplyr::filter((Phage == "ϕ 14-1" & Time == "4") |
                  (Phage == "ϕLUZ19" & Time == "2"))   # Remove T2 values from 14-1 dataset as only carried out for a single rep


phage_fitness_pfu_update_boxplots_normalised <- phage_fitness_pfu_update_boxplots %>%
  dplyr::group_by(Phage, temp) %>%
  dplyr::mutate(
    ancestor_mean = mean(pfu_ml[line == "Ancestor"], na.rm = TRUE),
    PFU_relative = pfu_ml / ancestor_mean
  ) %>%
  ungroup()

phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::mutate(Phage_temp = paste(Phage, temp, sep = "_"))

phage_fitness_pfu_update_boxplots_normalised$Phage_temp = factor(phage_fitness_pfu_update_boxplots_normalised$Phage_temp, levels = c("ϕ 14-1_37°C", "ϕ 14-1_42°C", "ϕLUZ19_37°C", "ϕLUZ19_42°C"), labels = c("ϕ14-1 (37°C)", "ϕ14-1 (42°C)", "ϕLUZ19 (37°C)", "ϕLUZ19 (42°C)"))


ancestor_se <- phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::filter(line == "Ancestor") %>%
  dplyr::group_by(Phage_temp) %>%
  dplyr::summarise(se = sd(PFU_relative, na.rm = TRUE) / sqrt(n()), .groups = "drop") %>%
  dplyr::mutate(
    ymin = 1 - se,
    ymax = 1 + se,
    xmin = -Inf,
    xmax = Inf
  )


phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  dplyr::filter(line != "Ancestor")


phage_fitness_pfu_update_boxplots_normalised$Legend_group <- dplyr::case_when(
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (mono.)" ~ "Static",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (mono.)" ~ "Static",
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (comp.)" ~ "Static\n(Co-culture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (comp.)" ~ "Static\n(Co-culture)",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (mono.)" ~ "Fluc.",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (comp.)" ~ "Fluc.\n(Co-culture)",
)


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

phage_fitness_pfu_update_boxplots_normalised$Legend_group <- factor(phage_fitness_pfu_update_boxplots_normalised$Legend_group, levels = c("Static", "Static\n(Co-culture)", "Fluc.", "Fluc.\n(Co-culture)"))

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
      (temp == "37°C" & evolution_treatment == "42_static" |
         (temp == "42°C" & evolution_treatment == "37_static"
         ))))

ggplot(phage_fitness_pfu_update_boxplots_normalised, aes(Legend_group, PFU_relative, color = as.factor(line))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_wrap2(~ Phage_temp, nrow = 1, strip.position = "top", strip = strip_themed(
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
  scale_color_manual(values = c("37°C (mono.)" = "#377EB8", 
                                "37°C (comp.)" = "black", 
                                "42°C (mono.)" = "#E41A1C", 
                                "42°C (comp.)" = "black",
                                "Fluctuating (mono.)" = "#FF7F00", 
                                "Fluctuating (comp.)" = "black")) +
  
  ylab("Phage growth rel. to ancestor") +
  xlab("Evolution treatment")+
  labs(colour="Evolution treatment")+
  scale_y_continuous(trans='log10', limits = c(0.1,5*10^5)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "grey", alpha = 0.3, inherit.aes = FALSE)+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )


ggsave("Figures/Fig2b.png", width = 8, height = 4, dpi = 300)




# Competition comparison

## Phage competitiveness - Figure 3

comp = read_csv("Data/Competition_assay_14_1_highMOI.csv")

comp2 = comp %>%
  #filter(Evolved_line %in% c("Ancestor", "37", "37C", "42", "42C")) %>%
  select(Biol_rep, Phage, Rep, Evolved_line, Temperature, Competition, PFU_ml)

pfu_avg <- comp2 %>%
  dplyr::group_by(Phage, Rep, Evolved_line, Temperature, Competition) %>%
  dplyr::summarise(mean_pfu = mean(PFU_ml), .groups = "drop")

phage_y_ancestor_means <- pfu_avg %>%
  filter(Phage == "LUZ19", Evolved_line == "Ancestor", Competition == "No_comp") %>%
  dplyr::group_by(Temperature) %>%
  dplyr::summarise(mean_pfu = mean(mean_pfu, na.rm = TRUE), .groups = "drop")

# Step 2: Copy phage_y Comp rows (non-Ancestor), but set competition to "No_comp" and pfu_ml to ancestor mean
phage_y_no_comp_estimates <- pfu_avg %>%
  dplyr::filter(Phage == "LUZ19", Competition == "Comp", Evolved_line != "Ancestor") %>%
  dplyr::select(-mean_pfu, -Competition) %>%
  left_join(phage_y_ancestor_means, by = "Temperature") %>%
  mutate(Competition = "No_comp")

# Step 3: Add these rows to your original dataframe
updated_df <- bind_rows(pfu_avg, phage_y_no_comp_estimates)


#updated_df <- updated_df %>%
 # mutate(initial_pfu = ifelse(Competition == "Comp", 5*10^3, 10^4))

updated_df <- updated_df %>%
  mutate(initial_pfu = ifelse(Competition == "Comp", 5*10^7, 10^8))

df_ratio <- updated_df %>%
  dplyr::mutate(fold_change = mean_pfu / initial_pfu) %>%
  dplyr::select(Phage, Rep, Evolved_line, Temperature, Competition, fold_change) %>%
  pivot_wider(names_from = Competition, values_from = fold_change) %>%
  dplyr::mutate(comp_no_comp_ratio = Comp / No_comp) %>%
  dplyr::select(Phage, Rep, Evolved_line, Temperature, comp_no_comp_ratio)

df_ratio2 = df_ratio %>%
  pivot_wider(names_from = "Phage", values_from = "comp_no_comp_ratio")

df_ratio2_comp = df_ratio2 %>%
  mutate(phage14_log = log10(`14_1` + min(`14_1`[`14_1` > 0]/100, na.rm = TRUE)),
         luz19_log = log10(LUZ19 + min(LUZ19[LUZ19 > 0]/100, na.rm = TRUE)),
         competitiveness = phage14_log - luz19_log) %>%
  mutate(
    regime = case_when(
      Evolved_line %in% c("37", "42")           ~ "Static",
      Evolved_line %in% c("37C", "42C")         ~ "Static (Co-culture)",
      Evolved_line == "F"                       ~ "Fluctuating",
      Evolved_line == "FC"                      ~ "Fluctuating (Co-culture)",
      TRUE                                       ~ NA_character_
    )
  ) %>%
  mutate(regime = factor(regime, levels = c("Static", "Static (Co-culture)", "Fluctuating", "Fluctuating (Co-culture)"), labels = c("Static", "Static\n(Co-culture)", "Fluctuating", "Fluctuating\n(Co-culture)")))


summary_df <- df_ratio2_comp %>%
  filter(!grepl("Ancestor", Evolved_line, ignore.case = TRUE)) %>%
  filter(
    !(
      (Temperature == "42" & Evolved_line %in% c("37", "37C")) |
        (Temperature == "37" & Evolved_line %in% c("42", "42C"))
    )) %>%
  group_by(regime, Temperature) %>%
  summarise(
    mean_comp = mean(competitiveness, na.rm = TRUE),
    n = sum(!is.na(competitiveness)),
    sd = sd(competitiveness, na.rm = TRUE),
    sem = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  na.omit(.)



df_ratio_low_14 = na.omit(summary_df) %>%
  mutate(label = "\U03D5 14-1 vs \U03D5LUZ19 ancestor")


### STATS

df_stats <- df_ratio2_comp %>%
  filter(!grepl("Ancestor", Evolved_line, ignore.case = TRUE)) %>%
  # drop mismatched static-temperature combos
  filter(
    !(
      (Temperature == "42" & Evolved_line %in% c("37", "37C")) |
        (Temperature == "37" & Evolved_line %in% c("42", "42C"))
    ))


df_stats$regime = as.factor(df_stats$regime)


mod <- lm(competitiveness ~ regime*Temperature, df_stats)
anova(mod)
emm = emmeans(mod, ~ regime | Temperature)
contrast(emm, method = "pairwise", adjust = "tukey")




## LUZ19

comp = read_csv("Data/Competition_assay_LUZ19_highMOI.csv")

comp2 = comp %>%
  #dplyr::filter(Evolved_line %in% c("Ancestor", "37", "37C", "42", "42C")) %>%
  dplyr::select(Biol_rep, Phage, Rep, Evolved_line, Temperature, Competition, PFU_ml)

pfu_avg <- comp2 %>%
  dplyr::group_by(Phage, Rep, Evolved_line, Temperature, Competition) %>%
  dplyr::summarise(mean_pfu = mean(PFU_ml), .groups = "drop")

phage_y_ancestor_means <- pfu_avg %>%
  filter(Phage == "14_1", Evolved_line == "Ancestor", Competition == "No_comp") %>%
  dplyr::group_by(Temperature) %>%
  dplyr::summarise(mean_pfu = mean(mean_pfu, na.rm = TRUE), .groups = "drop")

# Step 2: Copy phage_y Comp rows (non-Ancestor), but set competition to "No_comp" and pfu_ml to ancestor mean
phage_y_no_comp_estimates <- pfu_avg %>%
  dplyr::filter(Phage == "14_1", Competition == "Comp", Evolved_line != "Ancestor") %>%
  dplyr::select(-mean_pfu, -Competition) %>%
  left_join(phage_y_ancestor_means, by = "Temperature") %>%
  dplyr::mutate(Competition = "No_comp")

# Step 3: Add these rows to your original dataframe
updated_df <- bind_rows(pfu_avg, phage_y_no_comp_estimates)


#updated_df <- updated_df %>%
  #mutate(initial_pfu = ifelse(Competition == "Comp", 5*10^3, 10^4))

updated_df <- updated_df %>%
  mutate(initial_pfu = ifelse(Competition == "Comp", 5*10^7, 10^8))

df_ratio <- updated_df %>%
  dplyr::mutate(fold_change = mean_pfu / initial_pfu) %>%
  dplyr::select(Phage, Rep, Evolved_line, Temperature, Competition, fold_change) %>%
  pivot_wider(names_from = Competition, values_from = fold_change) %>%
  dplyr::mutate(comp_no_comp_ratio = Comp / No_comp) %>%
  dplyr::select(Phage, Rep, Evolved_line, Temperature, comp_no_comp_ratio)

df_ratio2 = df_ratio %>%
  pivot_wider(names_from = "Phage", values_from = "comp_no_comp_ratio")

df_ratio2_comp = df_ratio2 %>%
  mutate(phage14_log = log10(`14_1` + min(`14_1`[`14_1` > 0]/100, na.rm = TRUE)),
         luz19_log = log10(LUZ19 + min(LUZ19[LUZ19 > 0]/100, na.rm = TRUE)),
         competitiveness = luz19_log - phage14_log) %>%
  mutate(
    regime = case_when(
      Evolved_line %in% c("37", "42")           ~ "Static",
      Evolved_line %in% c("37C", "42C")         ~ "Static (Co-culture)",
      Evolved_line == "F"                       ~ "Fluctuating",
      Evolved_line == "FC"                      ~ "Fluctuating (Co-culture)",
      TRUE                                       ~ NA_character_
    )
  ) %>%
  mutate(regime = factor(regime, levels = c("Static", "Static (Co-culture)", "Fluctuating", "Fluctuating (Co-culture)"), labels = c("Static", "Static\n(Co-culture)", "Fluctuating", "Fluctuating\n(Co-culture)")))
  

summary_df <- df_ratio2_comp %>%
  filter(!grepl("Ancestor", Evolved_line, ignore.case = TRUE)) %>%
  filter(
    !(
      (Temperature == "42" & Evolved_line %in% c("37", "37C")) |
        (Temperature == "37" & Evolved_line %in% c("42", "42C"))
    )) %>%
  group_by(regime, Temperature) %>%
  summarise(
    mean_comp = mean(competitiveness, na.rm = TRUE),
    n = sum(!is.na(competitiveness)),
    sd = sd(competitiveness, na.rm = TRUE),
    sem = sd / sqrt(n),
    .groups = "drop"
  ) %>%
  na.omit(.)


df_ratio_low_luz19 = na.omit(summary_df) %>%
  mutate(label = "\U03D5LUZ19 vs \U03D5 14-1 ancestor")


combined_df = rbind(df_ratio_low_14, df_ratio_low_luz19)

combined_df = combined_df %>%
  mutate(
    Temperature_label = case_when(
      Temperature == "37"          ~ "37°C",
      Temperature == "42"          ~ "42°C")
  )
      

label_map <- c(
  "37°C" = "37°C",
  "42°C" = "42°C"
)

my_labeller <- as_labeller(label_map)

ggplot(na.omit(combined_df), aes(x = regime, y = mean_comp))+
  geom_col(fill = "grey70", col = "black") +
  geom_errorbar(
    aes(ymin = mean_comp - sem, ymax = mean_comp + sem),
    width = 0.2
  )+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  labs(y = "Focal phage competitiveness", x = "Evolution treatment")+
  ylim(-4,6)+
  facet_grid2(label~Temperature_label, strip = strip_themed(
      background_x = list(
        A = element_rect(fill = "lightblue"),
        B = element_rect(fill = "tomato")
      ),
      
      text_x = list(
        A = element_text(color = "black", face = "bold"),
        B = element_text(color = "black", face = "bold")
      )
    ),
    labeller = labeller(Temperature_label = my_labeller) 
    )+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )

ggsave("Figures/Fig2c.png", width = 8, height = 5, dpi = 300)



### STATS

df_stats <- df_ratio2_comp %>%
  filter(!grepl("Ancestor", Evolved_line, ignore.case = TRUE)) %>%
  # drop mismatched static-temperature combos
  filter(
    !(
      (Temperature == "42" & Evolved_line %in% c("37", "37C")) |
        (Temperature == "37" & Evolved_line %in% c("42", "42C"))
    ))
    


df_stats$regime = as.factor(df_stats$regime)
df_stats$Temperature = as.factor(df_stats$Temperature)

# Category stats
mod <- lm(competitiveness ~ regime*Temperature, df_stats)
anova(mod)
emm = emmeans(mod, ~ regime | Temperature)
contrast(emm, method = "pairwise", adjust = "tukey")

emtrends(mod, ~ Temperature, var = "regime_num")

# Numeric stats

mod <- lm(competitiveness ~ regime_num*Temperature, df_stats)
anova(mod)

emtrends(mod, ~ Temperature, var = "regime_num")

ggplot(df_stats, aes(x= regime, y = competitiveness))+
  geom_point()+
geom_line()+
  facet_grid(~Temperature)

### 




## Mutations and pop gen - monoculture (Figure 3)

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

cbPalette4 <- c("grey", "#377EB8", "#FF7F00", "#E41A1C")

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
cbPalette4 <- c("grey", "#377EB8", "#FF7F00", "#E41A1C")

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

ggsave("Figures/Fig3a.png", width = 8, height = 4, dpi = 300)





#14-1

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
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))


groups <- c("14-1-F", "14-1-37", "14-1-42")

# make a mutation ID and filter to those groups
mut_df <- breseq_14_filtered %>%
  filter(evol_group %in% groups) %>%
  mutate(mut_id = paste(gene, sep = "_")) %>%
  distinct(evol_group, mut_id)

# pivot to wide presence/absence
# pivot to wide presence/absence
presence <- mut_df %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = evol_group,
    values_from = present,
    values_fill = 0
  )

# classify into the categories you want
summary_presence <- presence %>%
  transmute(
    mut_id,
    category = case_when(
      `14-1-F` == 0 & `14-1-37` == 1 & `14-1-42` == 0 ~ "37°C\nonly",
      `14-1-F` == 0 & `14-1-37` == 0 & `14-1-42` == 1 ~ "42°C\nonly",
      `14-1-F` == 0 & `14-1-37` == 1 & `14-1-42` == 1 ~ "37°C\n+ 42°C",
      `14-1-F` == 1 & `14-1-37` == 0 & `14-1-42` == 0 ~ "Fluc.\nonly",
      `14-1-F` == 1 & `14-1-37` == 1 & `14-1-42` == 0 ~ "Fluc.\n+ 37°C",
      `14-1-F` == 1 & `14-1-37` == 0 & `14-1-42` == 1 ~ "Fluc.\n+ 42°C",
      `14-1-F` == 1 & `14-1-37` == 1 & `14-1-42` == 1 ~ "All three",
      TRUE ~ NA_character_  # drop other patterns
    )
  ) %>%
  filter(!is.na(category)) %>%
  count(category)


summary_presence$category <- factor(summary_presence$category,
                                    levels = c("37°C\nonly",
                                               "42°C\nonly",
                                               "37°C\n+ 42°C",
                                               "Fluc.\nonly",
                                               "Fluc.\n+ 37°C",
                                               "Fluc.\n+ 42°C",
                                               "All three"))



summary_presence = summary_presence %>%
  complete(category, fill = list(value = NA_real_))


phage14_muts = ggplot(summary_presence, aes(x = category, y = n, fill = category)) +
  geom_col(fill = "lightgrey", col = "black", show.legend = FALSE) +
  labs(
    x = "Evolved line",
    y = "Number of mutations") +
  ylim(0,15)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
    theme(
      axis.title.x = element_text(size = 10, face = "bold"),
      axis.title.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 8))
    



df_unique_mutations <- breseq_14_filtered %>%
  dplyr::filter(evol_group %in% c("14-1-F", "14-1-F-C")) %>%
  dplyr::group_by(type, position, annotation, description, evol_group) %>%
  dplyr::summarise(
    num_isolates = n_distinct(phage),
    .groups = "drop"
  ) %>%
  filter(num_isolates > 1) %>%
  pivot_wider(names_from = "evol_group", values_from = "num_isolates", values_fill = 0)


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
  )  # Scale size for visualization

phage_data_F <- phage_data %>%
  filter(evol_group == "14-1-F")

phage_data_F$phage <- factor(phage_data_F$phage, levels = c("14-1-F-1","14-1-F-2","14-1-F-3","14-1-F-4","14-1-F-5","14-1-F-6"))

phage_data_F = phage_data_F %>%
  complete(phage, fill = list(value = NA_real_))

phage_data_F_C <- phage_data %>%
  filter(evol_group == "14-1-F-C")

phage_data_F_C$phage <- factor(phage_data_F_C$phage, levels = c("14-1-F-C-1","14-1-F-C-2","14-1-F-C-3","14-1-F-C-4","14-1-F-C-5","14-1-F-C-6"))

phage_data_F_C = phage_data_F_C %>%
  complete(phage, fill = list(value = NA_real_))

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
  theme(axis.title=element_text(size=10,face="bold"))+
  theme(legend.position = "none")+ 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


snp_F_C <-  ggplot() +
  geom_segment(data = phage_data_F_C, 
               aes(x = 0, xend = 66298, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F_C %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F_C %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Insertion points
  geom_point(data = phage_data_F_C %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F_C %>% filter(type == "SNP"),
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
  theme(axis.title=element_text(size=10,face="bold"))+
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
  theme(axis.title=element_text(size=10,face="bold"))


#snp_F / genes  + plot_layout(guides = "collect", heights = c(0.8, 0.1))
snp_F_C / genes  + plot_layout(guides = "collect", heights = c(0.8, 0.1))

ggsave("Figures/Fig4B1.png", width = 4, height = 2.5, dpi = 300)




## LUZ19 

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


groups <- c("LUZ19-F", "LUZ19-37", "LUZ19-42")

# make a mutation ID and filter to those groups
mut_df <- breseq_luz19_filtered %>%
  filter(evol_group %in% groups) %>%
  mutate(mut_id = paste(position, gene, sep = "_")) %>%
  distinct(evol_group, mut_id)

# pivot to wide presence/absence
presence <- mut_df %>%
  mutate(present = 1) %>%
  pivot_wider(
    names_from = evol_group,
    values_from = present,
    values_fill = 0
  )

# classify into the categories you want
summary_presence <- presence %>%
  transmute(
    mut_id,
    category = case_when(
      `LUZ19-F` == 0 & `LUZ19-37` == 1 & `LUZ19-42` == 0 ~ "37°C\nonly",
      `LUZ19-F` == 0 & `LUZ19-37` == 0 & `LUZ19-42` == 1 ~ "42°C\nonly",
      `LUZ19-F` == 0 & `LUZ19-37` == 1 & `LUZ19-42` == 1 ~ "37°C\n+ 42°C",
      `LUZ19-F` == 1 & `LUZ19-37` == 0 & `LUZ19-42` == 0 ~ "Fluc.\nonly",
      `LUZ19-F` == 1 & `LUZ19-37` == 1 & `LUZ19-42` == 0 ~ "Fluc.\n+ 37°C",
      `LUZ19-F` == 1 & `LUZ19-37` == 0 & `LUZ19-42` == 1 ~ "Fluc.\n+ 42°C",
      `LUZ19-F` == 1 & `LUZ19-37` == 1 & `LUZ19-42` == 1 ~ "All three",
      TRUE ~ NA_character_  # drop other patterns
    )
  ) %>%
  filter(!is.na(category)) %>%
  count(category)

summary_presence$category <- factor(summary_presence$category,
                                    levels = c("37°C\nonly",
                                               "42°C\nonly",
                                               "37°C\n+ 42°C",
                                               "Fluc.\nonly",
                                               "Fluc.\n+ 37°C",
                                               "Fluc.\n+ 42°C",
                                               "All three"))



summary_presence = summary_presence %>%
  complete(category, fill = list(value = NA_real_))


luz19_muts = ggplot(summary_presence, aes(x = category, y = n)) +
  geom_col(fill = "lightgrey", show.legend = FALSE, col = "black") +
  labs(
    x = "Evolution line",
    y = NULL) +
  ylim(0,15)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )+
  theme(
    axis.title.x = element_text(size = 10, face = "bold"),
    axis.title.y = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(size = 8))

phage14_muts + luz19_muts

ggsave("Figures/Fig3b.png", width = 8, height = 2, dpi = 300)


df_unique_mutations <- breseq_luz19_filtered %>%
  dplyr::filter(!evol_group %in% c("LUZ19-F", "LUZ19-F-C")) %>%
  dplyr::group_by(type, position, annotation, description, evol_group) %>%
  dplyr::summarise(
    num_isolates = n_distinct(phage),
    .groups = "drop"
  ) %>%
  filter(num_isolates > 1) %>%
  pivot_wider(names_from = "evol_group", values_from = "num_isolates", values_fill = 0)



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
  filter(evol_group == "LUZ19-F")

phage_data_F$phage <- factor(phage_data_F$phage, levels = c("LUZ19-F-1","LUZ19-F-2","LUZ19-F-3","LUZ19-F-4","LUZ19-F-5","LUZ19-F-6"))

phage_data_F = phage_data_F %>%
  complete(phage, fill = list(value = NA_real_))

phage_data_F_C <- phage_data %>%
  filter(evol_group == "LUZ19-F-C")

phage_data_F_C$phage <- factor(phage_data_F_C$phage, levels = c("LUZ19-F-C-1","LUZ19-F-C-2","LUZ19-F-C-3","LUZ19-F-C-4","LUZ19-F-C-5","LUZ19-F-C-6"))

phage_data_F_C = phage_data_F_C %>%
  complete(phage, fill = list(value = NA_real_))

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
  theme(axis.title=element_text(size=10,face="bold"))+
  theme(legend.position = "none")+ 
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )


snp_F_C <-  ggplot() +
  geom_segment(data = phage_data_F_C, 
               aes(x = 0, xend = 43187, y = phage, yend = phage), 
               color = "black", size = 0.5) +
  
  # Deletions > 1bp
  geom_rect(data = phage_data_F_C %>% filter(type == "DEL" & mutation_length > 1), 
            aes(xmin = position, 
                xmax = position + mutation_length, 
                ymin = as.numeric(phage) - frequency / 5, 
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Small deletions (1 bp)
  geom_rect(data = phage_data_F_C %>% filter(type == "DEL" & mutation_length == 1), 
            aes(xmin = position - 250, 
                xmax = position + 250, 
                ymin = as.numeric(phage) - frequency / 5,  
                ymax = as.numeric(phage) + frequency / 5, 
                fill = "Deletion"), color = "black",
            alpha = 0.5) +
  
  # Insertion points
  geom_point(data = phage_data_F_C %>% filter(type == "INS"),
             aes(x = position, y = phage, size = frequency, shape = "Insertion"), 
             fill = "red", color = "black", alpha = 0.6) +
  
  # SNP points
  geom_point(data = phage_data_F_C %>% filter(type == "SNP"),
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
  theme(axis.title=element_text(size=10,face="bold"))+
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
  theme(axis.title=element_text(size=10,face="bold"))


snp_F / genes  + plot_layout(guides = "collect", heights = c(0.8, 0.1))
snp_F_C / genes  + plot_layout(guides = "collect", heights = c(0.8, 0.1))

ggsave("Figures/Fig4B2.png", width = 4, height = 2.5, dpi = 300)




### Population genetics - Figure 3

#14-1
    
breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.2)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

# No Ancestor

ancestor_mutations <- breseq_14_filtered %>%
  filter(grepl("anc", phage, ignore.case = TRUE)) %>%
  distinct(type, position, mutation)

# Step 2: remove rows that match any of those mutations
breseq_14_filtered <- breseq_14_filtered %>%
  anti_join(ancestor_mutations, by = c("type", "position", "mutation"))

df_wide <- breseq_14_filtered %>%
  unite(mutation_id, position, mutation, remove = FALSE) %>%  # Ignore seq_id
  select(phage, mutation_id, frequency) %>%
  spread(key = mutation_id, value = frequency, fill = 0)

# Compute Euclidean distance matrix
distance_matrix <- dist(df_wide[, -1], method = "euclidean")

# Convert to matrix format
distance_matrix <- as.matrix(distance_matrix)
rownames(distance_matrix) <- df_wide$phage
colnames(distance_matrix) <- df_wide$phage 

remove_labels <- grepl("anc", rownames(distance_matrix)) | grepl("anc", colnames(distance_matrix)) | grepl("C", rownames(distance_matrix)) | grepl("C", colnames(distance_matrix))
phage14_distance_matrix <- distance_matrix[!remove_labels, !remove_labels]



# LUZ19

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.2)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

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

remove_labels <- grepl("anc", rownames(distance_matrix)) | grepl("anc", colnames(distance_matrix)) | grepl("C", rownames(distance_matrix)) | grepl("C", colnames(distance_matrix))
LUZ19_distance_matrix_nc <- distance_matrix[!remove_labels, !remove_labels]



#14-1

PcoA_distance <- pcoa(phage14_distance_matrix_nc)

eigenvalues <- PcoA_distance$values$Eigenvalues
variance_explained <- (eigenvalues / sum(eigenvalues)) * 100
print(variance_explained)


PcoA_data <- data.frame(PcoA_distance$vectors)
Best_axes <- PcoA_data[,c("Axis.1","Axis.2")]

Best_axes <- Best_axes %>%
  mutate(evol_group = sub("-[^-]+$", "", rownames(.)),
         temp_group = sub("-C", "", evol_group))

cbPalette3 <- c("#1f78b4", "#FF7F00", "#e31a1c")
Best_axes$evol_group <- factor(Best_axes$evol_group, levels = c("14-1-37", "14-1-F", "14-1-42"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))

phage14_pcoa = ggplot(data=Best_axes, aes(x=Axis.1,y=Axis.2, colour = evol_group))+ 
  geom_point(size=2)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  stat_ellipse(size=1)+
  labs(x="PCo1 (49.6%)",y="PCo2 (10.9%)")+
  scale_colour_manual(values=cbPalette3)+
  theme(axis.title=element_text(size=10,face="bold"))+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=10))+
  labs(color="Evolved line")+
  theme(strip.text = element_text(face="bold", size=10))

#evol_group <- str_remove(rownames(phage14_distance_matrix_nc_noanc), "-[^-]+$")
#evol_group
#ano = anosim(phage14_distance_matrix_nc_noanc, evol_group, distance = "bray", permutations = 9999)
#ano # p = 0.0023, R = 1





PcoA_distance <- pcoa(LUZ19_distance_matrix_nc)

eigenvalues <- PcoA_distance$values$Eigenvalues
variance_explained <- (eigenvalues / sum(eigenvalues)) * 100
print(variance_explained)

PcoA_data <- data.frame(PcoA_distance$vectors)
Best_axes <- PcoA_data[,c("Axis.1","Axis.2")]
Best_axes <- Best_axes %>%
  mutate(evol_group = sub("-[^-]+$", "", rownames(.)),
         temp_group = sub("-C", "", evol_group))

Best_axes$evol_group <- factor(Best_axes$evol_group, levels = c("LUZ19-37", "LUZ19-F", "LUZ19-42"), labels = c("37\u00B0C (static)", "Fluctuating", "42\u00B0C (static)"))

LUZ19_pcoa = ggplot(data=Best_axes, aes(x=Axis.1,y=Axis.2, colour = evol_group))+ 
  geom_point(size=2)+
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  stat_ellipse(size=1)+
  labs(x="PCo1 (42.3%)",y="PCo2 (11.3%)")+
  scale_colour_manual(values=cbPalette3)+
  theme(axis.title=element_text(size=10,face="bold"))+
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=10))+
  labs(color="Evolved line")

#evol_group <- str_remove(rownames(luz19_distance_matrix_nc_noanc), "-[^-]+$")
#evol_group
#ano = anosim(luz19_distance_matrix_nc_noanc, evol_group, distance = "bray", permutations = 9999)
#ano # p = 0.003, R = 0.9093


phage14_pcoa + LUZ19_pcoa + plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave("Figures/Fig3c.png", width = 6, height = 2.5, dpi = 300)





# Fig 5b - Distance to ancestor

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.1)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

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

distance_long$evol_type = factor(distance_long$evol_type, levels = c("Static\n37\u00B0C", "Static\n42\u00B0C", "Fluctuating", "Static\n37\u00B0C\n(Co-culture)", "Static\n42\u00B0C\n(Co-culture)", "Fluctuating\n(Co-culture)"))

cbPalette3 <- c("#1f78b4", "#FF7F00", "#e31a1c")

phage14_dist_from_anc = ggplot(distance_long, aes(evol_type, distance)) +
  geom_boxplot(color = "black") +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("Genetic distance from ancestor") +
  xlab("Evolved line")+
  ylim(0,2.7)+
  labs(colour="Evolved line")+
  #scale_colour_manual(values=cbPalette3) +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"),
        axis.text.x = element_text(size = 7))

mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)
summary(mod)

emmeans(mod, pairwise ~ evol_type)



# LUZ19 

breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.1)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

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

distance_long$evol_type = factor(distance_long$evol_type, levels = c("Static\n37\u00B0C", "Static\n42\u00B0C", "Fluctuating", "Static\n37\u00B0C\n(Co-culture)", "Static\n42\u00B0C\n(Co-culture)", "Fluctuating\n(Co-culture)"))


luz19_dist_from_anc = ggplot(distance_long, aes(evol_type, distance)) +
  geom_boxplot(colour = "black") +
  geom_jitter() + 
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  ylab("Genetic distance from ancestor") +
  xlab("Evolved line")+
  ylim(0,2.7)+
  labs(colour="Evolved line")+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10),legend.text.align = 0)+
  theme(axis.title=element_text(size=10,face="bold"),
        axis.text.x = element_text(size = 7))

mod = lm(distance ~ evol_type, distance_long)
anova(mod)
plot(mod)

emmeans(mod, pairwise ~ evol_type)

phage14_dist_from_anc + luz19_dist_from_anc

ggsave("Figures/Fig4a.png", width = 8, height = 3, dpi = 300)



# Within-group diversity

# 14-1

breseq_14 <- read_csv("Data/14_1_combined_downsampled_breseq.csv")

breseq_14_filtered <- breseq_14 %>%
  filter(frequency > 0.1)

breseq_14_filtered$position <- as.factor(breseq_14_filtered$position)
breseq_14_filtered$mutation <- as.factor(breseq_14_filtered$mutation)
breseq_14_filtered$phage <- as.factor(breseq_14_filtered$phage)

# Anc comparison

df_wide <- breseq_14_filtered %>%
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
  filter(phage1 != phage2) %>%
  filter(evol_group_first == evol_group_second) %>%
  filter(evol_group_first %in% c("14-1-37", "14-1-42", "14-1-F"))


ggplot(distance_long, aes(evol_group_first, distance, colour = as.factor(evol_group_first))) +
  geom_boxplot() +
  geom_jitter() + 
  theme_bw()+
  ylab("Fst") +
  xlab("Evolved line")+
  labs(colour="Evolved line")+
  #scale_colour_manual(values=cbPalette7) +
  #ylim(0,0.1) +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10),legend.text.align = 0)+
  theme(axis.title=element_text(size=11,face="bold"), axis.text.x = element_text(angle = 45, hjust = 1))





breseq_luz19 <- read_csv("Data/LUZ19_combined_downsampled_breseq.csv")

breseq_luz19_filtered <- breseq_luz19 %>%
  filter(frequency > 0.1)

breseq_luz19_filtered$position <- as.factor(breseq_luz19_filtered$position)
breseq_luz19_filtered$mutation <- as.factor(breseq_luz19_filtered$mutation)
breseq_luz19_filtered$phage <- as.factor(breseq_luz19_filtered$phage)

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
  filter(phage1 != phage2) %>%
  filter(evol_group_first == evol_group_second) %>%
  filter(evol_group_first %in% c("LUZ19-37", "LUZ19-42", "LUZ19-F"))


ggplot(distance_long, aes(evol_group_first, distance, colour = as.factor(evol_group_first))) +
  geom_boxplot() +
  geom_jitter() + 
  theme_bw()+
  ylab("Fst") +
  xlab("Evolved line")+
  labs(colour="Evolved line")+
  #scale_colour_manual(values=cbPalette7) +
  #ylim(0,0.1) +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=10),legend.text.align = 0)+
  theme(axis.title=element_text(size=11,face="bold"), axis.text.x = element_text(angle = 45, hjust = 1))


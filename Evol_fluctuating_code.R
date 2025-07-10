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


# Population dynamics Figure 1

pfu <- read_csv("Data/PFU_values_each_passage.csv")
pfu$Temperature_treatment <-  factor(pfu$Temperature_treatment, levels = c("37C_static","42C_static","Fluctuating"), labels = c("37\u00B0C", "42\u00B0C", "Fluctuating"))
pfu$Phage_measured <- factor(pfu$Phage_measured, levels = c("phage14_1", "LUZ19"), labels = c("\U03D5 14-1", "\U03D5LUZ19"))
pfu$Phage_treatment <- factor(pfu$Phage_treatment, levels = c("Single_phage", "Competition"), labels = c("Monoculture", "Competition"))

pfu <- pfu %>%
  filter(Temperature_treatment == "Fluctuating") %>%
  filter(Passage %% 1 == 0)


ggplot(pfu, aes(Passage, PFU_ML, group = interaction(Phage_treatment, Rep), 
                color = Phage_treatment)) +
  geom_line(linewidth = 0.2, linetype = "dashed") +
  geom_line(aes(group = Phage_treatment), 
            stat = "summary", fun = mean, linewidth = 1.2, alpha = 0.5) +
  facet_wrap(~Phage_measured) +  
  scale_y_continuous(trans = 'log10', limits = c(1e7, 1e11)) +
  annotation_logticks(sides = "l") +
  scale_color_manual(values = c("Monoculture" = "#377EB8", "Competition" = "#E41A1C")) +
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

ggsave("Figures/Fig1a.png", width = 8, height = 4, dpi = 300)


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

ggplot(phage_fitness_pfu_update_boxplots_single_normalised, aes(line, PFU_relative, colour = as.factor(line))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(~Phage_temp, scales = "free_y") +
  ylab("Phage growth rel. to ancestor") +
  xlab("Evolved line")+
  labs(colour="Evolved line")+
  scale_y_continuous(trans='log10', limits = c(0.01,10^6)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "#4DAF4A")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "#4DAF4A", alpha = 0.3, inherit.aes = FALSE)+
  scale_colour_manual(values=cbPalette2) +
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )

ggsave("Figures/Fig1b.png", width = 10, height = 3, dpi = 300)



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
  phage_fitness_pfu_update_boxplots_normalised$line == "37°C (mono.)" ~ "37°C (mono.)",
  phage_fitness_pfu_update_boxplots_normalised$line == "42°C (mono.)" ~ "42°C (mono.)",
  phage_fitness_pfu_update_boxplots_normalised$line == "Fluctuating (mono.)" ~ "Monoculture",
  grepl("Competition", phage_fitness_pfu_update_boxplots_normalised$competition) ~ "Competition"
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

phage_fitness_pfu_update_boxplots_normalised = phage_fitness_pfu_update_boxplots_normalised %>%
  filter(evolution_treatment == "Fluctuating")

phage_fitness_pfu_update_boxplots_normalised$line <- factor(phage_fitness_pfu_update_boxplots_normalised$line, levels = c("anc", "Fluctuating (mono.)", "Fluctuating (comp.)"), labels = c("Ancestor","Monoculture", "Competition"))

ggplot(phage_fitness_pfu_update_boxplots_normalised, aes(line, PFU_relative, colour = as.factor(Legend_group))) +
  geom_boxplot() +
  geom_jitter(width = 0.2, height = 0) +
  theme_bw()+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(~Phage_temp, scales = "free_y") +
  scale_color_manual(values = c(
    "Monoculture" = "#FF7F00",
    "Competition" = "grey10"
  )) +
  ylab("Phage growth rel. to ancestor") +
  xlab("Evolution treatment")+
  labs(colour="Evolution treatment")+
  scale_y_continuous(trans='log10', limits = c(0.1,5*10^4)) +
  annotation_logticks(sides="l")+
  geom_hline(yintercept = 1, linetype = "dashed", color = "#4DAF4A")+
  geom_rect(data = ancestor_se,
            aes(ymin = ymin, ymax = ymax, xmin = xmin, xmax = xmax),
            fill = "#4DAF4A", alpha = 0.3, inherit.aes = FALSE)+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(
    legend.title = element_text(colour = "black", size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.position = "bottom",
    axis.title = element_text(size = 10, face = "bold"),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey95", colour = "black", size = 0.3)
  )


ggsave("Figures/Fig1c.png", width = 10, height = 3, dpi = 300)

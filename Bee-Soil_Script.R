### READ THIS ###

# To work, this script must be placed in the same folder with:
# Covercrops_Dataset.xlsx


# 0 - PREPARATION ----

.start_time <- Sys.time() #Set the clock!


packages <- c(
  #Operate with data
  "openxlsx", "readxl", "writexl", "dplyr", "tidyr", "tibble",
  "stats", "devtools", "broom.mixed", "BiocManager",
  
  #Models
  "vegan", "lme4", "lmerTest", "car", "emmeans", "Maaslin2",
  
  #Plots
  "ggplot2", "viridis", "VennDiagram",
  
  #Wake up!
  "beepr", "here" 
)


# Install and load missing packages
install.packages(setdiff(packages, rownames(installed.packages())))
invisible(lapply(packages, library, character.only = T))

#BiocManager::install("DESeq2")
library(DESeq2)

# Create folders to save results and figures
directory <- c("figures", "results")
for(name in directory) {
  if (!dir.exists(name)){
    dir.create(name)
  }else{
    print("dir already exists")
  }
}


# Clean the environment
rm(list = ls())


# Is the directory correct?
here()


### YOU SHOULD NOW BE ABLE TO RUN THE SCRIPT ###
beep()
cat("Required packages installed, now the script should run.\n\n\n")


## FOR BIOCMANAGER ##

# Installa BiocManager se manca
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")





#LOAD DATASETS

#Load metadata
metadata <- read.xlsx("Dataset_Lmal2023.xlsx", sheet = "Metadata")
metadata <- metadata[,-4]
metadata$Site_code <- as.factor(metadata$Site_code)
metadata$Type <- as.factor(metadata$Type)
metadata$Type <- relevel(metadata$Type, ref = "Control") #Control soil is our reference level
metadata$Type <- factor(metadata$Type, levels = c("Bee", "Nest", "Control"))
metadata$Pair <- as.factor(metadata$Pair)  #Paired study design

#Chemistry dataset
chem <- read.xlsx("Dataset_Lmal2023.xlsx", sheet = "Chemistry")
chem$Site_code <- as.factor(chem$Site_code)
chem$Type <- as.factor(chem$Type)
chem$Type <- relevel(chem$Type, ref = "Control")
chem$Type <- factor(chem$Type, levels = c("Nest", "Control"))
chem <- merge(chem, metadata[, c("ID", "Pair")], by = "ID", all.x = T)
chem$Pair <- as.factor(chem$Pair)

str(chem)
head(chem)


#Load OTU table 16s
bact <- read.xlsx("Dataset_Lmal2023.xlsx", sheet = "16s")
bact$Phylum <- as.factor(bact$Phylum)
bact$Class <- as.factor(bact$Class)
bact$Order <- as.factor(bact$Order)
bact$Family <- as.factor(bact$Family)
bact$Genus <- as.factor(bact$Genus)

str(bact)
head(bact)


#Load OTU table ITS
fung <- read.xlsx("Dataset_Lmal2023.xlsx", sheet = "ITS")
fung$Phylum <- as.factor(fung$Phylum)
fung$Class <- as.factor(fung$Class)
fung$Order <- as.factor(fung$Order)
fung$Family <- as.factor(fung$Family)
fung$Genus <- as.factor(fung$Genus)

str(fung)
head(fung)


#Custom aspects
my_palette <- c(
  "Bee"     = "purple4",
  "Control" = "skyblue3",
  "Nest"    = "gold3"
)

my_points <- c(
  "A" = 21,
  "C" = 22,
  "R" = 23,
  "S" = 24
)



#Finish
beep()
cat("Tables prepared.\n\n\n")





#---------------------------------------------------------------------------#





# 1 - CHEMISTRY ----

chem_vars <- colnames(chem[,5:23])

## 1.1 - Descriptive statistics ----
means_clean <- chem %>%
  pivot_longer(cols = all_of(chem_vars), names_to = "variable", values_to = "value") %>%
  group_by(variable, Type) %>%
  summarise(
    mean = mean(value, na.rm = T),
    sd   = sd(value, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.3f ± %.3f", mean, sd)) %>%
  select(variable, Type, mean_sd) %>%
  pivot_wider(names_from = Type, values_from = c(mean_sd))



## 1.2 - Mixed models ----
chem_results <- lapply(chem_vars, function(v) {
  f <- as.formula(paste0(v, " ~ Type + (1|Site_code/Pair)"))
  mod <- lmer(f, data = chem)
  
  # Type-II ANOVA
  anova_res <- Anova(mod, type = "II") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term == "Type") %>%
    transmute(
      variable  = v,
      Chisq     = round(Chisq, 3),
      df        = Df,
      p.value   = round(`Pr(>Chisq)`, 3)
    )
  
  # Post-hoc pairwise comparisons
  posthoc <- emmeans(mod, pairwise ~ Type, adjust = "tukey")$contrasts %>%
    as.data.frame() %>%
    transmute(
      variable  = v,
      contrast  = contrast,
      estimate  = round(estimate, 3),
      std.error = round(SE, 3),
      df        = round(df, 3),
      statistic = round(t.ratio, 3),
      p.value   = round(p.value, 3)
    )
  
  bind_rows(anova_res, posthoc)
}) %>% 
  bind_rows() %>%
  mutate(p_fdr = round(p.adjust(p.value, method = "fdr"), 3))



## 1.3 - Boxplots ----

plots <- lapply(chem_vars, function(varname){
  
  #Prepare the data
  data <- chem %>%
    dplyr::mutate(
      Type = Type,
      Site = Site_code
    )
  
  #Calculate medians
  site_medians <- data %>%
    dplyr::group_by(Site, Type) %>%
    dplyr::summarise(
      value = median(.data[[varname]], na.rm = T),
      .groups = "drop"
    )
  
  lab <- site_medians %>%
    dplyr::filter(Type == "Control")
  
  #Make the plot
  p <- ggplot(data, aes(x = Type, y = .data[[varname]])) +
    geom_boxplot(aes(color = Type), outlier.shape = NA) +
    geom_point(aes(color = Type, shape = Site),
               position = position_jitter(width = 0.1),
               alpha = 0.7, size = 2) +
    geom_point(data = site_medians,
               aes(x = Type, y = value, shape = Site),
               size = 3, color = "black") +
    geom_line(data = site_medians,
              aes(x = Type, y = value, group = Site),
              inherit.aes = FALSE, alpha = 0.6, color = "grey35") +
    scale_color_manual(values = my_palette) +
    scale_shape_manual(values = my_points) +
    theme_bw() +
    theme(text = element_text(size = 20)
          , axis.text = element_text(size = 15, color = "black")
          , legend.position = "none") +
    labs(x = NULL, y = varname)
  
  #Save the plot
  ggsave(filename = paste0("figures/CHEM-", varname, ".pdf"),
         width = 10, height = 20, units = "cm", plot = p)
  
  return(p)
  
})



## 1.4 - PCA + PERMANOVA ----

# Create a matrix
chem_matrix <- chem %>%
  select(ID, all_of(chem_vars[1:7])) %>%
  column_to_rownames("ID") %>%
  scale(center = T, scale = T)


# PCA
pca <- prcomp(chem_matrix, center = F, scale. = F)
pca_var <- round(summary(pca)$importance[2, 1:2] * 100, 1)

# Extract scores
pca_scores <- as.data.frame(pca$x[, 1:2]) %>%
  rownames_to_column("ID") %>%
  left_join(metadata[, c("ID", "Type", "Site_code", "Pair")], by = "ID")

# Convex hulls
hull_df <- pca_scores %>%
  group_by(Site_code, Type) %>%
  dplyr::slice(chull(PC1, PC2))

# Plot
ggplot(pca_scores, aes(PC1, PC2)) +
  geom_polygon(data = hull_df,
               aes(color = Type, group = interaction(Site_code, Type)),
               fill = NA, linewidth = 0.7) +
  geom_point(aes(fill = Type, shape = Site_code), size = 2, colour = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  scale_shape_manual(values = my_points) +
  labs(x = paste0("PC1 (", pca_var[1], "%)"),
       y = paste0("PC2 (", pca_var[2], "%)"))

ggsave("figures/CHEM-PCA.pdf", height = 15, width = 15, units = "cm")


# PERMANOVA
permanova <- adonis2(chem_matrix ~ Type * Site_code, data = chem,
                     method = "euclidean", permutations = 9999, by = "terms")
disp <- betadisper(vegdist(chem_matrix, method = "euclidean"), chem$Type)


# Results
sink("results/CHEM-PERMANOVA.txt")
anova(disp) #Ok
permutest(disp) #OK
tapply(disp$distances, chem$Type, mean)
anova(disp)
TukeyHSD(disp)
permanova
sink()



## 1.5 - Coefficient of Variation ----

# Calculate CV
cv <- function(x) sd(x, na.rm = T) / mean(x, na.rm = t) * 100

cv_results <- chem %>%
  filter(Type %in% c("Nest", "Control")) %>%
  pivot_longer(cols = all_of(chem_vars), names_to = "variable", values_to = "value") %>%
  group_by(Site_code, Type, variable) %>%
  summarise(CV = round(cv(value), 2), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = CV, names_prefix = "CV_")

cv_long <- cv_results %>%
  pivot_longer(cols = c(CV_Nest, CV_Control), 
               names_to = "Type", 
               values_to = "CV") %>%
  mutate(Type = gsub("CV_", "", Type))


# t-Test on CV
cv_ttest <- lapply(chem_vars, function(v) {
  nest    <- cv_long %>% filter(variable == v, Type == "Nest")    %>% pull(CV)
  control <- cv_long %>% filter(variable == v, Type == "Control") %>% pull(CV)
  
  t.test(nest, control, paired = FALSE) %>%
    broom::tidy() %>%
    mutate(variable = v) %>%
    select(variable, estimate, statistic, p.value, conf.low, conf.high)
  }) %>%
  bind_rows() %>%
  mutate(p_fdr = round(p.adjust(p.value, method = "fdr"), 3),
         across(where(is.numeric), ~ round(.x, 3)))


# Plot

cv_mean <- cv_long %>%
  group_by(variable, Type) %>%
  summarise(
    CV_mean = mean(CV, na.rm = TRUE),
    CV_se   = sd(CV, na.rm = TRUE) / sqrt(n()),
    .groups = "drop"
  )


ggplot(cv_mean, aes(x = variable, y = CV_mean, fill = Type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8),
           width = 0.7, color = "black", linewidth = 0.3, alpha = 0.7) +
  geom_errorbar(aes(ymin = CV_mean - CV_se, ymax = CV_mean + CV_se),
                position = position_dodge(width = 0.8),
                width = 0.2, linewidth = 0.2) +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        text = element_text(size = 14),
        legend.position = "none") +
  labs(x = NULL, y = "Coefficient of Variation (%)", fill = "Type", color = "Type")

ggsave("figures/CHEM-CV.pdf", width = 20, height = 15, units = "cm")



## 1.6 - Save results ----
write.xlsx(chem_results, "results/CHEM-LMM.xlsx")
write.xlsx(means_clean, "results/CHEM-SumStat.xlsx")
write.xlsx(cv_results, "results/CHEM-CV.xlsx")
write.xlsx(cv_ttest, "results/CHEM-CV-Test.xlsx")

#Finish
rm(list = setdiff(ls(), c("bact", "chem", "fung", "metadata", "my_palette", "my_points")))
dev.off()
beep()
cat("Chemistry analysed. \n\n\n")





#---------------------------------------------------------------------------#





# 2 - BACTERIA ----

## 2.1 - Stacked barcharts ----

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

plots_16S <- lapply(tax_levels, function(tax_level) {
  
  tax_sym <- sym(tax_level)
  
  # Prepare the data
  bact_long <- bact %>%
    pivot_longer(cols = matches("^[A-Z]{2}[0-9]$"),
                 names_to = "ID", values_to = "Abundance") %>%
    left_join(metadata[, c("ID", "Type")], by = "ID") %>%
    mutate(!!tax_sym := as.character(!!tax_sym))

  # Top 10 taxa plotted
  top_taxa <- bact_long %>%
    dplyr::group_by(!!tax_sym) %>%
    summarise(Total = sum(Abundance)) %>%
    slice_max(Total, n = 9) %>%
    pull(!!tax_sym)
  
  df <- bact_long %>%
    mutate(!!tax_sym := ifelse(!!tax_sym %in% top_taxa, !!tax_sym, "Other")) %>%
    group_by(ID, Type, !!tax_sym) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(ID) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>% 
    ungroup()
  
  df %>% dplyr::count(ID) %>% print(n = Inf)
  
  # Colour palette
  taxa_names <- c(sort(unique(df[[tax_level]][df[[tax_level]] != "Other"])), "Other")
  my_colors <- setNames(
    c(viridis::viridis(length(taxa_names) - 1), "grey60"),
    taxa_names
  )
  
  # Plot
  p <- ggplot(df, aes(x = ID, y = RelAbund, fill = !!tax_sym)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, alpha = 0.9) +
    facet_wrap(~ Type, scales = "free_x") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    labs(y = "Relative Abundance", x = NULL,
         title = paste("Bacterial community composition by", tax_level, "(only top-10 taxa showed)")) +
    scale_fill_manual(values = my_colors, name = tax_level)
  
  ggsave(filename = paste0("figures/16S-", tax_level, ".pdf"),
         width = 25, height = 20, units = "cm", plot = p)
  
  return(p)
  
})



# Heatmap

heat_bact <- bact %>%
  pivot_longer(cols = matches("^[A-Z]{2}[0-9]$"), names_to = "ID", values_to = "Abundance") %>%
  left_join(metadata[, c("ID", "Type")], by = "ID") %>%
  mutate(Genus = as.character(Genus)) %>%
  group_by(Type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  ungroup()

top_genera_bact <- heat_bact %>%
  filter(Type == "Bee") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  slice_max(Total, n = 10) %>%
  pull(Genus)

heat_bact %>%
  filter(Genus %in% top_genera_bact) %>%
  mutate(log10_abundance = log10(Abundance + 1)) %>%
  ggplot(aes(x = Type, y = Genus, fill = log10_abundance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "white", mid = "pink2", high = "red4",
                       midpoint = 3,
                       name = "log10\n(Reads + 1)") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic", size = 10),
        axis.text.x = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL, title = "Top-10 Bacterial genera in bees")

ggsave("figures/16S-Heatmap.pdf", width = 15, height = 15, units = "cm")



## 2.2 - A-Dversity models ----

# OTU matrix
bact_matrix <- bact %>%
  dplyr::select(matches("^[A-Z]{2}[0-9]$")) %>%
  t()

# Calculate indices
alpha_div <- data.frame(
  ID       = rownames(bact_matrix),
  Richness = specnumber(bact_matrix),
  Shannon  = diversity(bact_matrix, index = "shannon"),
  Pielou_Evenness = diversity(bact_matrix, index = "shannon") / log(specnumber(bact_matrix))
)


# Build dataset for models
df_model <- alpha_div %>%
  left_join(metadata[, c("ID", "Type", "Pair", "Site_code")], by = "ID")

bio_vars <- colnames(df_model[,2:4])


# Calculate and export means
means_clean <- df_model %>%
  pivot_longer(cols = all_of(bio_vars), names_to = "variable", values_to = "value") %>%
  group_by(variable, Type) %>%
  summarise(
    mean = mean(value, na.rm = T),
    sd   = sd(value, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.3f ± %.3f", mean, sd)) %>%
  select(variable, Type, mean_sd) %>%
  pivot_wider(names_from = Type, values_from = c(mean_sd))


# Models
bio_results <- lapply(bio_vars, function(v) {
  f <- as.formula(paste0(v, " ~ Type + (1|Site_code/Pair)"))
  mod <- lmer(f, data = df_model)
  
  # Type-II ANOVA
  anova_res <- Anova(mod, type = "II") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term == "Type") %>%
    transmute(
      variable  = v,
      Chisq     = round(Chisq, 3),
      df        = Df,
      p.value   = round(`Pr(>Chisq)`, 3)
    )
  
  # Post-hoc pairwise comparisons
  posthoc <- emmeans(mod, pairwise ~ Type, adjust = "tukey")$contrasts %>%
    as.data.frame() %>%
    transmute(
      variable  = v,
      contrast  = contrast,
      estimate  = round(estimate, 3),
      std.error = round(SE, 3),
      df        = round(df, 3),
      statistic = round(t.ratio, 3),
      p.value   = round(p.value, 3)
    )
  
  bind_rows(anova_res, posthoc)
}) %>% 
  bind_rows() %>%
  mutate(p_fdr = round(p.adjust(p.value, method = "fdr"), 3))



## 2.3 - A-Dversity boxplots ----

plots <- lapply(bio_vars, function(varname){
  
  #Prepare the data
  data <- df_model %>%
    dplyr::mutate(
      Type = Type,
      Site = Site_code
    )
  
  #Calculate medians
  site_medians <- data %>%
    dplyr::group_by(Site, Type) %>%
    dplyr::summarise(
      value = median(.data[[varname]], na.rm = T),
      .groups = "drop"
    )
  
  lab <- site_medians %>%
    dplyr::filter(Type == "Control")
  
  #Make the plot
  p <- ggplot(data, aes(x = Type, y = .data[[varname]])) +
    geom_boxplot(aes(color = Type), outlier.shape = NA) +
    geom_point(aes(color = Type, shape = Site),
               position = position_jitter(width = 0.1),
               alpha = 0.7, size = 1.5) +
    geom_point(data = site_medians,
               aes(x = Type, y = value, shape = Site),
               size = 3, color = "black") +
    geom_line(data = site_medians,
              aes(x = Type, y = value, group = Site),
              inherit.aes = FALSE, alpha = 0.6, color = "grey35") +
    scale_color_manual(values = my_palette) +
    scale_shape_manual(values = my_points) +
    theme_bw() +
    theme(text = element_text(size = 20)
          , axis.text = element_text(size = 15, color = "black")
          , legend.position = "none") +
    labs(x = NULL, y = varname)
  
  #Save the plot
  ggsave(filename = paste0("figures/16S-", varname, ".pdf"),
         width = 10, height = 20, units = "cm", plot = p)
  
  return(p)
  
})


# Save results
write.xlsx(bio_results, "results/16S-LMM.xlsx")
write.xlsx(means_clean, "results/16S-SumStat.xlsx")



## 2.4 - MAASLIN ----

# Get the OTUs
get_otus <- function(type) {
  samples <- metadata %>% dplyr::filter(Type == type) %>% dplyr::pull(ID)
  bact %>%
    select(OTU, all_of(samples)) %>%
    filter(rowSums(dplyr::select(., -OTU)) > 0) %>%
    pull(OTU)
}

otus_control <- get_otus("Control")
otus_bee     <- get_otus("Bee")
otus_nest    <- get_otus("Nest")


# Venn Diagram
venn.diagram(
  x = list(Control = otus_control, Bee = otus_bee, Nest = otus_nest),
  filename = "figures/16S-Venn.tiff",
  output = T,
  fill = c(my_palette["Control"], my_palette["Bee"], my_palette["Nest"]),
  alpha = 0.4,
  cex = 1.5,
  fontface = "bold",
  cat.cex = 0,
  rotation.degree = 240)


# OTU matrix
otu_mat <- bact %>%
  column_to_rownames("OTU") %>%
  select(matches("^[A-Z]{2}[0-9]$")) %>%
  as.matrix()


# Prepare the data
otu_table <- bact[, !(names(bact) %in% c("Phylum","Class","Order","Family","Genus","Species"))]
otu_table$Genus <- bact$Genus

genus_table <- otu_table %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
  as.data.frame() %>%
  column_to_rownames("Genus")

genus_table_t <- as.data.frame(t(genus_table))


# Metadata
metadata2 <- metadata[metadata$Type %in% c("Nest","Control"), ]
rownames(metadata2) <- metadata2$ID
genus_table_t <- genus_table_t[rownames(genus_table_t) %in% rownames(metadata2), ]


# Relative abundance
genus_rel <- genus_table_t / rowSums(genus_table_t)
genus_rel[is.na(genus_rel)] <- 0
prev <- colSums(genus_rel > 0) / nrow(genus_rel)
abund <- colMeans(genus_rel)


# Prevalence (>=50%) + Abundance (>=0.1%)
keep <- prev >= 0.5 & abund >= 0.001
genus_table_t <- genus_table_t[, keep]
metadata2 <- metadata2[match(rownames(genus_table_t), metadata2$ID), ]
metadata2$Type <- relevel(as.factor(metadata2$Type), ref = "Control")


# Run the analysis
maaslin_genus <- Maaslin2(
  input_data = genus_table_t,
  input_metadata = metadata2,
  output = "Genus-Bact",
  fixed_effects = "Type",
  random_effects = "Pair",
  normalization = "CLR",
  transform = "NONE",
  correction = "BH"
)


# Save the results
res <- maaslin_genus$results
write.xlsx(res, "results/Maaslin-Bact.xlsx")


# Plot the results
bee_samples <- metadata$ID[metadata$Type == "Bee"]
bee_genus <- rownames(genus_table)[
  rowSums(genus_table[, bee_samples] > 0) > 0]

plot_df <- res %>%
  dplyr::arrange(coef) %>%
  dplyr::mutate(
    feature = factor(feature, levels = feature),
    signif = ifelse(qval < 0.05, "significant", "ns"),
    direction = ifelse(coef > 0, "Nest", "Control"),
    color_group = ifelse(signif == "ns", "ns", direction),
    in_bee = ifelse(feature %in% bee_genus, "yes", "no")
  )

ggplot(plot_df, aes(x = coef, y = feature, col = color_group, shape = in_bee)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  scale_colour_manual(values = c("Nest" = "gold3",
                                 "Control" = "skyblue3",
                                 "ns" = "lightgrey")) +
  scale_shape_manual(values = c("yes" = 16, "no"  = 4)) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none"
  ) +
  labs(x = "Effect size", y = NULL)

ggsave("figures/Maaslin_Bact.pdf", width = 15, height = 30, units = "cm")



## 2.5 - B-Diversity ----

# NMDS
samples_nc <- metadata %>%
  filter(Type %in% c("Nest", "Control")) %>%
  pull(ID)

bact_matrix2 <- bact %>%
  select(matches("^[A-Z]{2}[0-9]$")) %>%
  select(all_of(samples_nc)) %>%
  t()

set.seed(123)

nmds <- metaMDS(bact_matrix2, distance = "bray", k = 2, try = 999, trymax = 9999)


pdf("figures/CHEM-NMDS-Stressplot.pdf")
stressplot(nmds)
dev.off()
nmds$stress # 0.114 Good


# Extract the scores to plot
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("ID") %>%
  left_join(metadata[, c("ID", "Type", "Site_code", "Pair")], by = "ID")


# Plot
hull_df <- nmds_scores %>%
  group_by(Site_code, Type) %>%
  dplyr::slice(chull(NMDS1, NMDS2))


ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_polygon(data = hull_df,
               aes(color = Type, group = interaction(Site_code, Type)),
               fill = NA, linewidth = 0.7) +
  geom_point(aes(fill = Type, shape = Site_code), size = 3, colour = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  scale_shape_manual(values = my_points)

ggsave("figures/16S-NMDS.pdf", height = 15, width = 15, units = "cm")



# Centroid distances
centroids <- nmds_scores %>%
  group_by(Site_code, Type) %>%
  summarise(c1 = mean(NMDS1), c2 = mean(NMDS2), .groups = "drop")

dist_df <- nmds_scores %>%
  left_join(centroids, by = c("Site_code", "Type")) %>%
  mutate(dist = sqrt((NMDS1 - c1)^2 + (NMDS2 - c2)^2))

# Mean per group
dist_summary <- dist_df %>%
  group_by(Site_code, Type) %>%
  summarise(mean_dist = mean(dist), .groups = "drop")

# T-test Nest vs Control
nest    <- dist_summary %>% filter(Type == "Nest")    %>% pull(mean_dist)
control <- dist_summary %>% filter(Type == "Control") %>% pull(mean_dist)
ttest_dist <- t.test(nest, control, paired = TRUE) %>% broom::tidy()

# Plot
dist_plot <- dist_summary %>%
  group_by(Type) %>%
  summarise(
    mean = mean(mean_dist),
    se   = sd(mean_dist) / sqrt(n()),
    .groups = "drop"
  )


ggplot(dist_plot, aes(x = Type, y = mean, fill = Type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") +
  labs(x = NULL, y = "Mean distance to centroid")

ggsave("figures/16S-Dispersion.pdf", width = 10, height = 15, units = "cm")


write.xlsx(ttest_dist, "results/16S-Dispersion-Test.xlsx")
write.xlsx(dist_summary, "results/16S-Dispersion.xlsx")


# PERMANOVA
metadata2 <- metadata %>%
  filter(Type %in% c("Nest", "Control"))

permanova <- adonis2(bact_matrix2 ~ Type * Site_code, data = chem,
                     method = "euclidean", permutations = 9999, by = "terms")

disp <- betadisper(vegdist(bact_matrix2, method = "bray"), metadata2$Type)


# Results
sink("results/16S-PERMANOVA.txt")
anova(disp) #Ok
permutest(disp) #OK
tapply(disp$distances, metadata2$Type, mean)
anova(disp)
TukeyHSD(disp)
permanova
sink()



#Finish
rm(list = setdiff(ls(), c("bact", "chem", "fung", "metadata", "my_palette", "my_points")))
dev.off()
beep()
cat("Bacteria analysed. \n\n\n")





#---------------------------------------------------------------------------#





# 3 - FUNGI ----

## 3.1 - Stacked barcharts ----

tax_levels <- c("Phylum", "Class", "Order", "Family", "Genus")

plots_ITS <- lapply(tax_levels, function(tax_level) {
  
  tax_sym <- sym(tax_level)
  
  # Prepare the data
  fung_long <- fung %>%
    pivot_longer(cols = matches("^[A-Z]{2}[0-9]$"),
                 names_to = "ID", values_to = "Abundance") %>%
    left_join(metadata[, c("ID", "Type")], by = "ID") %>%
    mutate(!!tax_sym := as.character(!!tax_sym))
  
  # Top 10 taxa plotted
  top_taxa <- fung_long %>%
    dplyr::group_by(!!tax_sym) %>%
    summarise(Total = sum(Abundance)) %>%
    slice_max(Total, n = 9) %>%
    pull(!!tax_sym)
  
  df <- fung_long %>%
    mutate(!!tax_sym := ifelse(!!tax_sym %in% top_taxa, !!tax_sym, "Other")) %>%
    group_by(ID, Type, !!tax_sym) %>%
    summarise(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(ID) %>%
    mutate(RelAbund = Abundance / sum(Abundance)) %>% 
    ungroup()
  
  df %>% dplyr::count(ID) %>% print(n = Inf)
  
  # Colour palette
  taxa_names <- c(sort(unique(df[[tax_level]][df[[tax_level]] != "Other"])), "Other")
  my_colors <- setNames(
    c(viridis::viridis(length(taxa_names) - 1), "grey60"),
    taxa_names
  )
  
  # Plot
  p <- ggplot(df, aes(x = ID, y = RelAbund, fill = !!tax_sym)) +
    geom_bar(stat = "identity", color = "black", linewidth = 0.3, alpha = 0.9) +
    facet_wrap(~ Type, scales = "free_x") +
    scale_y_continuous(labels = scales::percent) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "bottom") +
    labs(y = "Relative Abundance", x = NULL,
         title = paste("fungerial community composition by", tax_level, "(only top-10 taxa showed)")) +
    scale_fill_manual(values = my_colors, name = tax_level)
  
  ggsave(filename = paste0("figures/ITS-", tax_level, ".pdf"),
         width = 25, height = 20, units = "cm", plot = p)
  
  return(p)
  
})


# Heatmap

heat_fung <- fung %>%
  pivot_longer(cols = matches("^[A-Z]{2}[0-9]$"), names_to = "ID", values_to = "Abundance") %>%
  left_join(metadata[, c("ID", "Type")], by = "ID") %>%
  mutate(Genus = as.character(Genus)) %>%
  group_by(Type, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  ungroup()

top_genera_fung <- heat_fung %>%
  filter(Type == "Bee") %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  slice_max(Total, n = 10) %>%
  pull(Genus)

heat_fung %>%
  filter(Genus %in% top_genera_fung) %>%
  mutate(log10_abundance = log10(Abundance + 1)) %>%
  ggplot(aes(x = Type, y = Genus, fill = log10_abundance)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient2(low = "white", mid = "pink2", high = "red4",
                       midpoint = 1.8,
                       name = "log10\n(Reads + 1)") +
  theme_bw() +
  theme(axis.text.y = element_text(face = "italic", size = 10),
        axis.text.x = element_text(size = 12),
        legend.position = "right") +
  labs(x = NULL, y = NULL, title = "Top-10 Fungal genera in bees")

ggsave("figures/ITS-Heatmap.pdf", width = 15, height = 15, units = "cm")



## 3.2 - A-Dversity models ----

# OTU matrix
fung_matrix <- fung %>%
  select(matches("^[A-Z]{2}[0-9]$")) %>%
  t()

# Calculate indices
alpha_div <- data.frame(
  ID       = rownames(fung_matrix),
  Richness = specnumber(fung_matrix),
  Shannon  = diversity(fung_matrix, index = "shannon"),
  Pielou_Evenness = diversity(fung_matrix, index = "shannon") / log(specnumber(fung_matrix))
)


# Build dataset for models
df_model <- alpha_div %>%
  left_join(metadata[, c("ID", "Type", "Pair", "Site_code")], by = "ID")

bio_vars <- colnames(df_model[,2:4])


# Calculate and export means
means_clean <- df_model %>%
  pivot_longer(cols = all_of(bio_vars), names_to = "variable", values_to = "value") %>%
  group_by(variable, Type) %>%
  summarise(
    mean = mean(value, na.rm = T),
    sd   = sd(value, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(mean_sd = sprintf("%.3f ± %.3f", mean, sd)) %>%
  select(variable, Type, mean_sd) %>%
  pivot_wider(names_from = Type, values_from = c(mean_sd))


# Models
bio_results <- lapply(bio_vars, function(v) {
  f <- as.formula(paste0(v, " ~ Type + (1|Site_code/Pair)"))
  mod <- lmer(f, data = df_model)
  
  # Type-II ANOVA
  anova_res <- Anova(mod, type = "II") %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term == "Type") %>%
    transmute(
      variable  = v,
      Chisq     = round(Chisq, 3),
      df        = Df,
      p.value   = round(`Pr(>Chisq)`, 3)
    )
  
  # Post-hoc pairwise comparisons
  posthoc <- emmeans(mod, pairwise ~ Type, adjust = "tukey")$contrasts %>%
    as.data.frame() %>%
    transmute(
      variable  = v,
      contrast  = contrast,
      estimate  = round(estimate, 3),
      std.error = round(SE, 3),
      df        = round(df, 3),
      statistic = round(t.ratio, 3),
      p.value   = round(p.value, 3)
    )
  
  bind_rows(anova_res, posthoc)
}) %>% 
  bind_rows() %>%
  mutate(p_fdr = round(p.adjust(p.value, method = "fdr"), 3))



## 3.3 - A-Dversity boxplots ----

plots <- lapply(bio_vars, function(varname){
  
  #Prepare the data
  data <- df_model %>%
    dplyr::mutate(
      Type = Type,
      Site = Site_code
    )
  
  #Calculate medians
  site_medians <- data %>%
    dplyr::group_by(Site, Type) %>%
    dplyr::summarise(
      value = median(.data[[varname]], na.rm = T),
      .groups = "drop"
    )
  
  lab <- site_medians %>%
    dplyr::filter(Type == "Control")
  
  #Make the plot
  p <- ggplot(data, aes(x = Type, y = .data[[varname]])) +
    geom_boxplot(aes(color = Type), outlier.shape = NA) +
    geom_point(aes(color = Type, shape = Site),
               position = position_jitter(width = 0.1),
               alpha = 0.7, size = 1.5) +
    geom_point(data = site_medians,
               aes(x = Type, y = value, shape = Site),
               size = 3, color = "black") +
    geom_line(data = site_medians,
              aes(x = Type, y = value, group = Site),
              inherit.aes = FALSE, alpha = 0.6, color = "grey35") +
    scale_color_manual(values = my_palette) +
    scale_shape_manual(values = my_points) +
    theme_bw() +
    theme(text = element_text(size = 20)
          , axis.text = element_text(size = 15, color = "black")
          , legend.position = "none") +
    labs(x = NULL, y = varname)
  
  #Save the plot
  ggsave(filename = paste0("figures/ITS-", varname, ".pdf"),
         width = 10, height = 20, units = "cm", plot = p)
  
  return(p)
  
})


# Save results
write.xlsx(bio_results, "results/ITS-LMM.xlsx")
write.xlsx(means_clean, "results/ITS-SumStat.xlsx")



## 1.4 MAASLIN ----

# Get the OTUs
get_otus <- function(type) {
  samples <- metadata %>% 
    dplyr::filter(Type == type) %>% 
    dplyr::pull(ID) %>%
    intersect(colnames(fung))  # <- solo campioni presenti in fung
  fung %>%
    select(OTU, all_of(samples)) %>%
    filter(rowSums(select(., -OTU)) > 0) %>%
    pull(OTU)
}

otus_control <- get_otus("Control")
otus_bee     <- get_otus("Bee")
otus_nest    <- get_otus("Nest")


# Venn Diagram
venn.diagram(
  x = list(Control = otus_control, Bee = otus_bee, Nest = otus_nest),
  filename = "figures/ITS-Venn.tiff",
  output = T,
  fill = c(my_palette["Control"], my_palette["Bee"], my_palette["Nest"]),
  alpha = 0.4,
  cex = 1.5,
  fontface = "bold",
  cat.cex = 0,
  rotation.degree = 240)



# OTU matrix
otu_mat <- fung %>%
  tibble::column_to_rownames("OTU") %>%
  select(matches("^[A-Z]{2}[0-9]$")) %>%
  as.matrix()


# Prepare the data
otu_table <- fung[, !(names(fung) %in% c("Phylum","Class","Order","Family","Genus","Species"))]
otu_table$Genus <- fung$Genus

genus_table <- otu_table %>%
  dplyr::group_by(Genus) %>%
  dplyr::summarise(dplyr::across(where(is.numeric), sum)) %>%
  as.data.frame() %>%
  column_to_rownames("Genus")

genus_table_t <- as.data.frame(t(genus_table))


# Metadata
metadata2 <- metadata[metadata$Type %in% c("Nest","Control"), ]
rownames(metadata2) <- metadata2$ID
genus_table_t <- genus_table_t[rownames(genus_table_t) %in% rownames(metadata2), ]


# Relative abundance
genus_rel <- genus_table_t / rowSums(genus_table_t)
genus_rel[is.na(genus_rel)] <- 0
prev <- colSums(genus_rel > 0) / nrow(genus_rel)
abund <- colMeans(genus_rel)


# Prevalence (>=50%) + Abundance (>=0.1%)
keep <- prev >= 0.5 & abund >= 0.001
genus_table_t <- genus_table_t[, keep]
metadata2 <- metadata2[match(rownames(genus_table_t), metadata2$ID), ]
metadata2$Type <- relevel(as.factor(metadata2$Type), ref = "Control")


# Run the analysis
maaslin_genus <- Maaslin2(
  input_data = genus_table_t,
  input_metadata = metadata2,
  output = "Genus-Fung",
  fixed_effects = "Type",
  random_effects = "Pair",
  normalization = "CLR",
  transform = "NONE",
  correction = "BH"
)


# Save the results
res <- maaslin_genus$results
write.xlsx(res, "results/Maaslin-Fung.xlsx")


# Plot the results
bee_samples <- metadata$ID[metadata$Type == "Bee"]
bee_samples <- intersect(bee_samples, colnames(genus_table))
bee_genus <- rownames(genus_table)[
  rowSums(genus_table[, bee_samples] > 0) > 0]

plot_df <- res %>%
  dplyr::arrange(coef) %>%
  dplyr::mutate(
    feature = factor(feature, levels = feature),
    signif = ifelse(qval < 0.05, "significant", "ns"),
    direction = ifelse(coef > 0, "Nest", "Control"),
    color_group = ifelse(signif == "ns", "ns", direction),
    in_bee = ifelse(feature %in% bee_genus, "yes", "no")
  )

ggplot(plot_df, aes(x = coef, y = feature, col = color_group, shape = in_bee)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = coef - stderr, xmax = coef + stderr), height = 0.2) +
  scale_colour_manual(values = c("Nest" = "gold3",
                                 "Control" = "skyblue3",
                                 "ns" = "lightgrey")) +
  scale_shape_manual(values = c("yes" = 16, "no"  = 4)) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    axis.text = element_text(size = 10, color = "black"),
    axis.text.y = element_text(face = "italic"),
    legend.position = "none"
  ) +
  labs(x = "Effect size", y = NULL)

ggsave("figures/Maaslin_Fung.pdf", width = 15, height = 30, units = "cm")



## 1.5 - B-Diversity ----

# NMDS
samples_nc <- metadata %>%
  filter(Type %in% c("Nest", "Control")) %>%
  pull(ID)

fung_matrix2 <- fung %>%
  select(matches("^[A-Z]{2}[0-9]$")) %>%
  select(all_of(samples_nc)) %>%
  t()

set.seed(123)

nmds <- metaMDS(fung_matrix2, distance = "bray", k = 2, try = 999, trymax = 9999)


pdf("figures/ITS-NMDS-Stressplot.pdf")
stressplot(nmds)
dev.off()
nmds$stress # 0.171 Good


# Extract the scores to plot
nmds_scores <- as.data.frame(scores(nmds, display = "sites")) %>%
  rownames_to_column("ID") %>%
  left_join(metadata[, c("ID", "Type", "Site_code", "Pair")], by = "ID")


# Plot
hull_df <- nmds_scores %>%
  group_by(Site_code, Type) %>%
  dplyr::slice(chull(NMDS1, NMDS2))


ggplot(nmds_scores, aes(NMDS1, NMDS2)) +
  geom_polygon(data = hull_df,
               aes(color = Type, group = interaction(Site_code, Type)),
               fill = NA, linewidth = 0.7) +
  geom_point(aes(fill = Type, shape = Site_code), size = 3, colour = "black") +
  theme_bw() +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 15, color = "black"),
        legend.position = "none") +
  scale_fill_manual(values = my_palette) +
  scale_color_manual(values = my_palette) +
  scale_shape_manual(values = my_points)

ggsave("figures/ITS-NMDS.pdf", height = 15, width = 15, units = "cm")


# Centroid distances
centroids <- nmds_scores %>%
  group_by(Site_code, Type) %>%
  summarise(c1 = mean(NMDS1), c2 = mean(NMDS2), .groups = "drop")

dist_df <- nmds_scores %>%
  left_join(centroids, by = c("Site_code", "Type")) %>%
  mutate(dist = sqrt((NMDS1 - c1)^2 + (NMDS2 - c2)^2))

# Mean per group
dist_summary <- dist_df %>%
  group_by(Site_code, Type) %>%
  summarise(mean_dist = mean(dist), .groups = "drop")

# T-test Nest vs Control
nest    <- dist_summary %>% filter(Type == "Nest")    %>% pull(mean_dist)
control <- dist_summary %>% filter(Type == "Control") %>% pull(mean_dist)
ttest_dist <- t.test(nest, control, paired = TRUE) %>% broom::tidy()

# Plot
dist_plot <- dist_summary %>%
  group_by(Type) %>%
  summarise(
    mean = mean(mean_dist),
    se   = sd(mean_dist) / sqrt(n()),
    .groups = "drop"
  )


ggplot(dist_plot, aes(x = Type, y = mean, fill = Type)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3, width = 0.6, alpha = 0.7) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.15) +
  scale_fill_manual(values = my_palette) +
  theme_bw() +
  theme(text = element_text(size = 15),
        axis.text = element_text(size = 12, color = "black"),
        legend.position = "none") +
  labs(x = NULL, y = "Mean distance to centroid")

ggsave("figures/ITS-Dispersion.pdf", width = 10, height = 15, units = "cm")


write.xlsx(ttest_dist, "results/ITS-Dispersion-Test.xlsx")
write.xlsx(dist_summary, "results/ITS-Dispersion.xlsx")


# PERMANOVA

metadata2 <- metadata %>%
  filter(Type %in% c("Nest", "Control"))

permanova <- adonis2(fung_matrix2 ~ Type * Site_code, data = chem,
                     method = "euclidean", permutations = 9999, by = "terms")

disp <- betadisper(vegdist(fung_matrix2, method = "bray"), metadata2$Type)


# Results
sink("results/ITS-PERMANOVA.txt")
anova(disp) #Ok
permutest(disp) #OK
tapply(disp$distances, metadata2$Type, mean)
anova(disp)
TukeyHSD(disp)
permanova
sink()



#Finish
rm(list = setdiff(ls(), c("bact", "chem", "fung", "metadata", "my_palette", "my_points")))
dev.off()
beep()
cat("Fungi analysed. \n\n\n")





#---------------------------------------------------------------------------#




# 4 - EXTRA ----

# Table for shared bacterial and fungal genera
bee_samples  <- metadata$ID[metadata$Type == "Bee"]
nest_samples <- metadata$ID[metadata$Type == "Nest"]
ctrl_samples <- metadata$ID[metadata$Type == "Control"]

make_presence_table <- function(df, kingdom_name = NA) {
  
  if (!"Kingdom" %in% colnames(df)) {
    df$Kingdom <- kingdom_name
  }
  
  bee_cols  <- intersect(bee_samples,  colnames(df))
  nest_cols <- intersect(nest_samples, colnames(df))
  ctrl_cols <- intersect(ctrl_samples, colnames(df))
  
  df %>%
    dplyr::group_by(Genus) %>%
    dplyr::summarise(
      Kingdom = Kingdom[1],
      Phylum  = Phylum[1],
      Class   = Class[1],
      Order   = Order[1],
      Family  = Family[1],
      Bee     = any(rowSums(dplyr::across(dplyr::all_of(bee_cols))  > 0) > 0),
      Nest    = any(rowSums(dplyr::across(dplyr::all_of(nest_cols)) > 0) > 0),
      Control = any(rowSums(dplyr::across(dplyr::all_of(ctrl_cols)) > 0) > 0),
      .groups = "drop"
    ) %>%
    dplyr::select(Kingdom, Phylum, Class, Order, Family, Genus, Bee, Nest, Control)
}

# Make the tables
bact_genus_presence  <- make_presence_table(bact,  kingdom_name = "Bacteria")
fungi_genus_presence <- make_presence_table(fung, kingdom_name = "Fungi")

openxlsx::write.xlsx(
  list(
    Bacteria = bact_genus_presence,
    Fungi    = fungi_genus_presence
  ),
  file = "results/Taxonomy-Genus-Shared.xlsx",
  overwrite = TRUE
)





#---------------------------------------------------------------------------#








#FINISH
gc()
beep(3)
.end_time <- Sys.time()
elapsed <- difftime(.end_time, .start_time, units = "mins")
cat("Script run in", round(elapsed, 2), "minutes\n\n\n")
rm(list = ls())
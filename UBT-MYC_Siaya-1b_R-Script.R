#****************************************************************************************************************************************************************#
#***********************************************************************Siaya-1b R-Script************************************************************************#
#********************************************************************Gkoutselis et al. 2023**********************************************************************#
#****************************************************************************************************************************************************************#

#**********************************************************************General information***********************************************************************#
# metacommunity: entire OTUs (community/mycobiome) found on MP or in the soil
# each metacommunity is equal-depth rarefied as implemented in the NCM calculations
# entire dataset: unrarefied: 2441 OTUs; rarefied 2162 OTUs
# plastisphere metacommunity: unrarefied: 1174 OTUs; rarefied: 927 OTUs
# soil metacommunity: rarefied: 2194; unrarefied: 1991 OTUs
# both rarefied metacommunities combined: 2208 OTUs

#**********************************************************Multiplot 1: Nested and idiosyncratic diversity*******************************************************#

# relevant data sets (C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunites)
# OTU_p_sites: rarefied, site-collapsed plastisphere metacommunity
# OTU_s_sites: rarefied, site-collapsed soil metacommunity
# OTU_p_specific_sites: rarefied, plastic-specific & site-collapsed plastisphere  metacommunity 
# OTU_s_specific_sites: rarefied, plastic-specific & site-collapsed soil  metacommunity 
# OTU_ps_dikarya: rarefied, plastic-specific & site-collapsed plastisphere community only including dikarya
# OTU_ss_dikarya: rarefied, plastic-specific & site-collapsed soil community only including dikarya
# shared: OTUs shared between soil and plastisphere metacommunity
# fac_sites: minimalistic, site-collapsed factor table (no other factors included; same for plastic, soil & shared)
# tax_p_part-fract: fully annotated taxonomy of the plastisphere metacommunity including partitions and fractions
# tax_s_part-fract: fully annotated taxonomy of the soil metacommunity including partitions and fractions

# Shannon diversity boxplots (BXPs) of classes, guilds & biosphere fractions (abundance/prevalence fractions)
# multiple boxplots in one figure
# alpha diversity according to phyla and guilds
# library calls
library(reshape2)
library(ggplot2)
library(viridis)
library(ggpubr)
library(phyloseq)
library(microbiome)

# data import
# plastisphere metacommunity
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Partitions/OTU_p_partition.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/fac_p.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Partitions/tax_p_partition.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

S1b_p = phyloseq(OTU, TAX, FAC)
S1b_p

rank_names(S1b_p)
sample_variables(S1b_p)
S1b_p
summarize_phyloseq(S1b_p)

Site <- factor(sample_data(S1b_p)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# calculate the Shannon index first and export the file
# guilds
# pathogens
S1b_p_path = subset_taxa(S1b_p, Primary_lifestyle == "pathogen")
H_p_path <- estimate_richness(S1b_p_path, measures = c("Shannon"))
write.csv(H_p_path, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_path.txt"))

# saprotrophs
S1b_p_sap = subset_taxa(S1b_p, Primary_lifestyle == "saprotroph")
H_p_sap <- estimate_richness(S1b_p_sap, measures = c("Shannon"))
write.csv(H_s_mut, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_sap.txt"))

# mutualists
S1b_p_mut = subset_taxa(S1b_p, Primary_lifestyle == "mutualist")
H_p_mut <- estimate_richness(S1b_p_mut, measures = c("Shannon"))
write.csv(H_s_mut, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_mut.txt"))

# unidentified
S1b_p_unid = subset_taxa(S1b_p, Primary_lifestyle == "unidentified")
H_p_unid <- estimate_richness(S1b_p_mut, measures = c("Shannon"))
write.csv(H_s_mut, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_unid.txt"))

# classes
# Agaricomycetes
S1b_p_agar = subset_taxa(S1b_p, Class == "Agaricomycetes")
H_p_agar <- estimate_richness(S1b_p_agar, measures = c("Shannon"))
write.csv(H_p_agar, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_agar.txt"))

# Tremellomycetes
S1b_p_trem = subset_taxa(S1b_p, Class == "Tremellomycetes")
H_p_trem <- estimate_richness(S1b_p_trem, measures = c("Shannon"))
write.csv(H_p_trem, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_trem.txt"))

# Eurotiomycetes
S1b_p_euro = subset_taxa(S1b_p, Class == "Eurotiomycetes")
H_p_euro <- estimate_richness(S1b_p_euro, measures = c("Shannon"))
write.csv(H_p_euro, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_euro.txt"))

# Sordariomycetes
S1b_p_sord = subset_taxa(S1b_p, Class == "Sordariomycetes")
H_p_sord <- estimate_richness(S1b_p_sord, measures = c("Shannon"))
write.csv(H_p_sord, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_sord.txt"))

# Dothideomycetes
S1b_p_doth = subset_taxa(S1b_p, Class == "Dothideomycetes")
H_p_doth <- estimate_richness(S1b_p_doth, measures = c("Shannon"))
write.csv(H_p_doth, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_doth.txt"))

# other (unidentified + remaining classes)
S1b_p_other = subset_taxa(S1b_p, Class_top == "other")
H_p_other <- estimate_richness(S1b_p_other, measures = c("Shannon"))
write.csv(H_p_other, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_other.txt"))

# phyla
# Ascomycota
S1b_p_asco = subset_taxa(S1b_p, Phylum == "Ascomycota")
H_p_asco <- estimate_richness(S1b_p_asco, measures = c("Shannon"))
write.csv(H_p_asco, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_asco.txt"))

# Basidiomycota
S1b_p_basi = subset_taxa(S1b_p, Phylum == "Basidiomycota")
H_p_basi <- estimate_richness(S1b_p_basi, measures = c("Shannon"))
write.csv(H_p_basi, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_basi.txt"))

# Chytridiomycota
S1b_p_chyt = subset_taxa(S1b_p, Phylum == "Chytridiomycota")
H_p_chyt <- estimate_richness(S1b_p_chyt, measures = c("Shannon"))
write.csv(H_p_chyt, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_chyt.txt"))

# Glomeromycota
S1b_p_glom = subset_taxa(S1b_p, Phylum == "Glomeromycota")
H_p_glom <- estimate_richness(S1b_p_glom, measures = c("Shannon"))
write.csv(H_p_glom, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_glom.txt"))

# Mucoromycota
S1b_p_muco = subset_taxa(S1b_p, Phylum == "Mucoromycota")
H_p_muco <- estimate_richness(S1b_p_muco, measures = c("Shannon"))
write.csv(H_p_muco, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_muco.txt"))

# Unidentified
S1b_p_uni = subset_taxa(S1b_p, Phylum == "unidentified")
H_p_uni <- estimate_richness(S1b_p_uni, measures = c("Shannon"))
write.csv(H_p_uni, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_p_uni.txt"))


# render the boxplots
# guilds
df_p_guild_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_p_guilds_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

p_guild <- ggboxplot(df_p_guild_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("firebrick2", "gold", "darkolivegreen1", "azure3")) + ylim(0, 4)
p_guild

# phyla
df_p_phyla_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_p_phyla_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

p.phylum.2 <- ggboxplot(df_p_phyla_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("royalblue4", "darkcyan", "chartreuse4", "palegreen2", "goldenrod3","azure3"))
p.phylum.2 + ylim(0, 4)

# classes
df_p_classes_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_p_classes_long.txt", 
                                header = TRUE, check.names = FALSE, row.names = 1)

p_classes <- ggboxplot(df_p_classes_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("deepskyblue3", "goldenrod1", "springgreen2", "seagreen3", "palegreen4", "khaki2"))
p_classes + ylim(0, 4)

# fractions
df_p_fract_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_p_fractions_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

p_fract <- ggplot(df_p_fract_long, aes(x = Variable, y = Value, fill = Variable, legend = "right")) + geom_boxplot() +
  scale_fill_viridis(discrete=TRUE, guide="none", option = "D") + theme_bw() + theme_classic() + ylim(0, 4)
p_fract


# soil metacommunity
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Partitions/OTU_s_partition.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/fac_s.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Partitions/tax_s_partition.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

S1b_s = phyloseq(OTU, TAX, FAC)
S1b_s

rank_names(S1b_s)
sample_variables(S1b_s)
S1b_s
summarize_phyloseq(S1b_s)

Site <- factor(sample_data(S1b_s)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# calculate the Shannon index first and export the file
# guilds
# pathogens
S1b_s_path = subset_taxa(S1b_s, Primary_lifestyle == "pathogen")
H_s_path <- estimate_richness(S1b_s_path, measures = c("Shannon"))
write.csv(H_s_path, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_path.txt"))

# saprotrophs
S1b_s_sap = subset_taxa(S1b_s, Primary_lifestyle == "saprotroph")
H_s_sap <- estimate_richness(S1b_s_sap, measures = c("Shannon"))
write.csv(H_s_sap, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_sap.txt"))

# mutualists
S1b_s_mut = subset_taxa(S1b_s, Primary_lifestyle == "mutualist")
H_s_mut <- estimate_richness(S1b_s_mut, measures = c("Shannon"))
write.csv(H_s_mut, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_mut.txt"))

# unidentified
S1b_s_unid = subset_taxa(S1b_s, Primary_lifestyle == "unidentified")
H_s_unid <- estimate_richness(S1b_s_mut, measures = c("Shannon"))
write.csv(H_s_unid, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_unid.txt"))

# classes
# Agaricomycetes
S1b_s_agar = subset_taxa(S1b_s, Class == "Agaricomycetes")
H_s_agar <- estimate_richness(S1b_s_agar, measures = c("Shannon"))
write.csv(H_s_agar, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_agar.txt"))

# Tremellomycetes
S1b_s_trem = subset_taxa(S1b_s, Class == "Tremellomycetes")
H_s_trem <- estimate_richness(S1b_s_trem, measures = c("Shannon"))
write.csv(H_s_trem, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_trem.txt"))

# Eurotiomycetes
S1b_s_euro = subset_taxa(S1b_s, Class == "Eurotiomycetes")
H_s_euro <- estimate_richness(S1b_s_euro, measures = c("Shannon"))
write.csv(H_s_euro, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_euro.txt"))

# Sordariomycetes
S1b_s_sord = subset_taxa(S1b_s, Class == "Sordariomycetes")
H_s_sord <- estimate_richness(S1b_s_sord, measures = c("Shannon"))
write.csv(H_s_sord, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_sord.txt"))

# Dothideomycetes
S1b_s_doth = subset_taxa(S1b_s, Class == "Dothideomycetes")
H_s_doth <- estimate_richness(S1b_s_doth, measures = c("Shannon"))
write.csv(H_s_doth, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_doth.txt"))

# other (unidentified + remaining classes)
S1b_s_other = subset_taxa(S1b_s, Class_top == "other")
H_s_other <- estimate_richness(S1b_s_other, measures = c("Shannon"))
write.csv(H_s_other, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_other.txt"))

# phyla
# Ascomycota
S1b_s_asco = subset_taxa(S1b_s, Phylum == "Ascomycota")
H_s_asco <- estimate_richness(S1b_s_asco, measures = c("Shannon"))
write.csv(H_s_asco, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_asco.txt"))

# Basidiomycota
S1b_s_basi = subset_taxa(S1b_s, Phylum == "Basidiomycota")
H_s_basi <- estimate_richness(S1b_s_basi, measures = c("Shannon"))
write.csv(H_s_basi, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_basi.txt"))

# Chytridiomycota
S1b_s_chyt = subset_taxa(S1b_s, Phylum == "Chytridiomycota")
H_s_chyt <- estimate_richness(S1b_s_chyt, measures = c("Shannon"))
write.csv(H_s_chyt, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_chyt.txt"))

# Glomeromycota
S1b_s_glom = subset_taxa(S1b_s, Phylum == "Glomeromycota")
H_s_glom <- estimate_richness(S1b_s_glom, measures = c("Shannon"))
write.csv(H_s_glom, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_glom.txt"))

# Mucoromycota
S1b_s_muco = subset_taxa(S1b_s, Phylum == "Mucoromycota")
H_s_muco <- estimate_richness(S1b_s_muco, measures = c("Shannon"))
write.csv(H_s_muco, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_muco.txt"))

# Unidentified
S1b_s_uni = subset_taxa(S1b_s, Phylum == "unidentified")
H_s_uni <- estimate_richness(S1b_s_uni, measures = c("Shannon"))
write.csv(H_s_uni, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/H_s_uni.txt"))

# render the boxplots
# guilds
df_s_guild_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_s_guilds_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

s_guild <- ggboxplot(df_s_guild_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("firebrick2", "gold", "darkolivegreen1", "azure3"))
s_guild + ylim(0, 4)


# phyla
df_s_phyla_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_s_phyla_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

s_phylum <- ggboxplot(df_s_phyla_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("royalblue4", "darkcyan", "chartreuse4", "palegreen2", "goldenrod3","azure3"))
s_phylum + ylim(0, 4)

# classes
df_s_classes_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_s_classes_long.txt", 
                                header = TRUE, check.names = FALSE, row.names = 1)

s_classes <- ggboxplot(df_s_classes_long, x = "Variable", y = "Value", fill = "Variable", scales = "free") + 
  scale_fill_manual(values = c("deepskyblue3", "goldenrod1", "springgreen2", "seagreen3", "palegreen4", "khaki2"))
s_classes + ylim(0, 4)

# fractions
df_s_fract_long <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/Shannon_s_fractions_long.txt", 
                              header = TRUE, check.names = FALSE, row.names = 1)

s_fract <- ggplot(df_s_fract_long, aes(x = Variable, y = Value, fill = Variable, legend = "right")) + geom_boxplot() +
  scale_fill_viridis(discrete=TRUE, guide="none", option = "D") + theme_classic()
s_fract + ylim(0, 4)
#_________________________________________________________________________________________________________________________________________________________________
# compartment specific mycobiomes (plastic, soil and shared)
# taxa barplots (TBPs) of classes and guilds as well as biosphere fractions (abundance/prevalence fractions)

# library calls
library(phyloseq)
library(vegan)
library(microbiome)
library(ggplot2)
library(viridis)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

# plastic-specific
# data import
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/OTU_p_specific_sites.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_p_part-fract.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

# remove samples with 0 OTUs
OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_ps = phyloseq(OTU, TAX, FAC)
S1b_ps

rank_names(S1b_ps)
sample_variables(S1b_ps)
S1b_ps
summarize_phyloseq(S1b_ps)

Site <- factor(sample_data(S1b_ps)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# transform to relative abundance
rel_ps = transform_sample_counts(S1b_ps, function(x) x / sum(x))
rel_ps
get_sample(rel_ps)

# render the taxa barplots (TBP; width: 550, height: 450, no zoom)
# guilds
S1b_ps_guild = tax_glom(rel_ps, taxrank = "Primary_lifestyle", NArm = FALSE)
S1b_ps_guild
ps_guild <- plot_bar(S1b_ps_guild, fill = "Primary_lifestyle") + theme_bw() + scale_fill_manual(values = c("firebrick2", "gold", "darkolivegreen1", "azure3"))
ps_guild + geom_bar(aes(), stat="identity", position="stack")

# fractions
S1b_ps_fraction = tax_glom(rel_ps, taxrank = "Fraction", NArm = FALSE)
S1b_ps_fraction
ps_fraction <- plot_bar(S1b_ps_fraction, fill = "Fraction") + theme_bw() + scale_fill_viridis(discrete = TRUE)
ps_fraction + geom_bar(aes(), stat="identity", position="stack")

# plastic-specific dikarya
# data import
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/OTU_ps_dikarya.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_ps_dikarya.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_ps_dikarya = phyloseq(OTU, TAX, FAC)
S1b_ps_dikarya

rank_names(S1b_ps_dikarya)
sample_variables(S1b_ps_dikarya)
S1b_ps_dikarya
summarize_phyloseq(S1b_ps_dikarya)

Site <- factor(sample_data(S1b_ps_dikarya)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# transform to rel. abund.
rel_ps_dikarya = transform_sample_counts(S1b_ps_dikarya, function(x) x / sum(x) )
rel_ps_dikarya
get_sample(rel_ps_dikarya)

# classes
# "Class_top" is the descriptor to delineate between the top fungal classes and "other"
S1b_ps_class_dikarya = tax_glom(rel_ps_dikarya, taxrank = "Class_top", NArm = FALSE)
S1b_ps_class_dikarya
ps_class_dikarya <- plot_bar(S1b_ps_class_dikarya, fill = "Class_top") + theme_bw() + 
  scale_fill_manual(values = c("deepskyblue3", "goldenrod1", "goldenrod3", "springgreen2", "seagreen3", "palegreen4", "khaki2"))
ps_class_dikarya + geom_bar(aes(), stat="identity", position="stack")


# soil-specific
# import data
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/OTU_s_specific_sites.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_s_part-fract.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

# remove samples containing 0 OTUs
OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_ss = phyloseq(OTU, TAX, FAC)
S1b_ss

rank_names(S1b_ss)
sample_variables(S1b_ss)
S1b_ss
summarize_phyloseq(S1b_ss)

Site <- factor(sample_data(S1b_ss)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# transform to rel. abund.
rel_ss = transform_sample_counts(S1b_ss, function(x) x / sum(x) )
rel_ss
get_sample(rel_ss)

# guilds
S1b_ss_guild = tax_glom(rel_ss, taxrank = "Primary_lifestyle", NArm = FALSE)
S1b_ss_guild
ss_guild <- plot_bar(S1b_ss_guild, fill = "Primary_lifestyle") + theme_bw() + 
  scale_fill_manual(values = c("firebrick2", "gold", "darkolivegreen1", "azure3")) 
ss_guild + geom_bar(aes(), stat="identity", position="stack")

# fractions
S1b_ss_fraction = tax_glom(rel_ss, taxrank = "Fraction", NArm = FALSE)
S1b_ss_fraction
ss_fraction <- plot_bar(S1b_ss_fraction, fill = "Fraction") + theme_bw() +
  scale_fill_manual(values = c("#287D8EFF", "#73D055FF", "#FDE725FF"))
ss_fraction + geom_bar(aes(), stat="identity", position="stack")

# soil-specific dikarya
# data import
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/OTU_ss_dikarya.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_ss_dikarya.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_ss_dikarya = phyloseq(OTU, TAX, FAC)
S1b_ss_dikarya

rank_names(S1b_ss_dikarya)
sample_variables(S1b_ss_dikarya)
S1b_ss_dikarya
summarize_phyloseq(S1b_ss_dikarya)

Site <- factor(sample_data(S1b_ss_dikarya)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

rel_ss_dikarya = transform_sample_counts(S1b_ss_dikarya, function(x) x / sum(x) ) #creates a transformed phyloseq document with relative abundance data
rel_ss_dikarya
get_sample(rel_ss_dikarya)

# classes
S1b_ss_class_dikarya = tax_glom(rel_ss_dikarya, taxrank = "Class_top", NArm = FALSE)
S1b_ss_class_dikarya
ss_class_dikarya <- plot_bar(S1b_ss_class_dikarya, fill = "Class_top") + theme_bw() + 
  scale_fill_manual(values = c("deepskyblue3", "goldenrod1", "goldenrod3", "springgreen2", "seagreen3", "palegreen4", "khaki2"))
ss_class_dikarya + geom_bar(aes(), stat="identity", position="stack")


# Shared
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/shared.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_p_part-fract.txt", 
          header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_shared = phyloseq(OTU, TAX, FAC)
S1b_shared

rank_names(S1b_shared)
sample_variables(S1b_shared)
S1b_shared
summarize_phyloseq(S1b_shared)

Site <- factor(sample_data(S1b_shared)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

rel_shared = transform_sample_counts(S1b_shared, function(x) x / sum(x) ) #creates a transformed phyloseq document with relative abundance data
rel_shared
get_sample(rel_shared)

# guilds
S1b_shared_guild = tax_glom(rel_shared, taxrank = "Primary_lifestyle", NArm = FALSE)
S1b_shared_guild
shared_guild <- plot_bar(S1b_shared_guild, fill = "Primary_lifestyle") + theme_bw() + 
  scale_fill_manual(values = c("firebrick2", "gold", "darkolivegreen1", "azure3")) 
shared_guild + geom_bar(aes(), stat="identity", position="stack")

# fractions
S1b_shared_fraction = tax_glom(rel_shared, taxrank = "Fraction", NArm = FALSE)
S1b_shared_fraction
shared_fraction <- plot_bar(S1b_shared_fraction, fill = "Fraction") + theme_bw() +  scale_fill_viridis(discrete = TRUE)
shared_fraction + geom_bar(aes(), stat="identity", position="stack")

# shared dikarya
# data import
OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/shared.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/fac_sites.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Metacommunities/tax_shared.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

S1b_shared_dikarya = phyloseq(OTU, TAX, FAC)
S1b_shared_dikarya

rank_names(S1b_shared_dikarya)
sample_variables(S1b_shared_dikarya)
S1b_shared_dikarya
summarize_phyloseq(S1b_shared_dikarya)

Site <- factor(sample_data(S1b_shared_dikarya)$Site, levels = c("Landfill 1", "Roadside", "Market", "Courtyard", "Landfill 2"))
levels(Site)

# transform to rel. abund.
rel_shared_dikarya = transform_sample_counts(S1b_shared_dikarya, function(x) x / sum(x))
rel_shared_dikarya
get_sample(rel_shared_dikarya)

# classes
S1b_shared_class_dikarya = tax_glom(rel_shared_dikarya, taxrank = "Class_top", NArm = FALSE)
S1b_shared_class_dikarya
shared_class_dikarya <- plot_bar(S1b_shared_class_dikarya, fill = "Class_top") + theme_bw() + 
  scale_fill_manual(values = c("deepskyblue3", "goldenrod1", "goldenrod3", "springgreen2", "seagreen3", "palegreen4", "khaki2"))
shared_class_dikarya + geom_bar(aes(), stat="identity", position="stack")



#*************************************************************Multiplot 2: Variation in mycobiome composition***********************************************************#
# fractions of OTUs according to relative abundances
# fractioning of the rarefied metacommunities according to preset abundance and prevalence (occurrence frequency) thresholds
# fractioning performed separately for soil and MP; OTU tables consolidated in the aftermath

# relevant data sets
# AAT: OTU table of always abundant taxa with abundances >= 1% in all samples
# ART: OTU table of always rare taxa, with abundance < 0.01% in all samples
# CAT: OTU table of conditionally abundant taxa, which show abundances >= 0.01% in all & >= 1% in some samples
# CRAT: OTU table of conditionally rare and abundant taxa with abundances between <0.01% to >=1%) in all samples
# CRT: OTU table of conditionally rare taxa with abundances <1% in all and <0.01% in some samples
# MT: OTU table of moderate taxa, which show abundances between 0.01% und 1% in all samples
# fac: factors of the (original) entire data set
# tax: fully annotated taxonomy of the entire data set

# library calls
library("phyloseq")
library("ggplot2")
library("scales")
library("grid")
library(viridis)
library(microbiome)
library(vegan)
theme_set(theme_bw())

# AAT
# import from TSV
AAT <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/AAT.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/fac.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

# transform the three tables to a format that is convertible into Phyloseq
OTU = otu_table(AAT, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)
AAT = phyloseq(OTU, TAX, FAC) #combining the three tables into a single Phyloseq Document
AAT

rank_names(AAT)
sample_variables(AAT)
AAT
summarize_phyloseq(AAT)

Substrate <- factor(sample_data(AAT)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(AAT)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity and export fiile
AAT.H <- estimate_richness(AAT, measures = ("Shannon"))
AAT.H
write.csv(AAT.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/AAT.H.txt"))

# PCO plot based on Bray-Curtis dissimilarity
ord.PCoA.AAT <- ordinate(AAT, method="PCoA",k=3, distance="bray",trymax=20, wascores = TRUE, zerodist = "add")
PCOA.AAT <- plot_ordination(AAT, ord.PCoA.AAT, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3"))  + 
  ggtitle("AAT") + theme_bw() + geom_point(size = 5) +
  theme(panel.background = element_blank()) + stat_ellipse() 
PCOA.AAT

# ANOSIM
AAT.anosim <- phyloseq::distance(AAT, method = "bray")
AAT.anosim <- anosim(AAT.anosim, Factors$Substrate)
AAT.anosim

# transform to relative abundance
rel = transform_sample_counts(AAT, function(x) x / sum(x) ) #creates a transformed phyloseq document with relative abundance data
rel
get_sample(rel)

# taxa barplots
AAT.genus = tax_glom(rel, taxrank = "Genus", NArm = FALSE)
AAT.genus
plot_bar(AAT.genus, fill = "Genus") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + 
  scale_fill_brewer(palette = "Paired") + ggtitle("AAT")

# ART 
ART <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/ART.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/ART.fac.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt",
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(ART, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

# check for samples with zero counts
OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

# check for OTUs with zero counts
any(taxa_sums(OTU) == 0)
sum(taxa_sums(OTU) == 0) #gives the number of taxa with 0 read counts
OTU <- prune_taxa(taxa_sums(OTU) > 0, OTU) #remove all taxa with 0 read counts
any(taxa_sums(OTU) == 0)
OTU

ART = phyloseq(OTU, TAX, FAC)
ART

rank_names(ART)
sample_variables(ART)
ART
summarize_phyloseq(ART)

Substrate <- factor(sample_data(ART)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(ART)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity
ART.H <- estimate_richness(ART, measures = ("Shannon"))
ART.H
write.csv(ART.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/ART.H.txt"))

ord.PCoA.ART <- ordinate(ART, method="PCoA",k=3, distance="bray",trymax=200, wascores = TRUE, zerodist = "add")
PCOA.ART <- plot_ordination(ART, ord.PCoA.ART, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3"))  + 
  ggtitle("ART") + theme_bw() + geom_point(size = 5) + stat_ellipse()
PCOA.ART

ART.anosim <- phyloseq::distance(ART, method = "bray")
ART.anosim <- anosim(ART.anosim, Factors$Substrate)
ART.anosim

plot_richness(ART, x="Substrate", color="Substrate", measures = c("Observed", "Shannon", "Chao1")) + 
  geom_boxplot() + scale_color_manual(values = c("cornflowerblue", "darkred")) + theme_bw() + ggtitle("ART")

rel.ART = transform_sample_counts(ART, function(x) x / sum(x) )
rel.ART
get_sample(rel.ART)

ART.phylum = tax_glom(rel.ART, taxrank = "Phylum", NArm = FALSE)
ART.phylum
plot_bar(ART.phylum, fill = "Phylum") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + 
  scale_fill_brewer(palette = "Paired") 

# CAT
CAT <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/CAT.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/fac.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt",
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(CAT, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

CAT = phyloseq(OTU, TAX, FAC)
CAT
rank_names(CAT)
sample_variables(CAT)
CAT
summarize_phyloseq(CAT)

Substrate <- factor(sample_data(CAT)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(CAT)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity
CAT.H <- estimate_richness(CAT, measures = ("Shannon"))
CAT.H
write.csv(CAT.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/CAT.H.txt"))

ord.PCoA.CAT <- ordinate(CAT, method="PCoA",k=3, distance="bray",trymax=20, wascores = TRUE, zerodist = "add")
PCOA.CAT <- plot_ordination(CAT, ord.PCoA.CAT, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3"))  + ggtitle("CAT") + 
  theme_bw() + geom_point(size = 5) + stat_ellipse()
PCOA.CAT

CAT.anosim <- phyloseq::distance(CAT, method = "bray")
CAT.anosim <- anosim(CAT.anosim, Factors$Substrate)
CAT.anosim

plot_richness(CAT, x="Substrate", color="Substrate", measures = c("Observed", "Shannon", "Chao1")) + 
  geom_boxplot() + scale_color_manual(values = c("cornflowerblue", "darkred")) + theme_bw()

rel.CAT = transform_sample_counts(ART, function(x) x / sum(x) )
rel.CAT
get_sample(rel.CAT)

CAT.phylum = tax_glom(rel.CAT, taxrank = "Phylum", NArm = FALSE)
CAT.phylum
plot_bar(CAT.phylum, fill = "Phylum") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + 
  scale_fill_brewer(palette = "Paired") 

# CRAT
CRAT <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/CRAT.txt",
                   header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/CRAT.fac.txt",
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(CRAT, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU

CRAT = phyloseq(OTU, TAX, FAC)
CRAT

rank_names(CRAT)
sample_variables(CRAT)
CRAT
summarize_phyloseq(CRAT)

Substrate <- factor(sample_data(CRAT)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(CRAT)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity
CRAT.H <- estimate_richness(CRAT, measures = ("Shannon"))
CRAT.H
write.csv(CRAT.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/CRAT.H.txt"))

ord.PCoA.CRAT <- ordinate(CRAT, method="PCoA",k=3, distance="bray",trymax=20, wascores = TRUE, zerodist = "add")
PCOA.CRAT <- plot_ordination(CRAT, ord.PCoA.CRAT, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3"))  + ggtitle("CRAT") + 
  theme_bw() + geom_point(size = 5) + stat_ellipse()
PCOA.CRAT

CRAT.anosim <- phyloseq::distance(CRAT, method = "bray")
CRAT.anosim <- anosim(CRAT.anosim, Factors$Substrate)
CRAT.anosim

plot_richness(CRAT, x="Substrate", color="Substrate", measures = c("Observed", "Shannon", "Chao1")) + 
  geom_boxplot() + scale_color_manual(values = c("cornflowerblue", "darkred")) + theme_bw()

rel.CRAT = transform_sample_counts(CRAT, function(x) x / sum(x) )
rel.CRAT
get_sample(rel.CRAT)

CRAT.phylum = tax_glom(rel.CRAT, taxrank = "Phylum", NArm = FALSE)
CRAT.phylum
plot_bar(CRAT.phylum, fill = "Phylum") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + 
  scale_fill_brewer(palette = "Paired") 

# CRT
CRT <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/CRT.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/fac.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(CRT, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

CRT = phyloseq(OTU, TAX, FAC)
CRT
rank_names(CRT)

sample_variables(CRT)
CRT
summarize_phyloseq(CRT)

Substrate <- factor(sample_data(CRT)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(CRT)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity
CRT.H <- estimate_richness(CRT, measures = ("Shannon"))
CRT.H
write.csv(CRT.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/CRT.H.txt"))

ord.PCoA.CRT <- ordinate(CRT, method="PCoA",k=3, distance="bray",trymax=20, wascores = TRUE, zerodist = "add")
PCOA.CRT <- plot_ordination(CRT, ord.PCoA.CRT, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3")) + 
  ggtitle("CRT") + theme_bw() + geom_point(size = 5) + stat_ellipse()
PCOA.CRT

CRT.anosim <- phyloseq::distance(CRT, method = "bray")
CRT.anosim <- anosim(CRT.anosim, Factors$Substrate)
CRT.anosim

plot_richness(CRT, x="Substrate", color="Substrate", measures = c("Observed", "Shannon", "Chao1")) + 
  geom_boxplot() + scale_color_manual(values = c("cornflowerblue", "darkred")) + theme_bw()

rel.CRT = transform_sample_counts(CRT, function(x) x / sum(x) )
rel.CRT
get_sample(rel.CRT)

CRT.phylum = tax_glom(rel.CRT, taxrank = "Phylum", NArm = FALSE)
CRT.phylum
plot_bar(CRT.phylum, fill = "Phylum") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + 
  scale_fill_brewer(palette = "Paired") 

# MT
MT <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/MT.txt", 
                 header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/MT.fac.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

OTU = otu_table(MT, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

OTU <- prune_samples(sample_sums(OTU) > 0, OTU)
OTU
MT = phyloseq(OTU, TAX, FAC)
MT

rank_names(MT)
sample_variables(MT)
MT
summarize_phyloseq(MT)

Substrate <- factor(sample_data(MT)$Substrate, levels = c("plastic", "soil"))
levels(Substrate)
Site <- factor(sample_data(MT)$Site, levels = c("S-1", "S-2", "S-3", "S-4", "S-5"))
levels(Site)

# estimate Shannon diversity
MT.H <- estimate_richness(MT, measures = ("Shannon"))
MT.H
write.csv(MT.H, file = ("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/MT.H.txt"))

ord.PCoA.MT <- ordinate(MT, method="PCoA",k=3, distance="bray",trymax=20, wascores = TRUE, zerodist = "add")
PCOA.MT <- plot_ordination(MT, ord.PCoA.MT, color = "Substrate") + 
  scale_color_manual(values = c("cornflowerblue", "coral3"))  + 
  ggtitle("MT") + theme_bw() + geom_point(size = 5) + stat_ellipse()
PCOA.MT

MT.anosim <- phyloseq::distance(MT, method = "bray")
MT.anosim <- anosim(MT.anosim, Factors$Substrate)
MT.anosim

plot_richness(MT, x="Substrate", color="Substrate", measures = c("Observed", "Shannon", "Chao1")) + 
  geom_boxplot() + scale_color_manual(values = c("cornflowerblue", "darkred")) + theme_bw()

rel = transform_sample_counts(MT, function(x) x / sum(x) ) #creates a transformed phyloseq document with relative abundance data
rel
get_sample(rel)

MT.class = tax_glom(rel, taxrank = "Class", NArm = FALSE)
MT.class
plot_bar(MT.class, fill = "Class") + facet_wrap(~Substrate, scales= "free_x", nrow=1) + scale_fill_brewer(palette = "Paired") 

#_________________________________________________________________________________________________________________________________________________________________
# niche width
# correlation between niche width and mean abundance
# scatter plots with regression lines and point colors according to fraction

library(spaa)
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
theme_set(
  theme_bw() +
    theme(legend.position = "top")
)

# plastic
niche.p <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/niche.p.width.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
MP.niche <- niche.width(niche.p, method = c("shannon"))
# MP.niche <- niche.width(niche.p, method = c("levins"))

# export and import data
write.csv(MP.niche,"C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/MP.niche-width.txt", row.names = FALSE)
niche.p <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/Plastisphere_niche-width.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)

# render the scatter plot with regression line(ScP; width: 600, height: 450)
ggplot(niche.p, aes(x=Mean_abundance,y=Niche_width, color=Fraction)) + geom_point(size=2) + 
  scale_color_viridis(discrete=TRUE, guide=FALSE) + labs(x="Mean abundance [%]", y= "Niche width") +
  scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') + geom_smooth(method=lm , color="black", se=FALSE) 

# Spearman correlation and significance testing
cor(niche.p$Niche_width, niche.p$Mean_abundance, method = "spearman")
cor.test(niche.p$Niche_width, niche.p$Mean_abundance, method = "spearman")


# soil
niche.s <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/niche.s.width.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)
soil.niche <- niche.width(niche.s, method = c("shannon"))
# soil.niche <- niche.width(niche.s, method = c("levins"))

# export und import data
# write.csv(Soil.niche,"C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/soil.niche-width.txt", row.names = FALSE)
niche.s <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Fractions/Soil_niche-width.txt", 
           header = TRUE, check.names = FALSE, row.names = 1)

# render the scatter plot (ScP; width: 600, height: 450)
ggplot(niche.s, aes(x=Mean_abundance,y=Niche_width, color=Fraction)) +
  geom_point(size=2) + scale_color_viridis(discrete=TRUE, guide=FALSE, option = "D") + labs(x="Mean abundance [%]", y= "Niche width") +  
  scale_y_continuous(trans='log10') + scale_x_continuous(trans='log10') + geom_smooth(method=lm , color="black", se=FALSE)

# Spearman correlation and significance testing
cor(niche.s$Niche_width, niche.s$Mean_abundance, method = "spearman")
cor.test(niche.s$Niche_width, niche.s$Mean_abundance, method = "spearman")


#*************************************************************Plot 3: Phylogeny of plastic-specific fungi***************************************************************#
# phylogenetic analysis has been conducted with MEGA-X
# phylogram has been constructed with iTOL
# no R code present


#******************************************************************Multiplot 4: Neutral Community Model*****************************************************************#
# NCM according to Sloan (2006) as implemented in Burns et al. (2015)
# R code from Chen et al. (2019) & Burns et al. (2015)

# library calls
library(Hmisc)
library(minpack.lm)
library(stats4)
library(phyloseq)

# relevant data files - community table with taxa as rows and samples as columns
# spp: OTU table of the entire dataset, that is soil + plastic (unfiltered, unrarified, 2441 OTUs, 95 samples)
# spp_p: OTU table plastic (unfiltered, unrarified, 1174 OTUs, 46 samples)
# spp_s: OTU table soil (unfiltered, unrarified, 2194 OTUs, 49 samples)

# run the NCM for each compartment separately and adjust legends, labels and output file names
# non-linear least squares (NLS) to calculate R2:

# read in metacommunity file (OTU table containing all plastisphere or soil mycobiome OTUs, respectively)
spp <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Neutral community model/spp_s.txt", 
       header = TRUE, check.names = FALSE, row.names = 1)
spp = otu_table(spp, taxa_are_rows = TRUE)

# rarefy metacommunity to even depth
# rarefaction removes several OTUs for each compartment, defining a new metacommunity
spp = rarefy_even_depth(spp, rngseed=1, sample.size=1*min(sample_sums(spp)), replace=F) 

# transpose the metacommunity file
spp <- t(spp)

# calculate the number of individuals per community
N <- mean(apply(spp, 1, sum))

# calculate the average relative abundance of each taxa across communities
p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

# calculate the average relative abundance of each taxa across communities using a pool table
#p.m <- apply(pool, 2, mean)
#p.m <- p.m[p.m != 0]
#p <- p.m/N

# calculate the occurrence frequency of each taxa across communities
spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]

# combine
C <- merge(p, freq, by=0)
C <- C[order(C[,2]),]
C <- as.data.frame(C)
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]

# calculate the limit of detection
d = 1/N

# fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.1))
m.fit # get the m value
m.ci <- confint(m.fit, 'm', level=0.95)
m.ci

# fit neutral model parameter m (or Nm) using Maximum likelihood estimation (MLE)
sncm.LL <- function(m, sigma){
  R = freq - pbeta(d, N*m*p, N*m*(1-p), lower.tail=FALSE)
  R = dnorm(R, 0, sigma)
  -sum(log(R))
}
m.mle <- mle(sncm.LL, start=list(m=0.1, sigma=0.1), nobs=length(p))
m.mle

# calculate goodness-of-fit (R-squared and Root Mean Squared Error)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr # get the R2 value
RMSE <- sqrt(sum((freq-freq.pred)^2)/(length(freq)-1))
RMSE

# calculate Akaike's Information Criterion (AIC)
aic.fit <- AIC(m.mle, k=2)
aic.fit

# calculate the Bayesian Information Criterion (BIC)
bic.fit <- BIC(m.mle)
bic.fit

# calculate AIC and BIC for binomial model
bino.LL <- function(mu, sigma){
  R = freq - pbinom(d, N, p, lower.tail=FALSE)
  R = dnorm(R, mu, sigma)
  -sum(log(R))
}
bino.mle <- mle(bino.LL, start=list(mu=0, sigma=0.1), nobs=length(p))
bino.mle

aic.bino <- AIC(bino.mle, k=2)
aic.bino

bic.bino <- BIC(bino.mle)
bic.bino

# goodness of fit for binomial model
bino.pred <- pbinom(d, N, p, lower.tail=FALSE)
Rsqr.bino <- 1 - (sum((freq - bino.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr.bino
RMSE.bino <- sqrt(sum((freq - bino.pred)^2)/(length(freq) - 1))
RMSE.bino

bino.pred.ci <- binconf(bino.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

# calculate AIC for Poisson model
pois.LL <- function(mu, sigma){
  R = freq - ppois(d, N*p, lower.tail=FALSE)
  R = dnorm(R, mu, sigma)
  -sum(log(R))
}
pois.mle <- mle(pois.LL, start=list(mu=0, sigma=0.1), nobs=length(p))

aic.pois <- AIC(pois.mle, k=2)
bic.pois <- BIC(pois.mle)

# goodness of fit for Poisson model
pois.pred <- ppois(d, N*p, lower.tail=FALSE)
Rsqr.pois <- 1 - (sum((freq - pois.pred)^2))/(sum((freq - mean(freq))^2))
RMSE.pois <- sqrt(sum((freq - pois.pred)^2)/(length(freq) - 1))

pois.pred.ci <- binconf(pois.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)

# drawing the figure using grid package:
# p is the mean relative abundance
# freq is occurrence frequency
# freq.pred is predicted occurrence frequency
bacnlsALL <-data.frame(p,freq,freq.pred,pred.ci[,2:3])
inter.col<-rep('black',nrow(bacnlsALL))
inter.col[bacnlsALL$freq <= bacnlsALL$Lower]<-'#A52A2A'#define the color of below points
inter.col[bacnlsALL$freq >= bacnlsALL$Upper]<-'#29A6A6'#define the color of above points
library(grid)
grid.newpage()
pushViewport(viewport(h=0.6,w=0.6))
pushViewport(dataViewport(xData=range(log10(bacnlsALL$p)), yData=c(0,1.02),extension=c(0.02,0)))
grid.rect()
grid.points(log10(bacnlsALL$p), bacnlsALL$freq,pch=16,gp=gpar(col=inter.col,cex=0.4))
grid.yaxis()
grid.xaxis()
grid.lines(log10(bacnlsALL$p),bacnlsALL$freq.pred,gp=gpar(col='royalblue3',lwd=3),default='native')
grid.lines(log10(bacnlsALL$p),bacnlsALL$Lower ,gp=gpar(col='royalblue3',lwd=3,lty=2),default='native') 
grid.lines(log10(bacnlsALL$p),bacnlsALL$Upper,gp=gpar(col='royalblue3',lwd=3,lty=2),default='native')  
grid.text(y=unit(0,'npc')-unit(2.5,'lines'),label='Mean Relative Abundance (log10)', gp=gpar(fontface=2)) 
grid.text(x=unit(0,'npc')-unit(3,'lines'),label='Frequency of Occurance',gp=gpar(fontface=2),rot=90) 
#grid.text(x=unit(0,'npc')-unit(-1,'lines'), y=unit(0,'npc')-unit(-15,'lines'),label='Mean Relative Abundance (log)', gp=gpar(fontface=2)) 
#grid.text(round(coef(m.fit)*N),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2)) 
#grid.text(label = "Nm=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-15,'lines'),gp=gpar(fontface=2))
#grid.text(round(Rsqr,2),x=unit(0,'npc')-unit(-5,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
#grid.text(label = "Rsqr=",x=unit(0,'npc')-unit(-3,'lines'), y=unit(0,'npc')-unit(-16,'lines'),gp=gpar(fontface=2))
draw.text <- function(just, i, j) {
  grid.text(paste("Rsqr=",round(Rsqr,3),"\n","Nm=",round(coef(m.fit)*N)), x=x[j], y=y[i], just=just)
  #grid.text(deparse(substitute(just)), x=x[j], y=y[i] + unit(2, "lines"),
  #          gp=gpar(col="grey", fontsize=8))
}
x <- unit(1:4/5, "npc")
y <- unit(1:4/5, "npc")
draw.text(c("centre", "bottom"), 4, 1)

# compile files listing calculated values
# compile text documents to import into Excel with the separator of choice
#write.csv(p, file = "C:/Users/makis/OneDrive - Universität Bayreuth/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b/02_Data/02_R Studio/p.txt")
#write.csv(freq, file = "C:/Users/makis/OneDrive - Universität Bayreuth/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b/02_Data/02_R Studio/freq.txt")
#write.csv(freq.pred, file = "C:/Users/makis/OneDrive - Universität Bayreuth/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b/02_Data/02_R Studio/freq.pred.txt")
#write.csv(pred.ci, file = "C:/Users/makis/OneDrive - Universität Bayreuth/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b/02_Data/02_R Studio/pred.ci.txt")
write.csv(bacnlsALL, file = "C:/Users/makis/OneDrive - Universität Bayreuth/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b/02_Data/02_R Studio/bac_soil.txt")

# NCM metacommunities contain 927 (plastic) and 1991 OTUs (soil), respectively


#**********************************************************Multiplot 5: Ecological variance between partitions******************************************************#

# fitting proportions in bar charts according to compartments (plastic & soil)
# values were calculated in the respective Excel document (UBT-MYC_Siaya-1b_Calculations)
library(ggplot2)

# abundance proportions
df <- data.frame(samp=c("p_above", "p_neutral", "p_below", "s_above", "s_neutral", "s_below"),
                 prop=c(0.133, 0.597, 0.270, 0.055, 0.929, 0.016), 
                 part=c("above", "neutral", "below", "above", "neutral", "below"), 
                 sub=c("plastic", "plastic", "plastic", "soil", "soil", "soil"))
head(df)

ggplot(data=df, aes(x=part, y=prop, fill=part))+
  geom_bar(stat="identity", width=0.75)+
  theme_bw() + scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  facet_wrap(~sub, scales= "free_x", nrow=1) + 
  ggtitle("Abundance") + ylab("Proportion [%]") + ylim(0, 1)

# richness proportions
df2 <- data.frame(samp=c("p_above", "p_neutral", "p_below", "s_above", "s_neutral", "s_below"),
                  prop=c(0.110, 0.854, 0.036, 0.095, 0.876, 0.029), 
                  part=c("above", "neutral", "below", "above", "neutral", "below"), 
                  sub=c("plastic", "plastic", "plastic", "soil", "soil", "soil"))
head(df2)

ggplot(data=df2, aes(x=part, y=prop, fill=part))+
  geom_bar(stat="identity", width=0.75)+
  theme_bw() + scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  facet_wrap(~sub, scales= "free_x", nrow=1) + 
  ggtitle("Richness") + ylab("Proportion [%]") + ylim(0, 1)


# fitting proportions for the most abundant classes on plastic (excluding "unidentified" and "Incertae sedis"; threshold: 0.1%)
# plastic
df3 <- data.frame(class=c(rep("Agaricomycetes", 3), rep("Cystobasidiomycetes", 3), rep("Dothideomycetes",3), 
                          rep("Eurotiomycetes", 3), rep("Leotiomycetes", 3), rep("Malasseziomycetes", 3), 
                          rep("Microbotryomycetes", 3), rep("Mortierellomycetes", 3), rep("Pezizomycetes", 3), 
                          rep("Rhizophlyctidomycetes", 3), rep("Sordariomycetes",3), rep("Spizellomycetes", 3), 
                          rep("Tremellomycetes", 3)), 
                  part=c("above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow"), 
                  val=c(0.015, 0.208, 0.777, 0.172, 0.828, 0, 0.139, 0.525, 0.335, 0.291, 0.675, 0.034, 
                        0.456, 0.544, 0, 0, 0.318, 0.682, 0.062, 0.925, 0.013, 0, 0.534, 0.466, 0, 0.532, 0.468, 
                        0, 1, 0, 0.161, 0.732, 0.106, 0, 1, 0, 0.170, 0.830, 0)) 

ggplot(df3, aes(fill=part, y=val, x=class)) + theme_bw() +
  geom_bar(position="stack", stat="identity", width = 0.9) + 
  scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  ggtitle("Abundance plastic") + ylab("Proportion [%]")


# soil
df4 <- data.frame(class=c(rep("Agaricomycetes", 3), rep("Cystobasidiomycetes", 3), rep("Dothideomycetes",3), 
                          rep("Eurotiomycetes", 3), rep("Leotiomycetes", 3), rep("Malasseziomycetes", 3), 
                          rep("Microbotryomycetes", 3), rep("Mortierellomycetes", 3), rep("Pezizomycetes", 3),
                          rep("Rhizophlyctidomycetes", 3), rep("Sordariomycetes",3), rep("Spizellomycetes", 3), 
                          rep("Tremellomycetes", 3)), 
                  part=c("above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow",
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow"), 
                  val=c(0.095, 0.557, 0.347, 0.417, 0.583, 0, 0.017, 0.973, 0.010, 0.110, 0.871, 0.019, 
                        0.384, 0.616, 0, 0, 1, 0, 0.002, 0.997, 0.001, 0.022, 0.978, 0, 0.328, 0.672, 0, 
                        0.467, 0.533, 0, 0.059, 0.933, 0.007, 0.020, 0.980, 0, 0.002, 0.998, 0 )) 

ggplot(df4, aes(fill=part, y=val, x=class)) + theme_bw() +
  geom_bar(position="stack", stat="identity", width = 0.9) + 
  scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  ggtitle("Abundance soil") + ylab("Proportion [%]")


# richness proportions of afore-mentioned most abundant classes
# plastic 
df5 <- data.frame(class=c(rep("Agaricomycetes", 3), rep("Cystobasidiomycetes", 3), rep("Dothideomycetes",3), 
                          rep("Eurotiomycetes", 3), rep("Leotiomycetes", 3), rep("Malasseziomycetes", 3), 
                          rep("Microbotryomycetes", 3), rep("Mortierellomycetes", 3), rep("Pezizomycetes", 3), 
                          rep("Rhizophlyctidomycetes", 3), rep("Sordariomycetes",3), rep("Spizellomycetes", 3), 
                          rep("Tremellomycetes", 3)), 
                  part=c("above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow"), 
                  val=c(0.012, 0.904, 0.084, 0.2, 0.8, 0, 0.245, 0.726, 0.029, 0.068, 0.915, 0.017, 0.143, 0.857, 0, 
                        0, 0.5, 0.5, 0.091, 0.818, 0.091, 0, 0.947, 0.053, 0, 0.923, 0.077, 0, 1, 0, 0.094, 0.865, 0.041, 
                        0, 1, 0, 0.273, 0.727, 0)) 

ggplot(df5, aes(fill=part, y=val, x=class)) + theme_bw() +
  geom_bar(position="stack", stat="identity", width = 0.9) + 
  scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  ggtitle("Richness plastic") + ylab("Proportion [%]")


# Soil
df6 <- data.frame(class=c(rep("Agaricomycetes", 3), rep("Cystobasidiomycetes", 3), rep("Dothideomycetes",3), 
                          rep("Eurotiomycetes", 3), rep("Leotiomycetes", 3), rep("Malasseziomycetes", 3), 
                          rep("Microbotryomycetes", 3), rep("Mortierellomycetes", 3), rep("Pezizomycetes", 3), 
                          rep("Rhizophlyctidomycetes", 3), rep("Sordariomycetes",3), rep("Spizellomycetes", 3), 
                          rep("Tremellomycetes", 3)), 
                  part=c("above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow", "above", "neutral", "xelow", "above", "neutral", "xelow", 
                         "above", "neutral", "xelow"), 
                  val=c(0.062, 0.907, 0.031, 0.1, 0.9, 0, 0.096, 0.874, 0.030, 0.098, 0.865, 0.038, 0.077, 0.923, 0, 
                        0, 1, 0, 0.118, 0.824, 0.059, 0.108, 0.892, 0, 0.095, 0.905, 0, 0.2, 0.8, 0, 0.112, 0.866, 0.022, 
                        0.083, 0.917, 0, 0.036, 0.945, 0.018)) 

ggplot(df6, aes(fill=part, y=val, x=class)) + theme_bw() +
  geom_bar(position="stack", stat="identity", width = 0.9) + 
  scale_fill_manual(values=c("darkcyan", "black", "darkred")) + 
  ggtitle("Richness soil") + ylab("Proportion [%]")


#***********************************************************************Alpha Diversity************************************************************************#
# alpha diversity between partitions and substrates
# OTU_super-partition: OTU table with partitions for each sample and compartment
# fac_partition: factors adujusted to partitions for each sample and compartment
# tax: taxonomy table of the entire (unfiltered) dataset

library(phyloseq)
library(ggplot2)
library(microbiome)
install.packages("ggpubr")
library(ggpubr)

OTU <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/OTU_super-partition.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
Factors <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/fac_super-partition.txt", 
                      header = TRUE, check.names = FALSE, row.names = 1)
Tax <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Original data/tax.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)
tax_matrix = as.matrix(Tax)

# Phyloseq specifications
# transforming the three tables to a format that is convertible into Phyloseq
OTU = otu_table(OTU, taxa_are_rows = TRUE)
TAX = tax_table(tax_matrix)
FAC = sample_data(Factors)

# combining the three tables into a single Phyloseq Document
S1b = phyloseq(OTU, TAX, FAC)
S1b

# inspecting the phyloseq document
rank_names(S1b)
sample_variables(S1b)
S1b

# summarizing the contents of a phyloseq object
summarize_phyloseq(S1b)

# reorganizing the variable order for the entire dataset
Partition <- factor(sample_data(S1b)$Partition, levels = c("above", "neutral", "xelow"))
levels(Partition)
Substrate <- factor(sample_data(S1b)$Substrate, levels = c("all", "soil", "plastic"))
levels(Substrate)

# transformation to relative abundance
rel = transform_sample_counts(S1b, function(x) x / sum(x) ) #creates a transformed phyloseq document with relative abundance data
rel
get_sample(rel)

# alpha diversity calculations
a <- estimate_richness(S1b, measure = c("Chao1", "Observed", "Shannon"))
a

H <- a$Shannon # makes Shannon diversity a vector
S1 <- a$Observed # makes species richness a vector
S <- log(S1) # vectoring the logarithm of species richness
J <- H/S # calculating and vectoring Pielou's species evenness

a$Evenness = J # defining vector J as evenness within the alpha_diversity dataframe
a # returning all contained diversity indices (including evenness)

# combining alpha diversity boxplots
df <- data.frame(a, sample_data(S1b))
df2 <- tidyr::gather(df, key = "Measure", 1, value = "Value", Shannon, Evenness)

# dot plots
#ggplot(data = df2, aes(x = Partition, y = Value, color = Partition, shape = Substrate)) + 
#facet_wrap(~Measure, scale = "free") + geom_point(size = 3) + 
#theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
#scale_color_manual(values = c("darkcyan", "darkred", "black")) + theme_bw() +
#geom_boxplot()

# boxplots
p <- ggplot(data = df2, aes(x = Substrate, y = Value, fill = Partition)) +
  facet_wrap(~Measure, scale = "free") +
  geom_boxplot() +
  scale_fill_manual(values = c("darkcyan", "black", "darkred")) + theme_bw()
p

write.csv(df, file = "C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/08_Results/03_R Studio/alpha-diversity.txt")

# no statistical test amenable due to single-sample conditions (implicit mean alpha diversity metrics)



#*******************************************************************Multiplot-6: Co-occurrence networks*******************************************************************#

### co-occurrence networks with igraph and psych

# load the packages required
library(igraph) 
library(psych)

## plastisphere 
# load the data
# p10 = plastisphere OTUs with >10% occurrence frequency
p10 <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Networks/plastisphere_full.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)

# calculate spearman correlations into an adjacency matrix from the OTU table
occur_p10 = corr.test(t(p10),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occur_p10_r = occur_p10$r # make the matrix of correlation a vector
occur_p10_p = occur_p10$p # make p values of the correlations a vector

# remove correlation from the dataset below a set threshold
occur_p10_r[occur_p10_p>0.01|abs(occur_p10_r)<0.6] = 0 #remove those r<0.6, p>0.01.

# create igraph graph from adjacency matrix
igraph.p10 = graph_from_adjacency_matrix(occur_p10_r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.p10.weight = E(igraph.p10)$weight
E(igraph.p10)$weight = NA

# explore the igraph data
# igraph summary
igraph.p10
gsize(igraph.p10)
gorder(igraph.p10)

# node list
V(igraph.p10)

# edge list
E(igraph.p10)

# convert graph to an edge list
edge_p10 <- as_edgelist(igraph.p10, names = TRUE)

# export adjacency matrices including r and p values
#write.csv(occur_p10[["r"]], file = ("C:/Users/geras/Desktop/p10_r.txt"))
#write.csv(occur_p10[["p"]], file = ("C:/Users/geras/Desktop/p10_p.txt"))

# export 
write.csv(edge_p10, file = ("C:/Users/geras/Desktop/edge_p_all_sel.txt"))

# set edge color postive correlation pink color, negative blue.
E.color.p10 = igraph.p10.weight
E.color.p10 = ifelse(E.color.p10>0, "pink",ifelse(E.color.p10<0, "blue","grey")) 
E(igraph.p10)$color = as.character(E.color.p10)

#change edge width
E(igraph.p10)$width = abs(igraph.p10.weight)

# set vertices color, modularity
fc.p10 = cluster_fast_greedy(igraph.p10,weights =NULL)
modularity.p10 = modularity(igraph.p10,membership(fc.p10))
comps.p10 = membership(fc.p10)
colbar.p10 = rainbow(max(comps.p10))
V(igraph.p10)$color = colbar.p10[comps.p10]

# plot the network
set.seed(123)
plot(igraph.p10,vertex.frame.color=NA,vertex.label=NA,edge.width=0.2,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


### soil 

# s10 = soil OTUs wit >10% occurrence frequency
s10 <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Networks/soil_0.1.txt", 
                  header = TRUE, check.names = FALSE, row.names = 1)

# calculate spearman correlations into an adjacency matrix from the OTU table
occur_s10 = corr.test(t(s10),use="pairwise",method="spearman",adjust="fdr",alpha=.05)
occur_s10_r = occur_s10$r # make the matrix of correlation a vector
occur_s10_p = occur_s10$p # make p values of the correlations a vector

# remove correlation from the dataset below a set threshold
occur_s10_r[occur_s10_p>0.01|abs(occur_s10_r)<0.6] = 0 #remove those r<0.6, p>0.01.

# create igraph graph from adjacency matrix
igraph.s10 = graph_from_adjacency_matrix(occur_s10_r,mode="undirected",weighted=TRUE,diag=FALSE)
igraph.s10.weight = E(igraph.s10)$weight
E(igraph.s10)$weight = NA

# convert graph to an edge list
edge_s10 <- as_edgelist(igraph.s10, names = TRUE)

# export adjacency matrices including r and p values
write.csv(occur_s10[["r"]], file = ("C:/Users/geras/Desktop/s10_r.txt"))
#write.csv(occur_s10[["p"]], file = ("C:/Users/geras/Desktop/s10_p.txt"))

# export 
write.csv(edge_s10, file = ("C:/Users/geras/Desktop/edge_s10.txt"))

# set edge color postive correlation pink color, negative blue.
E.color.s10 = igraph.s10.weight
E.color.s10 = ifelse(E.color.s10>0, "pink",ifelse(E.color.s10<0, "blue","grey")) 
E(igraph.s10)$color = as.character(E.color.s10)

#change edge width
E(igraph.s10)$width = abs(igraph.s10.weight)

# set vertices color, modularity
fc.s10 = cluster_fast_greedy(igraph.s10,weights =NULL)
modularity.s10 = modularity(igraph.s10,membership(fc.s10))
comps.s10 = membership(fc.s10)
colbar.s10 = rainbow(max(comps.s10))
V(igraph.s10)$color = colbar.s10[comps.s10]

# plot the network
set.seed(123)
plot(igraph.s10,vertex.frame.color=NA, edge.width=0.2, vertex.label=NA,
     vertex.size=5, edge.lty=1,edge.curved=TRUE,margin=c(0,0,0,0))


#************************************************************************Multiplot 7: Correlations************************************************************************#

# VIR score boxplots (virulence scoring according to PLS and compartment)
# including mutualists, pathogens, saprotrophs as well as opportunists and non-opportunists

#cor_p: plastisphere MC GVT data
#cor_s: soil MC GVT data

# library calls
library(ggplot2)
library(ggpubr)
library(psych)

# import data
cor_p <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/cor_p.txt",
                          header = TRUE, check.names = FALSE, row.names = 1)
cor_s <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/cor_s.txt",
                    header = TRUE, check.names = FALSE, row.names = 1)


# VIR boxplots according to PLS
VIR_p_box_PLS <- ggboxplot(cor_p, x = "Primary_lifestyle", y = "VIR", fill = "Primary_lifestyle", scales = "free") + 
  scale_fill_manual(values = c("darkolivegreen1", "firebrick2", "gold"))
VIR_p_box_PLS

VIR_s_box_PLS <- ggboxplot(cor_s, x = "Primary_lifestyle", y = "VIR", fill = "Primary_lifestyle", scales = "free") + 
  scale_fill_manual(values = c("darkolivegreen1", "firebrick2", "gold"))
VIR_s_box_PLS

# ANOVA
aov_p_VIR_PLS <- aov(cor_p$Primary_lifestyle,cor_p$VIR)

# VIR boxplots according to opportunism
VIR_p_box_opp <- ggboxplot(cor_p, x = "Opportunism", y = "VIR", fill = "Opportunism", scales = "free") + 
  scale_fill_manual(values = c( "azure", "brown4"))
VIR_p_box_opp

VIR_s_box_opp <- ggboxplot(cor_s, x = "Opportunism", y = "VIR", fill = "Opportunism", scales = "free") + 
  scale_fill_manual(values = c( "azure", "brown4"))
VIR_s_box_opp


# import data
df <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df.txt",
      header = TRUE, check.names = FALSE, row.names = 1)
df_p <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_p.txt",
        header = TRUE, check.names = FALSE, row.names = 1)
df_s <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_s.txt",
        header = TRUE, check.names = FALSE, row.names = 1)
df_path <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_path.txt",
           header = TRUE, check.names = FALSE, row.names = 1)
df_p_path <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_p_path.txt",
             header = TRUE, check.names = FALSE, row.names = 1)
df_s_path <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_s_path.txt",
             header = TRUE, check.names = FALSE, row.names = 1)
df_CAO <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_CAO.txt",
          header = TRUE, check.names = FALSE, row.names = 1)
df_p_CAO <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_p_CAO.txt",
            header = TRUE, check.names = FALSE, row.names = 1)
df_s_CAO <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_s_CAO.txt",
            header = TRUE, check.names = FALSE, row.names = 1)
df_p_class <- read.delim("C:/Users/geras/Desktop/PhD/01_Projects/SFB MP/02_Siaya-1b (S1b)/02_Data/02_R Studio/Correlations/df_p_class.txt",
                       header = TRUE, check.names = FALSE, row.names = 1)

# render the boxplots (width: 550, height: 450)
# number of selected CAO (cosmopolitan adaptable opportunist) species between MP and soil
CAO_S_BXP <- ggboxplot(df_CAO, x = "sub", y = "VIR", fill = "sub", scales = "free") + 
  scale_fill_manual(values = c("cornflowerblue", "darkred"))
CAO_S_BXP

# TAS of CAO species between MP and soil
CAO_SI_BXP <- ggboxplot(df_CAO, x = "sub", y = "TAS", fill = "sub", scales = "free") + 
  scale_fill_manual(values = c("cornflowerblue", "darkred"))
CAO_SI_BXP

# SI of CAO species between MP and soil
CAO_GR_BP <- ggboxplot(df_CAO, x = "sub", y = "SI", fill = "sub", scales = "free") + 
  scale_fill_manual(values = c("cornflowerblue", "darkred")) 
CAO_GR_BP


# Correlation (scatter) plots

# render the correlation scatter plots (width: 550, height: 450)
# VIR x WGS (adjust data frame)
ggplot(df, aes(x=vir_score, y=GS, color=Sub, shape=Sub)) +
  geom_point() + geom_smooth(method=lm, linetype="solid", aes(fill=Sub)) + theme_classic() + 
  scale_fill_manual(values = c("cornflowerblue", "darkred"))+ # fill color of the 95% confidence intervals
  scale_color_manual(values = c("cornflowerblue", "darkred")) # regression line color

# WGS x AMS (adjust data frame)
ggplot(df_CAO, aes(x=AMS, y=WGS, color=sub, shape=sub)) +
  geom_point(size=.1) + geom_smooth(method=lm, linetype="solid", aes(fill=sub)) + theme_classic() + 
  scale_fill_manual(values = c("cornflowerblue", "darkred"))+ # fill color of the 95% confidence intervals
  scale_color_manual(values = c("cornflowerblue", "darkred")) # regression line color

# VIR x TAS (adjust data frame)
ggplot(df_path, aes(x=VIR, y=TAS, color=sub, shape=sub)) +
  geom_point(size=.1) + geom_smooth(method=lm, linetype="solid", aes(fill=sub)) + theme_classic() + 
  scale_fill_manual(values = c("cornflowerblue", "darkred"))+ # fill color of the 95% confidence intervals
  scale_color_manual(values = c("cornflowerblue", "darkred")) # regression line color


# correlation testing
# VIR x WGS for path groups
# plastic
cor(df_p_path$vir_score, df_p_path$GS, method = "spearman")
cor.test(df_p_path$vir_score, df_p_path$GS, method = "spearman")

# soil
cor(df_s_path$vir_score, df_s_path$GS, method = "spearman")
cor.test(df_s_path$vir_score, df_s_path$GS, method = "spearman")


# AMS x WGS for path groups
# plastic
cor(df_p_path$selected_from_total, df_p_path$GS, method = "spearman")
cor.test(df_p_path$selected_from_total, df_p_path$GS, method = "spearman")

# soil
cor(df_s_path$selected_from_total, df_s_path$GS, method = "spearman")
cor.test(df_s_path$selected_from_total, df_s_path$GS, method = "spearman")


# VIR x WGS for all selected groups
# plastic
cor(df_p$vir_score, df_p$GS, method = "spearman")
cor.test(df_p$vir_score, df_p$GS, method = "spearman")

# soil
cor(df_s$vir_score, df_s$GS, method = "spearman")
cor.test(df_s$vir_score, df_s$GS, method = "spearman")


# GS x GR for path groups
# plastic
cor(df_p_CAO$GR, df_p_CAO$GS, method = "spearman")
cor.test(df_p_CAO$GR, df_p_CAO$GS, method = "spearman")

# soil
cor(df_s_CAO$GR, df_s_CAO$GS, method = "spearman")
cor.test(df_s_CAO$GR, df_s_CAO$GS, method = "spearman")


# SFT x GS for CAO groups
# plastic
cor(df_p_CAO$GS, df_p_CAO$selected_from_total, method = "spearman")
cor.test(df_p_CAO$GS, df_p_CAO$selected_from_total, method = "spearman")

# soil
cor(df_s_CAO$GS, df_s_CAO$selected_from_total, method = "spearman")
cor.test(df_s_CAO$GS, df_s_CAO$selected_from_total, method = "spearman")


# SFT x GS for path groups
# plastic
cor(df_p_CAO$GS, df_p_CAO$selected_from_total, method = "spearman")
cor.test(df_p_CAO$GS, df_p_CAO$selected_from_total, method = "spearman")

# soil
cor(df_s_CAO$GS, df_s_CAO$selected_from_total, method = "spearman")
cor.test(df_s_CAO$GS, df_s_CAO$selected_from_total, method = "spearman")


# Bubble plots (width: 550, height: 450)
# CAO genera / species
ggplot(df_CAO, aes(x=AMS, y=WGS, size=SI, color = sub)) + scale_color_manual(values = c("cornflowerblue", "darkred")) + 
  geom_point(alpha=0.9) + theme_classic() + geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(.1, 8), name="Selection")

# path groups
ggplot(df_path, aes(x=selected_from_total, y=GS, size=GR, color=Sub)) + scale_color_manual(values = c("cornflowerblue", "darkred")) + 
  geom_point(alpha=0.9) + theme_classic() + geom_point(alpha=0.5, shape=21, color="black") +
  scale_size(range = c(.1, 10), name="Selection")


# Correlation Testing
# testing correlations between relevant ranks / descriptors / traits
cor_p <- read.delim("C:/Users/geras/Desktop/cor_p.txt", header = TRUE, check.names = FALSE, row.names = 1)
cor_s <- read.delim("C:/Users/geras/Desktop/cor_s.txt", header = TRUE, check.names = FALSE, row.names = 1)

# PLS x opp
cor.test(cor_p$PLS, cor_p$opp, method = "spearman")
cor.test(cor_s$PLS, cor_s$opp, method = "spearman")

# PLS x ex
cor.test(cor_p$PLS, cor_p$ex, method = "spearman")
cor.test(cor_s$PLS, cor_s$ex, method = "spearman")

# PLS x oli
cor.test(cor_p$PLS, cor_p$oli, method = "spearman")
cor.test(cor_s$PLS, cor_s$oli, method = "spearman")

# PLS x pig
cor.test(cor_p$PLS, cor_p$pig, method = "spearman")
cor.test(cor_s$PLS, cor_s$pig, method = "spearman")

# PLS x pol
cor.test(cor_p$PLS, cor_p$pol, method = "spearman")
cor.test(cor_s$PLS, cor_s$pol, method = "spearman")

# PLS x VIR
cor.test(cor_p$PLS, cor_p$VIR, method = "spearman")
cor.test(cor_s$PLS, cor_s$VIR, method = "spearman")

# ex x opp
cor.test(cor_p$ex, cor_p$opp, method = "spearman")
cor.test(cor_s$ex, cor_s$opp, method = "spearman")

# ex x oli
cor.test(cor_p$ex, cor_p$oli, method = "spearman")
cor.test(cor_s$ex, cor_s$oli, method = "spearman")

# ex x pig
cor.test(cor_p$ex, cor_p$pig, method = "spearman")
cor.test(cor_s$ex, cor_s$pig, method = "spearman")

# ex x pol
cor.test(cor_p$ex, cor_p$pol, method = "spearman")
cor.test(cor_s$ex, cor_s$pol, method = "spearman")

# ex x VIR
cor.test(cor_p$ex, cor_p$VIR, method = "spearman")
cor.test(cor_s$ex, cor_s$VIR, method = "spearman")

# opp x oli
cor.test(cor_p$opp, cor_p$oli, method = "spearman")
cor.test(cor_s$opp, cor_s$oli, method = "spearman")

# opp x pig
cor.test(cor_p$opp, cor_p$pig, method = "spearman")
cor.test(cor_s$opp, cor_s$pig, method = "spearman")

# opp x pol
cor.test(cor_p$opp, cor_p$pol, method = "spearman")
cor.test(cor_s$opp, cor_s$pol, method = "spearman")

# opp x VIR
cor.test(cor_p$opp, cor_p$VIR, method = "spearman")
cor.test(cor_s$opp, cor_s$VIR, method = "spearman")

# oli x pig
cor.test(cor_p$oli, cor_p$pig, method = "spearman")
cor.test(cor_s$oli, cor_s$pig, method = "spearman")

# oli x pol
cor.test(cor_p$oli, cor_p$pol, method = "spearman")
cor.test(cor_s$oli, cor_s$pol, method = "spearman")

# oli x VIR
cor.test(cor_p$oli, cor_p$VIR, method = "spearman")
cor.test(cor_s$oli, cor_s$VIR, method = "spearman")

# pig x pol
cor.test(cor_p$pig, cor_p$pol, method = "spearman")
cor.test(cor_s$pig, cor_s$pol, method = "spearman")

# pig x VIR
cor.test(cor_p$pig, cor_p$VIR, method = "spearman")
cor.test(cor_s$pig, cor_s$VIR, method = "spearman")

# pol x VIR
cor.test(cor_p$pol, cor_p$VIR, method = "spearman")
cor.test(cor_s$pol, cor_s$VIR, method = "spearman")

#************************************************************************End of Analyses************************************************************************#


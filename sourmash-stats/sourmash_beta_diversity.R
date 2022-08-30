### Load libraries
library(tidyverse)
library(phyloseq)
library(qiime2R)
library(here)
library(vegan)
library(ggpubr)
library(egg)
library(rstatix)
library(grid)

### Load data and make phyloseq object

#### ASV table, taxonomy table, and metadata table
asv_table_path <- here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "sourmash_table.csv")
tax_table_path <- here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "sourmash_gtdb_mmetsp_taxonomy-domain.csv")
metadata_path <- here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "sourmash_metadata_phyloseq.csv")
sourmash_phyloseq <- merge_phyloseq(otu_table(read.csv(asv_table_path, row.names = 1), taxa_are_rows = TRUE),
                                       tax_table(as.matrix(read.csv(tax_table_path, row.names = 1))),
                                       sample_data(read.csv(metadata_path, row.names = 1)))

### Run PERMANOVA

# Load metadata
pca_metadata <- read.csv(here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats","permanova_metadata.csv"))

## Load DEICODE distance matrix

distance <- read_qza(here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "deicode", "Lasker2019_Illumina_PacBio_sourmash_GTDB_MMETSP_deicode-distance.qza"))
# Extract distance matrix
distance_matrix <- distance$data
# convert DEICODE matrix to "dist" class object
PCA_dist <- as.dist(distance_matrix)

##----Run individual PERMANOVAs-------------------------------------------------------------
# PERMANOVA for PacBio vs. Illumina
permanova <- adonis2(PCA_dist ~ sequencing_platform, data = pca_metadata, permutations=999)

# PERMANOVA for depth groups
permanova <- adonis2(PCA_dist ~ DepthGroup, data = pca_metadata, permutations=999)

### Plot 16S DEICODE RPCA - pairs labeled

pco <- read_qza(here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "deicode", "Lasker2019_Illumina_PacBio_sourmash_GTDB_MMETSP_deicode-ordination.qza"))
label.PC1 <- paste0("PC1 (", round(pco$data$ProportionExplained$PC1, 3)*100,"%)")
label.PC1
label.PC2 <- paste0("PC2 (", round(pco$data$ProportionExplained$PC2, 3)*100,"%)")
label.PC2
label.PC3 <- paste0("PC3 (", round(pco$data$ProportionExplained$PC3, 3)*100,"%)")
label.PC3
## Prepare PCA data for ggplot
pca_metadata %>% dplyr::rename(., "SampleID" = "Sample") -> pca_metadata
pca_data <- pco$data$Vectors
pca_data <- right_join(pca_data,pca_metadata,on = "SampleID")
pca_data <- subset(pca_data, !(is.na(pca_data$PC1)))
# Extract depth data
pca_data$depth <- parse_number(pca_data$DepthGroup)

# Shape guide for depths
depth_shapes <- c("5-10m" = 21, "11-28m" = 22, "45-50m" = 24)
sequencing_colors <- c("Illumina" = "#ff7f00", "PacBio" = "#1f78b4")

##----Make plot-----------------------------------------------------------------
PCA_pairs_labeled <- ggplot(pca_data,aes(x=PC1,y=PC2,color=sequencing_platform, 
                                         fill=sequencing_platform, shape=DepthGroup)) + 
  geom_point(size = 5, stroke = 1) +
  geom_line(lty = 2, aes(group = eDNA_sample)) +
  xlab(print(label.PC1))+
  ylab(print(label.PC2))+
  scale_color_manual(values = c("black", "black"))+
  scale_fill_manual(values = sequencing_colors)+
  scale_shape_manual(values = depth_shapes)+
  theme(panel.background = element_rect(fill = "white",size = 1),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.title = element_text(size = 20),
        axis.text = element_blank(),
        axis.ticks.length=unit(0.25, "cm"),
        axis.ticks=element_blank(),
        plot.margin = margin(0.25, 0.25, 0.25, 0.25, "cm"),
        axis.line.x.bottom = element_line(color = "black", size = 1),
        axis.line.y.left = element_line(color = "black", size = 1),
        legend.position = "none")+
  ylim(min(pca_data$PC2)-abs(max(pca_data$PC2)-min(pca_data$PC2))/30, max(pca_data$PC2)+abs(max(pca_data$PC2)-min(pca_data$PC2))/12)

PCA_pairs_labeled

ggsave(path = "/Users/nastassia.patin/GitHub/PacBio-Marine-Metagenomes/sourmash-stats/",
       filename = "PacBio_Illumina_sourmash_GTDB-MMETSP_PCA_pairs_labeled.png", PCA_pairs_labeled,width = 7,height = 5)

# Create a legend manually
PCA_legend <- ggplot(pca_data) +
  # PacBio vs. Illumina
  annotate("text",label = "Sequencing  Platform", x = 0, y = 2.8,size = 8,adj = 0)+ 
  annotate("point", x = 0, y = 2.6, shape = 21, colour = "black", fill = "#ff7f00", size = 7, stroke = 2)+ 
  annotate("text",label = "PacBio", x = 0.1, y = 2.6,size = 7,adj = 5.45)+
  annotate("point", x = 0, y = 2.7, shape = 21, colour = "black", fill = "#1f78b4", size = 7, stroke = 2)+
  annotate("text",label = "Illumina", x = 0.1, y = 2.7,size = 7,adj = 5)+
  # Depth
  annotate("text",label = "Depth Group (m)", x = 0, y = 2.4,size = 8,adj = 0)+ # Title 
  annotate("point", x = 0, y = 2.3, shape = 21, colour = "black", fill = "black", size = 7, stroke = 3)+ # circle point
  annotate("text",label = "5-10 m", x = 0.1, y = 2.3,size = 7,adj = 5.4)+
  annotate("point", x = 0, y = 2.2, shape = 22, colour = "black", fill = "black", size = 7, stroke = 3)+ # triangle point
  annotate("text",label = "11-28 m", x = 0.1, y = 2.2,size = 7,adj = 4.6)+
  annotate("point", x = 0, y = 2.1, shape = 24, colour = "black", fill = "black", size = 7, stroke = 3)+ # square point
  annotate("text",label = "45-50 m", x = 0.1, y = 2.1,size = 7,adj = 4.6)+
  theme(panel.background = element_rect(fill="white"), 
        panel.border = element_rect(colour = "black", fill=NA, size=2),
        legend.position = "none",
        plot.margin = unit(c(0,10,0,10),"cm"))

# xlim(0,10)+
#  ylim(0.2,3.2)+

PCA_legend

combined_PCA_gg <- ggpubr::ggarrange(ggpubr::ggarrange(PCA_pairs_labeled, ncol = 1, nrow = 1,
                                                       labels = c("sourmash PCA"), label.x = 0.05, label.y = 0.985,
                                                       font.label = list(size = 35, color = "black", face = "bold"), hjust = -0.12),
                                     PCA_legend, nrow = 1, ncol = 1, heights = c(7,1))

combined_PCA_gg

#=============================
# Pairwise comparisons
#=============================
## Write functions that will be used for each marker

### Data prep/reformatting function
# This function returns a dataframe that contains distances for each pair of samples as well as the relationship between the two samples in terms of their cruise, depth, and whether or not they were paired samples

distance_data_prep <- function(metadata, distance) {
  # Extract metadata
  pca_metadata <- metadata
  pca_metadata %>% dplyr::rename(., "SampleID" = "Sample") -> pca_metadata
  
  # Extract distance matrix
  distance_matrix <- distance$data
  
  # Convert distance matrix into a df of each pair
  distance_matrix %>%
    as.matrix() %>%
    as.data.frame() %>%
    rownames_to_column("Sample1") %>%
    pivot_longer(., cols = colnames(.[2:ncol(.)]), names_to = "Sample2") %>%
    dplyr::rename(distance = value) -> distances_long
  
  ggplot(distances_long, aes(x = distance)) +
    geom_histogram()
  
  # Reformat metadata
  pca_metadata %>%
    dplyr::rename(Sample1 = SampleID) %>%
    dplyr::select(Sample1, sequencing_platform, DepthGroup, eDNA_sample) %>%
    dplyr::rename(
      sequencing_platform_1 = sequencing_platform,
      DepthGroup_1 = DepthGroup,
      eDNA_sample_1 = eDNA_sample
    ) -> metadata_sample1
  
  pca_metadata %>%
    dplyr::rename(Sample2 = SampleID) %>%
    dplyr::select(Sample2, sequencing_platform, DepthGroup, eDNA_sample) %>%
    dplyr::rename(
      sequencing_platform_2 = sequencing_platform,
      DepthGroup_2 = DepthGroup,
      eDNA_sample_2 = eDNA_sample
    ) -> metadata_sample2
  
  
  # Add metadata
  distances_long %>%
    left_join(., metadata_sample1, by = "Sample1") %>%
    left_join(., metadata_sample2, by = "Sample2") %>%
    # Create new fields to show within vs. among group differences
    mutate(
      same_sequencing = ifelse(sequencing_platform_1 == sequencing_platform_2, "Same sequencing", "Different sequencing"),
      paired = ifelse(eDNA_sample_1 == eDNA_sample_2, "Paired", "Not paired"),
      same_depth = ifelse(DepthGroup_1 == DepthGroup_2, "Same depth", "Different depth"),
      # Paired and sequencing combined
      same_sequencing_paired = ifelse(
        same_sequencing == "Same sequencing" &
          paired == "Paired",
        "Same sequencing, paired",
        ifelse(
          same_sequencing == "Same sequencing" &
            paired == "Not paired",
          "Same sequencing, not paired",
          ifelse(
            same_sequencing == "Different sequencing" &
              paired == "Paired",
            "Different sequencing, paired",
            ifelse(
              same_sequencing == "Different sequencing" &
                paired == "Not paired",
              "Different sequencing, not paired",
              NA
            )
          )
        )
      )
    ) %>%
    # Remove all zero distances (comparing sample to itself) %>%
    subset(distance != 0) %>%
    # Remove all duplicate distances (i.e., sample 1 and sample 2 vs. sample 2 and sample 1)
    filter(duplicated(distance) == FALSE) -> distances_metadata
  
  return(distances_metadata)
  
}

### Plotting function
# This function plots histograms of pairwise distances by sequencing platform, depth group, and paired vs. not paired.

# Create a function to plot histogram based on input data and variable of interest
distance_histogram <- function(input_data, variable){
  # Store colors
  if (variable == "same_sequencing_paired"){
    color_values = c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
  }
  else {
    color_values = c("gray80", "gray20")
  }
  
  # Store legend title
  if (variable == "same_sequencing_paired"){
    legend_title = "Sequencing & Paired"
  }
  else if (variable == "paired"){
    legend_title = "Paired"
  }
  else if (variable == "same_sequencing") {
    legend_title = "Sequencing"
  }
  
  
  histogram <- ggplot(data = input_data, aes_string(x = "distance", fill = variable))+
    geom_histogram(bins = 100) +
    scale_fill_manual(values = color_values, name = legend_title) +
    theme(panel.background = element_rect(fill ="white", color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0.5,0.25,0,0), "cm"),
          ) +
    scale_x_continuous(limits = c(0,max(input_data$distance+0.1)), oob = scales::oob_keep, expand = c(0,0)) +
    xlab("Aitchison Distance") +
    # scale_y_continuous(limits = c(0, 70), expand = c(0,0))
    scale_y_continuous(expand = expansion(mult = c(0, .05)))
  return(histogram)
}

### Reformat data

distance <- read_qza(here::here("GitHub", "PacBio-Marine-Metagenomes", "sourmash-stats", "deicode", "Lasker2019_Illumina_PacBio_sourmash_GTDB_MMETSP_deicode-distance.qza"))
distances_metadata <- distance_data_prep(metadata = pca_metadata, distance = distance)

### Create figures

# Make a single column figure, multiple panels. Same x-axis on each panel.
# Split the data so that you have all paired samples, samples from the same depth, samples from the same cruise.

ggplot(distances_metadata, aes(x = distance))+
  geom_histogram(bins = 100)
paired_histogram <- distance_histogram(input_data = distances_metadata, variable = "paired")
sequencing_histogram <- distance_histogram(input_data = distances_metadata, variable = "same_sequencing")
depth_histogram <- distance_histogram(input_data = distances_metadata, variable = "same_depth")
sequencing_paired_histogram <- distance_histogram(input_data = distances_metadata, variable = "same_sequencing_paired")
histograms_16S <- egg::ggarrange(paired_histogram, sequencing_histogram, depth_histogram, sequencing_depth_histogram,
                                 ncol = 2)

paired_hist_legend <- ggpubr::get_legend(paired_histogram)
paired_hist_legend <- as_ggplot(paired_hist_legend)
sequencing_paired_hist_legend <- ggpubr::get_legend(sequencing_paired_histogram)
sequencing_paired_hist_legend <- as_ggplot(sequencing_paired_hist_legend)

combined_dist_hist <- egg::ggarrange(paired_histogram, sequencing_histogram,
                                     sequencing_paired_histogram,
                                     nrow = 3, ncol = 1
                                     ) 
# labels = c("(a) paired/unpaired", "(b) sequencing platform", 
# "(c) sequencing platform and paired/unpaired"),
# label.args = list(gp=gpar(font=1, cex = 0.8), x=unit(0.25,"cm"), y=unit(3.4,"cm"), hjust = 0)

ggsave(here("resubmission", "combined_distance_histogram_v2.pdf"), combined_dist_hist,
       height = 11, width = 8)

### Run ANOVA

## Here we will be running an ANOVA on pairwise distances, rather than a PERMANOVA.

## In this ANOVA we will test three variables: 
# 1) paired vs. not paired
# 2) same sequencing vs. different sequencing
# 3) same depth vs different depth

# Order matters: control for sequencing and depth first
# Use same distance matrix as generated above for histograms
anova <- aov(distance ~ same_sequencing + paired,  data = distances_metadata)
summary(anova)


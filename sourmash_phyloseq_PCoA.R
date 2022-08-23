library(phyloseq)
library(ggplot2)
library(vegan)
library(data.table)
library(textshape)
library(magrittr) # to use %>% as a pipe
library(RColorBrewer)
library(Polychrome)
library(DESeq2)
library(ggfortify)
library(LDM)

# import ASV table (eg "table_IDs_sterivex") and ID-to-taxonomy table

microbe <- column_to_rownames(Lasker2019_Illumina_PacBio_sourmash_species_bpnumbers, 'X1')
microbe[is.na(microbe)] <- 0
tax_table <- column_to_rownames(sourmash_taxonomy_phyloseq, 'ID')
metadata <- column_to_rownames(sourmash_metadata_phyloseq, 'sample_name')

# convert all to matrices
microbe <- as.matrix(microbe)
tax_table <- as.matrix(tax_table)
metadata <- as.data.frame(metadata)

# combine into phyloseq object
ASV = otu_table(microbe, taxa_are_rows = TRUE)
TAX = tax_table(tax_table)
MET = sample_data(metadata)
physeq = phyloseq(ASV, TAX, MET)

# Calculate distances
DistBC = distance(physeq, method = "bray")

# Ordinate
ordBC = ordinate(physeq, method = "PCoA", distance = DistBC)

# Create color palette
mypal <- light.colors(14)
names(mypal) <- NULL

# Plot Bray-Curtis PCoA
plot_ordination(physeq, ordBC, color="eDNA.sample", shape="sequencing.platform") + 
  geom_point(size = 4) + theme_classic() + 
  scale_color_manual(values=mypal) # + ggtitle("PCoA: Bray-Curtis")

### PERMANOVA ###
asv = microbe[, colnames(microbe) != "lineage"]

# Import metadata
metadata <- as.data.frame(permanova_metadata)

# Run vst transformation
# add pseudocount
asv_pc <- asv + 1
asv_pc <- as.matrix(asv)
asv_vst <- varianceStabilizingTransformation(asv_pc, blind = TRUE, fitType = "mean")

# convert to data frame, transform, and convert again
asv_vst <- as.data.frame.matrix(asv_vst)
asv_t <- t(asv_vst)
asv_t <- as.data.frame.matrix(asv_t)

# PERMANOVAs for sequencing type and depth groups, euclidean distance
permanova <- adonis2(asv_t ~ Sequencing, data = metadata, permutations=999, method = "euclid")

permanova <- adonis2(asv_t ~ DepthGroup, data = metadata, permutations=999, method = "euclid")

permanova

# LDM PERMANOVAS for paired-sample design
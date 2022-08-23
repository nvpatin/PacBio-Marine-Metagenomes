#load
library(ggplot2)
library(ggfortify)
library(RColorBrewer)
library(vegan)
library(reshape)
library(wesanderson)

# load the mash output as distance matrix
df <- data.frame(PacBio_Illumina_sourmash_compare_distance_v3, row.names=1)
colnames(df) <- rownames(df)
dist_mat <- as.matrix(df)

# load the metadata
## SAMPLES MUST BE IN THE EXACT SAME ORDER AS THE DISTANCE MATRIX!!! ##
metadata <- data.frame(sourmash_compare_metadata)

# calculate NMDS with user-supplied distances, 2 dimensions
nmds <- metaMDS(dist_mat, dist, k=2, trace=FALSE, engine="isoMDS", autotransform=FALSE)
nmds <- metaMDS(dist_mat, dist, k=3, trace=FALSE, engine="isoMDS", autotransform=FALSE)
nmds <- metaMDS(dist_mat_subset, dist, k=3, trace=FALSE, engine="isoMDS", autotransform=FALSE)
# nmds <- metaMDS(dist_mat, dist, k=2, trace=FALSE) 

MDS1 = nmds$points[,1]
MDS2 = nmds$points[,2]
MDS3 = nmds$points[,3]

#NMDS
NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, MDS3=MDS3,
                  Seq=metadata$Sequencing.Platform) 

#plot
g <- ggplot(NMDS, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(col=Seq), size=4) + 
  theme_classic() + 
  scale_color_manual(values=c("#00AFBB", "#E7B800")) +
  labs(title="sourmash compare") # + geom_text(aes(label=metadata$ID), size=2)

g
# scale_color_manual(values=wes_palette("Zissou1", n=5)) +
# scale_color_manual(values=c("#00AFBB", "#E7B800")) +
# scale_color_brewer(palette="Blues")
# values=wes_palette("GrandBudapest1", n=5)


#2018 depth group subset
NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, 
                  Arthropod=factor(depth_select_metadata$Arthropod.level), 
                  Cnidaria=factor(depth_select_metadata$Cnidaria.level),
                  Diatoms=factor(depth_select_metadata$Diatom.level),
                  Dinoflagellata=factor(depth_select_metadata$Dinoflagellate.level),
                  Depth=factor(depth_select_metadata$Depth_group), 
                  Season=depth_select_metadata$Season) # MDS3=MDS3, 

#2019
NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, MDS3=MDS3, Depth=factor(metadata$Depth_group_m),
                  Bottom_depth=factor(metadata$bottom_depth_group_m),
                  Date=factor(metadata$collection_date))

#2018-2019
NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, MDS3=MDS3,
                  NO3=metadata$NO3,
                  NO2=metadata$NO2,
                  Si=metadata$SIO4,
                  Chl=metadata$CHLA,
                  PO4=metadata$PO4,
                  Depth=factor(metadata$Depth_group_m),
                  Year=factor(metadata$Year))

#2018-2019 depth group subset
NMDS = data.frame(MDS1=MDS1, MDS2=MDS2, 
                  Depth=factor(depth_select_metadata$Depth_m),
                  Year=factor(depth_select_metadata$Year))

#2019
g <- ggplot(NMDS, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(col=Depth), size=3) + 
  theme_classic() + 
  labs(title="eCruise 2019 metagenomes: Mash Distances (k=25, MDS3)", size=2) +
  geom_text(aes(label=metadata$ID), size=1.5)

g

#2018-2019
g <- ggplot(NMDS, aes(x=MDS1, y=MDS2)) + 
  geom_point(aes(col=Chl, shape=Year), size=3) + 
  theme_classic() + 
  labs(title="eCruise 2018-2019 metagenomes: Mash Distances (k=25, MDS3)", size=2) #+geom_text(aes(label=metadata$ID))

g

ggsave("name")


# To select depth group(s)
depth_select_metadata <- metadata[metadata$Depth_group_m != "195-203",]
depth_select_IDs <- as.list(depth_select_metadata$ID)
dist_mat_subset <- dist_mat[rownames(dist_mat)%in%depth_select_IDs,
                            colnames(dist_mat)%in%depth_select_IDs]
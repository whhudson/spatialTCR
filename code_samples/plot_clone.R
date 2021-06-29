#This is an example script to plot localization of antigen receptors using seurat
#Note that the clones generated here are from downsampled sequencing read data 
#and do not include all clones or UMIs identified from the full data set

library("Seurat")
library("ggplot2")

#Open cellranger output with seurat
data_dir <- "path/to/visiumTCR/cellranger_output/"
brain16 <- Load10X_Spatial(data_dir)

#Open TCR data from combine_counts.py
clone_counts16 <- read.table("path/to/visiumTCR/output_counts.txt")
rownames(clone_counts16) <- paste(rownames(clone_counts16), "-1", sep="") #match rownames with 10X format by adding "-1"
clone_counts16 <- clone_counts16[colnames(brain16),] #ensure rownames match

#Add TCR data as metadata to Seurat object; identifier will be cloneN
for(i in 1:length(colnames(clone_counts16))){
  brain16 <- AddMetaData(brain16, clone_counts16[,i], colnames(clone_counts16)[i])
}

#Read in clone data from MiXCR; note that cloneCount here is NOT UMI-corrected
clone_info16 <- read.table("path/to/visiumTCR/mixcr_output/clones.clonotypes.ALL.txt",
                           sep = "\t", header = T, row.names = 1)

#Plot a sample clone; for example clone21, a TCR with CDR3 AA sequence CASSFRPGQPDNEQFF, TRBV7-8, and TRBJ2-1
#In scRNA-seq data, this clone is 82% in clusters B/C (less-exhausted)
SpatialPlot(brain16, features = "clone21") &
  scale_fill_gradient2(low = "white", mid = "blue", high = "red", midpoint = 0.5) & 
  scale_alpha_continuous(range=c(0,1))

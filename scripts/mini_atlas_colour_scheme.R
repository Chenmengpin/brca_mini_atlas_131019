# ==================================
# molecular subtype colour codes
# ==================================

# data-frame can be used in ComplexHeatmap
dfColourCodes <- data.frame(name = c("Basal", "Her2", "LumA", "LumB", "Normal", 
                                     "LumB_Her2", "HER2_ER", "ER", "TNBC", "HER2", "metaplastic",
                                     "Her2_Amp", "Her2_Non_Amp", "TBC",
                                     "LuminalA", "LuminalB", "LuminalB_HER2",
                                     "Normal_breast",
                                     "LumA_SC50", "LumB_SC50", "Her2_SC50", "Basal_SC50",
                                     "cancer", "unassigned", "normal"),
                            value = c("red", "purple", "darkblue", "cyan", "green", 
                                      "pink", "pink", "blue", "red", "purple", "yellow", 
                                      "black", "lightgrey", "yellow",
                                      "darkblue", "cyan", "pink",
                                      "grey50",
                                      "darkblue", "cyan", "pink", "red",
                                      "black", "grey50", "green")
)

# For ggplot
SC50_colours <- c("Basal_SC50" = "red", "Her2_SC50" = "pink", "LumA_SC50" = "darkblue", "LumB_SC50" = "cyan")
PAM50_colours <- c("Basal" = "red", "Her2" = "pink", "LumA" = "darkblue", "LumB" = "cyan", "Normal" = "green", "NA" = "grey90")
Clinical_IHC_colours <- c("HER2_ER" = "purple", "ER" = "blue", "TNBC" = "red", "HER2" = "pink")
inferCNV_colours <- c("cancer" = "black", "unassigned" = "lightgrey", "normal" = "green")


# =========================================
# SampleIDs - colours fixed for each sample
# =========================================

library("Polychrome")

# define sample colours
set.seed(9641)
sample_colours <- data.frame(
  sample_id = c("CID3586", "CID44041", "CID4461", "CID4463", "CID4465", "CID4471",
    "CID4495", "CID44971", "CID44991", "CID4513", "CID4515", "CID3921",
    "CID45171", "CID4523", "CID4530N", "CID4535", "CID3941", "CID3948",
    "CID3963", "CID4066", "CID4067", "CID4290A", "CID4398", "CID3838", "CID3946",
    "CID4040"),
  colour = createPalette(
    26, c("#2A95E8", "#E5629C"), range = c(10, 60), M = 100000
  ),
  stringsAsFactors = F
)
rownames(sample_colours) <- sample_colours$sample_id

# input sample IDs vector and sample colours df to following function
# to get corresponding colours:
get_sample_colours <- function(sample_ids, sample_colours) {
  return(sample_colours[sample_ids,]$colour)
}

# e.g:
sample_ids <- c("CID4066", "CID44041", "CID4290A", "CID3963", "CID3838", "CID3946", "CID4040")
col_eg <- get_sample_colours(sample_ids, sample_colours)
# check:
barplot(seq(1:length(col_eg)), col=col_eg)


####################################

#library("Polychrome")
#
#set.seed(9641)
#num_samples <- length(unique(seurat_10X@meta.data$orig.ident))
#
#colour_pallete_mini_atlas <- createPalette(num_samples, c("#2A95E8", "#E5629C"), range = c(10, 60), M = 100000)
#dd <- computeDistances(colour_pallete_mini_atlas)
#names(colour_pallete_mini_atlas) <- unique(seurat_10X@meta.data$orig.ident)


### you will need to install all the packages listed below
library(stats)
library(ade4)
library(ape)
library(msa)
library(seqinr)
library(phylogram)
library(ggdendro)
library(ggtree)
library(dendextend)
library(circlize)
library(randomcoloR)

## loading GISAID sample information
gisaid_samples <- read.delim(file = "hcov_global.tsv", sep = "\t", stringsAsFactors = F)

### function to build phylogenetic tree from msa file
## msa_file - path to multiple sequnce alignemnt file
## gisaid_annotation - Table of GISAID samples
## out_file_name - name of the pdf file which will contain phylogenetic tree

phylogenetic_tree_from_fasta <- function(msa_file, gisaid_annotation, out_file_name) {
  
  msa_results <- readAAMultipleAlignment(filepath = msa_file)
  
  myAln2 <- msaConvert(msa_results, type="seqinr::alignment") # this object is a list object with 4
  
  d <- dist.alignment(myAln2, "identity")
  
  myTree <- nj(d)
  
  gisaid_annotation <- gisaid_annotation[which(gisaid_annotation$strain %in% myTree$tip.label),c("strain", "country", "region")]
  
  arm_sample_list <- data.frame(strain = myTree$tip.label[which(!(myTree$tip.label %in% gisaid_annotation$strain))],
                                country = c("Reference", rep("Armenia", (length(myTree$tip.label[which(!(myTree$tip.label %in% gisaid_annotation$strain))]) - 1) )),
                                region = c("Reference", rep("Armenia", (length(myTree$tip.label[which(!(myTree$tip.label %in% gisaid_annotation$strain))]) - 1) )),
                                stringsAsFactors = F
  )
  
  gisaid_annotation <- rbind(arm_sample_list, gisaid_annotation)
  
  ## automated coloring which is not so good  
  set.seed(1)
  label_colors <- setNames(object = randomcoloR::distinctColorPalette(k = length(unique(gisaid_annotation$region))),
                           nm = unique(gisaid_annotation$region)
                           )
  
  ## manual coloring which is better
  # label_colors <- setNames(object = c("#191970", "#006400", "#ff0000", "#ffd700", "#00ff00", "#00ffff", "#ff00ff", "#ffb6c1"),
  #                          nm = unique(gisaid_annotation$region)
  # )
  
  gisaid_annotation <- cbind(gisaid_annotation, label_colors[gisaid_annotation$region])
  
  colnames(gisaid_annotation)[ncol(gisaid_annotation)] <- "region_color"
  
  gisaid_annotation$region_color <- as.character(gisaid_annotation$region_color)
  
  rownames(gisaid_annotation) <- gisaid_annotation$strain
  
  gisaid_annotation$tree_label <- paste0(1:length(gisaid_annotation$region), "_", gisaid_annotation$region)
  
  node_names_and_colors <- setNames(nm = gisaid_annotation$tree_label, object = gisaid_annotation$region_color)
  
  myTree$tip.label <- gisaid_annotation[myTree$tip.label,"tree_label"]
  
  
  corona_dentro <- phylogram::as.dendrogram.phylo(myTree)
  
  corona_dendro_data <- ggdendro::dendro_data(corona_dentro)
  
  node_names <- as.character(corona_dendro_data$labels$label)
  
  node_names_and_colors <- node_names_and_colors[node_names]
  
  
  modified_dent <- corona_dentro %>% 
    # set("nodes_pch", 19)  %>% 
    # set("nodes_cex", 0.7) %>%
    set("labels_cex", 0.5) %>%
    set("nodes_col", "orange") %>%
    set("labels_col", labels = names(node_names_and_colors), value = unname(node_names_and_colors))
  
  
  pdf(file = out_file_name, paper='A4r', height = 8.27, width = 11.69, onefile = T)
  
  par(mar=c(1,1,1,7))
  
  circlize_dendrogram(modified_dent, facing = "outside", dend_track_height	= 0.7)
  
  legend("topright", unique(gisaid_annotation[,3:4])$region, fill = unique(gisaid_annotation[,3:4])$region_color, cex = 0.7, y.intersp = 0.6,border = NA, box.lty=0, bty = "o")
  
  
  dev.off()
}

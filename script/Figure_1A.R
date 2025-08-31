#script to generate Figure 1A - Pairwise genetic distance heatmap

#load libraries
library(ape)
library(ggplot2)
library(pheatmap)


#SET PATH DIRECTORY
path='../data/'


#FUNCTION TO READ FASTA AND CALCULATE PAIRWISE DISTANCES
read_and_calculate_distances <- function(gene_name, file) {
  sequences <- read.dna(file.path(paste0(path, file)), format = 'fasta')
  # Calculate genetic distances
  distances <- dist.dna(sequences, model = "N", pairwise.deletion = TRUE,as.matrix=TRUE)
  diag(distances) <- NA
  distances[upper.tri(distances)] <- NA
  return (distances)
}


all = read_and_calculate_distances('ALL','concatenated.fasta')
all_masked <- all
all_masked[is.na(all_masked)] <- -999

# Create label matrix and hide labels in masked positions
labels <- matrix(sprintf("%.0f", all_masked), nrow = nrow(all_masked))
labels[all_masked == -999] <- ""

# Color palette: white for dummy values (-999)
real_col <- colorRampPalette(c( "#ffffffff"))(255)
full_palette <- c("#FFFFFF", real_col)  # white for first value (dummy)
breaks <- c(-1000, seq(0, max(all_masked[all_masked != -999], na.rm = TRUE), length.out = 256))

# Read annotations
annotation <- read.table(paste0(path, 'sample_annotation.txt'), header = TRUE, row.names = 'sample')


#GENERATE FIGURE PLOT
figure <-   pheatmap( all_masked, display_numbers = labels, cellwidth = 30, cellheight = 30,
                        legend = FALSE,cluster_rows = FALSE,cluster_cols = FALSE,
                        number_format = "",color = full_palette,
                        breaks = breaks,annotation_row = annotation,
                        border_color = '#d9d9d9',
                        annotation_colors = list(host = c('swine' = '#d95f02', 'human' = '#073ae1ff')))


#export figure
ggsave('../figure/Figure_1A.pdf', dpi=600,
       height=110, width = 110, units = "mm", plot=figure )


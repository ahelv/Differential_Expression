prep_Fisher <- function(RNA, export_name1, peak_file, export_name2, 
                        peak_sort_name, nearest_gene_name, nearest_no_diff_name, 
                        distance = 100000){
  # downregulated genes 
  genes_down <- RNA[which(RNA$LogFC < 0) ,]
  # upregulated genes 
  genes_up <- RNA[which(RNA$LogFC > 0) ,]
  
  ## keep only one instance of each gene
  # down regulated
  genes_down_sort <- genes_down[order(genes_down[, 4]) ,]
  genes_down <- subset(genes_down_sort, !duplicated(ext_gene))
  names(genes_down)[1] <- "name"
  # up regulated
  genes_up_sort <- genes_up[order(genes_up[, 4]) ,]
  genes_up <- subset(genes_up_sort, !duplicated(ext_gene))
  names(genes_up)[1] <- "name"
  
  # find down regulated genes not in gene list 
  unlisted_down <- genes_down[which(!(genes_down$name %in% gene_bed$name)) ,]
  unlisted_up <- genes_up[which(!(genes_up$name %in% gene_bed$name)) ,]
  
  ## find set of non-dysregulated genes
  # remove down-regulated genes
  no_down <- gene_bed[-which(gene_bed$name %in% genes_down$name) ,]
  # remove up-regulated genes
  no_diff <- no_down[-which(no_down$name %in% genes_up$name) ,]
  
  # export
  write.table(no_diff, export_name1, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # Find Peaks that are in the gene list 
  if(ncol(peak_file) > 3) {
    peak_bed <- peak_file[c(2,3,4)]
  }else{
    peak_bed <- peak_file
  }
  full_down <- merge(gene_list_unique, genes_down, by = "name", all.y = TRUE)
  down_bed <- full_down[c(2,3,4,1,5,6,7)]
  full_up <- merge(gene_list_unique, genes_up, by = "name", all.y = TRUE)
  up_bed <- full_up[c(2,3,4,1,5,6,7)]
  
  # export
  write.table(peak_bed, export_name2, quote = FALSE, 
              row.names = FALSE, col.names = FALSE, sep = "\t")
  
  # use bedtools to perform sorting and intersection
  # NOTE: change the directory to where your bedtools is stored
  sort_command <- paste("/bedtools2/bin/bedtools sort -i", export_name2, ">", peak_sort_name)
  system(sort_command)
  int_cmd <- paste("/bedtools2/bin/bedtools closest -a", peak_sort_name,
                   "-b gene_bed_sort.bed -D b -t first >", nearest_gene_name)
  system(int_cmd)
  int_cmd_no_diff <- paste("/bedtools2/bin/bedtools closest -a", peak_sort_name,
                           "-b genes_no_diff_sort.bed -D b -t first >", nearest_no_diff_name)
  system(int_cmd_no_diff)
  
  file_path <- paste0("/Directory/", nearest_gene_name)
  nearest_gene <- read.delim(file_path, header = FALSE)
  file_path_no_diff <- paste0("/Directory/", nearest_no_diff_name)
  nearest_gene_no_diff <- read.delim(file_path_no_diff, header = FALSE)
  
  # dysregulated
  names(nearest_gene) <- c("chr_peak", "start_peak", "stop_peak",
                           "chr_gene", "start_gene", "stop_gene", "name",
                           "distance")
  gene_nearest_1 <- (nearest_gene[which(nearest_gene[, 8] <= 0), ])
  gene_nearest_full <- gene_nearest_1[which(gene_nearest_1[, 8] > -distance), ]
  gene_nearest <- unique(gene_nearest_1[which(gene_nearest_1[, 8] > -distance), 7])
  (perc_overlap_nearest <- length(gene_nearest) / nrow(gene_bed_flank))
  
  # non-dysregulated 
  names(nearest_gene_no_diff) <- c("chr_peak", "start_peak", "stop_peak",
                                   "chr_gene", "start_gene", "stop_gene", "name",
                                   "distance")
  gene_nearest_no_diff_1 <- (nearest_gene_no_diff[which(nearest_gene_no_diff[, 8] <= 0), ])
  gene_nearest_no_diff_full <- gene_nearest_no_diff_1[which(gene_nearest_no_diff_1[, 8] > -distance), ]
  gene_nearest_no_diff <- unique(gene_nearest_no_diff_1[which(gene_nearest_no_diff_1[, 8] > -distance), 7])
  (perc_overlap_nearest_no_diff <- length(gene_nearest_no_diff) / nrow(no_diff))
  
  # Find dysregulated genes near target of interest (originally MITF)
  down_nearest_target <- merge(gene_nearest_full, down_bed, by = "name")
  num_down_near_target_tot <- length(unique(down_nearest_target$name))
  
  up_nearest_target <- merge(gene_nearest_full, up_bed, by = "name")
  num_up_near_target_tot <- length(unique(up_nearest_target$name))
  
  # setup
  DR <- num_down_near_target_tot
  UR <- num_up_near_target_tot
  N <- length(gene_nearest_no_diff)
  tot <- DR + UR + N
  tot_DR <- nrow(genes_down)
  tot_UR <- nrow(genes_up)
  tot_N <- nrow(gene_bed) - (tot_DR + tot_UR)
  tot_tot <- nrow(gene_bed)
  NO_DR <- tot_DR - DR
  NO_UR <- tot_UR - UR
  NO_N <- tot_N - N
  NO_tot <- NO_DR + NO_UR + NO_N
  
  # full contingency table 
  full_table <- matrix(c(DR, NO_DR, tot_DR, 
                         UR, NO_UR, tot_UR, 
                         N, NO_N, tot_N,
                         tot, NO_tot, tot_tot), nrow = 3, ncol = 4)
  rownames(full_table) <- c("Nearest target", "Not Nearest target", "Total")
  colnames(full_table) <- c("Down-Regulated", "Up-Regulated", "Neither", "Total")
  
  # Fisher's Exact Test
  (down_reg <- matrix(c(DR, (tot_DR - DR),
                        N, (tot_N - N)), nrow = 2, ncol = 2))
  
  (up_reg <- matrix(c(UR, (tot_UR - UR),
                      N, (tot_N - N)), nrow = 2, ncol = 2))
  fisher.down <- fisher.test(down_reg)
  fisher.up <- fisher.test(up_reg)
  
  return(list(full_table = full_table, fisher.down = fisher.down, fisher.up = fisher.up))
}
library(data.table)
library(dplyr)
library(IRanges)
library(stringdist)
start <- Sys.time()
#setwd("../../../../Desktop/snpdata/")
sample_names <- list.files(pattern = ".snps")
ref_name     <- list.files(pattern = ".ptt")

# reading tables
list_of_tables <- lapply(sample_names, fread)

s_ids <- sapply(sample_names, function(file_name) strsplit(file_name, split = "[.]")[[1]][1])
names(list_of_tables) <- s_ids

#renaming and deleting extra col
list_of_tables <- lapply(list_of_tables, function(df) {
  df <- df[ , 1:3]
  colnames(df) <- c("pos_ref", "ref", "alt")
  return(df)
})

nuc <- c("A", "C", "T", "G")

# filterin anything exect "A","C", "T", "G"
filtered_tables <- lapply(list_of_tables, function(df) {
  indexes <- intersect(which(df$alt %in% nuc), which(df$ref %in% nuc))
  print(names(df))
  return(df[indexes, ])
})

#renaming columns for future merging
filtered_tables <- lapply(s_ids, function(s_id) {
  temp <- filtered_tables[[s_id]]
  colnames(temp)[3] <- s_id
  return(temp)
  })

# reading ref
ann <- fread(ref_name, header = T)

# adding "start", "end" col
ann$start <- sapply(ann$Location, function(entry) 
  as.numeric(strsplit(entry, split = "[..]")[[1]][1]) ) 
ann$end   <- sapply(ann$Location, function(entry) 
  as.numeric(strsplit(entry, split = "[..]")[[1]][3]) ) 

# creating Iranges object for cheking positions
ir <- IRanges(start = ann$start, end = ann$end)
cover <- coverage(ir)

# merging al SNPS in one table
big_table <- Reduce(function(x, y) merge.data.frame(x, y, by = c("pos_ref", "ref"), all = T), filtered_tables)
# removing SNP that are not in range of annotation
big_table <- subset(big_table, pos_ref < max(ann$end))

indexes   <- cover[big_table$pos_ref] %>% as.vector() %>% as.logical()
big_table <- big_table[indexes, -1 ]
# filling gaps with references
big_table <- sapply(big_table, function(col) {
  col[which(is.na(col))] <- big_table$ref[which(is.na(col))]
  return(col)
})

# collapsing to strings 
SNPS <- apply(big_table, 2, paste, collapse = "")

MATRX <- as.matrix(stringdistmatrix(SNPS, method = 'hamming', useNames = "names"))
Sys.time() - start


phyl_tree <- hclust(as.dist(MATRX), method = "single")
plot(phyl_tree)

#growth rate predictions using gRodon
devtools::install_github("jlw-ecoevo/gRodon2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings", force = TRUE)
BiocManager::install("coRdon")
install.packages("matrixStats")

library(gRodon)
library(Biostrings)

#read in genes
path_prodigal <- "C:/Users/byersa/OneDrive - Lincoln University/Documents/Rutherford/Metagenomics/Growth potential/prodigal_files"
list.files(path_prodigal)
prodigal_names <- sort(list.files(path_prodigal, pattern="_prodigal_genes.fna", full.names = TRUE))
# Extract sample names
sample.names <- sapply(strsplit(basename(prodigal_names), "_prodigal_genes.fna"), `[`, 1)

#read in files
prodigal_files <- lapply(prodigal_names, function(x) readDNAStringSet(x))
gene_IDs <- lapply(prodigal_files, function(x) gsub(" .*","",names(x)))
names(prodigal_files) <- sample.names
names(gene_IDs) <- sample.names

#Search for genes annotated as ribosomal proteins
path_blast <- "C:/Users/byersa/OneDrive - Lincoln University/Documents/Rutherford/Metagenomics/Growth potential/blast_files"
list.files(path_blast)
blast_names <- sort(list.files(path_blast, pattern="_ribosomal_protiens.blast", full.names = TRUE))
# Extract sample names
sample.names2 <- sapply(strsplit(basename(blast_names), "_ribosomal_protiens.blast"), `[`, 1)

#read in files
blast_files <- lapply(blast_names, function(x) read.table(x, sep = '\t', header = FALSE))
names(blast_files) <- sample.names2
lapply(blast_files, function(x) unique(x$V1))

highly_expressed <- gene_IDs$La_1 %in% blast_files$La_1$V1
highly_expressed <- lapply(seq_along(gene_IDs), function(i) gene_IDs[[i]] %in% blast_files[[i]]$V1)
names(highly_expressed) <- sample.names
summary(highly_expressed$La_1)

#read in coverage depths
read_depths <- "C:/Users/byersa/OneDrive - Lincoln University/Documents/Rutherford/Metagenomics/Growth potential/coverage_files"
list.files(read_depths)
depth_names <- sort(list.files(read_depths, pattern="_gene_coverage.txt", full.names = TRUE))
# Extract sample names
sample.names3 <- sapply(strsplit(basename(depth_names), "_gene_coverage.txt"), `[`, 1)

#read in files
depth_files <- lapply(depth_names, function(x) read.delim(x, stringsAsFactors = FALSE))
names(depth_files) <- sample.names3
depths <- lapply(depth_files, function(x) x$meandepth)

#make sure in correct order
for (i in seq_along(depths)) {
names(depths[[i]]) <- depth_files[[i]]$X.rname  
}
depth_of_coverage <- lapply(seq_along(depths), function(i) depths[[i]][gsub(" .*", "", names(prodigal_files[[i]]))])
head(depth_of_coverage, 10)

##RUNNING GRODON, per sample
growthRate <- lapply(seq_along(prodigal_files), function(i) predictGrowth(prodigal_files[[1]], 
                    highly_expressed[[1]], depth_of_coverage = depth_of_coverage[[i]], 
                    bg="individual", mode = "meta_testing"))

# APPENDING STEREOSCOPE DATA
# Sunny Z. Wu
# 
#
# Run in screen 
# qrsh -pe smp 4 -l mem_requested=10G -P TumourProgression
# source activate r_spatial
# R
# 
#
# 01: SETUP---------------------------------------

# setup
library(zeallot)
library(Seurat, lib.loc = "/rd1/apps/R-3.5.1/library_ext")
library(STutility)
library(ggplot2)
library(magrittr)
library(magick)
library(dplyr)
library(RColorBrewer)
library(ggpubr)

# DIRECTORY
dir.create("append_stereoscope", showWarnings = F)

# 02: LOAD DATA -----------------------------------------------------------

# processed STUtility objects
se.list <- readRDS("data_processing/RDATA_Visium_brca_objects.Rdata")

infoTable <- read.csv(file = "data_processing/infoTable.csv", header = T, stringsAsFactors = F)

# 03: APPEND STEREOSCOPE RESULTS ----------------------------------------------

# load data from each tier and append to object as reduction data slot
for(sample in unique(infoTable$patientid)){
  print(sample)
  
  for(tier in c("major", "minor", "subset")){
    print(tier)
    if(tier == "major"){
      temp_file_string <- "W.2020-06-04043728.683825.tsv"
    }
    if(tier == "minor"){
      temp_file_string <- "W.2020-06-04043854.061760.tsv"
    }
    if(tier == "subset"){
      temp_file_string <- "W.2020-06-05030141.613517.tsv"
    }
    
    # load stereoscope
    sc.files <- c(paste0("to-sunny-zip/",tier,"/",temp_sampleid,"-count-matrix.tsv/",temp_file_string))
    sc.data.tnbc <- read.table(sc.files, sep = "\t", header = T, row.names = 1)
    
    # append
    sc <- sc.data.tnbc
    meta.data <- GetStaffli(se.list[[sample]])@meta.data
    meta.data$spotids <- rownames(meta.data)
    tissue_positions <- subset(read.csv(infoTable$spotfiles[infoTable$patientid == sample], 
                                        header = F, 
                                        stringsAsFactors = FALSE), 
                               V2 == 1)
    convert_ids <- 
      data.frame(xy = paste(tissue_positions$V5, 
                            tissue_positions$V6, 
                            sep = "x"), 
                 spotid = tissue_positions$V1, 
                 stringsAsFactors = FALSE)
    rownames(sc) <- paste0(convert_ids$spotid, "_", sample)
    sc.data.tnbc <- sc
    
    # Set sc.data as cell.embeddings
    cell.embeddings <- sc.data.tnbc
    colnames(cell.embeddings) <- paste0("tnbc_", 1:ncol(cell.embeddings))
    
    if(!length(rownames(cell.embeddings))==length(colnames(se.list[[sample]]))){
      print("barcodes dont match")
      cell.embeddings <- cell.embeddings[rownames(cell.embeddings) %in% colnames(se.list[[sample]]),,drop=F]
      cell.embeddings <- cell.embeddings[colnames(se.list[[sample]]),,drop=F]
      print("filtered")
    }
    
    if(all.equal(colnames(se.list[[sample]]), rownames(cell.embeddings))){
      print("barcodes match")
      
      reduction.data <- CreateDimReducObject(
        embeddings = as.matrix(cell.embeddings),
        loadings = matrix(),
        assay = "SCT",
        key = paste0("tnbc",tier,"_")
      )
      
      # Store reduction object in reductions slot
      se.list[[sample]][[paste0("stereoscope_tnbc",tier)]] <- reduction.data
      
    } else{
      print("barcodes dont match")
    }
    
    
  }
}

# SAVE RDS ----------------------------------------------------------------

saveRDS(se.list, "RDATA_Visium_brca_objects_stereoscope.Rdata")

#! /usr/bin/env Rscript

######################################################################################
####
#### Description: This script classifies sequences in a fasta file using the RDP 
####              classifier with the eHOMD training dataset. 
####
#### Usage:       RDPclassify.R file.fasta eHOMDtrainingSet.fa.gz
####               
#### Author:      Yanmei Huang
#### Version:     1.0
#### Date:        2017-12-01
####
######################################################################################

message(paste0("Analysis started at ", Sys.time()))
message()

#### {r get commandArgs}
# Make this backwards-compatible with "R --vanilla < foo.R" syntax
offset = 0
if (paste(commandArgs()[1], commandArgs()[2]) == "/usr/lib/R/bin/exec/R --vanilla")
{
  message("Running as R --vanilla ....")
  offset = 3
}

infile = as.character(commandArgs()[6 - offset])
trainingSet = as.character(commandArgs()[7 - offset])

####{r, warning=FALSE, message=FALSE, load required packages}
library(tidyverse); packageVersion("tidyverse")
library(Biostrings); packageVersion("Biostrings")
library(dada2); packageVersion("dada2")

####{r a function to reduce taxonomy level from 8 to 7}
# this function replace "NA"s in "Species" with values
# from "SuperSpecies" 
eightTOseven = function(taxa)
{
  if (length(names(taxa)) > 0)
  {
    taxa$boot[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$boot[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$tax[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax = cbind(taxa$tax[, 1:6], taxa$tax[, 8])
    taxa$boot = cbind(taxa$boot[, 1:6], taxa$boot[, 8])
  } else {
    taxa[is.na(taxa[, 'Species']), 'Species'] = taxa[is.na(taxa[, 'Species']), 'SuperSpecies']
    taxa = cbind(taxa[, 1:6], taxa[, 8])
  }
  dimnames(taxa$tax)[[2]][7] = "Species"
  dimnames(taxa$boot)[[2]][7] = "Species"
  return (taxa)
}

####{r read fasta file}

input = readDNAStringSet(infile)

####{r run RDP Classifier}

classified = assignTaxonomy(paste(input), 
                                   trainingSet, 
                                   outputBootstraps=TRUE, 
                                   taxLevels=c("Kingdom", "Phylum", 
                                               "Class", "Order", "Family", 
                                               "Genus", "SuperSpecies", 
                                               "Species"))

# # save workspace image
# outfile = paste0(gsub(".fasta", "", infile), ".classified.", dis, ".Rdata")
# save.image(file = outfile )

classified = eightTOseven(classified)
####{r format output table}
classified.tax = as.data.frame(classified$tax)
rownames(classified.tax) = names(input)

classified.boot = as.data.frame(classified$boot)
rownames(classified.boot) = names(input)

####{r write output file}
outfile = paste0(gsub(".fasta", "", infile), ".taxonomy.csv")
write.csv(classified.tax, file = outfile, quote = FALSE)

outfile = paste0(gsub(".fasta", "", infile), ".boot.csv")
write.csv(classified.boot, file = outfile, quote = FALSE)

message()
message(paste0("Analysis finished at ", Sys.time()))
message()

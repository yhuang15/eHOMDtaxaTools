#! /usr/bin/env Rscript

######################################################################################
####
#### Description: This script runs RDP classifier on MED output fasta using a
####              training data set
####
#### Usage:       classifyMEDout.R MEDoutfile trainingset_file
####               
#### Author:      Yanmei Huang
#### Version:     1.0
#### Date:        2017-11-17
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
Tfile = as.character(commandArgs()[7 - offset])
dis = gsub(".*Taxa_", "", Tfile)
dis = gsub(".seqTab.fa.gz", "", dis)



####{r, warning=FALSE, message=FALSE, load required packages}
library(tidyverse); packageVersion("tidyverse")
library(Biostrings); packageVersion("Biostrings")
library(readr); packageVersion("readr")
library(stringr); packageVersion("stringr")
library(dada2); packageVersion("dada2")

####{r read fasta file}

nodes = readDNAStringSet(infile)

####{r run RDP Classifier}

SuperTaxa_8levels = assignTaxonomy(paste(nodes), 
                                   Tfile, 
                                   outputBootstraps=TRUE, 
                                   taxLevels=c("Kingdom", "Phylum", 
                                               "Class", "Order", "Family", 
                                               "Genus", "SuperSpecies", 
                                               "Species"))

####{r format output table}
SuperTaxa_8levels = as.data.frame(cbind(select(as.data.frame(SuperTaxa_8levels$tax), 
                                               Kingdom, Phylum, Class, Order, Family, 
                                               Genus, SuperSpecies, Species), 
                                        select(as.data.frame(SuperTaxa_8levels$boot), 
                                               Kingdom, Phylum, Class, Order, Family, 
                                               Genus, SuperSpecies, Species)))
colnames(SuperTaxa_8levels) = c("Kingdom","Phylum","Class","Order","Family","Genus",
                                "SuperSpecies","Species", "Kingdom.B","Phylum.B",
                                "Class.B","Order.B","Family.B","Genus.B",
                                "SuperSpecies.B","Species.B")

rownames(SuperTaxa_8levels) = names(nodes)
SuperTaxa_8levels$Reads = gsub(".*size:", "", names(nodes))
SuperTaxa_8levels$Reads = as.integer(SuperTaxa_8levels$Reads)

####{r write output file}
outfile = paste0(gsub(".fasta", "", infile), ".classified.", dis, ".csv")
write.csv(SuperTaxa_8levels, file = outfile, quote = FALSE)


message()
message(paste0("Analysis finished at ", Sys.time()))
message()

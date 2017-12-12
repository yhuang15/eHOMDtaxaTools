#! /usr/bin/env Rscript

######################################################################################
####
#### Description: This script runs RDP classifier on MED output fasta using an
####              eHOMD training data set
####
#### Usage:       classifyMEDout.R MEDoutfile trainingset_file
####               
#### Author:      Yanmei Huang
#### Version:     1.11
#### Date:        2017-12-11
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
dis = gsub(".fa.gz", "", dis)

#### {r set minimal boot strap value. eHOMD research group decided it should be 90 }
minBoot_val = 90

####{r load required packages}
library(Biostrings); packageVersion("Biostrings")
library(dada2); packageVersion("dada2")

####{r, a function to merge SuperSpecies with Species}
eightTOseven = function(taxa)
{
  if (length(names(taxa)) > 0)
  {
    taxa$boot[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$boot[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$tax[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax = cbind(taxa$tax[, 1:6], taxa$tax[, 8])
    taxa$boot = cbind(taxa$boot[, 1:6], taxa$boot[, 8])
    dimnames(taxa$tax)[[2]][7] = "Species"
    dimnames(taxa$boot)[[2]][7] = "Species"
  } else {
    taxa[is.na(taxa[, 'Species']), 'Species'] = taxa[is.na(taxa[, 'Species']), 'SuperSpecies']
    taxa = cbind(taxa[, 1:6], taxa[, 8])
    dimnames(taxa)[[2]][7] = "Species"
  }
  return (taxa)
}

####{r read fasta file}

MEDnodes = readDNAStringSet(infile)

####{r run RDP Classifier}

classified = assignTaxonomy(paste(MEDnodes), 
                                   Tfile, 
                                   outputBootstraps=TRUE,
                                   minBoot = minBoot_val,
                                   taxLevels=c("Kingdom", "Phylum", 
                                               "Class", "Order", "Family", 
                                               "Genus", "SuperSpecies", 
                                               "Species"))

# save workspace image
outfile = paste0(gsub(".fasta", "", infile), ".classified.", dis, ".Rdata")
save.image(file = outfile )

# merge SuperSpecies with Species
classified = eightTOseven (classified)

####{r format output table}
classified = cbind(as.data.frame(classified$tax)[, c('Kingdom', 'Phylum', 'Class', 
                                                     'Order', 'Family', 
                                                     'Genus', 'Species')], 
                  as.data.frame(classified$boot)[, c('Kingdom', 'Phylum', 'Class', 
                                                     'Order', 'Family', 
                                                     'Genus', 'Species')],
                  row.names = names(MEDnodes))
                                 
colnames(classified) = c("Kingdom","Phylum","Class","Order","Family","Genus",
                                "Species", "Kingdom.B","Phylum.B",
                                "Class.B","Order.B","Family.B","Genus.B",
                                "Species.B")

classified$Reads = gsub(".*size:", "", names(MEDnodes))
classified$Reads = as.integer(classified$Reads)

####{r write output file}
outfile = paste0(gsub(".fasta", "", infile), ".classified.", dis, ".csv")
write.csv(classified, file = outfile, quote = FALSE)

message()
message(paste0("Analysis finished at ", Sys.time()))
message()

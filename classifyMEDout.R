#! /usr/bin/env Rscript

######################################################################################
####
#### Description: This script runs RDP classifier on MED output fasta using a
####              eHOMD training data set
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
dis = gsub(".fa.gz", "", dis)

####{r load required packages}
library(Biostrings); packageVersion("Biostrings")
library(dada2); packageVersion("dada2")

####{r read fasta file}

MEDnodes = readDNAStringSet(infile)

####{r run RDP Classifier}

classified = assignTaxonomy(paste(MEDnodes), 
                                   Tfile, 
                                   outputBootstraps=TRUE, 
                                   taxLevels=c("Kingdom", "Phylum", 
                                               "Class", "Order", "Family", 
                                               "Genus", "SuperSpecies", 
                                               "Species"))

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

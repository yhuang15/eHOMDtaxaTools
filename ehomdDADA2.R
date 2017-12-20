#! /usr/bin/env Rscript

######################################################################################
####
#### Description: This script performs DADA2 denoising for paired-end 16S sequencing 
####              reads and assign taxonomy to the denoised sequences using the RDP 
####              classifier with a training dataset provided by eHOMD specifically
####              for 16S data from microbiome samples of the aerodigestive track. 
####              * The script runs directly above the directory containing the raw 
####                fastq files. Successful completion of the script creates an 
####                output directory (specified in the "outdir" command argument)
####                containing quality-filtered reads and the read count table as  
####                well as various intermediate R objects created during the run so  
####                that if run terminated before finishing, the saved intermediate  
####                R objects can be used in two additional scirpts, "assignTax.R"  
####                and "fixBracketNcreateWB.R" to proceed to finish.
####              * Parameters for the analysis are set in the configFile
####              * The raw output is in "?_countTable_minBoot?.txt". 
####
#### Usage:       Dada2rdp.R fastq_dir trainingSet configFile output_dir
####
#### Author:      Yanmei Huang
#### Version:     1.0
#### Date:        2017-12-12
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
fastqDir = as.character(commandArgs()[6 - offset])
trainingSet = as.character(commandArgs()[7 - offset])
configFile = as.character(commandArgs()[8 - offset])
outdir = as.character(commandArgs()[9 - offset])

projectName = sapply(strsplit(outdir, "/"), `[`, length(strsplit(outdir, "/")[[1]]))

#### {r get parameter values }
parameters = read.table(configFile, sep="\t", header=TRUE, stringsAsFactors = FALSE)
Vregion = as.character(parameters$value[which(parameters$parameter == "Vregion")])
minBoot_val = as.numeric(parameters$value[which(parameters$parameter == "minBoot")])
max_consist = as.numeric(parameters$value[which(parameters$parameter == "max_consist")])
truncLen_L = as.integer(parameters$value[which(parameters$parameter == "truncLen_L")])
truncLen_R = as.integer(parameters$value[which(parameters$parameter == "truncLen_R")])
maxEE_L = parameters$value[which(parameters$parameter == "maxEE_L")]
maxEE_R = parameters$value[which(parameters$parameter == "maxEE_R")]
minQ_L = parameters$value[which(parameters$parameter == "minQ_L")]
minQ_R = parameters$value[which(parameters$parameter == "minQ_R")]

message("Parameters are set as:")
message(paste0("minBoot = ", minBoot_val))
message(paste0("MAX_CONSIST = ", max_consist))
message(paste0("truncLen = c(", truncLen_L,",", truncLen_R, ")"))
message(paste0("maxEE = c(", maxEE_L, ",", maxEE_R, ")"))
message(paste0("minQ = c(", minQ_L,",", minQ_R, ")"))
message()

#### {r load packages}
library(dada2); message(paste0("using Dada2 version ", packageVersion("dada2")))
library(ShortRead); message(paste0("using ShortRead version ", packageVersion("ShortRead")))
library(ggplot2); message(paste0("using ggplot2 version ", packageVersion("ggplot2")))
message()

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

#### {r creat output directory}
if(!file_test("-d", outdir)) dir.create(outdir)

#### {r set Dada2 parameters}
setDadaOpt(MAX_CONSIST = max_consist)

#### {r read files}
fnsFiles = list.files(fastqDir)
fns = file.path(fastqDir, fnsFiles)
fastqs = fns[grepl(".fastq.gz$", fns)]
fastqs = sort(fastqs) # Sort ensures forward/reverse reads are in same order
fnFs = fastqs[grepl("_R1.fastq.gz", fastqs)] # Just the forward read files
fnRs = fastqs[grepl("_R2.fastq.gz", fastqs)] # Just the reverse read files
# Get sample names, assuming files named as so: SAMPLENAME_XXX.fastq
fnsFiles = fnsFiles[grepl(".fastq.gz$", fnsFiles)]
fnsFiles = sort(fnsFiles)
sample.names = sapply(strsplit(fnsFiles, "_R"), `[`, 1)
sample.names = sample.names[!duplicated(sample.names)]

#### {r filter and trim}
# Make directory and filenames for the filtered fastqs
filt_path = file.path(outdir, "filtered_fastq")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs = file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs = file.path(filt_path, paste0(sample.names, "_R_filt.fastq.gz"))

# Filter
for(i in seq_along(file.path(fastqDir, fnFs))) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(truncLen_L, truncLen_R), 
#                    trimleft=c(5, 5),
                    maxEE=c(maxEE_L, maxEE_R),
                    minQ=c(minQ_L, minQ_R),
                    compress=TRUE,
                    verbose=TRUE)
}

#### {r dereplication}
derepFs = derepFastq(filtFs, verbose=TRUE)
derepRs = derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) = sample.names
names(derepRs) = sample.names

#### {r learn error rates}
dadaFs.lrn = dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE, pool=TRUE)
errF = dadaFs.lrn[[1]]$err_out
dadaRs.lrn = dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE, pool=TRUE)
errR = dadaRs.lrn[[1]]$err_out
pdf(file.path(outdir, paste0(outdir,"_errorPlots.pdf")), w=8, h=11)
plotErrors(dadaFs.lrn[[1]], nominalQ=TRUE)
plotErrors(dadaRs.lrn[[1]], nominalQ=TRUE)
temp=dev.off()

#### {r sequence inference}
dadaFs = dada(derepFs, err=errF, multithread=TRUE, pool=TRUE)
dadaFs[[1]]
dadaRs = dada(derepRs, err=errR, multithread=TRUE, pool=TRUE)
dadaRs[[1]]
# save workspace image
save.image(file=file.path(outdir, paste0(projectName, "_workspace.Rdata")))

#### {r merge paired reads}
if (Vregion == "v1v3")
{
  mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, justConcatenate = TRUE, verbose=TRUE)
} else
{
  mergers = mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
}
# Inspect the merger data.frame from the first sample
# head(mergers[[1]])
save(mergers, file=file.path(outdir, paste0(projectName, "_mergers.Rdata")))

#### {r Construct the sequence table}
seqtab = makeSequenceTable(mergers)
# dim(seqtab)
# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

#### {r remove chimera}
seqtab.nochim = removeBimeraDenovo(seqtab, method = "pooled", verbose=TRUE)
#dim(seqtab.nochim)
#sum(seqtab.nochim)/sum(seqtab)
save(seqtab.nochim, file=file.path(projectName, paste0(outdir, "_seqtab.nochim.Rdata")))

#### {r get stats for chimera removal}
before = rowSums(seqtab)
after = rowSums(seqtab.nochim)
chimStat = cbind(before, after)
write("chimeraStat starts here:", stdout())
chimStat
write("chimeraStat ends here.", stdout())

#### {r assign taxonomy}
taxa = assignTaxonomy(seqtab.nochim,
                      trainingSet,
                      outputBootstraps=TRUE,
                      minBoot=minBoot_val,
                      tryRC = TRUE,
                      taxLevels=c("Kingdom", "Phylum", 
                                        "Class", "Order", "Family", 
                                        "Genus", "SuperSpecies", 
                                        "Species"))

# merge "SuperSpecies" level with "Species"
taxa = eightTOseven(taxa)

# merge tax with counts
counts = as.data.frame(t(seqtab.nochim))
#countTable = cbind(taxa$tax, counts)
countTable = as.data.frame(cbind(taxa$tax, taxa$boot, counts), stringAsFacter = FALSE)
names(countTable)[8:14] = paste0(names(countTable)[8:14], ".B")
countTable$Total = rowSums(counts)

write.table(countTable, file=file.path(outdir, paste0(projectName,"_countTable_minBoot",minBoot_val,".txt")), 
            sep="\t", quote=FALSE, row.names=TRUE, col.names = NA)


total_readCount = sum(rowSums(counts))
kingdom_identified = sum(rowSums(counts[!(is.na(countTable$Kingdom)),]))
phylum_identified = sum(rowSums(counts[!(is.na(countTable$Phylum)),]))
class_identified = sum(rowSums(counts[!(is.na(countTable$Class)),]))
order_identified = sum(rowSums(counts[!(is.na(countTable$Order)),]))
family_identified = sum(rowSums(counts[!(is.na(countTable$Family)),]))
genus_identified = sum(rowSums(counts[!(is.na(countTable$Genus)),]))
species_identified = sum(rowSums(counts[!(is.na(countTable$Species)),]))

percent_kingdom_identified = round((kingdom_identified/total_readCount)*100,2)
percent_phylum_identified = round((phylum_identified/total_readCount)*100,2)
percent_class_identified = round((class_identified/total_readCount)*100,2)
percent_order_identified = round((order_identified/total_readCount)*100,2)
percent_family_identified = round((family_identified/total_readCount)*100,2)
percent_genus_identified = round((genus_identified/total_readCount)*100,2)
percent_species_identified = round((species_identified/total_readCount)*100,2)

classiStat = data.frame(total_readCount,
                        kingdom_identified,
                        phylum_identified,
                        class_identified,
                        order_identified,
                        family_identified,
                        genus_identified,
                        species_identified,
                        percent_kingdom_identified,
                        percent_phylum_identified,
                        percent_class_identified,
                        percent_order_identified,
                        percent_family_identified,
                        percent_genus_identified,
                        percent_species_identified)
classiStat= t(classiStat)
classiStat = data.frame(Stat = rownames(classiStat), Value = classiStat[,1])
classiStat[1:8,2]= as.integer(classiStat[1:8,2])

write.table(classiStat, file=file.path(outdir, paste0(projectName,"_classifyStat_minBoot", minBoot_val, ".txt")), 
            sep="\t", quote=FALSE, row.names=FALSE)

message(paste0("Analysis finished at ", Sys.time()))

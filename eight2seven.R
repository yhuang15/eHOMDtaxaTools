######################################################################################
####
#### Description: This R function is needed to format the output object of the DADA2 
####              package's "assignTaxonomy()" function when the eHOMD RDP classifier
####              training dataset is used for assigning Taxonomy.  Without using this
####              function, the object produced by "assignTaxonomy()" using eHOMD
####              training dataset contains 8 levels of taxonomy, i.e. there is an 
####              additional "SuperSpecies" level above the "Species" level.  This
####              function merge "SuperSpecies" with "Species", replacing cases where
####              "NA"s are assigned in "Species" with the assignment and corresponding
####              bootstrap value at "SuperSpecies" level.
####
#### Usage:       eightTOseven(taxa)
####
#### Arguments:   taxa  (Required). The output object of the "assignTaxonomy()"
####                    function when the eHOMD training dataset was used to train
####                    the RDP classifier. Note that taxLevels needs to be set as
####                    taxLevels=c("Kingdom", "Phylum", "Class", "Order", "Family", 
####                    "Genus", "SuperSpecies", "Species") in "assignTaxonomy()"
####               
#### Author:      Yanmei Huang
#### Version:     1.1
#### Date:        2017-11-30
####
######################################################################################

eightTOseven = function(taxa)
{
  if (length(names(classified)) > 0)
  {
    taxa$boot[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$boot[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$tax[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
    taxa$tax = cbind(taxa$tax[, 1:6], taxa$tax[, 8])
    taxa$boot = cbind(taxa$boot[, 1:6], taxa$boot[, 8])
  } else {
    taxa[is.na(taxa[, 'Species']), 'Species'] = taxa[is.na(taxa[, 'Species']), 'SuperSpecies']
    taxa = cbind(taxa[, 1:6], taxa[, 8])
  }
  return (taxa)
}

######################################################################################
####
#### Description: This R function is needed to format the output object of the DADA2 
####              package's "assignTaxonomy" function when the eHOMD RDP classifier
####              training dataset is used for assigning Taxonomy.  Without using this
####              function, the object produced by "assignTaxonomy" using eHOMD
####              training dataset contains 8 levels of taxonomy, i.e. there is an 
####              additional "SuperSpecies" level above the "Species" level.  This
####              function merge "SuperSpecies" with "Species", replacing cases where
####              "NA"s are assigned in "Species" with the assignment and corresponding
####              bootstrap value at "SuperSpecies" level.
####              
####               
#### Author:      Yanmei Huang
#### Version:     1.0
#### Date:        2017-11-29
####
######################################################################################

eightTOseven = function(taxa)
{
  taxa$boot[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$boot[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
  taxa$tax[is.na(taxa$tax[, 'Species']), 'Species'] = taxa$tax[is.na(taxa$tax[, 'Species']), 'SuperSpecies']
  taxa$tax = cbind(taxa$tax[, 1:6], taxa$tax[, 8])
  taxa$boot = cbind(taxa$boot[, 1:6], taxa$boot[, 8])
  return (taxa)
}


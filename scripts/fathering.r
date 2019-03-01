library(sjPlot)

# Load data
datapath <- '~/git/amplexus_success/input'
outpath <- '~/git/amplexus_success/output/fathering'

dir.create(outpath)

father_data <- read.table(file.path(datapath, 'fathering_data.txt'), header=T)

str(father_data)

### All data summarised by genotype
summary(father_data$geno)

# Define Y differentiation 
## define by whole chromosome haplotype
father_data$genoDiff

father_data$genodiff <- as.character(father_data$geno)
father_data$genodiff[father_data$geno=="YB1a" | father_data$geno=="YB2a" | father_data$geno=="YB2b" | father_data$geno=="YA1a"] <- "diff"
father_data$genodiff[father_data$geno=="YB10" | father_data$geno=="YB20" | father_data$geno=="YB30" | father_data$geno=="YB40" | father_data$geno=="YB50"] <- "semidiff"
father_data$genodiff[father_data$geno=="XX"] <- "undiff"
father_data$genodiff <- as.factor(father_data$genodiff)
summary(father_data$genodiff)

### define by major dmrt regions
father_data$genodmrt <- as.character(father_data$geno)
father_data$genodmrt[father_data$geno=="YB1a" | father_data$geno=="YB10"] <- "YB1"
father_data$genodmrt[father_data$geno=="YB2a" | father_data$geno=="YB2b" | father_data$geno=="YB20"] <- "YB2"
father_data$genodmrt[father_data$geno=="YB30" | father_data$geno=="YB40" | father_data$geno=="YB50"] <- "YBother"
father_data$genodmrt <- as.factor(father_data$genodmrt)
summary(father_data$genodmrt)
father_data$amplexus <- factor(father_data$amplexus)

# Tables of fathering success compared to all swimming frogs in the population by haplotype differentiation or Dmrt allele

sjt.xtab(father_data$genodiff, father_data$mating, title='Fathering success (differentiation) All years', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genotype_year_fathering_genodiff_Mating.html'))

sjt.xtab(father_data$genodmrt, father_data$mating, title='Fathering success (dmrt type) All years', show.row.prc=T, show.legend=T, show.summary=T, file=file.path(outpath, 'contingency_genotype_year_fathering_genodmrt_Mating.html'))


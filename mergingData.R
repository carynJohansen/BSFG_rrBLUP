#purpose of this scriot is to load the expression data, the phenotype data, and the genotype data
# and to merge the data based on one list of maize lines.

setwd("~/Dropbox/Projects/BSFG_rrBLUP/")

#phenotype data
phenotype.raw <- read.table("~/Dropbox/UCD/rotations/RRI/maize_genomics_project/data/processed/WIDIV_2010_paper_phenotypic_data_cleaned.csv", header = T, sep=",")

pheno_lines <- as.character(phenotype[phenotype.raw$Year == "2008" & phenotype.raw$Rep == 1,]$Entry)

#imputed phenotype data (I think I need to re-do this. I can replicate what Rongkui did)
phenotype.imputed <- read.table("~/Dropbox/UCD/rotations/RRI/maize_genomics_project/data/traits/trait_impute.csv", header =T, sep=",")

#expression data
expression.raw <- read.table("~/Dropbox/UCD/rotations/RRI/maize_genomics_project/data/processed/maize_503genotypes_raw_expression_values_FPKM_cleaned.txt", header = T, sep="\t")
  
#master list of names, correcting for differences in spelling and use of marks
namesMaster <- read.table("~/Dropbox/UCD/rotations/RRI/maize_genomics_project/maizeLines_meta.csv", header=T, sep=",", na.strings = "")
namesMaster$Phenotype.Lines <- as.character(namesMaster$Phenotype.Lines)
namesMaster$Expression.Lines <- as.character(namesMaster$Expression.Lines)
namesMaster$MetaName <- as.character(namesMaster$MetaName)
namesMaster$FullName <- as.character(namesMaster$FullName)
namesMaster$GBS_NAME <- as.character(namesMaster$GBS_NAME)

#maize lines with Genotype, Phenotype, and Expression Data
lines_inAll <- namesMaster[complete.cases(namesMaster),]
lines_inAll$Phenotype.Lines <- as.character(lines_inAll$Phenotype.Lines)
lines_inAll$Expression.Lines <- as.character(lines_inAll$Expression.Lines)
lines_inAll$MetaName <- as.character(lines_inAll$MetaName)
lines_inAll$FullName <- as.character(lines_inAll$FullName)
lines_inAll$GBS_NAME <- as.character(lines_inAll$GBS_NAME)

#subset of imputed phenotypes with common maize liens
pheno_imputed_common <- phenotype.imputed[phenotype.imputed$Entry %in% lines_inAll$Phenotype.Lines,]
pheno_imputed_common$Entry <- as.character(pheno_imputed_common$Entry)
pheno_imputed_common <- pheno_imputed_common[order(pheno_imputed_common$Entry),]

#replace the entry column with the re-formated maize lines
pheno_Entry <- as.character(pheno_imputed_common$Entry)
pheno_entry_formated <- lines_inAll[lines_inAll$Phenotype.Lines %in% pheno_Entry,]$MetaName
pheno_entry_formated <- sort(pheno_entry_formated)
compare_formating <- data.frame(cbind(pheno_Entry, pheno_entry_formated))
pheno_imputed_common$Entry <- pheno_entry_formated

#subset of Expression data for just the common maize lines.
exp_lines <- colnames(expression.raw)[5:507]
exp_lines <- sort(exp_lines, decreasing = F)
exp_lines_common <- exp_lines[lines_inAll$Phenotype.Lines %in% exp_lines]
exp_lines_common_formated <- lines_inAll[lines_inAll$Expression.Lines %in% exp_lines_common,]$MetaName
exp_lines_common_formated <- sort(exp_lines_common_formated)

#lines in both expression and in phenotype data...
pheno_exp_lines <- pheno_entry_formated[pheno_entry_formated %in% exp_lines_common_formated]

#subset imputed data
pheno_imputed_subset <- pheno_imputed_common[pheno_imputed_common$Entry %in% pheno_exp_lines,]
write.table(pheno_imputed_subset, "Dropbox/UCD/rotations/RRI/maize_genomics_project/data/processed/imputedPhenotypes_commonLines.txt", sep="\t", quote = F, row.names = F)


#subset expression data
exp_info <- expression.raw[,1:4]
exp <- expression.raw[,5:507]
exp <- exp[,order(colnames(exp))]
names_exp <- namesMaster[complete.cases(namesMaster$Expression.Lines),]
exp_names_formated <- names_exp$MetaName
exp2 <- exp
colnames(exp2) <- exp_names_formated

#subset expression data with the pheno_exp_lines vector
exp_commonLines <- exp2[,colnames(exp2) %in% pheno_exp_lines]

expression_commonLines <- as.data.frame(cbind(exp_info, exp_commonLines))

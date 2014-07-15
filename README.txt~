


#### Here I store comments on what is what



#########expression_data folder
one file for each "condition", which can be a cell type, a time point, an activation
for format of each file must be
expression_[condition].RData
where [condition] will vary depending on the data.
Each of these RData file must contain 2 objects:
[condition] and support.[condition]
Condition is a matrix of numbers, with row.names that match the ProbeID column of the support file and colnames that match the different individuals in the dataset.


######### genotypes folder
One file per chromosome, stored in snpStats format, called
chr1
chr2
...
chr22
chrX
chrY
(no RData suffix, for no good reason)


The name of the object must be genotypes in each file. 
genotypes is a list that contains:
Mandatory: genotypes$genotypes which is a snpStats object, snps ideally in rsid format
Mandatory: genotypes$map file, with columns: SNP, allele.1, allele.2, position
Optional: a ped file, genotypes$fam

######### phenotypes folder
Two potential files:
binary.RData contain a single object named binary.pheno
continuous.RData contains a single object named continuous.pheno

numeric matrix, columns are measurements and rows are sample IDs that MUST match the genotypes
Don't forget the row and column names.

######## covariates folder
basic is a file called
covariates.tab
tab delimited, must contain the colum "id" that matches the sample names





########
Summary stats file, choice of headers:
SNPID, CHR, POS, then BETA.biom, SE.biom, PVAL.biom, N.biom in biom_chrXX.RData where biom is a code (ht), df
SNPID, CHR, POS, then BETA.cc, SE.cc, PVAL.cc, N.cc, Ncases.cc in cc_chrXX.RData where cc is a code (SCZ...).df
name of dataframes: biomarker or cc
automated generation of these files from : genotypes, phenotypes, and covariates




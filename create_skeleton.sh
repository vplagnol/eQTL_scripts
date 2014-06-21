code=$1

if [[ "$code" != "" ]]; then

    mkdir data/$code
    mkdir data/${code}/scripts
    mkdir data/${code}/genotypes
    mkdir data/${code}/phenotypes
    mkdir data/${code}/expression_data
    mkdir data/${code}/covariates
    mkdir data/${code}/support
    mkdir data/${code}/rawData
    mkdir data/${code}/summaryStats
    mkdir data/${code}/eQTLs
    mkdir data/${code}/eQTLs/matrixEQTL
    mkdir data/${code}/eQTLs/figs_matrixEQTL
    mkdir data/${code}/eQTLs/fgwas
    mkdir data/${code}/eQTLs/fgwas/fgwas_individual_files
    mkdir data/${code}/eQTLs/fgwas/fgwas_individual_figs
    mkdir data/${code}/figs
else
    echo "Need to specify a code"
fi




#### summary stats file: SNPID, CHR, POS, then BETA.biom, SE.biom, PVAL.biom, N.biom in biom_chrXX.RData where biom is a code (ht), df
#### summary stats file: SNPID, CHR, POS, then BETA.cc, SE.cc, PVAL.cc, N.cc, Ncases.cc in cc_chrXX.RData where cc is a code (SCZ...).df
##name of dataframes: biomarker or cc
##automated from : genotypes, phenotypes, and covariates



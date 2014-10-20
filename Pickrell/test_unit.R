
source('scripts/Pickrell/create_Pickrell_input.R')

prepare <- FALSE

if (prepare) {
  load('data/liver_Schadt/genotypes/chr1')
  load('data/liver_Schadt/expression_data/expression_Liver.RData')
}

dataset <- 'liver_Schadt'
ProbeID <- '10025913027'
snp.name <- 'rs12143028'

condition <- 'Liver'
expression <- Liver

print(genotypes$map[ 'rs12143028',])

te <- create.Pickrell.input.file (dataset = dataset,
                                  condition = condition,
                                  genotypes = genotypes , expression = expression,
                                  ProbeID = probeID, chromosome = '1', snp.name = snp.name, min.MAF = 0.03, min.certain.calls = 0.1,
                                  base.folder = '/cluster/project8/vyp/eQTL_integration')

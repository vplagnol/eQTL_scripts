source('scripts/eQTLs_scripts/find_all_eQTLs.R')

dataset <- 'WB_dexamethasone_DiRienzo'
condition <- 'logFC'

test <- run.eQTL ( dataset = dataset, condition = condition, chromosome = 22, start = 1, end = 300*10^6, pvOutputThreshold = 1e-5, force = FALSE)

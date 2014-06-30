source('scripts/plotting_functions/plot_eQTLs_per_SNP.R')

LCL.WB.logFC <- list( LCL_dexamethasone_DiRienzo = c('logFC'), WB_dexamethasone_DiRienzo = c('logFC'))  ###key argument
LCL.logFC <- list( LCL_dexamethasone_DiRienzo = c('logFC'))  ###key argument
WB.logFC <- list( WB_dexamethasone_DiRienzo = c('logFC'))  ###key argument

DC.logFC <- list ( DC_MTB_Barreiro = c('logFC'))
monocytes.LPS <- list ( monocytes_Knight = c('normal', 'LPS2', 'LPS24'))
monocytes.LPS24logFC <- list ( monocytes_Knight = c('LPS24logFC'))


GTex.brain <- list (GTex = c('Brain'))
GTex.muscle <- list (GTex = c('Muscle'))


#LCL.logFC.3 <- list( LCL_dexamethasone_DiRienzo = c('untreated', 'treated', 'logFC'))
#plot.eQTL (choice.sets = LCL.logFC.3, snp.name = 'rs12280343', chromosome = '11', gene.names = 'BIRC3', output.pdf = 'BIRC3.pdf', free.scale = TRUE); stop()




#DC.WB.LCL.logFC <- list ( DC_MTB_Barreiro = c('logFC'), WB_dexamethasone_DiRienzo = c('logFC'), LCL_dexamethasone_DiRienzo = c('logFC'))



######## the trans eQTL module from Fairfax et al
#test <- plot.eQTL ( choice.sets = monocytes.LPS24logFC, snp.name = 'rs2275888', gene.names = c('IRF7', 'CCL8', 'EPSTI1', 'IRF9'), output.pdf = 'IFNB_trans_monocytes.pdf', chromosome = '9', ylim = c(-2, 4), hline = 0); stop()
#test <- plot.eQTL ( choice.sets = DC.logFC, snp.name = 'rs2275888', gene.names = c('IRF7', 'CCL8', 'EPSTI1', 'IRF9'), output.pdf = 'IFNB_trans_Barreiro.pdf', chromosome = '9', ylim = c(-2, 4), hline = 0); stop()

test <- plot.eQTL ( choice.sets = monocytes.LPS, snp.name = 'rs2275888', gene.names = 'IFNB1', output.pdf = 'LPS_IFNB1.pdf', chromosome = '9'); stop()






######################
#test <- plot.eQTL ( choice.sets = GTex.brain, snp.name = 'rs13178576', gene.names = 'ENSG00000227355', output.pdf = 'check_error.pdf', chromosome = '5')
test <- plot.eQTL ( choice.sets = GTex.muscle, snp.name = 'rs7101840', gene.names = 'LDHAP4', output.pdf = 'check_error.pdf', chromosome = '11')
test <- plot.eQTL ( choice.sets = GTex.muscle, snp.name = 'rs6649', gene.names = 'RNF26', output.pdf = 'check_real.pdf', chromosome = '11')



#plot.eQTL ( choice.sets = DC.logFC, snp.name = '22:27200942', gene.names = 'CRYBA4', output.pdf = 'CRYBA4.pdf', chromosome = '22')



#plot.eQTL (choice.sets = LCL.WB.logFC, snp.name = 'rs812956', chromosome = '5', gene = 'SCGB3A1', output.pdf = 'data/WB_dexamethasone_DiRienzo/figs/SCGB3A1_steroids.pdf')  ##WB only but not LCL

#plot.eQTL (choice.sets = LCL.WB.logFC, snp.name = 'rs11211317', chromosome = '1', gene.names = 'MKNK1', output.pdf = 'data/WB_dexamethasone_DiRienzo/figs/MKNK1_steroids.pdf')   ## also WB only

#plot.eQTL (choice.sets = LCL.logFC, snp.name = 'rs11211317', chromosome = '1', gene.names = c('MKNK1', 'IL2RA'), output.pdf = 'test_multi_genes.pdf')   ## also WB only


#plot.eQTL (choice.sets = LCL.logFC, snp.name = 'rs79894583', chromosome = '2', gene.names = 'PRUNE2', output.pdf = 'data/LCL_dexamethasone_DiRienzo/figs/PRUNE2_transeQTL.pdf')   ## also WB only

#plot.eQTL (choice.sets = LCL.logFC, snp.name = 'rs113220072', chromosome = '11', gene.names = 'HS.126513', output.pdf = 'outlier.pdf')

#plot.eQTL (choice.sets = LCL.logFC, snp.name = 'rs4571654', chromosome = '7', gene.names = 'TRIB1', output.pdf = 'data/LCL_dexamethasone_DiRienzo/figs/TRIB1_transeQTL.pdf')   ## also WB only




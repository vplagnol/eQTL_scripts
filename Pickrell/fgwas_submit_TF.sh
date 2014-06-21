#step1=/cluster/project8/vyp/vincent/Software/pipeline/fgwas/fgwas_step1.sh
#sh $step1 --listBinAnnot $binPheno --inputFrame data/eQTLs_no_annotations.tab --outputFrame $annotationFile  --TF TRUE --ENCODE FALSE   ##should be part of the initial R scripts, so no need to call that


step2=/cluster/project8/vyp/vincent/Software/pipeline/fgwas/fgwas_step2.sh
step3=/cluster/project8/vyp/vincent/Software/pipeline/fgwas/fgwas_step3.sh

#root=data/monocytes_Knight/eQTLs/fgwas/fgwas_LPS2
root=data/WB_dexamethasone_DiRienzo/eQTLs/fgwas/fgwas_logFC
#root=data/brain_UKBEC/eQTLs/fgwas/fgwas_lncRNA_cerebellar_cortex
#root=data/DC_MTB_Barreiro/eQTLs/fgwas/fgwas_logFC
#root=data/WB_dexamethasone_DiRienzo/eQTLs/fgwas/fgwas_logFC
#root=data/LCL_dexamethasone_DiRienzo/eQTLs/fgwas/fgwas_logFC

binPheno=${root}_list_bin_annot.tab
annotationFile=${root}_fine_mapping_input_with_annotation.tab.gz
outputFolder=${root}_results
output=${outputFolder}/fgwas

if [ ! -e $outputFolder ]; then mkdir $outputFolder; fi

for file in $binPheno $annotationFile; do
    if [ ! -e $file ]; then echo "File $file does not exist"; exit; fi
done


script=cluster/submission/fgwasTF.sh
#sh $step2 --listPhenotypes_bin $binPheno --annotationFile $annotationFile --output $output --script $script
sh $step3 --listPhenotypes_bin $binPheno --output $output 


echo $script

#script=cluster/submission/fgwasTF_RELA.sh
#### with RELA as baseline now
#sh $step2 --listPhenotypes_bin $binPheno --annotationFile $annotationFile --output ${output}_bRELA --script $script --baseline RELA
#sh $step3 --listPhenotypes_bin $binPheno --output ${output}_bRELA 

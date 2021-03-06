######## important to set this before submission, maybe link here via a bash script from project folder.
project="WHII"
########
shopt -s expand_aliases
source ~/.bashrc
export PATH=${PATH}:/share/apps/R-3.0.2/bin
alias runR='sh /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/runRonCluster.sh'
alias R=/share/apps/R-3.0.2/bin/R
output_directory="/cluster/project8/vyp/eQTL_integration/data/"$project"/summaryStats/submissionScripts"
rm -r $output_directory
mkdir $output_directory
cd $output_directory
for chr in {1..22}
do
scriptname="script_"$chr"_sumstat.R"
cp /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/scripts/sumstattemplate.R $scriptname
rOutputFileName="biom_chr"$chr".Rdata"
rInput='oFile <- '"'$rOutputFileName'"';chr <- "'$chr'"; project='"'$project'"''
echo $rInput | cat - $scriptname > temp && mv temp $scriptname
runR $scriptname
done


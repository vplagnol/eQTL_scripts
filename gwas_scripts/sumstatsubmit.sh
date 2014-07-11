shopt -s expand_aliases
source ~/.bashrc
export PATH=${PATH}:/share/apps/R-3.0.2/bin
alias runR='sh /cluster/project8/vyp/eQTL_integration/summaryStats/scripts/runRonCluster.sh'
alias R=/share/apps/R-3.0.2/bin/R
output_directory="/cluster/project8/vyp/eQTL_integration/data/WHII/summaryStats/submissionScripts"
rm -r $output_directory
mkdir $output_directory
cd $output_directory
for chr in {1..22}
do
scriptname="script_"$chr"_sumstat.R"
cp /cluster/project8/vyp/eQTL_integration/summaryStats/scripts/sumstattemplate.R $scriptname
rOutputFileName="biom_chr"$chr".Rdata"
rInput='oFile <- '"'$rOutputFileName'"';chr <- "'$chr'"'
echo $rInput | cat - $scriptname > temp && mv temp $scriptname
runR $scriptname
done


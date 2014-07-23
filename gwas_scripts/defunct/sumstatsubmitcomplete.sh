######## important to set this before submission, maybe link here via a bash script from project folder.
project="WHII"
########
shopt -s expand_aliases
source ~/.bashrc
export PATH=${PATH}:/share/apps/R-3.0.2/bin
alias R=/share/apps/R-3.0.2/bin/R
output_directory="/cluster/project8/jonathan/testscripts/"$project"/summaryStats/submissionScripts"
rm -r $output_directory
mkdir $output_directory
cd $output_directory
for chr in {1..22}
do
Rscriptname="script_"$chr"_sumstat.R"
scriptname="script_"$chr"_sumstat.sh"
cp /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/scripts/sumstattemplate.R $Rscriptname
rOutputFileName="biom_chr"$chr".Rdata"
rInput='oFile <- '"'$rOutputFileName'"';chr <- "'$chr'"; project='"'$project'"''
echo $rInput | cat - $Rscriptname > temp && mv temp $Rscriptname

f=$scriptname
y=${f%R}
scriptname=$y"sh"
cp '
#$ -S /bin/sh
#$ -l h_vmem=8G
#$ -l tmem=8G
#$ -l h_rt=24:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/cluster/output/
#$ -e /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/cluster/error/ $script' $scriptname
echo R CMD BATCH --no-save $f >> $scriptname
echo "Running" $f "on cluster as" $scriptname
#qsub $scriptname
done


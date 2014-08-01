######## important to set this before submission
project="MRC1946"
repoFolder=/cluster/project8/vyp/jonathan/git/eQTL_integration/
outputfolder=/cluster/project8/vyp/jonathan/cluster/
########
shopt -s expand_aliases
source ~/.bashrc
export PATH=${PATH}:/share/apps/R-3.0.2/bin
alias R=/share/apps/R-3.0.2/bin/R
script_output_directory=$outputfolder'submission/summaryStats'$project
rm -r $script_output_directory
mkdir $script_output_directory
cd $script_output_directory
for chr in {1..22}
do
Rscriptname="script_"$chr"_sumstat.R"
scriptname="script_"$chr"_sumstat.sh"
cp $repofolder/gwas_scripts/sumstattemplate.R $Rscriptname
rOutputFileName="biom_chr"$chr".RData"
rInput='oFile <- '"'$rOutputFileName'"';chr <- "'$chr'"; project='"'$project'"''
echo $rInput | cat - $Rscriptname > temp && mv temp $Rscriptname
f=$scriptname
y=${f%R}
scriptname=$y
echo '
#$ -S /bin/sh
#$ -l h_vmem=20G
#$ -l tmem=20G
#$ -l h_rt=24:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o '$output_folder'output/summaryStats'$project'
#$ -e '$output_folder'error/summaryStats'$project > $scriptname
echo R CMD BATCH --no-save $f >> $scriptname
echo "Running" $f "on cluster as" $scriptname
qsub $scriptname
done


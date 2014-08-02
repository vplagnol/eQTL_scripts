######## important to set this before submission
project="WHII"
repoFolder='/cluster/project8/vyp/jonathan/git/eQTL_scripts/'
outputfolder='/cluster/project8/vyp/jonathan/cluster/'
########
shopt -s expand_aliases
source ~/.bashrc
export PATH=${PATH}:/share/apps/R-3.0.2/bin
alias R=/share/apps/R-3.0.2/bin/R
script_output_directory=$outputfolder"submission/sumstathits"$project
rm -r $script_output_directory
mkdir $script_output_directory
cd $script_output_directory
for chr in {1..22}
do
Rscriptname="script_"$chr"_hits.R"
scriptname="script_"$chr"_hits.sh"
cp $repoFolder"gwas_scripts/get_hits.R" $Rscriptname
rOutputFileName="hits"$chr".RData"
rInput='oFile <- '"'$rOutputFileName'"';chr <- "'$chr'"; project='"'$project'"''
echo $rInput | cat - $Rscriptname > temp && mv temp $Rscriptname
f=$Rscriptname
echo '
#$ -S /bin/sh
#$ -l h_vmem=3G
#$ -l tmem=3G
#$ -l h_rt=1:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o '$outputfolder'output/summaryStatHits'$project'
#$ -e '$outputfolder'error/summaryStatHits'$project > $scriptname
echo R CMD BATCH --no-save $f >> $scriptname
echo "Running" $Rscriptname "on cluster as" $scriptname
qsub $scriptname
done


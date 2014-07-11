
#$ -S /bin/sh
#$ -l h_vmem=20G
#$ -l tmem=20G
#$ -l h_rt=48:00:0
#$ -V
#$ -R y
#$ -pe smp 1
#$ -cwd
#$ -o /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/cluster/output/
#$ -e /cluster/project8/vyp/eQTL_integration/scripts/gwas/summaryStats/cluster/error/

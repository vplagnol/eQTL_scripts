f=$1
y=${f%R}
script=$y"sh"
cp /cluster/project8/vyp/eQTL_integration/summaryStats/runRonCluster_input.sh $script
echo R CMD BATCH --no-save $f >> $script
echo "Running" $f "on cluster as" $script
qsub $script

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
replicates=$METADATA/replicates_1.txt;
sample_info=$METADATA/sample_info_1.txt;
popmap=$METADATA/popmap_1.txt;
order=$METADATA/order_1.txt;

### Create single column list with individuals and populations to remove to remove 
nano $METADATA/remove_inds.1a.txt;
nano $METADATA/remove_pops.1.txt;

### Create new files, excluding listed individuals
cat $METADATA/remove_inds.1a.txt | cut -f 1 \
> $METADATA/remove_inds.1b.txt;
remove=$METADATA/remove_inds.1b.txt;
cp $replicates $METADATA/replicates_2.txt;
cp $sample_info $METADATA/sample_info_2.txt;
cp $popmap $METADATA/popmap_2.txt
replicates=$METADATA/replicates_2.txt;
sample_info=$METADATA/sample_info_2.txt;
popmap=$METADATA/popmap_2.txt;
patterns=$(tr '\n' '|' < $remove);
patterns="${patterns::-1}"
sed -ri "/$patterns/d" $replicates;
sed -ri "/$patterns/d" $sample_info;
sed -ri "/$patterns/d" $popmap;

### Create new files, excluding listed populations
remove=$METADATA/remove_pops.1.txt;
cp $order $METADATA/order_2.txt;
order=$METADATA/order_2.txt;
patterns=$(tr '\n' '|' < $remove);
patterns="${patterns::-1}";
sed -ri "/$patterns/d" $order;

#### For some reason 236 disappears here even if I haven't specified it. So I put it back manually

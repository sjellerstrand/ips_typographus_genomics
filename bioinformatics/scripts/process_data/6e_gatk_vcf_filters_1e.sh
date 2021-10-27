### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
sample_info=$METADATA/sample_info_2.txt;
popmap=$METADATA/popmap_2.txt;

### Create single column list with individuals to remove 
nano $METADATA/remove_inds.2a.txt;

### Create new files, excluding listed individuals
cat $METADATA/remove_inds.2a.txt | cut -f 1 \
> $METADATA/remove_inds.2b.txt;
remove=$METADATA/remove_inds.2b.txt;
cp $sample_info $METADATA/sample_info_all.txt;
cp $popmap $METADATA/popmap_all.txt
sample_info=$METADATA/sample_info_all.txt;
popmap=$METADATA/popmap_all.txt;
patterns=$(tr '\n' '|' < $remove);
patterns="${patterns::-1}";
sed -ri "/$patterns/d" $sample_info;
sed -ri "/$patterns/d" $popmap;

#### 236 might still need some checking if it is still there and with proper tab delimitation
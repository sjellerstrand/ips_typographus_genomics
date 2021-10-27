### Calulate genotype differences based on likelihoods
bcftools gtcheck $VCF_OUT.vcf.gz  | grep "ERR" > $VCF_OUT.err;

### Plot replicate comparisons
Rscript $MAINDIR/scripts/radseq2020/quality_control/replicate_comp.r \
--args MAINDIR=$MAINDIR OUTDIR=$OUTDIR VCF_OUT=$VCF_OUT \
sample_info=$sample_info replicates=$replicates;

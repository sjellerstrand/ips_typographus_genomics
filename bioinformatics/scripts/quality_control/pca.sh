### Do a PCA
## Prune variants
plink --vcf $VCF_OUT.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --indep-pairwise 50 10 0.1 --out $VCF_OUT;
## Calculate Pricipal components
plink --vcf $VCF_OUT.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# --extract $VCF_OUT.prune.in --pca 20 --out $VCF_OUT;

### Plot Pricipal component 1 and 2
Rscript $MAINDIR/scripts/radseq2020/quality_control/pca.r \
--args VCF_OUT=$VCF_OUT DATA=$DATA METADATA=$METADATA;

### Plot annotations for evaluation of GATK hard cut-off settings

## Retrieve statistics
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 13 | cut -d'=' -f2 > $VCF_OUT.QD;
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 7 | cut -d'=' -f2 > $VCF_OUT.FS;
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 11 | cut -d'=' -f2 > $VCF_OUT.MQ;
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 12 | cut -d'=' -f2 > $VCF_OUT.MQRankSum;
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 14 | cut -d'=' -f2 > $VCF_OUT.ReadPosRankSum;
bcftools view -H $VCF_OUT.vcf.gz | cut -f8 | cut -d';' -f 6 | cut -d'=' -f2 > $VCF_OUT.ExcessHet;
paste $VCF_OUT.QD $VCF_OUT.FS $VCF_OUT.MQ $VCF_OUT.MQRankSum $VCF_OUT.ReadPosRankSum \
$VCF_OUT.ExcessHet > $VCF_OUT.gatk.annotations.temp;
echo -e "QD\tFS\tMQ\tMQRankSum\tReadPosRankSum\tExcessHet" | \
cat - $VCF_OUT.gatk.annotations.temp > $VCF_OUT.gatk.annotations;
rm $VCF_OUT.gatk.annotations.temp;

## Plot statistics
Rscript $MAINDIR/scripts/radseq2020/quality_control/gatk_filters.r \
--args VCF_OUT=$VCF_OUT;

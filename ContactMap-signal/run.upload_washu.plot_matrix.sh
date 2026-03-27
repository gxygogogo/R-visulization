#####################################################################################################################################################################
# for genotype in WT SA1KO SA2KO
# do
# DIR=/data/tang/cohesin_project/ChIA_PET_FINAL/Final.for.Downstream.Analysis/feature_analysis_of_signal_stripe_denseregion/composition.of.stripes
# cat ${DIR}/${genotype}.stripe.for.visualization/*sliding_anchor_CTCF.bedpe | awk 'BEGIN{OFS="\t";}{print $1,$2+1,$3,$4,$5+1,$6,"CTCF";}' > tmp.txt
# cat ${DIR}/${genotype}.stripe.for.visualization/*sliding_anchor_CTCF_NIPBL.bedpe | awk 'BEGIN{OFS="\t";}{print $1,$2+1,$3,$4,$5+1,$6,"CTCF_NIPBL";}' >> tmp.txt
# cat ${DIR}/${genotype}.stripe.for.visualization/*sliding_anchor_NIPBL.bedpe | awk 'BEGIN{OFS="\t";}{print $1,$2+1,$3,$4,$5+1,$6,"NIPBL";}' >> tmp.txt
# cat ${DIR}/${genotype}.stripe.for.visualization/*sliding_anchor_None.bedpe | awk 'BEGIN{OFS="\t";}{print $1,$2+1,$3,$4,$5+1,$6,"None";}' >> tmp.txt
# cat tmp.txt | sort -k1,1 -k2,2n -k5,5n | uniq > ${genotype}.stripe.signal_type.bedpe
# done

#####################################################################################################################################################################
chr=chr5
start=73735988
end=74945987

mkdir ${chr}_${start}_${end}
cd ${chr}_${start}_${end}
DIR=/data/tang/cohesin_project/ChIA_PET_FINAL/Final.for.Downstream.Analysis/run.call.signals/summary_of_signal/finalFile
for genotype in WT SA1KO SA2KO
do
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a ../categorized.signals.for.visualization/${genotype}.stripe.signals.bedpe -b - -type both | cut -f1,2,5,7 | uniq | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Stripe"}'> ${genotype}_SMC1A.Stripe.signal.${chr}_${start}_${end}.bedpe
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a ../categorized.signals.for.visualization/${genotype}.reseting.signals.bedpe -b - -type both | cut -f1,2,5,7 | uniq | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"Reseting"}'> ${genotype}_SMC1A.reseting.signal.${chr}_${start}_${end}.bedpe
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a ../categorized.signals.for.visualization/${genotype}.CTCF_surrounding.signals.bedpe -b - -type both | cut -f1,2,5,7 | uniq | awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"CTCF_surrounding"}'> ${genotype}_SMC1A.CTCF_surrounding.signal.${chr}_${start}_${end}.bedpe
cat ${genotype}_SMC1A.Stripe.signal.${chr}_${start}_${end}.bedpe ${genotype}_SMC1A.reseting.signal.${chr}_${start}_${end}.bedpe ${genotype}_SMC1A.CTCF_surrounding.signal.${chr}_${start}_${end}.bedpe > ${genotype}_SMC1A.signal.${chr}_${start}_${end}.bedpe
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a ${DIR}/${genotype}_ChIAPET_SMC1A.iPET.5kbBin.rawMatrix.bedpe -b - -type both | cut -f1,2,5,7 | uniq > ${genotype}_SMC1A.rawMatrix.${chr}_${start}_${end}.bedpe
done


## 绘图
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr1_165013116_166618115 0.6
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr1_65593352_66613351 0.9
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr2_186751998_187598710 1.1
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr5_111442517_112497516 0.9
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr9_17771901_18356900 1.5
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr13_21518404_23963403 0.4
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr20_39176835_40009300 1.2
/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript run.plot.matrix.stripe_signal.R chr5_73735988_74945987 0.8


# reg <- "chr2:146975600-147987867"
# res <- 5000
# reg=gsub(",", "", reg)
# chr=strsplit(reg, split="[_:-]")[[1]][1]
# start=as.numeric(strsplit(reg, split="[_:-]")[[1]][2])
# end=as.numeric(strsplit(reg, split="[_:-]")[[1]][3])
# start = start - start%%res
# end = end + res - end%%res
# paste0(chr,":",start,"-",end+5000)

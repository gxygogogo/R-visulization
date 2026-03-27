#!/bin/bash

alias hicConvertFormat=/home/xiongdan/software/HiCExplorer-3.7.2/bin/hicConvertFormat
##### matrix #################################################################################
##### 10kb  ###############################################
##### rheMac8 ##################

path=/data2/yuhan/Project/3D_evolution/HiCDCPlus/macaque
res=5000
chr=18
r=$((${res}/1000))

name=LCL8664_ChIAPET_CTCF_merged_B1_2
name=RheNeu_d12_ChIAPET_CTCF_merged_B1_2
name=RheNeu_ChIAPET_RNAPII_merge_d15_B1_2_d18_B1_2_3
hicConvertFormat -m ${path}/${name}.qvalue.nb_hurdle.${r}kb.hic --inputFormat hic -o ${name}.chr${chr}.qvalue.nb_hurdle.${r}kb.cool --outputFormat cool --chromosome ${chr} --resolutions ${res}
hicConvertFormat -m ${name}.chr${chr}.qvalue.nb_hurdle.${r}kb_${res}.cool --inputFormat cool -o ${name}.chr${chr}.qvalue.nb_hurdle.${r}kb.h5 --outputFormat h5 --resolutions ${res}

## 将h5格式转换为interactions格式
hicConvertFormat -m /data3/xiongdan/plot_hicdc_out_matrix/DLGAP1/9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.h5 --inputFormat h5 --outputFormat ginteractions -o 9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions --resolutions 10000
hicConvertFormat -m /data3/xiongdan/plot_hicdc_out_matrix/DLGAP1/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.h5 --inputFormat h5 --outputFormat ginteractions -o RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions --resolutions 10000

## 提取相应区域信号
chr=18
start=3200000
end=5100000
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a 9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > 9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

chr=18
start=4400000
end=6250000
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

chr=18
start=4389691
end=6612342
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

chr=18
start=2875138
end=5128616
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a 9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > 9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

chr=18
start=4300000
end=6450000
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

chr=18
start=4400000
end=6290000
awk -v CHR=${chr} -v START=${start} -v ENDcor=${end} 'BEGIN{OFS="\t";print CHR,START,ENDcor;}' | bedtools pairtobed -a RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.tsv -b - -type both | cut -f1,2,5,7 | uniq > RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.${chr}_${start}_${end}.bedpe

## 使用绘图脚本进行热图绘制
# 参数1: 输入文件的路径
# 参数2: 绘制范围的标签
# 参数3: 绘制热图的分辨率
# 参数4: 热图展示的最大截断值
# 参数5: 热图展示的最小截断值
# 参数5: 输出热图的路径

Rscript=/public1/xinyu/Software/Miniconda3/envs/R4.2/bin/Rscript
plotHeatmap=/public1/xinyu/3D.ContactHeatmap/script/run.Contact_map_plot.R

${Rscript} ${plotHeatmap} "/public1/xinyu/3D.ContactHeatmap/9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.18_3200000_5100000.bedpe" \
                          "18_3200000_5100000" \
                          10000 \
                          25 \
                          "/public1/xinyu/3D.ContactHeatmap/9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.18_3200000_5100000.Clip25.pdf"

${Rscript} ${plotHeatmap} "/public1/xinyu/3D.ContactHeatmap/9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.18_2875138_5128616.bedpe" \
                          "18_2875138_5128616" \
                          10000 \
                          25 \
                          "/public1/xinyu/3D.ContactHeatmap/9_GW11B3_d20_ChIAPET_CTCF_merged_B1_2_3.matrix.10kb.interactions.18_2875138_5128616.new.pdf"

${Rscript} ${plotHeatmap} "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4400000_6250000.bedpe" \
                          "18_4400000_6250000" \
                          10000 \
                          10 \
                          "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4400000_6250000.Clip10.pdf"

${Rscript} ${plotHeatmap} "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4389691_6612342.bedpe" \
                          "18_4400000_6200000" \
                          10000 \
                          10 \
                          "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4389691_6612342.V3.pdf"

${Rscript} ${plotHeatmap} "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4400000_6290000.bedpe" \
                          "18_4400000_6290000" \
                          10000 \
                          12 \
                          2 \
                          "/public1/xinyu/3D.ContactHeatmap/RheNeu_d12_ChIAPET_CTCF_merged_B1_2.matrix.10kb.interactions.18_4400000_6290000.Clip12_2.pdf"

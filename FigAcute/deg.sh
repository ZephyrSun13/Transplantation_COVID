
set -e 

#cp ../FactorAll2Batch.txt Factor.xls

#Rscript statClin.r

#Rscript stat.r ../Expr.xls

#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/limmaVoomNum.r ../Expr.match.xls Factor.xls unpaired Group Score 0 0.05 Ensembl Group No voom 1

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/limmaVoomNum.r ./Expr.match.xls Factor.xls unpaired Group Score 0 0.05 Ensembl Group "Age;Sex;Batch;Lymph" voom 1

head -1 Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls > Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.selected.xls

awk '$5<=0.01' Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls >> Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.selected.xls

#awk -F "\t" '(($5 <= 0.01) && ($2 > 0)){printf("%s\t%s\n", $8, $2)}' Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls > PosGene

#awk -F "\t" '(($5 <= 0.01) && ($2 < 0)){printf("%s\t%s\n", $8, $2)}' Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls > NegGene

awk -F "\t" '($5 <= 0.01){printf("%s\t%s\n", $8, $2)}' Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls | sort | uniq > AllGene


#Rscript volcano_plot.r Condition_Case-Control_Case-Control_diff.ann.xls Diff

#ml R/4.0.3
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r PosGene Symbol
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r NegGene Symbol
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r AllGene Symbol
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r PosGene Symbol
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r NegGene Symbol
#
#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r AllGene Symbol


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichR.r AllGene

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/plotEnrichNP.r EnrichR.sig.xls Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls EnrichR.plot 

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/plotEnrichNP.r EnrichR.sig.xls EnrichR.plot


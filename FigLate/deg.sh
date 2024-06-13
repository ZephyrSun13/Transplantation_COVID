
set -e 

ml R/3.5.3

#cp ../FactorAll2Batch.txt Factor.xls

#Rscript stat.r ../Expr.xls

Rscript statClin.r

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/limmaVoom.r ../Expr.match.xls Factor.xls unpaired Group Chronic_Acute 0 0.01 Ensembl ChronicVSAcute "Age;Sex;Batch;Lymph;SeverityScore" voom 1

head -1 ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.xls > ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.selected.xls 

awk '$5<=0.01' ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.xls >> ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.selected.xls 

awk -F "\t" '(($5 <= 0.01) && ($2 > 0)){printf("%s\t%s\n", $8, $2)}' ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.xls > PosGene

awk -F "\t" '(($5 <= 0.01) && ($2 < 0)){printf("%s\t%s\n", $8, $2)}' ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.xls > NegGene

awk -F "\t" '($5 <= 0.01){printf("%s\t%s\n", $8, $2)}' ChronicVSAcute_Group_Chronic_Acute_Age_Sex_Batch_Lymph_SeverityScore_Chronic_Acute_diff.ann.xls | sort | uniq > AllGene


#Rscript volcano_plot.r Condition_Case-Control_Case-Control_diff.ann.xls Diff

ml R/4.0.3

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r PosGene Symbol

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r NegGene Symbol

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichGO.r AllGene Symbol

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r PosGene Symbol

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r NegGene Symbol

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichKEGG.r AllGene Symbol


Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/enrichR.r AllGene

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/plotEnrichNP.r EnrichR.sig.xls EnrichR.plot



head -1 /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/AcuteNum.lymph/Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls > CheckGene.xls

awk -F "\t" '$8=="TNF"' /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/AcuteNum.lymph/Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls >> CheckGene.xls

awk -F "\t" '$8=="IL1"' /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/AcuteNum.lymph/Group_Group_Score_Age_Sex_Batch_Lymph_Group_diff.ann.xls >> CheckGene.xls


awk -F "\t" '$8=="TNF"' /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/Public/GSE152418/Severe//Group_Group_Score_gender_Group_diff.ann.xls >> CheckGene.xls

awk -F "\t" '$8=="IL6"' /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/Public/GSE152418/Severe//Group_Group_Score_gender_Group_diff.ann.xls >> CheckGene.xls

awk -F "\t" '$8=="IL1A"' /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/Public/GSE152418/Severe//Group_Group_Score_gender_Group_diff.ann.xls


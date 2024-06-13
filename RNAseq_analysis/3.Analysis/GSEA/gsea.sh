
module load R/3.5.3

Rscript gsea_pre.r 

grep GOBP /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c5.all.v7.4.symbols.gmt.txt | grep NEUTROPHIL > GOBP_NEUTROPHIL.gmt

grep NEUTROPHIL /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt | grep REACTOME >> GOBP_NEUTROPHIL.gmt

grep KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt >> GOBP_NEUTROPHIL.gmt

grep Neutrophil_degranulation /sc/arion/projects/GOCAR/Sun/5.DRWork/8.COVID/3.AnaAll/All.lymph.SeveAdj/GSEA/GOBP_CELLCYCLE.gmt >> GOBP_NEUTROPHIL.gmt

#awk '$1~/TH1/' /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt > GOBP_NEUTROPHIL.gmt

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.r Ranks GOBP_NEUTROPHIL.gmt 0.1

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.plot.r fgseaRes.rds GSEA.plot Ranks GOBP_NEUTROPHIL.gmt

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.r Ranks /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.kegg.gmt 0.1



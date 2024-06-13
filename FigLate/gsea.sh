
module load R/3.5.3

Rscript gsea_pre.r 

grep GOBP /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c5.all.v7.4.symbols.gmt.txt | grep GOBP_T_CELL_RECEPTOR_SIGNALING_PATHWAY > GOBP_CELLCYCLE.gmt

#grep NEUTROPHIL /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt | grep REACTOME >> GOBP_NEUTROPHIL.gmt

grep KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt >> GOBP_CELLCYCLE.gmt

grep BIOCARTA_TH1TH2_PATHWAY /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt >> GOBP_CELLCYCLE.gmt

grep KEGG_T_CELL_RECEPTOR_SIGNALING_PATHWAY /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt >> GOBP_CELLCYCLE.gmt

grep GOBP /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c5.all.v7.4.symbols.gmt.txt | grep NEUTROPHIL >> GOBP_CELLCYCLE.gmt

grep NEUTROPHIL /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt | grep REACTOME >> GOBP_CELLCYCLE.gmt

cut -f 3 QuickGO_Neutrophil_degranulation | sed '1d' |  tr '[a-z]' '[A-Z]' | sort | uniq | awk 'BEGIN{printf("Neutrophil_degranulation\tNeutrophil_degranulation")}{if(NR>1){printf("\t%s", $1)}}END{printf("\n")}' >> GOBP_CELLCYCLE.gmt

#awk '$1~/TH1/' /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.gmt.txt > GOBP_NEUTROPHIL.gmt

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.r Ranks GOBP_CELLCYCLE.gmt 0.1

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.plot.r fgseaRes.rds GSEA.plot Ranks GOBP_CELLCYCLE.gmt

#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/fgsea.r Ranks /sc/arion/projects/zhangw09a/PANDA/db_ZS/GSEA/c2.all.v7.4.symbols.kegg.gmt 0.1



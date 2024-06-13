
ls ../2.Map3/Tmp/Map/*star/*count.txt | awk '{split($0, tt, "/"); gsub("_star", "", tt[5]); printf("\t%s", tt[5])} END {printf("\n")}' > head

ls ../2.Map3/Tmp/Map/*star/*count.txt | head -1 | awk '{system("cut -f 1 "$1)}' > Expr.xls

#cut -f 1 ../2.Map3/Tmp/Map/Monte-23_star/Monte-23_count.txt > Expr.xls 

for expr in ../2.Map3/Tmp/Map/*star/*count.txt
do

	cut -f 2 $expr > tt

	paste Expr.xls tt > ttt

	mv ttt Expr.xls

done

cp head Expr.xls.stat.xls

cat head Expr.xls | grep '__' >> Expr.xls.stat.xls

cat head Expr.xls | grep -v '__' > ttt

mv ttt Expr.xls

Rscript statMerge.r Expr.xls ../3.Ana2/Expr.xls 

Rscript stat.r Expr.xls

#awk -F "\t" '(NR > 1){printf "Sample_%sN_norm\t%s\t%s\tN\nSample_%sT_norm\t%s\t%s\tT\n",$1,$2,$3,$1,$2,$3}' Sample.info > Factor.xls

#perl -ne 'chomp; @tt = split("\t"); $ID = shift @tt; $content = join("\t", @tt); printf "Sample_%sN_norm\t%s\tN\nSample_%sT_norm\t%s\tT\n", $ID, $content, $ID, $content' Sample.info > Factor.xls

module load R/3.4.3

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/norm.r -e Expr.xls -l TRUE -q TRUE -c 100

Rscript pca.r

Rscript metaExpr.r

awk 'BEGIN{printf("\tGroup\n")} (NR>1){printf("%s\tCase\n", $2)}' ../../1.L1/sample.txt > Factor_meta.xls

Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/norm.r -e Expr.meta.xls -l TRUE -q TRUE -c 100 -b Factor_meta.xls

Rscript pca2.r

#Rscript norm.r -e Expr_L1FGGYOver.xls -l TRUE -q TRUE -c 100

#Rscript norm.r -e Expr_L1FGGYKnock.xls -l TRUE -q TRUE -c 100

#Rscript norm.r -e Expr_NVR_DMSO.xls -l TRUE -q TRUE -c 100


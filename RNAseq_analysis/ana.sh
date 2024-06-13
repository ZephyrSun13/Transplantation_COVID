
ls /sc/arion/scratch/sunz04/Luca/*star/*count.txt | awk '{split($0, tt, "/"); gsub("_star", "", tt[7]); printf("\t%s", tt[7])} END {printf("\n")}' > head

ls /sc/arion/scratch/sunz04/Luca/*star/*count.txt | head -1 | awk '{system("cut -f 1 "$1)}' > Expr.xls

#cut -f 1 ../2.Map3/Tmp/Map/Monte-23_star/Monte-23_count.txt > Expr.xls 

for expr in /sc/arion/scratch/sunz04/Luca/*star/*count.txt
do

	cut -f 2 $expr > tt

	paste Expr.xls tt > ttt

	mv ttt Expr.xls

done

cp head Expr.xls.stat.xls

cat head Expr.xls | grep '__' >> Expr.xls.stat.xls

cat head Expr.xls | grep -v '__' > ttt

mv ttt Expr.xls

Rscript pre.r

module load R/3.4.3

#Rscript /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/code/norm.r -e Expr.xls -l TRUE -q TRUE -c 100

Rscript pca.r Expr.corr.xls 1 "Cell;Stiffness"

Rscript metaExpr.r

Rscript idMap.r


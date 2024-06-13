
#ls Tmp/Map/*/Log.final.out | head -1 | awk '{system("cut -f 1 "$1)}' > map.stat

ls /sc/arion/scratch/sunz04/LRG1/*/*Log.final.out | head -1 | awk '{system("cut -f 1 "$1)}' > map.stat

for sam in /sc/arion/scratch/sunz04/LRG1/*/*Log.final.out
do

	key=${sam%/*}

	awk -F "\t" '{print $2}' $sam > tem.txt

	#echo $key >> tem.txt

        paste map.stat tem.txt > union.txt

        mv union.txt map.stat

done

rm -f tem.txt

rm -f union.txt

ls /sc/arion/scratch/sunz04/LRG1/*/*Log.final.out | awk '{split($1, tt, "/"); printf "\t%s", tt[7]} END{printf("\n")}' > head

ls /sc/arion/scratch/sunz04/LRG1/*/*Log.final.out | awk '{split($1, tt, "/"); printf "\t%s", tt[6]} END{printf("\n")}' >> head

cat head map.stat > tt && mv tt map.stat

module load R/3.5.3 && Rscript mapStat.r


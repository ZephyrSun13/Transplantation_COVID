
if [ -f InsertSize ]
then

	rm -f InsertSize

fi

module load samtools/1.1

for fd in $(cat Stars)
do

	samtools view $fd/Aligned.sortedByCoord.out.bam | awk -v fd="$fd" '{if ($9 > 0) {S+=$9; T+=1}}END{printf "%s\t%s\n", fd, S/T}' >> InsertSize 

done

head -10000 mappings.sam | awk 'M=YOUR_MEAN{if ($9 > 0) {S+=($9-M)*($9-M); T+=1}}END{print "StdDev: " sqrt(S/T)}'


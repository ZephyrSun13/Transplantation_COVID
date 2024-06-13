#!/usr/bin/bash

if [ -f StandardDeviation ]
then

	rm -f StandardDeviation

fi

filename="$1"

module load samtools/1.1

while IFS=$'\t' read -r -a line; do

    	echo "${line[0]}"
	echo "${line[1]}"

	samtools view ${line[0]}/Aligned.sortedByCoord.out.bam | awk -v M="${line[1]}" -v F="${line[0]}" '{if ($9 > 0) {S+=($9-M)*($9-M); T+=1}}END{printf("%s\t%s\n"), F, sqrt(S/T)}' >> StandardDeviation

done < "$filename"


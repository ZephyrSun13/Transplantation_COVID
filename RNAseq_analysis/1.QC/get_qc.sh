
#ls /sc/arion/projects/weijiaTemp/TMP/MadavCOVID3/Data/*.fastq.gz | awk '{if(NR%2==1){split($1, tt, "/"); gsub("_R1_001.fastq.gz", "", tt[9]);printf("%s\t%s", tt[9], $1)}else{printf("\t%s\n", $1)}}' > Sam

ls /sc/arion/projects/zhangw09a/Data/Zeguo_Sun/project/6.Transplant/8.LILR/7.Target/0.Data/Target/*.fastq.gz | awk '{if(NR%2==1){split($1, tt, "/"); gsub("_R1_001.fastq.gz", "", tt[14]);printf("%s\t%s", tt[14], $1)}else{printf("\t%s\n", $1)}}' > Sam

mkdir -p QC

rm -f qc.sh

mkdir -p Clean

DIR=/sc/arion/projects/zhangw09a/Data/Zeguo_Sun/project/6.Transplant/8.LILR/7.Target/1.QC

awk -v DIR="$DIR" '{

	printf("module load fastqc/0.11.8");

	printf(" && fastqc -o %s/QC/ %s %s", DIR, $2, $3);

	#printf("&& /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/SOAPnuke filter -l 10 -q 0.1 -n 0.01 -Q 2 -1 %s -2 %s -o %s/Clean/%s/ -C %s.clean.fq.gz -D %s.clean.fq.gz", DIR, $2, $3, $1, $2, $3);

	#printf(" && module load python/3.7.3 && cutadapt -m 50 -j 5 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA -o %s/%s.clean.R1.fq.gz -p %s/%s.clean.R2.fq.gz %s %s", DIR, $1, DIR, $1, $2, $3);

	#printf(" && module load fastqc/0.11.8 && fastqc -o %s/QC/ %s/Clean/%s.clean.R1.fq.gz %s/Clean/%s.clean.R2.fq.gz", DIR, DIR, $1, DIR, $1);

        #printf(" && fastqc -o /sc/arion/projects/weijiaTemp/TMP/MadavCOVID3/QC/ /sc/arion/projects/weijiaTemp/TMP/MadavCOVID3/Clean/%s/%s_R1_001.fastq.gz.clean.fq.gz /sc/arion/projects/weijiaTemp/TMP/MadavCOVID3/Clean/%s/%s_R2_001.fastq.gz.clean.fq.gz", $1, $1, $1, $1);

	printf("\n");

}' Sam > qc.sh

#nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --pro_code acc_zhangw09a --reqsub --jobprefix qc --resource 2000 --time 360 --convert no --verbose qc.sh &

nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --pro_code acc_zhangw09a --reqsub --jobprefix qc --resource 2000 --time 360 --convert no --verbose qc.sh &


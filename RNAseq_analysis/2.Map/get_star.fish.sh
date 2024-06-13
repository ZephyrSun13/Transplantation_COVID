
#ls /sc/arion/projects/zhangw09a/Data/Zeguo_Sun/project/3.Kidney/Fish/RNAseq-NPHS1/0.Data/*.gz | awk '{if(NR%2==1){printf("%s\t", $0)}else{printf("%s\n", $0)}}' > Sample.lst

rm -f run_star.sh

DIREC=/sc/arion/scratch/sunz04/Fish/

mkdir -p $DIREC

awk -v DIREC="$DIREC" '{

	printf("module load python/3.7.3");

	printf(" && rm -rf %s/%s_star && mkdir %s/%s_star && module load star/2.7.5b && STAR --runThreadN 5 --genomeDir /sc/arion/projects/weijiaTemp/TMP/DB/Ensembl/Zebrafish/GRCz11_ensembl_default/ --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within KeepPairs --readFilesIn %s %s --readFilesCommand zcat --outFileNamePrefix %s/%s_star/ ", DIREC, $1, DIREC, $1, $2, $3, DIREC, $1);

	printf(" && /sc/arion/projects/zhangw09a/PANDA/ext/samtools/samtools index %s/%s_star/Aligned.sortedByCoord.out.bam", DIREC, $1);

	printf(" && htseq-count -f bam -r pos -s no -m union --nonunique none %s/%s_star/Aligned.sortedByCoord.out.bam /sc/arion/projects/weijiaTemp/TMP/DB/Ensembl/Zebrafish/Danio_rerio.GRCz11.109.corr.gff3 > %s/%s_star/%s_count.txt", DIREC, $1, DIREC, $1, $1);

        printf(" && htseq-count -f bam -r pos -s no -m union --nonunique all %s/%s_star/Aligned.sortedByCoord.out.bam /sc/arion/projects/weijiaTemp/TMP/DB/Ensembl/Zebrafish/Danio_rerio.GRCz11.109.corr.gff3 > %s/%s_star/%s_count.multi.txt", DIREC, $1, DIREC, $1, $1);

	printf("\n");

}' Sample.lst > run_star.sh

#nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --convert no --pro_code acc_zhangw09a --jobprefix Star --resource 40000 --time 1440 --verbose run_star.sh &

#nohup /sc/arion/projects/zhangw09a/PANDA/ext_ZS/bin/common/common/qsub-sge.pl --queue premium --convert no --pro_code acc_zhangw09a --reqsub --jobprefix htseq --resource 20000 --time 600 --verbose run_star.sh &


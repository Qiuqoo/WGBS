find /data3/users/Qiuqoo/C17_genome/WGBS/rawdata/ -name "*_R1.fq.gz" | \
    sed 's/_R1.fq.gz//g' | \
    xargs -I {} -P 16 bash -c 'trim_galore --paired --gzip --phred33 --fastqc --fastqc_args "-t 16" -o /data3/users/Qiuqoo/C17_genome/WGBS/trimgalore {}_R1.fq.gz {}_R2.fq.gz'

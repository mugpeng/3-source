[TOC]



# 0. 环境配置

创建环境，激活：

```
conda create -y -n wes_py3 python=3 && conda activate wes_py3
```

安装软件：

```
conda install -y multiqc qualimap samtools trim-galore fastqc bwa 
```



下载gatk：

```
mkdir -p gatk4 &&  cd gatk4
wget  https://github.com/broadinstitute/gatk/releases/download/4.1.6.0/gatk-4.1.6.0.zip
unzip gatk-4.1.6.0.zip
```

版本：gatk-4.1.6.0



# 1. 搭建分析流程

```
mkdir -p 0.log 1.raw_fq 2.clean_fq 3.qc/{raw_qc,clean_qc} 4.align/qualimap 5.gatk/{gatk_marked,gatk_bqsr} 6.snv/{gatk_hc,gatk_clean_hc}
```



创建bin 文件，并添加到环境：

```
echo 'export PATH=~/0.bin:$PATH' >> ~/.bashrc

# 激活一下
$source ~/.bashrc

# 软连接到bin 文件
ln -s ~/3.biosoft/01_wes/gatk-4.1.6.0/gatk .

# 改名一下
mv gatk gatk4.0
```



# 2. 参考数据

文章描述是 ``279 primary tumor tissue and paired blood samples`` ，  本来想尝试一下使用配对样本数据走一遍流程，可是无奈样本中只有279个外显子建库数据，并没有配对的数据：

```
$ grep 'TT_WES' 'filereport_read_run_PRJNA486023_tsv (1).txt' | wc -l
     279
```



测序数据的位置：

```
/home/tnbc/public_data/FUSCCTNBC
```



测序用到的参考数据：

```
/home/data/server/
```



设置名称文件：

```
grep 'TT_WES' 'filereport_read_run_PRJNA486023_tsv (1).txt' | awk '{print $1}' > name.txt

 head name.txt 
SRR7696207
SRR8517853
SRR8517854
SRR8517855
SRR8517856
SRR8517857
SRR8517858
SRR8517859
SRR8517860
SRR8517861
```



文件已经补全了：

```
bio2@bio2-2288H-V5:/home/tnbc/public_data/FUSCCTNBC/WES$ ls *_1.fastq.gz|wc -l
279
bio2@bio2-2288H-V5:/home/tnbc/public_data/FUSCCTNBC/WES$ ls *_2.fastq.gz|wc -l 
279
bio2@bio2-2288H-V5:/home/tnbc/public_data/FUSCCTNBC/WES$
```



# 3. 质控

```
nohup fastqc -t 4 ./SRR* --outdir ../3.qc/raw_qc/ > ../0.log/raw_fq.log 2>&1 &
```



用multiqc 合并结果：

```
multiqc ./*zip -o ./multiqc
```



![image.png](https://upload-images.jianshu.io/upload_images/19725743-df81d033eb8c7f35.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)

# 4. 过滤及过滤后质控

批量处理：

```
dir=~/1.pipeline/01_wes/03_PRC_TNBC_ERBB2
cat ../name.txt | while read id
do
fq1=$dir/1.raw_fq/${id:0:10}_1.fastq.gz
fq2=$dir/1.raw_fq/${id:0:10}_2.fastq.gz
echo "trim_galore --paired --phred33 -q 28 --length 36 --stringency 3 --cores 4 --max_n 3 -o $dir/2.clean_fq $fq1 $fq2"
done > trim_galore.sh
```



```
dir=~/1.pipeline/01_wes/03_PRC_TNBC_ERBB2
cat ../2_name.txt | while read id
do
fq1=$dir/1.raw_fq/${id:0:10}_1.fastq.gz
fq2=$dir/1.raw_fq/${id:0:10}_2.fastq.gz
echo "trim_galore --paired --phred33 -q 28 --length 36 --stringency 3 --cores 4 --max_n 3 -o $dir/2.clean_fq $fq1 $fq2"
done > trim_galore2.sh
```



后台运行：

```
nohup bash ./trim_galore.sh > ../0.log/trim_galore.log 2>&1 &
```



过滤后数据：

```
$ ls -lh *.fq.gz | cut -d' ' -f 6-
1月  18 17:41 SRR8517889_1_val_1.fq.gz
1月  18 17:41 SRR8517889_2_val_2.fq.gz
1月  18 20:57 SRR8517909_1_val_1.fq.gz
1月  18 20:57 SRR8517909_2_val_2.fq.gz
1月  18 17:04 SRR8517919_1_val_1.fq.gz
1月  18 17:04 SRR8517919_2_val_2.fq.gz
1月  18 21:21 SRR8518016_1_val_1.fq.gz
1月  18 21:21 SRR8518016_2_val_2.fq.gz
1月  18 17:23 SRR8518103_1_val_1.fq.gz
1月  18 17:23 SRR8518103_2_val_2.fq.gz
1月  18 21:14 SRRzhaoyua_1_val_1.fq.gz
1月  18 21:14 SRRzhaoyua_2_val_2.fq.gz

```



再质控一波：

```
nohup fastqc -t 8 ./SRR* --outdir ../3.qc/clean_qc/ > ../0.log/clean_fq.log 2>&1 &
```



合并质控结果：

```
multiqc ./*zip -o ./multiqc
```

![image.png](https://upload-images.jianshu.io/upload_images/19725743-7bd99060726ac94d.png?imageMogr2/auto-orient/strip%7CimageView2/2/w/1240)



# 5. 比对及结果质控

```
dir=~/1.pipeline/01_wes/03_PRC_TNBC_ERBB2
cat ../name.txt | while read id
do
  echo "start bwa $id" $(date)
  bwa mem -M -t 8 -R "@RG\tID:${id}\tSM:${id}\tLB:WXS\tPL:Illumina" /home/data/server/reference/index/bwa/hg38 $dir/2.clean_fq/${id}_1_val_1.fq.gz $dir/2.clean_fq/${id}_2_val_2.fq.gz | samtools sort -@ 8 -m 8G -o $dir/4.align/${id}.bam -
  echo "finish bwa $id" $(date)
done

# 后台运行
nohup sh bwa.sh > ../0.log/bwa.log 2>&1 &
```



比对文件：

```
$ ls -lh *.bam | cut -d' ' -f 6-

```



qualimap质控一下：

```
dir=~/1.pipeline/01_wes/01_PRC_TNBC
cat $dir/name.txt | while read id
do
	qualimap bamqc --java-mem-size=10G -gff ~/2.data/UCSC/annotation/exon-regions/hg38.exon.bed -nr 100000 -nw 500 -nt 16 -bam $dir/4.align/${id}.bam -outdir ./${id}
done

nohup sh qualimap.sh 1>~/1.pipeline/01_wes/01_PRC_TNBC/0.log/qualimap.log 2>&1 &
```



合并质控结果：

```
multiqc SRR*
```



# 6. 质控结果

整合fq raw/clean qc 及bam qc结果：

范例：

| Sample  Name | %  GC | Ins.  size | ≥ 30X  | Median  cov | Mean  cov | %  Aligned | M Aligned | % Dups | % GC | length | M Seqs | %raw Dups | %raw GC | raw Mseqs |
| ------------ | ----- | ---------- | ------ | ----------- | --------- | ---------- | --------- | ------ | ---- | ------ | ------ | --------- | ------- | --------- |
| SRR8517860   | 52%   | 200        | 91.80% | 169.0X      | 199.7X    | 100.00%    | 123.1     | 30.30% | 54%  | 143 bp | 62     | 27.60%    | 50%     | 63.2      |
| SRR8517861   | 56%   | 176        | 87.80% | 145.0X      | 217.3X    | 100.00%    | 124       | 15.10% | 54%  | 141 bp | 37.4   | 30.50%    | 54%     | 63.6      |
| SRR8517862   | 56%   | 196        | 76.00% | 58.0X       | 87.6X     | 100.00%    | 85.3      | 14.30% | 54%  | 141 bp | 36.7   | 12.80%    | 52%     | 43.9      |
| SRR8517863   | 53%   | 176        | 93.40% | 229.0X      | 269.6X    | 100.00%    | 153.6     | 31.60% | 53%  | 142 bp | 72.6   | 30.90%    | 51%     | 78.8      |
| SRR8517864   | 51%   | 172        | 93.30% | 228.0X      | 253.9X    | 100.00%    | 153.2     | 12.70% | 52%  | 144 bp | 42.7   | 28.40%    | 50%     | 78.2      |
| SRR8517865   | 53%   | 203        | 90.00% | 124.0X      | 144.8X    | 100.00%    | 91        | 30.50% | 51%  | 144 bp | 76.8   | 14.00%    | 51%     | 51.1      |
| SRR8517866   | 56%   | 184        | 77.60% | 70.0X       | 108.9X    | 100.00%    | 65.3      | 14.20% | 51%  | 142 bp | 48.9   | 15.40%    | 54%     | 38.8      |
| SRR8517867   | 55%   | 180        | 90.50% | 126.0X      | 175.6X    | 100.00%    | 102       | 27.30% | 50%  | 144 bp | 61.6   | 32.00%    | 53%     | 74.7      |
| SRR8517868   | 56%   | 194        | 79.40% | 67.0X       | 112.8X    | 100.00%    | 70.6      | 28.20% | 50%  | 141 bp | 76.6   | 14.30%    | 54%     | 38.5      |
| SRR8517869   | 52%   | 200        | 89.10% | 120.0X      | 135.6X    | 100.00%    | 82.2      | 24.80% | 50%  | 144 bp | 56.4   | 25.00%    | 50%     | 57.4      |



# 7. 比对后gatk 处理

- 标注和去重复

```
dir=~/1.pipeline/01_wes/01_PRC_TNBC
cat $dir/name.txt | while read id
do
	BAM=$dir/4.align/${id}.bam
	if [ ! -f $dir/5.gatk/gatk_marked/status/ok.${id}_marked.status ]
	then
		echo "start MarkDuplicates for ${id}" `date`
		gatk4.0 --java-options "-Xmx20G -Djava.io.tmpdir=./" MarkDuplicates \
		-I ${BAM} \
		--REMOVE_DUPLICATES=true \
		-O $dir/5.gatk/gatk_marked/${id}_marked.bam \
		-M $dir/5.gatk/gatk_marked/${id}.metrics \
		1>~/1.pipeline/01_wes/01_PRC_TNBC/0.log/gatk_marked/${id}.log 2>&1 
		
		if [ $? -eq 0 ]
		then
			touch $dir/5.gatk/gatk_marked/status/ok.${id}_marked.status
		fi
		echo "end MarkDuplicates for ${id}" `date`
		samtools index -@ 16 -m 4G -b $dir/5.gatk/gatk_marked/${id}_marked.bam $dir/5.gatk/gatk_marked/${id}_marked.bai 1>~/1.pipeline/01_wes/01_PRC_TNBC/0.log/gatk_marked/${id}_index.log 2>&1 
	fi
done

nohup bash gatk_marked.sh 1>./marked-time.log 2>&1 &
```



status 文件：

```
$ ls -lh | cut -d' ' -f6- 

1月  12 21:31 ok.SRR8517860_marked.status
1月  12 21:50 ok.SRR8517861_marked.status
1月  12 22:14 ok.SRR8517862_marked.status
1月  12 22:37 ok.SRR8517863_marked.status
1月  12 23:01 ok.SRR8517864_marked.status
1月  12 23:15 ok.SRR8517865_marked.status
1月  12 23:26 ok.SRR8517866_marked.status
1月  12 23:41 ok.SRR8517867_marked.status
1月  12 23:52 ok.SRR8517868_marked.status
1月  13 00:05 ok.SRR8517869_marked.status
```



- 碱基质量重校正

```
snp=~/2.data/GATK/gatk-bundle/dbsnp_146.hg38.vcf.gz
indel=~/2.data/GATK/gatk-bundle/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
ref=~/2.data/GATK/gatk-bundle/Homo_sapiens_assembly38.fasta
dir=~/1.pipeline/01_wes/01_PRC_TNBC
cat $dir/name.txt | while read id
do
	if [ ! -f $dir/5.gatk/gatk_bqsr/status/ok.${id}_bqsr.status ]
	then
		echo "start BQSR for ${id}" `date`
		gatk4.0 --java-options "-Xmx20G -Djava.io.tmpdir=./"  BaseRecalibrator \
		-R $ref  \
		-I $dir/5.gatk/gatk_marked/${id}_marked.bam  \
		--known-sites ${snp} \
		--known-sites ${indel} \
		-O $dir/5.gatk/gatk_bqsr/${id}_recal.table \
		1>~/1.pipeline/01_wes/01_PRC_TNBC/0.log/gatk_bqsr/${id}_BaseRecalibrator.log 2>&1 
		
		gatk4.0 --java-options "-Xmx20G -Djava.io.tmpdir=./"  ApplyBQSR \
		-R $ref  \
		-I $dir/5.gatk/gatk_marked/${id}_marked.bam  \
		-bqsr $dir/5.gatk/gatk_bqsr/${id}_recal.table \
		-O $dir/5.gatk/gatk_bqsr/${id}_bqsr.bam \
		1>~/1.pipeline/01_wes/01_PRC_TNBC/0.log/gatk_bqsr/${id}_ApplyBQSR.log  2>&1 
        
    	if [ $? -eq 0 ]
		then
			touch $dir/5.gatk/gatk_bqsr/status/ok.${id}_bqsr.status
		fi
        echo "end BQSR for ${id}" `date`
    fi
done

nohup bash gatk_bqsr.sh 1>./bqsr-time.log 2>&1 &
```

status 文件：

```
$ ls -lh | cut -d' ' -f6- 

1月  13 11:31 ok.SRR8517860_bqsr.status
1月  13 13:10 ok.SRR8517861_bqsr.status
1月  13 14:38 ok.SRR8517862_bqsr.status
1月  13 16:38 ok.SRR8517863_bqsr.status
1月  13 18:44 ok.SRR8517864_bqsr.status
1月  13 20:06 ok.SRR8517865_bqsr.status
1月  13 21:07 ok.SRR8517866_bqsr.status
1月  13 22:32 ok.SRR8517867_bqsr.status
1月  13 23:45 ok.SRR8517868_bqsr.status
1月  14 01:01 ok.SRR8517869_bqsr.status
```





# 8. gatk-HC 找germline 变异

使用并行命令加速运行。



生成文件：

```
mkdir status vcf_{log,time_log,sh}
ls
gatk_vcf.sh  status  vcf_log  vcf_sh  vcf_time_log
```



批量打印：

```
snp=~/2.data/GATK/gatk-bundle/dbsnp_146.hg38.vcf.gz
ref=~/2.data/GATK/gatk-bundle/Homo_sapiens_assembly38.fasta
dir=~/1.pipeline/01_wes/01_PRC_TNBC

cat $dir/name.txt | while read id
do
(echo "if [ ! -f ./status/ok.${id}.status ]
then echo "start operation ${id}"" '$(date)'
echo "gatk4.0 --java-options \"-Xmx20G -Djava.io.tmpdir=./\" HaplotypeCaller -I $dir/5.gatk/gatk_bqsr/${id}_bqsr.bam -R ${ref} --dbsnp ${snp} -O $dir/6.snv/gatk_hc/${id}.vcf > ./vcf_log/${id}.log 2>&1"
echo "if [ \$? -eq 0 ]; then touch ./status/ok.${id}.status; fi
fi
echo "end for ${id}"" '$(date)') > ./vcf_sh/vcf_${id}.sh
done

bash vcf_script.sh &
```



查看某个批量文件：

```
$ cat vcf_SRR8517860.sh
if [ ! -f ./status/ok.SRR8517860.status ]
then echo start operation SRR8517860 $(date)
gatk --java-options "-Xmx20G -Djava.io.tmpdir=./" HaplotypeCaller -I /home/yzpeng/1.pipeline/01_wes/01_PRC_TNBC/5.gatk/gatk_bqsr/SRR8517860_bqsr.bam -R /home/yzpeng/2.data/GATK/gatk-bundle/Homo_sapiens_assembly38.fasta --dbsnp /home/yzpeng/2.data/GATK/gatk-bundle/dbsnp_146.hg38.vcf.gz -O /home/yzpeng/1.pipeline/01_wes/01_PRC_TNBC/6.snv/gatk_hc/SRR8517860.vcf > ./vcf_log/SRR8517860.log 2>&1
if [ $? -eq 0 ]; then touch ./status/ok.SRR8517860.status; fi
fi
echo end for SRR8517860 $(date)
```



接着直接在环境下运行指令：

```
cat $dir/name.txt | tail -3 | while read id
do
nohup bash ./vcf_sh/vcf_${id}.sh 1>./vcf_time_log/${id}_time.log 2>&1 &
done
```



并行命令，还是很快的，本来估计要到明天才能完成呢：

```
$ ls -lh *.vcf | cut -d ' ' -f6-
1月  14 13:17 SRR8517860.vcf
1月  14 13:44 SRR8517861.vcf
1月  14 14:43 SRR8517862.vcf
1月  15 01:14 SRR8517863.vcf
1月  14 14:01 SRR8517864.vcf
1月  14 12:56 SRR8517865.vcf
1月  14 12:24 SRR8517866.vcf
1月  14 12:56 SRR8517867.vcf
1月  14 13:38 SRR8517868.vcf
1月  14 12:40 SRR8517869.vcf
```



# 9. gatk 变异质控-ApplyRecalibration | 硬过滤

先来看一下过滤前的germilne 突变数目：

```
$ ls *.vcf | while read id ; do wc -l $id ;done
636878 SRR8517860.vcf
589118 SRR8517861.vcf
1278243 SRR8517862.vcf
614053 SRR8517863.vcf
793482 SRR8517864.vcf
510152 SRR8517865.vcf
367110 SRR8517866.vcf
489222 SRR8517867.vcf
411711 SRR8517868.vcf
446910 SRR8517869.vcf
```



使用gatk 最佳实践的ApplyRecalibration（建模）、VariantRecalibrator（质控过滤）。（或是硬过滤，取决于样本量多少，以30为界，少则硬过滤）



# 10. ANNOVAR 注释

采用文章中的注释方式：

```
dir=~/1.pipeline/01_wes/03_PRC_TNBC_ERBB2
cat $dir/name.txt | while read id
do
	if [ ! -f $dir/7.annotation/annovar/status/ok.${id}_annovar.status ]
	then
		echo "start annovar for ${id}" `date`
      table_annovar.pl $dir/6.snv/gatk_clean_hc/${id}.clean.vcf ~/2.data/snv_annotation/annovar/ \
      -buildver hg38 \
      -out $dir/7.annotation/annovar/${id} \
      -remove \
      -protocol refGene,knownGene,clinvar_20170905 \
      -operation g,g,f \
      -nastring . \
      -vcfinput
    if [ $? -eq 0 ]
			then
				touch $dir/7.annotation/annovar/status/ok.${id}_annovar.status
			fi
        echo "end annovar for ${id}" `date`
    fi
done

nohup bash annovar.sh 1>./annovar-time.log 2>&1 &
```







# 其他

## 服务器使用scp下载和上传

- 下载

```
scp -oPort=6652 user@120.77.173.108:/path/filename /home/folder 
```

- 上传

```
scp -oPort=6652 /home/filename user@120.77.173.108:/path/folder
```



```
scp -oPort=6652 -r yzpeng@120.77.173.108:~/1.pipeline/01_wes/03_PRC_TNBC_ERBB2/3.qc/clean_qc/multiqc/multiqc_report.html ~/Downloads/

scp -oPort=6000 -r  yzpen@119.45.49.20:~/project/001_BRCA-KI/wes/2.align/qualimap/WGC_all ~/Downloads/
```



## 使用奶牛快传

https://www.yuque.com/mugpeng/pc_skills/kqduxr
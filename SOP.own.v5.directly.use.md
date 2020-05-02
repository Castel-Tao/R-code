# RNA-seq SOP

## Part three: Linux scripts and results

### A0. Preparation: Make directory

```shell

#!/bin/bash

#1. 存放原始fastq文件的目录
mkdir ./fastq
mkdir ./fastq/unzip

#2. 存放qc的目录
mkdir ./qc
mkdir ./qc/fastqc
mkdir ./qc/fastqScreen

#3. 存放trimmomatic结果的目录
mkdir ./trimmomatic
mkdir ./trimmomatic/trim_fastq
mkdir ./trimmomatic/trim_report
mkdir ./trimmomatic/trim_fastqc

#4. 存放bam文件的目录
mkdir ./bam
mkdir ./bam/sam
mkdir ./bam/sort_bam
mkdir ./bam/hisat2_report
mkdir ./bam/Samtools_report

#5. Stringtie结果目录
mkdir ./stringtie
mkdir ./stringtie/st_report
mkdir ./stringtie/st_output

#6. featureCounts结果目录
mkdir ./feature_counts
mkdir ./feature_counts/fc_report
mkdir ./feature_counts/fc_output

#后面这步不是必须的

#7. rseqc结果目录
mkdir ./rseqc
mkdir ./rseqc/rs_report
mkdir ./rseqc/rs_output


#在 /st_output 下：  
for i in $(cat sample.list);do mkdir $i;done
```



### A. fastqc-fastqScreen 对原始测序数据进行质量控制

```shell
#!/bin/bash                   
#PBS -l nodes=compute-node10:ppn=4        
#PBS -l walltime=72:00:00    
#PBS -l mem=12gb               
#PBS -l vmem=12gb              
#PBS -j oe                    
#PBS -N fastqc-Screen-1-10            

#fastqc
##project directory, you may change accordingly
project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/RenXiaojun
fastq_Dir=${project_dir}/fastq 
fastqc_output=${project_dir}/qc/fastqc
conf_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/fastqScreen
Screen_output=${project_dir}/qc/fastqScreen

R1=R1_001.fastq
R2=R2_001.fastq
nt=4

##load modules
module load fastqc-0.11.7
module load bowtie2-2.3.4
module load fastq-screen-0.11.3

##fastqc for raw data
fastqc -t ${nt} -o ${fastqc_output} ${fastq_Dir}/${name}_${R1}.gz
fastqc -t ${nt} -o ${fastqc_output} ${fastq_Dir}/${name}_${R2}.gz

gunzip -c ${fastq_Dir}/${name}_${R1}.gz > ${fastq_Dir}/${name}_${R1}
gunzip -c ${fastq_Dir}/${name}_${R2}.gz > ${fastq_Dir}/${name}_${R2}

fastq_screen --aligner bowtie2 --threads ${nt} --conf ${conf_Dir}/fastq_screen.conf --outdir ${Screen_output} ${fastq_Dir}/${name}_${R1} ${fastq_Dir}/${name}_${R2}

rm ${fastq_Dir}/${name}_${R1}
rm ${fastq_Dir}/${name}_${R2}

```





### B. Trimmomatic  去接头

```shell
#!/bin/bash
#PBS -l nodes=compute-node10:ppn=4
#PBS -l walltime=72:00:00
#PBS -l mem=12gb               
#PBS -l vmem=12gb 
#PBS -j oe
#PBS -N Trimmomatic

##load modules
module load Trimmomatic-0.36	

project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/RenXiaojun
fastq_Dir=${project_dir}/fastq
TrimmOutput=${project_dir}/trimmomatic/trim_fastq
adapter=/mnt/pgx_src_data_pool_4/home/taoyichen/software/trimmomatic_adapter/TruSeq3-PE.fa

R1=R1_001.fastq
R2=R2_001.fastq
nt=4

trimmomatic PE -phred33 -threads ${nt} ${fastq_Dir}/${name}_${R1}.gz ${fastq_Dir}/${name}_${R2}.gz -baseout ${TrimmOutput}/${name}_fastq.gz ILLUMINACLIP:${adapter}:2:30:10:1:true HEADCROP:15 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

##################################注 释#####################################
#PE: pair end sequence
#-phred33: 还可以是phred64，目前最常用33
#-baseout： 软件对自动为4个输出文件命名
#MINLEN:如果reads修剪之后低于这个长度，则丢弃
#ILLUMINACLIP：接头路径：第一步 seed 搜索时允许的错配碱基个数，例如2：PE的palindrome clip（D模式）模式下，R1、R2 之间至少比对分值，才会进行接头切除，例如30：指定切除接头序列的最低比对分值（A/B 模式），通常7-15之间：palindrome 模式下可以切除的接头序列最短长度，例如1bp：D 模式下，R1和R2在去除了接头序列之后剩余的部分是完全反向互补的，默认参数false。
#HEADCROP：从 reads 的起始开始直接切除部分碱基，一般是10bp,每次都要更改
#4g大小的fastq文件，8核，24g内存,15min;
#trimmomatic这一步内存一定要大一些；至少是文件大小的4倍吧
###########################################################################

```


### A2. fastqc-fastqScreen 对Trim数据进行质量控制
```shell

#!/bin/bash                   
#PBS -l nodes=compute-node09:ppn=4        
#PBS -l walltime=72:00:00    
#PBS -l mem=12gb               
#PBS -l vmem=12gb              
#PBS -j oe                    
#PBS -N fastqc-Screen-1-10            

#fastqc
##project directory, you may change accordingly
project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/RenXiaojun
fastq_Dir=${project_dir}/trimmomatic/trim_fastq 
fastqc_output=${project_dir}/trimmomatic/trim_fastqc
conf_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/fastqScreen
Screen_output=${project_dir}/qc/fastqScreen

R1=fastq_1P
R2=fastq_2P
nt=4

##load modules
module load fastqc-0.11.7
module load bowtie2-2.3.4
module load fastq-screen-0.11.3

##fastqc for raw data
fastqc -t ${nt} -o ${fastqc_output} ${fastq_Dir}/${name}_${R1}.gz
fastqc -t ${nt} -o ${fastqc_output} ${fastq_Dir}/${name}_${R2}.gz

gunzip -c ${fastq_Dir}/${name}_${R1}.gz > ${fastq_Dir}/${name}_${R1}
gunzip -c ${fastq_Dir}/${name}_${R2}.gz > ${fastq_Dir}/${name}_${R2}

fastq_screen --aligner bowtie2 --threads ${nt} --conf ${conf_Dir}/fastq_screen.conf --outdir ${Screen_output} ${fastq_Dir}/${name}_${R1} ${fastq_Dir}/${name}_${R2}

rm ${fastq_Dir}/${name}_${R1}
rm ${fastq_Dir}/${name}_${R2}
```



### C. Hisat2 and Samtools

```shell
#!/bin/bash                   
#PBS -l nodes=compute-node06:ppn=7       
#PBS -l walltime=72:00:00    
#PBS -l mem=21gb               
#PBS -l vmem=21gb              
#PBS -j oe                    
#PBS -N hisat_ren-1-5

#load modules
module load hisat2-2.1.0	
module load samtools-1.8

#index directory
project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/RenXiaojun
TrimmOutput=${project_dir}/trimmomatic/trim_fastq
index_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human/hisat2_download_index/grch38_snp_tran
sam_output=${project_dir}/bam/sam
bam_Output=${project_dir}/bam/sort_bam

nt=7

#mapping
hisat2 -p ${nt} --dta -x ${index_Dir}/genome_snp_tran -1 ${TrimmOutput}/${name}_fastq_1P.gz -2 ${TrimmOutput}/${name}_fastq_2P.gz -S ${sam_output}/${name}.sam 

#将sam文件转换为bam文件
samtools view -bS ${sam_output}/${name}.sam >${bam_Output}/${name}.bam
#对bam文件进行排序
samtools sort -@ ${nt} ${bam_Output}/${name}.bam -o ${bam_Output}/${name}.sort.bam
#对排序后的bam文件创建索引
samtools index ${bam_Output}/${name}.sort.bam
#得到比对率
samtools flagstat ${bam_Output}/${name}.sort.bam >${bam_Output}/${name}.txt


## 这一步基本是越舍得花线程越快
#################################注 释######################################
#--trim5 10  :从5'端，左端开始去除10个碱基
#-t:时间
#-x:index
#-1:左端文件
#-2:右端文件
#-S：sam文件输出路径
#-dta：后面如果用stringtie重构转录本的话这一步是必须的
#4g大小的fastq文件，10核，30g内存,20min；
# update 根据sam文件大小，差不多2min/5G。所以50G sam ~ 20min
###########################################################################
#################################注 释######################################
#4g大小的fastq文件，10核，30g内存,30min
###########################################################################


```



### E. StringTie 转录本重构

```shell
StringTie要先按照样本名在./Project/stringtie/st_output/Ballgown 下面新建文件夹
project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/Ren
for i in $(cat sample.list);do mkdir ${project_dir}/stringtie/st_output/$i;done

#!/bin/bash                   
#PBS -l nodes=compute-node01:ppn=4          
#PBS -l walltime=72:00:00    
#PBS -l mem=12gb               
#PBS -l vmem=12gb              
#PBS -j oe                    
#PBS -N stringtie

#这里建议按照样本名新建文件夹， 这样ballgown不会被覆盖
module load stringtie-1.3.4

project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/Ren
bam_Input=${project_dir}/bam/sort_bam
GTF_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human
ST_ouputdir=${project_dir}/stringtie/st_output/${name}

stringtie -p 4 ${bam_Input}/${name}.sort.bam -G ${GTF_Dir}/Homo_sapiens.GRCh38.93.gtf -o ${ST_ouputdir}/${name}.gtf -C ${ST_ouputdir}/${name}.cov.ref.gtf -e -A ${ST_ouputdir}/${name}.gene.abundance.txt -B

################################注 释#######################################
#-G:基因组注释文件
#-o:输出gtf文件
#-B:用于下游Ballgown分析,不要用-b,亲测只输出一个文件夹，其它样本的ballgown都被覆盖掉了
#-A:用于输出Gene abundance文件
#-e:表示只对参考基因组注释文件中的转录组进行定量
#4g大小的fastq文件，8核，24g内存,20min；有的时候只要10min
###########################################################################

```

​	

## F. featureCounts 从bam文件中提取gene feature

```SHELL
#!/bin/bash                   
#PBS -l nodes=compute-node01:ppn=4         
#PBS -l walltime=72:00:00    
#PBS -l mem=12gb               
#PBS -l vmem=12gb              
#PBS -j oe                    
#PBS -N featureCounts

module load Subread-1.6.2

project_dir=/mnt/pgx_src_data_pool_4/home/taoyichen/Ren
gtfdir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human
outputdir=${project_dir}/feature_counts/fc_output
bam_Input=${project_dir}/bam/sort_bam

featureCounts -B -p -T 4 -g gene_id -a ${gtfdir}/Homo_sapiens.GRCh38.93.gtf -o ${outputdir}/${name}.counts ${bam_Input}/${name}.sort.bam
#这里我去掉了-C我不确定是不是合适
################################注 释#######################################
#-a:gtf文件路径
#-B:基因组注释文件
#-C:嵌合fragments will not be count
#-p:双端测序文件
#-t:specify the feature type; 默认就是exon,所以不指定也可以；-t exon
#-g:从注释文件中提取Meta-feature信息用于read count，默认是gene_id
#-o:输出文件
#这一步很快的 10min；所以还是用8 cpu,24g 内存吧
# 但是我用4 cpu 12G内存就很慢，20min~50min不等；好像4，6，8，10，12节点运行featureCount快一点？
###########################################################################

```



### G. RseQC 对比对结果进行质量控制

```shell
#!/bin/bash                   
#PBS -l nodes=compute-node01:ppn=10          
#PBS -l walltime=72:00:00    
#PBS -l mem=30gb               
#PBS -l vmem=30gb              
#PBS -j oe                    
#PBS -N rseqc

module load rseqc-2.6.4
module load python-2.7

bam_Input= ./Project/bam/sorted_bam
outputdir= ./Project/rseqc/rs_output
Bedir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

#刚才已经说了这是一个可执行文件，所以不需要安装，直接使用./genePredToGtf;它可以用来转化gtf成为bed12文件；GenePred ？？

# 1.clipping_profile.py  这个模块用于评估RNA-seq的BAM或SAM文件中的有切除核苷酸的reads情况.
# 这个模块会生成.r格式的作图脚本以及.pdf格式的报告文件以及.xls的数据文件; slow (55 cpu,240 RAM,5g bam,about 21min.)
clipping_profile.py -i ${bam_Input}/${name}.sort.bam -s "PE" -o ${outputdir}/${name}

# 2.deletion_profile.py: 也就是reads deletion位点的分布;fast (55 cpu,240 RAM,5g bam,about 12min.)
#deletion_profile.py -i ${datadir}/${name}.sort.bam -l 150 -o ${outputdir}/${name}  #这个我在实际中没有运行成功

# 3. FPKM_count.py:  slow (55 cpu,240 RAM,5g bam,about 30min.)
# 根据read count和gene注释文件(bed12格式)计算每个基因的FPM (fragment per million)或者FPKM(fragment per million mapped reads per kilobase exon).
#FPKM_count.py -r ${Bamdir}/rn6.bed -i ${datadir}/$nam.sort.bam -o ${outputdir}/$nam

# 4. geneBody_coverage.py: 计算RNA-seq  reads在基因上的覆盖度 ; very slow (55 cpu,240 RAM,5g bam,about 50min.)
geneBody_coverage.py -r ${Bedir}/GRCh38.bed12 -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name}

# 5.inner_distance.py: 针对双端测序,计算read pairs的内部距离或者插入距离; very fast (55 cpu,240 RAM,5g bam,about 2min.)
inner_distance.py -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name} -r ${bam_Input}/GRCh38.bed12  

# 6.read_distribution.py: fast (55 cpu,240 RAM,5g bam,about 18min.)
# 这个模块根据提供的BAM/SAM文件和bed12格式的gene模型文件就按比对上去的reads在基因组上的分布情况,比如在CDS exon,5'UTR exon,intro,基因间区域的reads分布
#read_distribution.py  -i ${datadir}/$nam.sort.bam -r ${Bamdir}/rn6.bed

# 7.read_GC.py 计算reads的GC含量分布；fast (55 cpu,240 RAM,5g bam,about 10min.)
read_GC.py -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name}

# 8.RNA_fragment_size.py: 在map后计算每个gene上的fragment的大小,包括:每个gene上所有的fragment的均值,中位数,方差 
#very fast (55 cpu,240 RAM,5g bam,about 10min.)
#RNA_fragment_size.py -r ${Bamdir}/GRCh38.bed12 -i ${datadir}/${name}.sort.bam > ${outputdir}/${name}.fragSize

# 9.bam_stat.py;默认output到提交script的文件夹,或者我也可以想上面一样用>输出？  fast (55 cpu,240 RAM,5g bam,about 15min.)
bam_stat.py -i ${bam_Input}/${name}.sort.bam  

# 10.tin.py:这个模块用来在转录本级别计算RNA完整性TIN (transcript integrity number)值; very slow(55 cpu,240 RAM,5g bam,about 55min.)
# tin.py -r ${Bamdir}/GRCh38.bed12 -i ${datadir}/${name}.sort.bam  #这个我在实际中没有运行成功


#!/bin/bash                   
#PBS -l nodes=compute-node06:ppn=30        
#PBS -l walltime=72:00:00    
#PBS -l mem=100gb               
#PBS -l vmem=100gb              
#PBS -j oe                    
#PBS -N wuxi_genebody

module load rseqc-2.6.4
module load python-2.7

bam_Input=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/NovaSeq/wuxi_polyA/bam/sort_bam
outputdir=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/NovaSeq/wuxi_polyA/rseqc/rs_output
Bedir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

#bam_stat 很有用,像这种没有返回值的可以把#PBS -N 改成${name}
bam_stat.py -i ${bam_Input}/${name}.sort.bam

#CDS exon,5'UTR exon,intro,基因间区域的reads分布,没有-o选项
read_distribution.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam 

#insert size
inner_distance.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam -o 
${outputdir}/${name}

#GC content 
read_GC.py -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name}

# geneBody_coverage.py
geneBody_coverage.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input} -o ${outputdir}/WUX

#tin值，这是一个贼慢的过程。要10个小时+，所以最好一次提交上所有的样本运行。10G的fastq.gz要24h+
tin.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam  



#!/bin/bash
#PBS -l nodes=compute-node10:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=12gb
#PBS -l vmem=12gb
#PBS -j oe
#PBS -N TNBC_batches_genebody

module load rseqc-2.6.4
module load python-2.7

bam_Input=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/TNBC_batches/bam/sort_bam
outputdir=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/TNBC_batches/rseqc/rs_output
Bedir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

# geneBody_coverage.py
geneBody_coverage.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name}



#!/bin/bash                   
#PBS -l nodes=compute-node11:ppn=1          
#PBS -l walltime=72:00:00    
#PBS -l mem=5gb               
#PBS -l vmem=5gb              
#PBS -j oe     
#PBS -N wuxi_ribo_TIN

module load rseqc-2.6.4
module load python-2.7

bam_Input=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/NovaSeq/wuxi_ribo/bam/sort_bam
outputdir=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/NovaSeq/wuxi_ribo/rseqc/rs_output
Bedir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

#tin值
tin.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam  




#insert size
#inner_distance.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam -o #${outputdir}/${name}
#CDS exon,5'UTR exon,intro,基因间区域的reads分布,没有-o选项
#read_distribution.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam 
#GC content 
#read_GC.py -i ${bam_Input}/${name}.sort.bam -o ${outputdir}/${name}.GC

#tin值
tin.py -r ${Bedir}/GRCh38.2.bed12 -i ${bam_Input}/${name}.sort.bam  
```



### H. qualimap比对

```
#!/bin/bash                   
#PBS -l nodes=compute-node10:ppn=30        
#PBS -l walltime=72:00:00    
#PBS -l mem=100gb               
#PBS -l vmem=100gb              
#PBS -j oe                    
#PBS -N quali_test

module load qualimap-2.1.3

bam_Input=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/AGR/bam/sort_bam
outputdir=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/AGR/qualimap/quali_report
GTF_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

qualimap rnaseq -bam ${bam_Input}/SEQC_ILM_AGR_A_1.sort.bam -gtf  ${GTF_Dir}/Homo_sapiens.GRCh38.93.gtf -outdir ${outputdir}/rnaseq_qc_results --java-mem-size=8G

#!/bin/bash                   
#PBS -l nodes=compute-node10:ppn=30        
#PBS -l walltime=72:00:00    
#PBS -l mem=30gb               
#PBS -l vmem=30gb              
#PBS -j oe                    
#PBS -N quali_test

module load qualimap-2.1.3

bam_Input=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/AGR/bam/sort_bam
outputdir=/mnt/pgx_src_data_pool_4/home/taoyichen/sampleAB/AGR/qualimap/quali_report
GTF_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human

qualimap bamqc -bam ${bam_Input}/SEQC_ILM_AGR_A_1.sort.bam -outdir ${outputdir}/rnaseq_qc_results --java-mem-size=30G
```





### Full-automatic script

说明：	新服务器40核，120g内存每个节点。9号11号好像不能用

​		老服务器60核，240g内存每个节点

#### 先同时进行fastqc，fastqScreen，Trimmomatic三个步骤

```shell
# 在./Project/trimmomatic/trim_report 目录中提交

#!/bin/bash
#PBS -l nodes=compute-node03:ppn=1
#PBS -l walltime=72:00:00
#PBS -l mem=20gb               
#PBS -l vmem=20gb 
#PBS -j oe
#PBS -N Trimmomatic

##load modules
module load Trimmomatic-0.36
module load fastqc-0.11.7

fastq_Dir= ./Project/fastq
Trim_Output=./Project/trimmomatic/trim_fastq
adapter=/mnt/pgx_src_data_pool_4/home/taoyichen/software/trimmomatic_adapter/TruSeq3-PE.fa

trimmomatic PE -phred33 ${fastq_Dir}/${name}_R1.fastq.gz ${fastq_Dir}/${name}_R2.fastq.gz -baseout ${Trim_Output}/${name}_fastq.gz ILLUMINACLIP:${adapter}:2:30:10:1:true HEADCROP:10 LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:50

#trimmomatic 之后直接fastqc
fastqc_output= ./Project/trimmomatic/trim_fastqc

fastqc -o ${fastqc_output} ${Trim_Output}/${name}_R1.fastq.gz
fastqc -o ${fastqc_output} ${Trim_Output}/${name}_R2.fastq.gz
```

```shell
# 在./Project/qc 目录中提交

#!/bin/bash                   
#PBS -l nodes=compute-node02:ppn=1         
#PBS -l walltime=72:00:00    
#PBS -l mem=10gb               
#PBS -l vmem=10gb              
#PBS -j oe                    
#PBS -N fastqc-fastqScreen  

nt=1

module load bowtie2-2.3.4
module load fastq-screen-0.11.3

conf_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/fastqScreen
fastq_Dir=./Project/fastq
Screen_output=./Project/qc/fastqScreen

fastq_screen --aligner bowtie2 --conf ${conf_Dir}/fastq_screen.conf --threads ${nt} --outdir ${Screen_output} ${fastq_Dir}/${name}_R1.fastq ${fastq_Dir}/${name}_R2.fastq
```

```shell
# 在./Project/qc 目录中提交

#!/bin/bash                   
#PBS -l nodes=compute-node01:ppn=1         
#PBS -l walltime=72:00:00    
#PBS -l mem=10gb               
#PBS -l vmem=10gb              
#PBS -j oe                    
#PBS -N fastqc-fastqScreen            

fastq_Dir= ./Project/fastq
fastqc_output= ./Project/qc/fastqc

module load fastqc-0.11.7

fastqc -o ${fastqc_output} ${fastq_Dir}/${name}_R1.fastq.gz
fastqc -o ${fastqc_output} ${fastq_Dir}/${name}_R2.fastq.gz

module load bowtie2-2.3.4
module load fastq-screen-0.11.3
```



#### 有结果之后check QC，如果没问题做mapping

```shell
# 在./Project/bam/hisat2_report 目录中进行提交

#!/bin/bash                   
#PBS -l nodes=compute-node04:ppn=8        
#PBS -l walltime=72:00:00    
#PBS -l mem=24gb               
#PBS -l vmem=24gb              
#PBS -j oe                    
#PBS -N hisat2-samtools-stringtie-featurecounts

#load modules
module load hisat2-2.1.0	
module load samtools-1.8

#index directory
index_Dir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human/hisat2_download_index/grch38_snp_tran
Trim_Output=./Project/trimmomatic/trim_fastq
sam_output=./Project/bam/sam
bam_Output= ./Project/bam/sort_bamm

#mapping
hisat2 -t --dta -x ${index_Dir}/genome_snp_tran -1 ${Trim_Output}/${name}_fastq_1P.gz -2 ${Trim_Output}/${name}_fastq_2P.gz -S ${sam_output}/${name}.sam 

# samtools
#将sam文件转换为bam文件
samtools view -bS ${sam_output}/${name}.sam >${bam_Output}/${name}.bam
#对bam文件进行排序
samtools sort -m 5000000000 ${bam_Output}/${name}.bam -o ${bam_Output}/${name}.sort.bam
#对排序后的bam文件创建索引
samtools index ${bam_Output}/${name}.sort.bam
```



#### 然后同时用stringtie，feature counts对bam进行转录本定量

```shell
#!/bin/bash                   
#PBS -l nodes=compute-node05:ppn=10          
#PBS -l walltime=72:00:00    
#PBS -l mem=30gb               
#PBS -l vmem=30gb              
#PBS -j oe                    
#PBS -N stringtie-featurecount

# feature counts
module load Subread-1.6.2

gtfdir=/mnt/pgx_src_data_pool_4/home/taoyichen/software/gtf/human
outputdir=./Project/stringtie/st_output
bamdir=./Project/bam/sort_bamm

featureCounts -T 5 -p -B -C -t exon -g gene_id -a ${gtfdir}/Homo_sapiens.GRCh38.93.gtf -o ${outputdir}/counts.txt ${bamdir}/*.sort.bam


#stringtie
module load stringtie-1.3.4

mkdir ./Project/stringtie/st_output/${name}
ST_ouputdir=./Project/stringtie/st_output/${name}

stringtie ${bamdir}/${name}.sort.bam -G ${GTF_Dir}/Homo_sapiens.GRCh38.93.gtf -o ${ST_ouputdir}/${name}.gtf -C ${ST_ouputdir}/${name}.cov.ref.gtf -e -A ${ST_ouputdir}/${name}.gene.abundance.txt -B
```


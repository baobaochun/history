#2020/8/2
samtools view -bu -f 12 -F 256 /data2/home/lichunhui/lunfish_meta/Dormancy/03_bwa_map/1-1_BDME202032413-1a.bam > \
~/data2/lunfish_meta/Dormancy/04_select_unmap/1-1_BDME202032413-1a_unmap1.bam

#2020/8/3
python kraken2_annotation.py /data2/home/yanxiaoting/project/Lungfish-Mete/Dormancy/06.kraken2/ \
/public/home/lichunhui/project/lungfish_meta/kraken2/dormancy_result/

python kraken2_annotation.py /data2/home/yanxiaoting/project/Lungfish-Mete/No_dormancy/06.kraken2/ \
/public/home/lichunhui/project/lungfish_meta/kraken2/no_dormancy_result/

python kraken2_anno_num.py /data2/home/yanxiaoting/project/Lungfish-Mete/Dormancy/06.kraken2/ \
/public/home/lichunhui/project/lungfish_meta/kraken2/dormancy_result/

python kraken2_anno_num.py /data2/home/yanxiaoting/project/Lungfish-Mete/No_dormancy/06.kraken2/ \
/public/home/lichunhui/project/lungfish_meta/kraken2/no_dormancy_result/

#计算bray-curtis距离并进行PCoA排序分析
library(dplyr)
library(permute)
library(lattice)
library(vegan)

in_file1 <- read.csv('C:/Users/bbc/Desktop/dphylum_abun_num.csv')
in_file2 <- read.csv('C:/Users/bbc/Desktop/ndphylum_abun_num.csv')

colnames(in_file1) <- c('name','dsample2','dsample3','dsample1')
colnames(in_file2) <- c('name','ndsample1','ndsample2','ndsample3')

ddata <- full_join(in_file1,in_file2,by='name')
rownames(ddata) <- ddata[,1]
ddata <- ddata[,-1]
ddata <- t(ddata)

b_rdis <- vegdist(ddata,method = 'bray',na.rm = TRUE)
b_rdis_data <- as.matrix(b_rdis)


pcoa <- cmdscale(b_rdis_data,k=(nrow(ddata)-1),eig = TRUE)
pcoa$eig
point <- data.frame(pcoa$points)

#plot
pcoa_eig <- (pcoa$eig[1:2]/sum(pcoa$eig))

sample_site <- data.frame(pcoa$points)[1:2]
sample_name <- rownames(sample_site)
names(sample_site)[1:2] <- c('PCoA1','PCoA2')
sample_site['sample'] <- sample_name

library(ggplot2)
ggplot(sample_site,aes(x=sample_site$PCoA1,y=sample_site$PCoA2,color=sample_site$sample))+
  geom_point(size=5)

#2020/8/4
#下载ncbi-nr数据库
nohup wget -c ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz >log.log 2>&1 &

#解压nr.gz
time unpigz -k -p 16 nr.gz

#2020/8/5
#构建diamond参考
nohup time diamond makedb --in nr -d nr -p 12 > log.log 2>&1 &
#使用nr数据库进行比对注释
nohup time diamond blastx --db /data2/home/lichunhui/database/nr/nr.dmnd -t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ -p 20 -q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa --daa /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.daa > /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &

nohup time diamond blastx --db /data2/home/lichunhui/database/nr/nr.dmnd -t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ -p 20 -q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa -o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.nr.m8 > /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &

#KEGG注释
#构建diamond参考
nohup time diamond makedb --in /public/home/yanxiaoting/database/kegg_all_clean/kegg_all_clean.fa -d kegg -p 10 > log.log 2>&1 &
#使用kegg数据库进行比对注释
nohup time diamond blastx --db /data2/home/lichunhui/database/kegg/kegg.dmnd -t /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/ -p 20 -q /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa -o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.kegg.m8 > /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/log.log 2>&1 &

#安装下载megan
wget -c https://software-ab.informatik.uni-tuebingen.de/download/megan6/MEGAN_Community_unix_6_19_6.sh
#megan下的NCBI-nr编号与物种和功能注释
wget -c https://software-ab.informatik.uni-tuebingen.de/download/megan6/megan-map-Jul2020-2.db.zip
#解压
unzip megan-map-Jul2020-2.db.zip
#转化daa文件为MEGAN特有的rma文件
~/megan/tools/daa2rma -i /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.daa -top 50  -mdb /data2/home/lichunhui/database/megan/megan-map-Jul2020-2.db  -o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.rma
#提取物种注释信息
~/megan/tools/rma2info -i /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy.rma -r2c Taxonomy -c2c Taxonomy -n true -p true -r true -o /data2/home/lichunhui/lunfish_meta/Dormancy/08_diamond_nr/dormancy_taxon.txt


#安装eggnog-mapper，环境依赖python=2.7
conda create -n eggnog python=2.7
wget https://github.com/eggnogdb/eggnog-mapper/archive/1.0.3.tar.gz -O ~/software/eggnog-mapper-1.0.3.tar.gz
tar zxf ~/software/eggnog-mapper-1.0.3.tar.gz -C /data2/home/lichunhui/database/eggnog/
cd /data2/home/lichunhui/database/eggnog/eggnog-mapper-1.0.3
nohup python download_eggnog_data.py -y -f euk bact arch viruses > log.log 2>&1 & #下载不完整，出错
nohup ./download_eggnog_data.py euk bact arch viruses > log.log 2>&1 & #换成这个下载命令,还是下载不了，缺文件，气死人，不用1.0版本了

#使用eggnog-mapper进行注释,默认基于hmmer
nohup python /data2/home/lichunhui/database/eggnog/eggnog-mapper-1.0.3/emapper.py -d bact -i /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa --data_dir /data2/home/yanxiaoting/biosoft/eggnog-mapper/eggnog-mapper-1.0.3/data -o eggNOG_hmmer --cpu 20 --no_file_comments --translate > log.log 2>&1 &

#重新安装2.0版本的eggnog-mapper
git clone https://github.com/jhcepas/eggnog-mapper.git
#使用v2.0 eggnog-mapper进行注释
nohup python /data2/home/lichunhui/database/eggnog/eggnog-mapper/emapper.py -m diamond -i /data2/home/lichunhui/lunfish_meta/Dormancy/06_megahit_zz/dormancy.fa -o eggNOG_diamond  --cpu 24 --translate > log_diamond.log 2>&1 &


#安装metawrap进行分箱
conda create -n metawrap python=2.7
conda config --add channels ursky
conda install -y -c ursky metawrap-mg
#分箱
nohup metawrap binning -t 24 --metabat2 --maxbin2 --concoct \
-a /data2/home/yanxiaoting/project/Golden_Snub-nosed_Monkey/Meta/Faece-8/03.assembly-IDBA/11.assembly.fa \
-o /data2/home/lichunhui/test/metawrap_bin \
/data2/home/lichunhui/test/reads_1.fastq \
/data2/home/lichunhui/test/reads_2.fastq > bin.log 2>&1 &
#合并分箱，合并分箱过程出现了点问题
nohup metawrap bin_refinement -o bin_refinement -t 20 -c 50 -x 10 \
-A /data2/home/lichunhui/test/metawrap_bin/metabat2_bins/ \
-B /data2/home/lichunhui/test/metawrap_bin/maxbin2_bins/ \
-C /data2/home/lichunhui/test/metawrap_bin/concoct_bins/ > bin_refinement.log 2>&1 &
#要在合并分箱后使用checkM进行分箱评估，需要提前下载好checkM数据库并添加进metawrap的搜索路径
mkdir checkm_folder
checkm data setRoot
# CheckM will prompt to to chose your storage location... Give it the path to the folder you just made.

# Now manually download the database:
cd checkm_folder
wget https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz
tar -xvf *.tar.gz
rm *.gz

concoct的分箱出现了问题，尝试单独运行concoct来得到结果
concoct单独运行也出现了问题



#安装kraken2
git clone https://github.com/DerrickWood/kraken2.git
#kraken2数据库路径
--db /data2/home/yanxiaoting/biosoft/kraken2/KrakenDB/KrakenDB_201809/



#2020/8/9
#单独使用concoct进行分箱
cut_up_fasta.py /data2/home/yanxiaoting/project/Golden_Snub-nosed_Monkey/Meta/Faece-8/03.assembly-IDBA/11.assembly.fa -c 10000 -o 0 --merge_last -b contigs_10K.bed > contigs_10K.fa

concoct_coverage_table.py /data2/home/lichunhui/test/metawrap_bin/work_files/assembly_10K.bed /data2/home/lichunhui/test/metawrap_bin/work_files/reads.bam > coverage_table.tsv

concoct --composition_file /data2/home/lichunhui/test/metawrap_bin/work_files/assembly_10K.fa --coverage_file coverage_table.tsv -b concoct_output/

merge_cutup_clustering.py concoct_output/clustering_gt1000.csv > concoct_output/clustering_merged.csv

extract_fasta_bins.py /data2/home/yanxiaoting/project/Golden_Snub-nosed_Monkey/Meta/Faece-8/03.assembly-IDBA/11.assembly.fa concoct_output/clustering_merged.csv --output_path concoct_output/fasta_bins
#还是有问题，分箱到中途报错退出
Traceback (most recent call last):
  File "/public/home/lichunhui/software/miniconda3/envs/metawrap/bin/extract_fasta_bins.py", line 39, in <module>
    main(args)
  File "/public/home/lichunhui/software/miniconda3/envs/metawrap/bin/extract_fasta_bins.py", line 27, in main
    seqs = [all_seqs[contig_id] for contig_id in contig_ids] 
KeyError: 'scaffold_365_1.70484'



#还是继续跑流程。。。
/data2/home/yanxiaoting/project/Lungfish-Mete/No_dormancy/03.mapping

#在NCBI上下载数据
#先下载对应服务器版本的sratoolkit软件
cd /public/home/lichunhui/software/
wget -c https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/sratoolkit.2.10.8-centos_linux64.tar.gz
#解压
tar xzvf sratoolkit.2.10.8-centos_linux64.tar.gz
#下载到指定文件夹下
/public/home/lichunhui/software/sratoolkit.2.10.8-centos_linux64/bin/prefetch -p -O /data2/home/lichunhui/human_stool SRR8427257
#拆包，转换为fastq
cd /data2/home/lichunhui/human_stool/SRR8427257/
/public/home/lichunhui/software/sratoolkit.2.10.8-centos_linux64/bin/fastq-dump SRR8427257.sra


#先用fastqc查看序列质量，再找找有没有专门进行3代测序数据评估的软件
fastqc -o /data2/home/lichunhui/human_stool/00_fastqc/ -t 10 /data2/home/lichunhui/human_stool/SRR8427257/SRR8427257.fastq
#创建新的conda环境
conda create -n 'nanopore' python='3.7'
#安装nanofilt进行质控
NanoFilt /data2/home/lichunhui/human_stool/SRR8427257/SRR8427257.fastq -l 1000 -q 8 --headcrop 40 > /data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq
#查看质控后的序列质量
fastqc -o /data2/home/lichunhui/human_stool/01_nanofilt/ -t 10 /data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq
#使用minimap2比对到人类参考基因组上，去除宿主污染
nohup time sh /data2/home/lichunhui/human_stool/02_minimap2/minimap.sh > /data2/home/lichunhui/human_stool/02_minimap2/minimap2.log 2>&1 &

#提取出没有比对到参考基因组上的结果
samtools view -bu -f 12 -F 256 /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257.bam > /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_unmap.bam

#生成fastq文件
bedtools bamtofastq -i /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_sort.bam -fq SRR8427257_meta.fastq


#canu组装genomesize=4.8m
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh canu.sh 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/canu_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu/SRR8427257.contigs.fasta

#canu组装genomesize=50m
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh canu.sh 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/canu50m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_50m/SRR8427257.contigs.fasta

#canu组装genosize=100m
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh canu.sh 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/canu100m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_100m/SRR8427257.contigs.fasta


#flye组装genomesize=1.2g
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh flye.sh > flye.log 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/flye_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye/assembly.fasta

#flye组装genosize=100m
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh flye.sh > flye_100m.log 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/flye100m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_100m/assembly.fasta

#flye组装genomesize=250m
cd /data2/home/lichunhui/human_stool/04_assemble
nohup time sh flye.sh > flye_250m.log 2>&1 &
#质量评估
python /public/home/lichunhui/software/quast/quast.py -o /data2/home/lichunhui/human_stool/04_assemble/flye250m_quast /data2/home/lichunhui/human_stool/04_assemble/SRR8427257_flye_250m/assembly.fasta


#使用kraken2对质控后的reads进行注释
cd /data2/home/lichunhui/human_stool/03_read_kraken2
nohup time sh read_kraken2.sh > read_kraken2.log 2>&1 &
#提取注释结果
python /public/home/lichunhui/project/lungfish_meta/kraken2_annotation.py /data2/home/lichunhui/human_stool/03_read_kraken2/ /data2/home/lichunhui/human_stool/03_read_kraken2/

#使用kraken2对组装后的contigs进行注释
cd /data2/home/lichunhui/human_stool/05_congtig_kraken2
nohup time sh contig_kraken2.sh > contig_kraken2.log 2>&1 &


#查看统计信息
cd /data2/home/lichunhui/human_stool/02_minimap2
samtools flagstat SRR8427257.bam

#修改flag值再提取没有比对上的序列
samtools view -bu -f 4 SRR8427257.bam > SRR8427257_unmap.bam
#sort排序
samtools sort SRR8427257_unmap.bam -o SRR8427257_unmap_sort.bam
#得到fastq文件
bedtools bamtofastq -i /data2/home/lichunhui/human_stool/02_minimap2/SRR8427257_unmap_sort.bam -fq SRR8427257_meta1.fastq



#2020/08/11
#使用v2.0 eggnog-mapper对粪便数据进行注释
cd /data2/home/lichunhui/human_stool/06_eggnog
nohup eggnog.sh > eggnog.log 2>&1 &

#重新用minimap2进行比对生成sam文件
cd /data2/home/lichunhui/human_stool/02_minimap2
nohup minimap2 -ax map-ont /public/home/renqingmiao2018/project/yak.rumen.metagenome/00.ref/human.genome/GCF_000001405.38_GRCh38.p12_genomic.fna /data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq -t 24 > SRR8427257.sam &

#安装nextpolish
cd software/
git clone --recursive https://github.com/jts/nanopolish.git

#racon抛光
cd /data2/home/lichunhui/human_stool/07_racon
racon -m 8 -x -6 -g -8 -t 12 \
/data2/home/lichunhui/human_stool/01_nanofilt/SRR8427257_nanofilt.fastq \
/data2/home/lichunhui/human_stool/02_minimap2/SRR8427257.sam \
/data2/home/lichunhui/human_stool/04_assemble/SRR8427257_canu_100m/SRR8427257.contigs.fasta
2020/07/29
#!/bin/bash
for i in ~/data/Lungfish-Mete/Dormancy/02.cleandata/*_1.clean.trim.fq.gz
do
echo $i
base=$(basename $i _1.clean.trim.fq.gz)
echo $base
bwa mem -t 10 -M -R '@RG\tID:1-1_BDME202032413-1a\tLB:1-1_BDME202032413-1a\tPL:Illumina\tSM:1-1_BDME202032413-1a' \
/public/home/wangkun/projects/lungfish/06.bionano_scaffold/lungfish_pilon_bionano.fasta \
~/data/Lungfish-Mete/Dormancy/02.cleandata/${base}_1.clean.trim.fq.gz \
~/data/Lungfish-Mete/Dormancy/02.cleandata/${base}_2.clean.trim.fq.gz \
| samtools view -bS > ~/data2/lunfish_meta/Dormancy/03_bwa_map/${base}.bam
done


nohup sh /data2/home/lichunhui/lunfish_meta/Dormancy/bwa_map/bwa_map.sh > /data2/home/lichunhui/lunfish_meta/Dormancy/bwa_map/bwa_map.log 2>&1 &

IDBA安装
#网站https://github.com/loneknightpy/idba
git clone https://github.com/loneknightpy/idba.git
cd ~/lichunhui/software/idba/
./build.sh
./configure
make

IDBA_UD组装,要将paired-end保存在一个fq文件中
bin/fq2fa --merge --filter /data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/1-1_BDME202032413-1a_1.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/1-1_BDME202032413-1a_2.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/1-1_BDME202032413-1a.cl.fq

bin/fq2fa --merge --filter /data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/2-2_BDME202032413-1a_1.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/2-2_BDME202032413-1a_2.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/2-2_BDME202032413-1a.cl.fq

bin/fq2fa --merge --filter /data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/3-1_BDME202032413-1a_1.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/3-1_BDME202032413-1a_2.cl.fq \
/data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/3-1_BDME202032413-1a.cl.fq

idba_ud默认最长只支持reads长度为128的序列，
如果你的reads序列大于128，可修改 src/sequence/short_sequence.h文件中的kMaxShortSequence值
这里我进行了修改，因为我的reads的长度为150

nohup bin/idba_ud -r /data2/home/lichunhui/lunfish_meta/Dormancy/05_meta_data/1-1_BDME202032413-1a.cl.fq \
--maxk 90 --step 10 -o ~/data2/lunfish_meta/Dormancy/IDBA_UD_contigs --num_threads 4 --min_contig 500 >~/idba.log 2>&1 &

#prodigal预测ORF
cd /data2/home/lichunhui/lunfish_meta/Dormancy/07_prodigal/
vim prodigal.sh

#!/bin/bash
prodigal -a 1-1_BDME202032413-1a.orf.faa -i ~/data2/lunfish_meta/Dormancy/06_megahit_zz/1-1_BDME202032413-1a.megahit_contigs.fa \
-f gff -o 1-1_BDME202032413-1a.gff -p single -q -d 1-1_BDME202032413-1a.orf.ffn

#!/bin/bash
for i in ~/data2/lunfish_meta/Dormancy/06_megahit_zz/*.megahit_contigs.fa
do
echo $i
base=$(basename $i .megahit_contigs.fa)
prodigal -a ${base}.orf.faa \
-i ~/data2/lunfish_meta/Dormancy/06_megahit_zz/${base}.megahit_contigs.fa \
-f gff -o ${base}.gff -p single -q -d ${base}.orf.ffn
done

cat *.ffn > dormancy.ffn

cd-hit-est -i dormancy.ffn -o dormancy.geneSet.ffn -n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 1 -T 10


2020/07/30
#R画物种组成的相对丰度图
#加载包
library(ggplot2)
library(readr)
library(tidyr)

#导入数据
dor_data <- read_csv('C:/Users/bbc/Desktop/dor_phylum_abun.csv')
nodor_data <- read_csv('C:/Users/bbc/Desktop/nodor_phylum_abun.csv')

#修改列名
names(dor_data) <- c('rank_code','name','dor_sample1','dor_sample2','dor_sample3')
names(nodor_data) <- c('rank_code','name','nodor_sample1','nodor_sample2','nodor_sample3')

#把宽数据转换为长数据
dordata_melt <- gather(dor_data,key = 'sample',value = 'ratio',dor_sample1:dor_sample3)
nodordata_melt <- gather(nodor_data,key = 'sample',value = 'ratio',nodor_sample1:nodor_sample3)

#把两个表的数据结合
total_data <- rbind(dordata_melt,nodordata_melt)

#作图
data_plot <- ggplot(data = total_data,aes(x=total_data$sample,y=total_data$ratio,fill=total_data$name))+
  geom_bar(stat = 'identity',position = 'stack',color='white',width = 1)+
  theme(panel.grid =element_blank()) +   #删除网格线
  theme(axis.text = element_blank()) +   #删除所有刻度标签
  theme(axis.ticks = element_blank()) +   #删除所有刻度线
  theme(panel.border = element_blank()) #删除外层边框

data_plot

#在R中，如果在与表的合并过程中遇到有一列在两个表中同名，但是值不同，合并的时候又都想保留下来，就可以用suffix给每个表的重复列名增加后缀
#例子：
data_join7 <- dplyr::left_join(df1,df2,by='x',suffix=c('.1','.2'))

2020/7/31
#做秩和检验
import pandas as pd
import scipy.stats as stats
out_file=open(r'C:\Users\bbc\Desktop\wilcox_ranksum.txt','a+')
df=pd.read_table(r'C:\Users\bbc\Desktop\phylum_abun.txt')
out_file.write('phylum_name\tp_value\n')
for i in range(len(df)):
    name=df.iloc[i,0]
    result=stats.ranksums(df.iloc[i,1:4],df.iloc[i,4:7])
    out_file.write(name+'\t'+str(result[1])+'\n')
out_file.close()
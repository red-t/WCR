## LACHESIS 的安装流程
### 安装 samtools-0.1.18
去 SOURCEFOREGE 的 页面找到 [samtools-0.1.18](https://sourceforge.net/projects/samtools/files/samtools/0.1.18/) 的压缩包，下载下来。
```
tar xf samtools-0.1.18.tar.bz2
cd samtools-0.1.18
make -j 20
```

### 安装 boost 1.59.0
按照 [官方文档](https://www.boost.org/doc/libs/1_59_0/more/getting_started/unix-variants.html) 的指引进行安装：
```
wget http://sourceforge.net/projects/boost/files/boost/1.59.0/boost_1_59_0.tar.bz2
tar --bzip2 -xf boost_1_59_0.tar.bz2
cd boost_1_59_0
mkdir opt
./bootstrap.sh --prefix=/data/tusers/zhongrenhu/Software/boost_1_59_0/opt
./b2 install

# 创建 LACHESIS 所需目录
mkdir -p opt/stage
ln -s /data/tusers/zhongrenhu/Software/boost_1_59_0/opt/lib /data/tusers/zhongrenhu/Software/boost_1_59_0/opt/stage/lib
```

### 下载LACHESIS 并解压
```
curl -o LACHESIS.zip https://codeload.github.com/shendurelab/LACHESIS/legacy.zip/master
unzip LACHESIS.zip
mv shendurelab-LACHESIS-2e27abb LACHESIS
cd LACHESIS
```

### 修改文件
```
# LACHESIS/src/include/gtools/ 目录下的 SAMStepper.h & SAMStepper.cc
# 对两个文件进行相同的修改

#include <bam/sam.h> ----> #include </data/tusers/zhongrenhu/Software/samtools-0.1.18/sam.h>
```

### 声明环境变量 LD_LIBRARY_PATH
```
将 /data/tusers.ds/zhongrenhu/Software/boost_1_59_0/opt/lib/ 添加到环境变量 LD_LIBRARY_PATH 中，并且在 ~/.profile 文件中进行声明
```

### 声明环境变量 LACHESIS_BOOST_DIR 以及 LACHESIS_SAMTOOLS_DIR
```
# LACHESIS_BOOST_DIR 是 boost 的安装路径
export LACHESIS_BOOST_DIR=/data/tusers/zhongrenhu/Software/boost_1_59_0/opt/

# LACHESIS_SAMTOOLS_DIR 也是 boost 的安装路径(它直接在原本的目录下编译了)
export LACHESIS_SAMTOOLS_DIR=/data/tusers/zhongrenhu/Software/samtools-0.1.18
```
  
### 运行 configure 生成 Makefile
```
./configure --with-samtools=/data/tusers/zhongrenhu/Software/samtools-0.1.18 --with-boost=/data/tusers.ds/zhongrenhu/Software/boost_1_59_0/opt
```

### 修改 LACHESIS/src/include/gtools/Makefile
```
INCLUDES= ----> INCLUDES= -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)/include
```

### 修改 LACHESIS/src/include/markov/Makefile
```
INCLUDES= ----> INCLUDES= -I$(LACHESIS_SAMTOOLS_DIR) -I$(LACHESIS_BOOST_DIR)/include

BOOST_LIBS= -lboost_system -lboost_filesystem -lboost_regex ----> BOOST_LIBS = -L$(LACHESIS_BOOST_DIR)/stage/lib -lboost_system -lboost_filesystem -lboost_regex
```

### 运行 make 进行编译
```
make
```

### 使用测试数据尝试运行
```
Lachesis INIS/test_case.ini
```

### 使用 Container 技术进行安装
```
# 安装 aakashsur 的镜像
singularity pull --name lachesis_aakashur docker://aakashsur/lachesis

# 安装 zitsens 的镜像
singularity pull --name lachesis_zitsen docker://zitsen/lachesis
```

----

## LACHESIS 程序 & 参数 & 结果 相关

### GenomeLinkMatrix
这是 **一种数据结构**，包含一个大小为 `_N_bins x _N_bins`，对应 `Hi-C heatmap` 的 `2-D Boost UBLAS compressed matrix`，这个数据结构所包含的方法如下：
```
- Load in Hi-C data from a set of SAM files. In a de novo GLM, each contig becomes a bin.

- Using the SAM files, create a matrix of link density between each pair of bins.

- Cluster these bins based on their link density. This creates a ClusterVec object. This is the main algorithm!

- Analyze and validate the clusters.
```
简单来说，`GenomeLinkMatrix` 这个数据结构就是用于做 **“group clustering”** 的。而其中 **需要注意的点** 有：
```
Short and/or repetitive contigs can screw up the clustering algorithm. It's best to leave these contigs out of the clustering process, though they can optionally be assigned to clusters after the clusters have been made (see SetClusters).

Use SkipShortContigs() to mark contigs for skipping if they are below a given size threshold.

Use SkipContigsWithFewREs() to mark contigs for skipping if they don't have enough RE (restriction endonuclease) sites.

Use SkipRepeats() to mark contigs for skipping if they are repetitive - i.e., they have a normalized number of Hi-C links that is much greater than average.
```



### heatmap.txt 文件
由 `GenomeLinkMatrix.cc` 当中 `WriteFile` 函数生成，提供给 `heatmap.R` 脚本绘制热图，


### CLM 文件
LACHESIS 程序运行之后会生成 `CLM_file` , 也就是 `chromosome link matrix file`, 由 `ChromLinkMatrix.cc` 当中的 `WriteFile` 函数生成。其主要内容如下：
```
WriteFile: Write the data in this ChromLinkMatrix to file CLM_file. The output format is a long tall table with (4*_N_contigs^2) rows and three (or more) columns: X, Y, Z, [data]. X and Y are bin IDs (bin ID = 2 * contig ID + contig orientation). Z is an integer representing the amount of data in bin [X,Y]. If heatmap = false (default), then following Z is actually a tab-separated list of all the integers in bin [X,Y]. This format is designed for easy input to ChromLinkMatrix::ReadFile() (if heatmap = false), or to R (if heatmap = true). If this is a de novo CLM, also write the contig lengths to an auxiliary file.
```

CLM 文件名格式为 `groupN.CLM`, 内容示例如下：
```
# ChromLinkMatrix file - see ChromLinkMatrix.h for documentation of this object type
# Species = wcr
# De novo CLM? true
# N_contigs = 1559
# contig_size = 0 (ignored if a contig_lens_file is supplied)
# contig_lens_file = lachesis_scaffold/round1/cached_data/group0.CLM.lens
# RE_sites_file = lachesis_scaffold/round1/cached_data/group0.CLM.RE_sites
# heatmap = false
# SAM files used in generating this dataset: HiC/HiC-WCR-ov2-10cycle.REduced.paired_only.bam
X       Y       Z       link_lengths
2       2       2       152    76
10      7       1       2431    106
```
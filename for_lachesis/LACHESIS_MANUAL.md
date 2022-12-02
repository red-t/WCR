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
git@github.com:shendurelab/LACHESIS.git
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

# 安装 huzr14830 的镜像 (目前的方案)
singularity pull --name lachesis_huzr14830 docker://huzr14830/lachesis
```

----
## LACHESIS 输入文件需求

### Draft assembly fasta
对应参数 `DRAFT_ASSEMBLY_FASTA` 以及变量 `_draft_assembly_fasta`
1. 如果是第一次运行，则需要提供完整的 fasta 文件，例如 `draft.fasta`
2. 如果不是第一次运行，可以不需要 draft fasta，但需要保证相应文件的 `names file` & `RE_counts file` 存在，例如 `draft.fasta.names` & `draft.fasta.counts_AAGCTT.txt`

### Reference assembly fasta
对应参数 `REF_ASSEMBLY_FASTA` 以及变量 `_ref_assembly_fasta`
1. 如果是第一次运行，则需要提供完整的 fasta 文件，例如 `ref.fasta`
2. 如果不是第一次运行，可以不需要 ref fasta，但需要保证相应文件的 `names file` 存在，例如 `ref.fasta.names`

### Hi-C alignments
对应参数 `SAM_FILES` 以及变量 `_SAM_files`，用于构造数据结构 `GenomeLinkMatrix` & `ChromLinkMatrix`
1. 不论是第几次运行，都需要提供完整的 BAM 文件
2. 如果不是第一次运行，提供缓存文件 `all.GLM`，可以加快程序的运行速度

### BLAST output
对应参数 `BLAST_FILE_HEAD` 以及变量 `_BLAST_file_head`，用于构造数据结构 `TrueMapping`，主要包含四个向量: `_target, _start, _stop, _qual_alignability, _qual_specificity`
1. 如果是第一次运行，需要提供完整的 `BLAST output`，例如 `draft.NNN.blast.out`
2. 如果不是第一次运行，可以不需要 `BLAST output`，但需要保证缓存文件 `TrueMapping.assembly.txt` 存在
3. 数据结构 `TrueMapping` 的第一个用处是在 clustering 之前，根据 contigs 与 reference genome 的比对位置(中点位置)，将 contigs 进行排序，然后绘制 heatmap_preprocess.jpg
4. 第二个用处是在 clustering 完成之后，ordering 开始之前，用于统计 `total_aligned, mis-clustered, informative mis-clustered` contigs 的数量；并且筛选能够比对到 reference genome 上的 contigs，用于绘制 `clustering_dotplot`
5. 第三个用处是在 ordering 完成之后，对一个 ordered cluster 当中的 contigs，标记它们在 reference genome 上比对的位置(中点)，并且根据它们在 order 当中的位置绘制 `clm.N.dotplot.txt`
6. 第四个用处是 ordering 结束之后，评估 clustering 的效果，与4类似
7. 第五个用处是 ordering 结束之后，评估每个 trunck/full order 的效果，包括：根据比对结果标记contig是否“aligned”；根据前后两个contig比对的chrom是否相同标记两个contig是否为“mismatch”；根据前后两个contig比对结果的相对位置信息以及在order中的相对位置标记后一个contig是否有“orient error”；根据前后三个contig比对结果的相对位置信息以及在order中的相对位置标记中间的contig是否有“order error”；

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

1. 一开始通过 `GenomeLinkMatrix::GenomeLinkMatrix`  函数构建 GenomeLinkMatrix，写入到 `outdir/cached_data/all.GLM` 文件， 该函数只涉及 `input BAMs/SAMs`、`species`、`RE_sites_file`，不涉及其他的 `INI` 参数，因此往后如果要继续跑其他的 round，可以将`cached_data/all.GLM` 拷到新的 `outdir` 下
	1. `GenomeLinkMatrix::GenomeLinkMatrix` **从 SAM/BAM 文件的 header 部分来判断总共有多少个 bin —— Each target contig is a bin.** 
	2. 并且使用原本的顺序来排列每个 bin (即 header 中的顺序，并不进行重排)
	3. 记录每个 bin 的长度
	4. 通过 `GenomeLinkMatrix::LoadRESitesFile` 函数从 `draft_assembly/wcr_v0.3.0.fasta.counts_AAGCTT.txt` 读取每个 bin (contig) 当中的 RE site counts，并且每个值都加上 1 (防止作为分母时，有 0 值存在)
	5. 初始化一个长度为 `_N_bins` 的布尔值向量，用于记录每个 contig 是否该被 “skip”，初始值都是 false
	6. 初始化一个 `N*N` 的空白 matrix (值都为0) —— `_matrix`
	7. 遍历每个 BAM 文件的每一条 read，如果符合以下条件：
		- read 和 mate read 比对的位置相同 (但似乎没考虑是比对到相同的 contig 上？)
		- mate read 比对不上
		- read mapping quality == 0
		- read 和 mate read 比对到相同的 contig 上
		则跳过这条 read，否则 `_matrix` 相应位置的值 +1 （比如两条 read 分别比对到 contig i & j，则 `_matix(i,j) += 1; _matix(j,i) += 1;`）
	9. 调用 `GenomeLinkMatrix::WriteFile` 函数将 `_matrix` 写出为 `outdir/cached_data/all.GLM`:
		- 写成三列的格式（X Y Z），X & Y 为矩阵坐标，Z为值
		- X == Y 不写出（即对角线上，“contig 自己跟自己的 link”）
		- Z == 0 不写出
2. 如果已经存在 `all.GLM`，则调用 `GenomeLinkMatrix::ReadFile` 进行读取，构建 `_matix`：
	1. 从 `all.GLM` 的 header 部分读取一些 meta 信息
	2. 初始化 `_matrix`，和上一步的 step 6 相同
	3. 读取 RE site counts，与上一步的 step 4 相同
	4. 如果不是 header line，则根据 X & Y & Z 填充矩阵 `_matrix`
3. 调用 `GenomeLinkMatrix::NormalizeToDeNovoContigLengths` 对 `_matix` 的值进行 normalization：
	- 判断使用 RE site count 还是 contig length 进行 normalize （默认是前者）
	- 找出 longest_contig （最长 RE 或 最大 length），计算 `longest_squared = longest_contig * longest_contig`
	- 对矩阵的每一个 `_matix(i,j)` 都进行 normalize：`_matrix(i,j) = _matrix(i,j)*longest_squared/(len(i)*len(j))*`
4. 调用 `GenomeLinkMatrix::SkipContigsWithFewREs` 根据 RE sites 数量标记后续该 "skip" 的contig：
	- min_N_REs 可以在 INI 文件中定义
	- 如果某条 contig 的 RE sites 小于该值，第一步 step 5 初始化的 bool 向量 `_contig_skip` 相应位置则标记为 true
	- 计算所有被 skip 的 contig 的平均长度 & 平均 RE site counts
5. 调用 `GenomeLinkMatrix::SkipRepeats` 标记 “可能来自于重复区域” 的 contig：
	- 如果 contig 的 links 数是平均 links 数的 `repeat_multiplicity` 以上，则被认为是 "repeat" 的
	- `repeat_multiplicity` 对应 INI 文件当中的 `CLUSTER_MAX_LINK_DENSITY` 参数
	- 同样是通过 `_contig_skip` 进行标记
	- 标记完之后，计算被 skip 的 contig 的平均长度

### heatmap.txt 文件
由 `GenomeLinkMatrix.cc` 当中 `WriteFile` 函数生成，提供给 `heatmap.R` 脚本绘制热图，而 `heatmap.R` 的内容也十分简单，其中硬编码了输入文件为 `heatmap.txt (X、Y 为 matrix 坐标，log10(Z+1) 为相应位置的值)`

在整个 LACHESIS 流程中，会调用两次 `DrawHeatmap` 函数，所生成的文件貌似都叫 `heatmap.txt` ?

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

### MakeWholeAssemblyHeatmap 函数
这是 `Reporter.cc` 当中定义的一个函数，在最后绘制一张 heatmap，用于展示 scaffold 的 “intra-cluster” & “inter-cluster” Hi-C link density，从而评估 scaffolding 的质量。

1. 该函数首先根据长度筛选出 top 1000 contigs
2. 根据 contigs 在 cluster 中的 order，对这 top 1000 contigs 进行排序，每个 contig 占一个 bin。`heatmap.chrom_breaks.txt` 用于记录 cluster 的边界 ---- 即一个 cluster 当中包含 top 1000 contigs 的数量(累加)
3. 根据长度(或 RE sites 数量)对这 top 1000 contigs 之间的 Hi-C links 数量进行 normalization
4. 经过 normalization 之后的数据再除以 1000 ("Rescale the data to put it into a more intuitive range for visual clarity.")
5. 将 matrix 以 `"X\tY\tZ\n"` 的格式写出，并调用 `heatmap.MWAH.R` 绘制热图
# 10 ATAC-seq 染色质可及性分析：DAR、PLG 顺式元件与 RUNX2 TFBS

> **前置依赖**：01（SRA 数据）。**对应原文**：Fig. S5A–C、Fig. S11；Methods "Data acquisition"（PRJNA725028, 44 样本）与结果部分 ATAC 段落。
> **产出**：MASH vs MASLD 的 DAR（上 145,812 / 下 10,567）、PLG 的 7 个顺式调控元件、27 个 RUNX2 TFBS、TSS 距离分布、模型基因启动子 ±2 kb 可及性。

---

## 1. ATAC-seq 原理速览

ATAC-seq（Assay for Transposase-Accessible Chromatin）用高活性 Tn5 转座酶切入**开放染色质**并插入测序接头：开放区域（启动子/增强子）被切割概率高 → 测序富集 → peak。因此 peak 强度 ≈ 局部染色质可及性 ≈ 顺式调控活性。TF 结合位点（TFBS）处常表现为 peak 内的 **footprint**（被 TF 占据而免受切割的小凹陷），但 footprint 检测信噪比低，原文采用 **motif 富集 + peak 注释** 路线（更稳健）。

---

## 2. 全流水线（shell + Python）

```bash
# code/10_atac.sh  （每样本；44 样本建议并行/集群）
# 0) 质控与接头去除
fastp -i SRR_R1.fastq.gz -I SRR_R2.fastq.gz -o R1.clean.gz -O R2.clean.gz \
      --detect_adapter_for_pe -w 8

# 1) 比对（hg38）
bowtie2 -X 2000 --no-discordant --no-mixed -p 8 -x hg38_index \
        -1 R1.clean.gz -2 R2.clean.gz -S aln.sam
# 2) 过滤：MAPQ>=30、proper pair、去线粒体/黑名单
samtools view -b -q 30 -f 2 aln.sam | samtools sort -o srt.bam
samtools index srt.bam
bedtools intersect -v -abam srt.bam -bbg hg38.blacklist.bed > clean.bam
# 3) 去重
picard MarkDuplicates I=clean.bam O=dedup.bam M=dup.txt REMOVE_DUPLICATES=true
# 4) Tn5 偏移（+链 +4 / -链 -5）后 call peak
MACS2 callpeak -t shifted.bam -f BAMPE -g hs -n sample \
      --nomodel --shift -100 --extsize 200 --keep-dup all -q 0.01
```

> 参数解释：ATAC 的插入片段中位 ~180 bp（核小体游离），`--nomodel --shift -100 --extsize 200` 是 ATAC 的 MACS2 标准设置（把切割位点还原为转座事件中心）。BAMPE 利用双端插入片段长度。

### 2.1 组级 peak 合并与计数矩阵

```bash
# 合并所有样本 peak 为统一坐标集（union peaks）
cat */*.narrowPeak | sort -k1,1 -k2,2n | bedtools merge > union_peaks.bed
# 每样本在 union peaks 上计数
featureCounts -a union_peaks.bed -b -o counts.txt -F SAF ... # 或 bedtools multicov
```

### 2.2 差异可及区（DAR）

```r
# code/10_dar.R
library(DESeq2)
dds <- DESeqDataSetFromMatrix(cnt, meta, ~ disease)     # MASH vs MASLD
dds <- dds[rowSums(counts(dds) >= 5) >= 3, ]
dds <- DESeq(dds)
res <- results(dds, contrast = c("disease","MASH","MASLD"), alpha = 0.05)
sum(res$padj < 0.05 & res$log2FoldChange > 0)    # 期望 ≈145,812（上调）
sum(res$padj < 0.05 & res$log2FoldChange < 0)    # 期望 ≈10,567（下调）
fwrite(as.data.frame(res), "results/FigS5A_DAR.csv")
```

> **数量级核对**：14.6 万上调 DAR 属于"宽松阈值 + union peak 数很大"的常见结果；若复现得 5–10 万，检查①peak 合并策略②DESeq2 的独立过滤③阈值（原文未给 padj/FC 细节，复现时记录所用 cutoff 并以"上/下比 ≈ 13.8:1"作为方向性核对）。

---

## 3. PLG 顺式调控元件与 RUNX2 motif（Fig. S5B–C）

```r
# 1) 取 PLG 基因座 ±扩展窗口（如 TSS±100 kb 或基因体+50 kb）内的 DAR/peak
library(GenomicRanges); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
plg <- genes(TxDb)[ "PLG" ]  # 或用 orgDb 取坐标
win <- resize(plg, width = width(plg) + 2e5, fix = "center")
plg_peaks <- subsetByOverlaps(union_gr, win)

# 2) 这些元件是否落在"上调簇"：与 DAR 上调集取交集
#    原文：7 个 PLG 相关顺式元件全部映射到 up-cluster
sum(plg_peaks %in% up_dar_gr)   # 期望 7/7

# 3) motif 富集：对 7 个元件序列做 RUNX2 motif 扫描/富集
#    工具：HOMER findMotifsGenome.pl 或 MEME-ChIP / AME（JASPAR RUNX2 motif MA0684）
#    输出 27 个 RUNX2 TFBS（原文 Supplementary Table 10）
```

```bash
# HOMER 路线
fastaFromBed -bed plg7_elements.bed -fo plg7.fa -tab > /dev/null
findMotifsGenome.pl plg7_elements.bed hg38 homer_out/ -size 200 -mknown \
      -p 8
# 或用 FIMO(MEME) 逐碱基扫描得位点级 TFBS 坐标，计数=27
fimo --thresh 1e-4 JASPAR_RUNX2.meme plg7.fa
```

### 3.1 TFBS 功能注释与 TSS 距离分布（Fig. S5B–C）

```r
library(ChIPseeker)
ann <- annotatePeak(tfbs_gr, TxDb = TxDb, tssRegion = c(-3000, 3000))
# 原文：RUNX2 主要结合 PLG 的 promoter 与 intron，少量远端 intergenic
plotAnnoBar(ann)
d <- distanceToTSS(ann)
# 直方图：多数 TFBS 位于 TSS 上下游 10–100 kb（Fig.S5C）
hist(abs(d), breaks = c(0,1e3,1e4,1e5,1e6))
```

---

## 4. 模型基因启动子 ±2 kb 可及性（Fig. S11，IGV）

```r
# 对 8 个模型基因（ARHGAP22, DHRS9, DTNA, ITGAX, RUNX2, TRPM2, TESK2, AUTS2）
# 取 TSS±2 kb 窗口内 MASH vs MASLD 的平均可及性（标准化计数），wilcox 检验
wins <- promoters(genes_of_8, upstream = 2000, downstream = 2000)
acc  <- summarizeOverlaps(wins, bamlist, mode = "Union")
# MASH 显著高于 MASLD（原文）
# IGV 可视化：生成 bigWig（deeptools bamCoverage --normalizeUsing RPKM）后载入 IGV 截图
```

```bash
bamCoverage -b dedup_MASH_rep1.bam -o MASH.bw --binSize 10 --normalizeUsing RPKM -p 8
# IGV: 载入 hg38 + MASH.bw + MASLD.bw，定位 8 基因 TSS±2kb 截图 → Fig.S11
```

---

## 5. 陷阱与导师提示

1. **批次/供体效应**：44 样本若来自多中心，DESeq2 设计式应加 `~ center + disease`；否则 DAR 数会膨胀（这可能正是原文 14.6 万的来源之一）。复现报告中给出两种设计式的 DAR 数对比。
2. **union peak 的宽度异质性**：合并 peak 长度不一会影响计数；可用固定宽度（peak 中心 ±250 bp）重计数做敏感性分析。
3. **motif 数据库版本**：JASPAR 2024 vs HOMER 内置版本的 RUNX2 motif 略有差异，TFBS 计数（27）对阈值敏感；记录 p/thresh 与数据库版本。
4. **ATAC 与表达的衔接**：DAR 上调 ≠ 基因上调；把 PLG 元件的可及性变化与 09 的 PLG 表达变化并列报告，才是"表观→转录"的完整证据链。

---
title: Sinz_AC_2019
tags: cleavable cross-linker DSBU
notebook: DSSO
---
We present a cross-linking/mass spectrometry
workflow for performing proteome-wide cross-linking analyses
within 1 week.
文章开发一种新的cxms方法, 在一周内完成蛋白组水平的交联实验。
#  文章亮点：
1. MetoX 进化到2.0 版本，实现全自动的搜索过程，支持了蛋白组水平的搜索。
2. 在果蝇的胚胎lysate水平鉴定到7436交联对，其中1611对分子间的交联
3. 多种方法对交联结果的进行验证，包括1.陷阱库，搜索库混入E.coli蛋白；2. 用已有晶体结构验证距离约束，3. 与蛋白相互作用数据库string进行对照，看分子间的交联的蛋白是否在string中也有相互作用。

## 技术
**使用仪器**： QE-plus(速度略慢)
**分馏技术**：size-exclusion chromatography (SEC)（superdex 30 Increase 10/300 GL）
收集了18个组分
**液相方法**：peptides were eluted
and separated using a linear gradient from 3% to 42% B (with solvent B: 0.1% (v/v) formic acid and 85% (v/v) acetonitrile) with a constant flow rate of 300 nl/min over $\color{red}{360}$ min, 42% to 99% B (5 min) and 99% B (5 min).
6个小时的分离梯度！！！（也是拼了）
**质谱方法**：
1. 使用step hcd (27, 30, 33)
2. 前三个（1-3)馏分碎裂的母离子价态3-7， 后面馏分（4-18）价态为3-4.

**数据分析：**
peptides’ length: 5 to 30
modification： c+57, m+16
link site: Lys and N-termini (这时候不考虑STY了！推测原因是搜索空间大，软件速度慢。)
过滤条件：minimum peptide score: 10, precursor mass accuracy: 4 ppm, product ion mass accuracy: 8 ppm (performing mass recalibration, average of deviations), signal-to-noise ratio: 1.5, precursor mass correction activated, prescore cut-off at 10% intensity, FDR cut-off: 1%, and minimum score cut-off: 70.（Sinz实验室在母离子和碎片离子误差要求上确实很严格。）
## MeroX 2.0的新功能是什么？
1. proteome-wide mode （plink早就支持）
2. Peptide score based filtering
3. Decoy-sequences now calculated at the peptide level. decoy peptides are generated on the peptide level by shuffling the peptide sequence of an identified candidate. （不了解与主流的蛋白水平reverse，究竟有多大差别）
4. FDR estimated separately for interprotein, intraprotein, intrapeptide and dead-end cross-
links. (plink早就这样做了)
5. Scoring modified to include all recorded subscores (peptide scores, prescore, deltascore)
6. Precursor mass correction implemented （pparse的功能）

总的来说，这篇文章是sinz提升DSBU影响力，赶超DSSO的重要文章，对于dsso来说，前有刘凡的两篇nature 子刊文章，后有yu haiyuan的maxlinker，不断刷新鉴定记录至9,319。本篇文章的7436对交联，虽然不及maxlinker的记录，但分子间的1611交联比maxlinker的1268也多出了不少。我推测sinz对这篇文章的期待至少是nature 子刊级别的，可最终发在了ac上。下面我理一下这篇文章不太好的方面。
##文章糟点
1. 从创新性来讲，无论是分离，质谱方法和软件上都没有特别出彩的地方。
2. 实验可重复性差。  只有15.4%的结果在三次生物重复中都出现了。（fig3a）
3. 可疑的灵敏度。 在分别用整个果蝇的蛋白组（21973蛋白）和用MaxQuant鉴定到的蛋白（9535个蛋白）作为搜索数据库，使用21973个蛋白搜索得到的结果居然更多。如果把这个实验也当作陷阱库实验的话，利用MaxQuant鉴定到的蛋白作为数据库的搜索结果中有4.2%的交联对落在了陷阱库里。
4. 使用宽松的距离约束获得更好的结构吻合度。 DSBU的交联剂臂长12.5 Å，加上k侧链6x2 Å，是24.5Å，文中把最大约束放到了32 Å.

##我印象深刻的一点：
作者使用陷阱库：在搜索的fasta里加入了E.coli的蛋白，来评估搜索结果的准确度。结果是在鉴定到了29957谱图中有325张谱落在了E.coli的蛋白组里，其中103是与果蝇蛋白组共享的肽段，另外的222全是分子间的交联。说明分子间交联的假阳性还是偏高。

#交联鉴定的validation
最后想谈一下交联鉴定的validation, 我理了一下目前大家都用到的validation方法：
1. 陷阱库方法，搜索数据库加入别的物种的蛋白序列。这是一种最为简单的方法。可以衡量鉴定结果是否可疑。MaXlinker文章5000蛋白加入20000蛋白作为诱饵库，这里是10000蛋白加入5000蛋白作为诱饵库。问题是加入多大规模的诱饵库才是合适的？
2. 晶体结构。拿鉴定到的交联对，map到晶体结构里，看是否小于交联剂可伸展的最大距离。但是大家通常都在放大这个距离约束，让自己的结果看起来很好。
3. 与蛋白相互作用数据库进行对照。这种方法真的很无力，因为即便两个一致也不能说明鉴定可靠。那些不一致的，作者会说这是交联新发现的相互作用蛋白（无语）。
4. 同位素标记肽段或者交联剂。成对的母离子共洗脱，一级和二级都可以相互检验。
5. 人工检查谱图。这个主观因素大，没有严格的标准。

以上这些方法中，我认为只有4 同位素标记的方法是能够规模化验证谱图的方法。其他几种都不直接，不有效（还会人为放宽约束）。这也是正式MaXlinker文章和这篇文章欠缺的。这两篇文章虽然给了一个很惊人的鉴定记录，但没有发到高水平杂志的，缺少的就是有效的自证清白方法。MaXlinker文章倒是可以采用贺老师之前提到的SILAC方法来做validation。我其实有点小庆幸，KArGO文章没有被问到“如何证明鉴定到的交联是对的”这个问题？ 如果真被问到了，那么合成同位素标记的交联剂或者使用同位素标记的蛋白会是首要的选择。

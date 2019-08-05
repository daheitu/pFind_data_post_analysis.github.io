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
2. 1-3馏分碎裂的母离子价态3-7， 馏分4-18价态为3-4.

**数据分析：**
peptides’ length: 5 to 30
modification： c+57, m+16
link site: Lys and N-termini (这时候不考虑sty了！)
filter：minimum peptide score: 10, precursor mass accuracy: 4 ppm, product ion mass accuracy: 8 ppm (performing mass recalibration, average of deviations), signal-to-noise ratio: 1.5, precursor mass correction activated, prescore cut-off at 10% intensity, FDR cut-off: 1%, and minimum score cut-off: 70.

## MeroX 2.0的新功能是什么？
1. proteome-wide mode （plink早就支持）
2. Peptide score based filtering
3. Decoy-sequences now calculated at the peptide level. decoy peptides are generated on the peptide level by shuffling the peptide sequence of an identified candidate. （不了解与主流的蛋白水平reverse，究竟有多大差别）
4. FDR estimated separately for interprotein, intraprotein, intrapeptide and dead-end cross-
links. (plink早就这样做了)
5. Scoring modified to include all recorded subscores (peptide scores, prescore, deltascore)
6. Precursor mass correction implemented （pparse的功能）

##文章糟点







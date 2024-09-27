## 表观遗传学之 WGBS

在说 WGBS 之前，先来说说 DNA 甲基化。DNA 甲基化是最早被发现、也是研究最深入的表观遗传调控机制之一。广义上的 DNA 甲基化是指 DNA 序列上特定的碱基，在 DNA 甲基转移酶（DNA methyltransferase，DNMT）的催化作用下，以 S-腺苷甲硫氨酸（S-adenosyl methionine，SAM）作为甲基供体，通过共价键结合的方式获得一个甲基基团的化学修饰过程。这种 DNA 甲基化修饰可以发生在胞嘧啶的 C-5 位、腺嘌呤的 N-6 位以及鸟嘌呤的 G-7 位等位点。

<!--more-->

一般研究中涉及的 DNA 甲基化主要指发生在 CpG 二核苷酸中胞嘧啶上第 5 位碳原子的甲基化过程，其产物称为 5-甲基胞嘧啶（5-mC），是植物、动物等真核生物 DNA 甲基化的主要形式，也是哺乳动物 DNA 甲基化的唯一形式。由于 DNA 甲基化是一种相对稳定的修饰状态，可随 DNA 的复制过程遗传给新生的子代 DNA，因此是一种重要的表观遗传机制。

因此，5-甲基胞嘧啶（或称甲基组）在全基因组分布备受关注。而全基因组甲基化测序（Whole Genome Bisulfite Sequencing，WGBS）是通过 Bisulfite 处理将基因组中未发生甲基化的 C 碱基进行转换用以区分具有甲基化修饰的 C 碱基，并结合高通量测序技术判断 CpG/CHG/CHH 位点是否发生甲基化的一种方法。已成功应用于真核系统发育树各分支、多个物种的甲基组分析，以及人类胚胎干细胞、诱导多能干细胞、外周血单个核细胞、结肠癌细胞等甲基组的分析中。这些 WGBS 数据带来了许多其他方法无法获得的新发现。随着测序成本的降低，WGBS 越来越成为研究的首选方法。然而，传统的 WGBS 方法对低投入量（low input）样本存在较大困难。随着甲基化应用场景的不断增加，无论从胚胎发育研究到肿瘤早期筛查的临床应用，低投入量样本中甲基化建库的需求越来越大。

![](https://images.yuanj.top/202404051232054.png)

<div style="text-align:center; margin-bottom:20px;">
<em>Application of bisulfite treatment in whole genome bisulfite sequencing to convert unmethylated cytosine, not 5-methylcytosine, to uracil. During amplification by polymerase chain reaction, uracil is converted to thymine.</em>
</div>

## 传统的 WGBS

虽然加权基因组亚硫酸氢盐测序（WGBS）已被证实是一项强大的技术，但其传统方法在实际应用中面临一些限制：

- 文库制备标准方案需要大量 input DNA（最少 200-500 ng，最多 5μg），但许多甲基化研究只有有限数量的样本
- 亚硫酸盐转化过程中，大约 90%的 DNA 会被破坏，需要微克级别的 DNA
- 亚硫酸盐处理后，DNA 严重降解，导致无法有效扩增
- 未甲基化的胞嘧啶区域更易断裂，影响测序覆盖均匀性
- 断裂片段两端缺乏完整接头，造成部分 DNA 片段无法成为有效文库分子
- 接头连接与亚硫酸盐处理顺序影响效率
- 亚硫酸氢盐处理本身回收率并不低，但由其引起的 DNA 降解会导致文库制备效率低下和可用于后续分析的完整文库分子数量减少
- 不能区分 DNA 甲基化（5mC）和 DNA 羟甲基化（5hmC）

## 新型 WGBS

为了解决传统 WGBS 所带来的问题，许多的技术人员也在传统技术流程的基础上做出了很多优化。亚硫酸盐后接头标记（PBAT-WGBS）和标签 WGBS（T-WGBS）通过改善 WGBS 文库制备来提高 WGBS 的效率。

### PBAT-WGBS

Post-Bisulfite Adaptor Tagging (PBAT) 是一种先转化后建库的方法。在文章 [Amplification-free whole-genome bisulfite sequencing by post-bisulfite adaptor tagging](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458524/) 中首次提出。

研究者开发了一种基于 PBAT 策略的新型 WGBS 方法，该方法的效率极高，并能避免亚硫酸氢盐引起的测序模板断裂问题。他们认为，省略 DNA 连接和凝胶纯化步骤也有助于提升效率。该方法能够在仅使用微量 DNA 的情况下，对哺乳动物进行 WGBS 分析。值得一提的是，这一过程中无需进行任何全局扩增。研究者成功地利用仅 100 ng 星形胶质细胞 DNA 为起始材料，实现了对小鼠基因组的 21.1 倍无扩增覆盖。这与以往哺乳动物 WGBS 研究形成了显著差异，后者需要使用微克级 DNA 并进行全局 PCR 扩增。由于广泛的扩增往往会导致覆盖偏差，PBAT 策略有望成为传统方法的有效替代，尤其是在 DNA 量非常有限且需要更多 PCR 循环的情况下。

由于亚硫酸盐处理的 DNA 是单链的，PBAT 的实现需要一种有效的方法将接头序列连接到单链 DNA（ssDNA）上。研究者们需要在单链 DNA 上做研究。一种被称为随机启动介导的 PBAT，即 random priming-mediated PBAT（rPBAT），的方法被开发出来，通过两轮随机引物扩增将经过亚硫酸盐处理后的样本连接上完整的接头从而形成完整的文库。

![](https://images.yuanj.top/202404051257070.png)

<div style="text-align:center; margin-bottom:20px;">
<em>(A) Schematic of the conventional WGBS protocols. Bisulfite treatment follows adaptor tagging, thereby leading to bisulfite-induced fragmentation of adaptor-tagged template DNAs. (B) Schematic of PBAT strategy. Bisulfite treatment precedes adaptor tagging, thereby circumventing bisulfite-induced fragmentation of adaptor-tagged template DNAs. (C) Random priming-mediated PBAT method. Two rounds of random priming on bisulfite-treated DNA generate directionally adaptor-tagged template DNAs.</em>
</div>

随机启动引导的 PBAT（Primer Binding At Terminal）使得甲基化文库构建的效率大大提高，甚至可能从 125 pg 基因组 DNA 中制备全基因组甲基化测序（WGBS）文库。通过不断的优化，rPBAT 甚至被用于单细胞甲基组测序，尽管基因组覆盖率有所降低。然而，事实上，rPBAT 文库，特别是那些从极少量 DNA（如 250 pg）中制备的文库，包含越来越多的不能够进行 Mapping 的 reads。这表明 rPBAT 文库中存在很大一部分的嵌合 reads。这些缺点可能源于随机引物启动反应。

具体而言，随机启动通常不是从 3'端开始，而是从目标 ssDNA 片段内的一个位点开始，从而留下一个边缘序列，该边缘序列不包括在新合成的待测序的 reads 中。因此，每个库片段不可避免地比其目标 ssDNA 短，特别是当使用两轮随机引物反应时。其次，随机引物的效率取决于随机引物-基因组 DNA 双链的稳定性。因此，随机启动反应通常会导致对 GC 丰富区域的偏好，而对 AT 丰富区域的基因组覆盖效率较低。最后，随机引物不仅发生在引物和基因组 DNA 之间，也发生在两个引物之间，以及两个基因组 DNA 片段之间，从而分别产生无法 Mapping 和嵌合的 reads。

因此，为了解决这一问题，可以通过消除一个或两个随机引物扩增的步骤，来延长库片段的长度、减少 GC 偏向的覆盖、减少无法 Mapping 和嵌合率，进一步改善 rPBAT。

尽管 rPBAT 方法在制备文库的效率上远超传统的 WGBS 方法，但只有大约 10%的输入 DNA 被转化为文库片段，这表明仍存在优化的空间。是否有方法能够减少随机引物扩增步骤呢？研究者们探索了两种替代方法进行接头连接：一种是利用 RNA 连接酶，该酶能够将 5'端磷酸化的接头连接至 DNA 的 3'末端，但其对 DNA 寡核苷酸间的连接效率低于 RNA 间的连接效率。另一种方法是采用末端脱氧核糖核酸转移酶（TdT）介导的 3'端加尾法，通过特定的末端序列作为引物，合成单链 DNA 的互补链。有研究将这两种技术结合，形成了一种称为 Tdt 辅助的 PBAT（tPBAT）方法，并成功从 125 pg 的人类基因组 DNA 中构建了文库。

如前所述，由于随机引物难以与目标单链 DNA 的 3'端有效杂交，导致亚硫酸盐转化后的 DNA 片段全长覆盖困难。此外，两轮随机引物扩增导致获得的 DNA 片段更短。tPBAT 方法通过 TdT 连接替代了 rPBAT 中的第二轮随机引物扩增，从而产生了插入片段更长的文库分子。

研究显示，在不同输入 DNA 量下，tPBAT 与 rPBAT 具有相似的灵敏度，都能对低输入量样本进行甲基化建库。然而，当输入 DNA 量低于 1 ng 时，rPBAT 的文库分子产量超过了 tPBAT（见图 4）。尽管如此，输入与输出之间仍保持近似线性关系。表 1 中的测序结果表明，尽管在低输入量时 tPBAT 文库产量较低，但其 mapping 率显著高于 rPBAT。随着输入 DNA 量的减少，rPBAT 和 tPBAT 的 mapping 率均有所下降。分析认为，尽管 tPBAT 和 rPBAT 的净产出相当，但 tPBAT 在文库片段长度上具有优势，同时嵌合体的产生也影响了文库产量。

![](https://images.yuanj.top/202404051306501.png)

<div style="text-align:center; margin-bottom:20px;">
<em>Two PBAT methods map the read length distribution on the genome.</em>
</div>

### T-WGBS

Tagmentation-based WGBS (T-WGBS) 是一种利用 Epicenter Tn5 转座体和亚硫酸氢盐转换来研究 5mC 的方法，在文章 [Tagmentation-based whole-genome bisulfite sequencing](https://www.nature.com/articles/nprot.2013.118) 中首次被提出。在该方法中，DNA 与含有甲基化引物的 Tn5 转座体一起孵育，甲基化引物会片段化 DNA 并连接接头。标记的 DNA 首先进行寡核苷酸置换，然后进行甲基化寡核苷酸替换和间隙修复，确保将甲基化接头添加到标记的 DNA 中。然后使用亚硫酸氢钠处理 DNA，进行 PCR 扩增并测序。

![](https://images.yuanj.top/202404051315282.png)

<div style="text-align:center; margin-bottom:20px;">
<em>Overview and components of T-WGBS library preparation. For a detailed description, refer to Experimental design in the INTRODUCTION. Roman numerals indicate stages described therein.</em>
</div>

这种 Tn5mC-seq 方法相对于传统方案，减少了起始材料的使用量超过 100 倍。因此，研究者可以使用 10 ng 的输入 DNA 生成高度复杂的亚硫酸氢盐测序文库，以及从 1 ng 的输入中生成大量有用的序列脱氧核糖核酸。研究者通过对人类淋巴母细胞系的甲基化组进行测序，每条链的高质量覆盖达到 8.6 倍，并证明了 Tn5mC-seq 的有效性。

## WGBS 实验流程

![](https://images.yuanj.top/202404051331600.png)

<div style="text-align:center; margin-bottom:20px;">
<em>WGBS technical process.</em>
</div>

## WGBS 的应用

- 表观遗传学研究：WGBS 可以用于研究不同细胞类型、组织或发育阶段之间的 DNA 甲基化差异，帮助揭示表观遗传变化在生物学过程中的作用，从而有助于理解基因表达调控、细胞分化、发育和疾病等方面的机制
- 疾病研究：WGBS 在疾病研究中扮演着重要角色。研究人员可以比较正常和疾病状态下的 DNA 甲基化模式，以找出与疾病发生和进展相关的甲基化变化，这对癌症、神经系统疾病、心血管疾病等研究具有重要意义
- 个体差异和种群遗传学：WGBS 可以用于研究个体之间的 DNA 甲基化差异，进而帮助我们了解甲基化在人群中的遗传变异，这有助于解析甲基化的遗传基础，以及它如何影响个体的健康和易感性
- 环境影响研究：外部环境因素（如营养、毒素、药物等）可以影响 DNA 甲基化，WGBS 可以帮助研究人员了解环境因素如何通过甲基化改变基因表达，从而影响个体的生理和疾病风险
- 进化研究：WGBS 还可以用于比较不同物种之间的 DNA 甲基化模式，揭示甲基化在进化中的角色。这有助于理解甲基化在物种适应和多样性产生中的贡献

## 参考资料

- [Whole genome bisulfite sequencing](https://en.wikipedia.org/wiki/Whole_genome_bisulfite_sequencing)
- [全基因组甲基化测序技术（WGBS）](https://mp.weixin.qq.com/s/KZD9aOkbreb5TnNpBjsGGA)
- [带你揭开 low input 甲基化建库的神秘面纱](https://zhuanlan.zhihu.com/p/573668799)
- [Tagmentation-based whole-genome bisulfite sequencing](https://www.nature.com/articles/nprot.2013.118)
- [Amplification-free whole-genome bisulfite sequencing by post-bisulfite adaptor tagging](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3458524/)
- [Post-Bisulfite Adaptor Tagging for PCR-Free Whole-Genome Bisulfite Sequencing](https://pubmed.ncbi.nlm.nih.gov/29224142/)
- [Post-bisulfite Adaptor Tagging Based on an ssDNA Ligation Technique (tPBAT)](https://pubmed.ncbi.nlm.nih.gov/36173563/)
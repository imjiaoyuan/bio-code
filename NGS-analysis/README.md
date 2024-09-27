# Description

高通量测序数据分析

# Environment

测试环境：

- 上游分析：成都理工大学超算平台 Red Hat 4.8.5-36
- 下游分析：Debian GNU/Linux 12 (bookworm)

conda 环境搭建：

```bash
git clone https://github.com/imjiaoyuan/NGS-analysis.git
cd NGS-analysis/environment
cp .condarc ~/
pip config set global.index-url https://pypi.tuna.tsinghua.edu.cn/simple
conda env create --file rna-seq-env.yml   # RNA-seq
conda env create --file atac-seq-env.yml  # ATAC-seq
```

# Docs

- [转录组学分析基础——测序技术](./docs/转录组学分析基础——测序技术.md)
- [RNA-seq 上游分析学习](./docs/RNA-seq上游分析学习.md)
- [RNA-seq 上游分析实践](./docs/RNA-seq上游分析实践.md)
- [详解 ATAC-seq](./docs/详解ATAC-seq.md)
- [单细胞转录组技术综述](./docs/单细胞转录组技术综述.md)
- [表观遗传学之 Chip-Seq](./docs/表观遗传学之Chip-seq.md)
- [表观遗传学之 WGBS](./docs/表观遗传学之WGBS.md)
- [挖掘 TFBS 的利器——DAP-Seq](./docs/挖掘TFBS的利器——DAP-seq.md)
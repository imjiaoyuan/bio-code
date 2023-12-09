# Description

使用 GitHub 的 action 功能来进行自动化的 hmmer，主要功能有 hmmersearch、基因筛选、蛋白质序列提取及自定义的 hmmer 命令。

# Usage

你需要在 `config.py`中修改一些参数。

| name | example | description |
| -- | -- | -- |
| PF_number | PF00031 | The Gene family number of the [pfam](http://pfam-legacy.xfam.org/) , such as PF00031 |
| evaluation_threshold | 1e-5 | The threshold for filtering IDs from the results of hmmsearch |
| species | Oryza_sativa | Species of gene family |

仓库文件夹中只有这几个物种的蛋白序列：Arabidopsis_thaliana , Oryza_sativa , Lolium_perenne, and Zea_mays。

## Local usage

使用 conda 建立一个单独的环境来使用，本地使用时需要自己安装 hmmer 和对应的 python 库。

```bash
conda create -n hmmer python=3.7
conda activate hmmer
conda install wget biopython hmmer
```

环境搭建完成后修改 `config.py`然后运行`main.py`。

## Cloud usage

`fork`此仓库，clone 到本地：

```bash
git clone https://github.com/UncleCAT4/auto-hmmer.git
```

修改`config.py`中的参数后提交更改，action 会自动开始。

```bash
git add .
git commit -m "update"
git push
```

你也可以直接在 GitHub 网页版修改参数，但并不建议这样做。

# Results

action 运行完成后你可以在产物中下载以下文件：

> result.out ( Results of hmmsearch )

> id_list.txt ( Filtered ID List )

> protein_for_target_id.fasta ( Protein sequence extracted from ID list )

# 自定义 hmmer 命令

你只需要在`main.py`中定义一个新的 hmmer 命令，例如：

```python
hmmbuild_command = 'hmmbuild model.hmm PFseed.txt'
```

随后在`main.py`中调用命令：

```python
run_command(hmmbuild_command)
```

# Problem

- 没有去除基因的不同转录本。

# More

如果你想使用额外物种的蛋白序列，将序列放在`protein_seq`文件夹下，修改`config.py`中的物种名称即可，请注意，序列的文件拓展名应该是`fasta`，如果是其他拓展名，你可以在`main.py`中修改。

```python
target_species_proteins = "./protein_seq/{}.fasta".format(config.species)
target_gene_proteins = 'protein_for_target_id.fasta'
```

由于 GitHub 的仓库大小限制为 1GB，所以建议`clone`到本地后先清除其它物种的蛋白序列。
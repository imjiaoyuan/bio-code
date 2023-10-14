# Description

This study utilized GitLab CI/CD and GitHub Actions to automate the construction of the Hidden Markov Model ( HMM ) and conduct a search for Homeotic genes. It further screened out gene IDs that satisfy the user-defined threshold and extracted protein sequences from the obtained results files.

# Usage

You need to define the following variables in `config.py`.

| name | example | description |
| -- | -- | -- |
| PF_number | PF00031 | The Gene family number of the [pfam](http://pfam-legacy.xfam.org/) , such as PF00031 |
| evaluation_threshold | 1e-5 | The threshold for filtering IDs from the results of hmmsearch |
| species | Oryza_sativa | Species of gene family |

The repository only contains the genome protein sequences of three species: Arabidopsis_thaliana , Oryza_sativa , Lolium_perenne, and Zea_mays .

## Local usage

I suggest that you use Conda for installation.

```bash
conda create -n hmmer python=3.7
conda activate hmmer
conda install wget biopython hmmer
```

After setting up the environment, modify the information in `config.py` and run `main.py`.

## Cloud usage

You need to clone the warehouse locally first.

```bash
git clone https://github.com/UncleCAT4/auto-hmmer.git
```

Fork this warehouse to your own account, modify the information in `config.py`, and submit the changes to the repository. CI/CD will help you complete the remaining tasks.

```bash
git add .
git commit -m "update"
git push
```

You can also make modifications directly on the webpage of the code repository, and the running results can be downloaded from the product.

# Results

If you are running locally, the results are clearly visible.If you are using it in the cloud,you can download the completed files from the GitHub Action or Gitlab CI/CD running artifacts.

Usually contains the following three main artifacts:

> result.out ( Results of hmmsearch )

> id_list.txt ( Filtered ID List )

> protein_for_target_id.fasta ( Protein sequence extracted from ID list )

# Customize the hmmer command

You only need to define the format of the new command in `main.py`, such as the following command:

```python
hmmbuild_command = 'hmmbuild model.hmm PFseed.txt'
```

Then run the command by calling the function in `main.py`.

```python
run_command(hmmbuild_command)
```

# Problem

Different transcripts of the same gene were not removed.

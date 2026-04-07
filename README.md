# Merge and Correlation Analysis

An exercise using Python that looks at merging two bioinformatics datasets using Pandas and performing correlation analysis between gene expression and mutation data.

## Setup

1. Clone the repository:

```
git clone https://github.com/edecocke-uncc/merge_practice.git
```

2. Go into the project folder:

```
cd merge_practice
```
3. Create the environment:
```bash
conda env create -f environment.yml
```

4. Activate the environment:

```
conda activate merge_practice
```

## Run

```
python3 merge_practice.py
```

## What it does

* Loads gene expression data and mutation count data for 5 genes (BRCA1, TP53, EGFR, MYC, KRAS) across Heart, Liver, and Brain samples
* Merges the two datasets on the shared Gene column using pd.merge()
* Computes mean expression per gene across all three tissues
* Performs a Pearson correlation analysis to test whether highly expressed genes tend to accumulate more mutations
* Interprets the correlation results in a biological context

## Output

A `gene_expression_mutations.csv` file is written to the current directory after each run.

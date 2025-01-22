# strokeDTI

Official implementation of **strokeDTI**, a tool designed for graph-based analysis of RNA-seq data and prediction of drug-target interactions.

It was tested for cell death pathways, but in theory should work for most pathways with shared genes.

![Graphical abstract](/strokeDTI/img/graphical_abstract_fig_1.png)

---

## Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
   - [Target Identification Module](#1-target-identification-module)
   - [DTI Prediction](#2-dti-prediction)
3. [Cite Us](#cite-us)
4. [References](#references)

---

# Installation

To set up the environment, use the following commands:

```
git clone
cd path/to/strokeDTI
conda env create -f environment.yml

```

To remove the conda environment

```
conda deactivate
conda env remove -n stroke_dti
```

## Usage

**strokeDTI** consists of three modules, which can be used independently.

### 1. Target Identification Module

#### Input:

1. **Deseq2 Data**: Processed RNA-seq data in `.tabular` format. The file must include the following headers:

   | GeneID | Base mean | log2(FC) | StdErr | Wald-Stats | P-value | P-adj | Chromosome | Start | End | Strand | Feature | Gene name |

   Example RNA-seq data was obtained from public repositories: [GSE137482](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE137482)[1].

2. **KEGG Pathways**: A `.txt` file with terms separated by commas.
   - KEGG pathways can be obtained from the official [KEGG website](https://www.genome.jp/kegg/pathway.html).
   - Example file: `strokeDTI/data/kegg.txt`.

Target identification module usage example:

```
conda activate stroke_dti

identify_targets -s mouse -rna_seq data/galaxy_mouse.tabular -kegg data/kegg.txt -out output/

```

### Output

1. CSV files of the kegg terms merged with Deseq2 data
2. Merged network graph 1
   ![Merged Graph 1](/strokeDTI/output/graph1.png)

3. Merged network graph 2
   ![Merged Graph 1](/strokeDTI/output/graph2.png)

4. CSV files of node importance:

- degree centrality
- closeness centrality
- betweenness centrality
- eigenvector centrality
- k_shell

## 2. DTI prediction

Coming soon! (Details will be added in future updates)

## Cite us:

If you use StrokeDTI in your research, please cite our work:

```

```

## References

[1] Androvic P, Kirdajova D, Tureckova J, Zucha D, Rohlova E, Abaffy P, et al. Decoding the Transcriptional Response to Ischemic Stroke in Young and Aged Mouse Brain. Cell Rep. 2020;31:107777.

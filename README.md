# Drug2cell

This is a collection of utility functions for gene group activity evaluation in scanpy, both for per-cell scoring and marker-based enrichment/overrepresentation analyses. Drug2cell makes use of established methodology, and offers it in a convenient/efficient form for easy scanpy use. The package comes with a set of ChEMBL derived drug target sets, but has straightforward input formatting so it can be easily used with any gene groups of choice.

## Installation

```bash
pip install drug2cell
```

Drug2cell's GSEA makes use of blitzGSEA, which needs to be installed from GitHub:

```bash
pip install git+https://github.com/MaayanLab/blitzgsea.git
```

## Usage

```python3
import drug2cell as d2c
```

**Per-cell scoring** can be done via `d2c.score()`, and creates a new gene group feature space object in `.uns['drug2cell']` of the supplied object. The `.obs` and `.obsm` of the original object are copied over for ease of downstream use, like plotting or potential marker detection.

**Enrichment** is done with GSEA, contained within the function `d2c.gsea()`. **Overrepresentation** analysis can be done with the hypergeometric test of `d2c.hypergeometric()`. Both of those require having ran `sc.tl.rank_genes_groups()` on the input object, and return one data frame per evaluated cluster with enrichment/overrepresentation results.

It's possible to provide your own gene groups as a dictionary, with the names of the groups as keys and the corresponding gene lists as the values. Pass this dictionary as the `targets` argument of any of those three functions.

**Please refer to the [demo notebook](https://nbviewer.org/github/Teichlab/drug2cell/blob/main/notebooks/demo.ipynb), all of the input arguments are detailed in [ReadTheDocs](https://drug2cell.readthedocs.io/en/latest/).**

## ChEMBL parsing

Drug2cell also features instructions on how to parse the ChEMBL database into drugs and their targets. Some additional notebooks are included. Both refer to a pre-parsed data frame of ChEMBL human targets, which can be accessed at `ftp://ftp.sanger.ac.uk/pub/users/kp9/chembl_30_merged_genesymbols_humans.pkl`.
 - [Filtering](https://nbviewer.org/github/Teichlab/drug2cell/blob/main/notebooks/chembl/filtering.ipynb) shows how this data frame was turned into the drugs:targets dictionary shipped with the package by default. There are some helper functions included in `d2c.chembl` which can assist you shall you wish to filter it in a different way.
 - [Initial database parsing](https://nbviewer.org/github/Teichlab/drug2cell/blob/main/notebooks/chembl/initial_database_parsing.ipynb) shows how the data frame was created from online resources.
 
## Citation

If you use drug2cell in your work, please cite the [preprint](https://www.biorxiv.org/content/10.1101/2023.01.30.526202v1)

```
@article {Kanemaru2023.01.30.526202,
	author = {Kanemaru, Kazumasa and Cranley, James and Muraro, Daniele and Miranda, Antonio M.A. and Pett, Jan Patrick and Litvinukova, Monika and Kumasaka, Natsuhiko and Ho, Siew Yen and Polanski, Krzysztof and Richardson, Laura and Mach, Lukas and Dabrowska, Monika and Richoz, Nathan and Barnett, Sam N. and Perera, Shani and Wilbrey-Clark, Anna L and Talavera-L{\'o}pez, Carlos and Mulas, Ilaria and Mahbubani, Krishnaa T. and Bolt, Liam and Mamanova, Lira and Tuck, Liz and Wang, Lu and Huang, Margaret M. and Prete, Martin and Pritchard, Sophie and Dark, John and Saeb-Parsy, Kourosh and Patel, Minal and Clatworthy, Menna R. and Chowdhury, Rasheda A. and Noseda, Michela and Teichmann, Sarah A.},
	title = {Spatially resolved multiomics of human cardiac niches},
	elocation-id = {2023.01.30.526202},
	year = {2023},
	doi = {10.1101/2023.01.30.526202},
	publisher = {Cold Spring Harbor Laboratory},
	abstract = {A cell{\textquoteright}s function is defined by its intrinsic characteristics and its niche: the tissue microenvironment in which it dwells. Here, we combine single-cell and spatial transcriptomic data to discover cellular niches within eight regions of the human heart. We map cells to micro-anatomic locations and integrate knowledge-based and unsupervised structural annotations. For the first time, we profile the cells of the human cardiac conduction system, revealing their distinctive repertoire of ion channels, G-protein coupled receptors, and cell interactions using a custom CellPhoneDB.org module. We show that the sinoatrial node is compartmentalised, with a core of pacemaker cells, fibroblasts and glial cells supporting paracrine glutamatergic signalling. We introduce a druggable target prediction tool, drug2cell, which leverages single-cell profiles and drug-target interactions, providing unexpected mechanistic insights into the chronotropic effects of drugs, including GLP-1 analogues. In the epicardium, we show enrichment of both IgG+ and IgA+ plasma cells forming immune niches which may contribute to infection defence. We define a ventricular myocardial-stress niche enriched for activated fibroblasts and stressed cardiomyocytes, cell states that are expanded in cardiomyopathies. Overall, we provide new clarity to cardiac electro-anatomy and immunology, and our suite of computational approaches can be deployed to other tissues and organs.Competing Interest StatementIn the past three years, S.A.T. has consulted or been a member of scientific advisory boards at Roche, Genentech, Biogen, GlaxoSmithKline, Qiagen and ForeSite Labs and is an equity holder of Transition Bio. The remaining authors declare no competing interests.},
	URL = {https://www.biorxiv.org/content/early/2023/02/01/2023.01.30.526202},
	eprint = {https://www.biorxiv.org/content/early/2023/02/01/2023.01.30.526202.full.pdf},
	journal = {bioRxiv}
}
```

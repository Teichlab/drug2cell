# Drug2cell

This is a collection of utility functions for gene group activity evaluation in scanpy, both for per-cell scoring and marker-based enrichment/overrepresentation analyses. Drug2cell makes use of established methodology, and offers it in a convenient/efficient form for easy scanpy use. The package comes with a set of ChEMBL derived drug target sets, but has straightforward input formatting so it can be easily used with any gene groups of choice.

## Installation

Once public, this will be on `pip`. For now install straight from this repository:

```bash
pip install git+https://github.com/Teichlab/drug2cell.git
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

**Please refer to the [demo notebook](notebooks/demo.ipynb), all of the input arguments are detailed in [ReadTheDocs](https://drug2cell.readthedocs.io/en/latest/).**

## ChEMBL parsing

Drug2cell also features instructions on how to parse the ChEMBL database into drugs and their targets. Some additional notebooks are included. Both refer to a pre-parsed data frame of ChEMBL human targets, which can currently be accessed at `/nfs/team205/kk18/data/ChEMBL/pkls/chembl_30_merged_genesymbols_humans.pkl`. Once public, this will be put up on FTP.
 - [Filtering](notebooks/chembl/filtering.ipynb) shows how this data frame was turned into the drugs:targets dictionary shipped with the package by default. There are some helper functions included in `d2c.chembl` which can assist you shall you wish to filter it in a different way.
 - [Initial database parsing](notebooks/chembl/initial_database_parsing.ipynb) shows how the data frame was created from online resources.
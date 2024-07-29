import blitzgsea
import anndata
import pandas as pd
import numpy as np

#gives access to submodules
from . import chembl
from . import data
from . import util

from statsmodels.stats.multitest import multipletests
from scipy.sparse import issparse
from scipy.stats import hypergeom

def _sparse_nanmean(X, axis):
    # function from https://github.com/scverse/scanpy/blob/034ca2823804645e0d4874c9b16ba2eb8c13ac0f/scanpy/tools/_score_genes.py
    """
    np.nanmean equivalent for sparse matrices
    """
    if not issparse(X):
        raise TypeError("X must be a sparse matrix")

    # count the number of nan elements per row/column (dep. on axis)
    Z = X.copy()
    Z.data = np.isnan(Z.data)
    Z.eliminate_zeros()
    n_elements = Z.shape[axis] - Z.sum(axis)

    # set the nans to 0, so that a normal .sum() works
    Y = X.copy()
    Y.data[np.isnan(Y.data)] = 0
    Y.eliminate_zeros()

    # the average
    s = Y.sum(axis, dtype='float64')  # float64 for score_genes function compatibility)
    m = s / n_elements

    return m

def _mean(X,names,axis):
    '''
    Helper function to compute a mean of X across an axis, respecting names and possible nans.
    
    Derived from sc.tl.score_genes() logic.
    '''
    if issparse(X):
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(X, axis=axis)).flatten(),  
            index=names,
        )  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(X, axis=axis), index=names
        )  # average expression of genes
    return obs_avg

def score(adata, targets=None, nested=False, categories=None, method="mean", layer=None, use_raw=False, n_bins=25, ctrl_size=50, sep=","):
    '''
    Obtain per-cell scoring of gene groups of interest. Distributed with a set of 
    ChEMBL drug targets that can be used immediately.
    
    Please ensure that the gene nomenclature in your target sets is compatible with your 
    ``.var_names`` (or ``.raw.var_names``). The ChEMBL drug targets use HGNC (human gene 
    names in line with standard cell ranger mapping output).
    
    Adds ``.uns['drug2cell']`` to the input AnnData, a new AnnData object with the same 
    observation space but with the scored gene groups as the features. The gene group 
    members used to compute the scores will be listed in ``.var['genes']`` of the new 
    object.
    
    Input
    -----
    adata : ``AnnData``
        Using log-normalised data is recommended.
    targets : ``dict`` of lists of ``str``, optional (default: ``None``)
        The gene groups to evaluate. Can be targets of known drugs, GO terms, pathway 
        memberships, anything you can assign genes to. If ``None``, will load the 
        ChEMBL-derived drug target sets distributed with the package.
        
        Accepts two forms:
        
            - A dictionary with the names of the groups as keys, and the entries being the 
            corresponding gene lists.
            - A dictionary of dictionaries defined like above, with names of gene group 
            categories as keys. If passing one of those, please specify ``nested=True``.
        
    nested : ``bool``, optional (default: ``False``)
        Whether ``targets`` is a dictionary of dictionaries with group categories as keys.
    categories : ``str`` or list of ``str``, optional (default: ``None``)
        If ``targets=None`` or ``nested=True``, this argument  can be used to subset the 
        gene groups to one or more categories (keys of the original dictionary). In case 
        of the ChEMBL drug targets, these are ATC level 1/level 2 category codes.
    method : ``str``, optional (default: ``"mean"``)
        The method to use to score the gene groups. The default is ``"mean"``, which 
        computes the mean over all the genes. The other option is ``"seurat"``, which 
        generates an appropriate background profile for each target set and subtracts it 
        from the mean. This is inspired by ``sc.tl.rank_genes()`` logic, which in turn 
        was inspired by Seurat's gene group scoring algorithm.
    layer : ``str``, optional (default: ``None``)
        Which ``.layers`` of the input AnnData to use for the expression values. If 
        ``None``, will default to ``.X``.
    use_raw : ``bool``, optional (default: ``False``)
        Whether to use ``.raw.X`` for the expression values.
    n_bins : ``int``, optional (default: 25)
        Only used with ``method="seurat"``. The number of expression bins to partition the 
        feature space into.
    ctrl_size : ``int``, optional (default: 50)
        Only used with ``method="seurat"``. The number of genes to randomly sample from 
        each expression bin.
    sep : ``str``, optional (default: ``","``)
        What delimiter to use when storing the corresponding gene groups for each feature 
        in ``.uns['drug2cell'].var['genes']``
    '''
    #select expression and gene names to use
    if layer is not None:
        if use_raw:
            raise ValueError("Cannot specify `layer` and have `use_raw=True`.")
        X = adata.layers[layer]
        var_names = adata.var_names
    else:
        if use_raw and adata.raw is not None:
            X = adata.raw.X
            var_names = adata.raw.var_names
        else:
            X = adata.X
            var_names = adata.var_names
    #get {group:[targets]} form of gene groups to evaluate based on arguments
    #skip target reconstruction to overwrite any potential existing scoring
    targets = util.prepare_targets(
        adata,
        targets=targets,
        nested=nested,
        categories=categories,
        sep=sep,
        reconstruct=False
    )
    #store full list of targets to have them on tap for later
    full_targets = targets.copy()
    #turn the list of gene IDs to a boolean mask of var_names
    for drug in targets:
        targets[drug] = np.isin(var_names, targets[drug])
    #perform scoring
    #the scoring shall be done via matrix multiplication
    #of the original cell by gene matrix, by a new gene by drug matrix
    #with the entries in the new matrix being the weights of each gene for that drug
    #the first part, the mean across targets, is constant; prepare weights for that
    weights = pd.DataFrame(targets, index=var_names)
    #kick out drugs with no targets
    weights = weights.loc[:, weights.sum()>0]
    #scale to 1 sum for each column, weights for mean acquired. get mean
    weights = weights/weights.sum()
    if issparse(X):
        scores = X.dot(weights)
    else:
        scores = np.dot(X, weights)
    #the second part only happens for seurat scoring
    #logic inspired by sc.tl.score_genes()
    if method == "seurat":
        #obtain per-gene means
        obs_avg = _mean(X, names=var_names, axis=0)
        #bin the genes (score_genes() logic in full effect here)
        n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
        obs_cut = obs_avg.rank(method='min') // n_items
        #we'll be working on the array for ease of stuff going forward
        obs_cut = obs_cut.values
        #compute a cell by control group matrix
        #using matrix multiplication again
        #so build a gene by control group matrix to enable it
        control_groups = {}
        for cut in np.unique(obs_cut):
            #get the locations of the value
            mask = obs_cut == cut
            #get the nonzero values, i.e. indices of the locations
            r_genes = np.nonzero(mask)[0]
            #shuffle these indices, and only keep the top N in the mask
            np.random.shuffle(r_genes)
            mask[r_genes[ctrl_size:]] = False
            #store mask
            control_groups[cut] = mask
        #turn to weights matrix, like earlier
        control_gene_weights = pd.DataFrame(control_groups, index=var_names)
        control_gene_weights = control_gene_weights/control_gene_weights.sum()
        #compute control profiles
        if issparse(X):
            control_profiles = X.dot(control_gene_weights)
        else:
            control_profiles = np.dot(X, control_gene_weights)
        #identify the bins for each drug
        drug_bins = {}
        #use the subset form that's in weights
        for drug in weights.columns:
            #targets is still a handy boolean vector of membership
            #get the bins for this drug
            bins = np.unique(obs_cut[targets[drug]])
            #mask the bin order in the existing variables with this
            #control_gene_weights has column names
            #control_profiles is ordered the same way but is just values
            drug_bins[drug] = np.isin(control_gene_weights.columns, bins)
        #turn to weights matrix again
        drug_weights = pd.DataFrame(drug_bins, index=control_gene_weights.columns)
        drug_weights = drug_weights/drug_weights.sum()
        #can now get the seurat reference profile via matrix multiplication
        #subtract from the existing scores
        seurat = np.dot(control_profiles, drug_weights)
        scores = scores - seurat
    #we now have the final form of the scores
    #create a little helper adata thingy based on them
    #store existing .obsm in there for ease of plotting stuff
    adata.uns['drug2cell'] = anndata.AnnData(scores, obs=adata.obs)
    adata.uns['drug2cell'].var_names = weights.columns
    adata.uns['drug2cell'].obsm = adata.obsm
    #store gene group membership, going back to targets for it
    for drug in weights.columns:
        #mask the var_names with the membership, and join into a single delimited string
        adata.uns['drug2cell'].var.loc[drug, 'genes'] = sep.join(var_names[targets[drug]])
        #pull out the old full membership dict and store that too
        adata.uns['drug2cell'].var.loc[drug, 'all_genes'] = sep.join(full_targets[drug])

def gsea(adata, targets=None, nested=False, categories=None, absolute=False, plot_args=True, sep=",", **kwargs):
    '''
    Perform gene set enrichment analysis on the marker gene scores computed for the 
    original object. Uses blitzgsea.
    
    Returns:
    
    - a dictionary with clusters for which the original object markers were computed \
    as the keys, and data frames of test results sorted on q-value as the items
    
    - a helper variable with plotting arguments for ``d2c.plot_gsea()``, if \
    ``plot_args=True``. ``['scores']`` has the GSEA input, and ``['targets']`` is the \
    gene group dictionary that was used.
    
    Input
    -----
    adata : ``AnnData``
        With marker genes computed via ``sc.tl.rank_genes_groups()`` in the original 
        expression space.
    targets : ``dict`` of lists of ``str``, optional (default: ``None``)
        The gene groups to evaluate. Can be targets of known drugs, GO terms, pathway 
        memberships, anything you can assign genes to. If ``None``, will use 
        ``d2c.score()`` output if present, and if not present load the ChEMBL-derived 
        drug target sets distributed with the package.
        
        Accepts two forms:
        
            - A dictionary with the names of the groups as keys, and the entries being the 
            corresponding gene lists.
            - A dictionary of dictionaries defined like above, with names of gene group 
            categories as keys. If passing one of those, please specify ``nested=True``.
        
    nested : ``bool``, optional (default: ``False``)
        Whether ``targets`` is a dictionary of dictionaries with group categories as keys.
    categories : ``str`` or list of ``str``, optional (default: ``None``)
        If ``targets=None`` or ``nested=True``, this argument  can be used to subset the 
        gene groups to one or more categories (keys of the original dictionary). In case 
        of the ChEMBL drug targets, these are ATC level 1/level 2 category codes.
    absolute : ``bool``, optional (default: ``False``)
        If ``True``, pass the absolute values of scores to GSEA. Improves statistical 
        power.
    plot_args : ``bool``, optional (default: ``True``)
        Whether to return the second piece of output that holds pre-compiled information 
        for ``d2c.plot_gsea()``.
    sep : ``str``, optional (default: ``","``)
        The delimiter that was used with ``d2c.score()`` for gene group storage.
    kwargs
        Any additional arguments to pass to ``blitzgsea.gsea()``.
    '''
    #get {group:[targets]} form of gene groups to evaluate based on arguments
    #allow for target reconstruction for when this is ran after scoring
    targets = util.prepare_targets(
        adata,
        targets=targets,
        nested=nested,
        categories=categories,
        sep=sep,
        reconstruct=True
    )
    #store the GSEA results in a dictionary, with groups as the keys
    enrichment = {}
    #the plotting-minded output can already store its targets
    #and will keep track of scores as they get made during the loop
    plot_gsea_args = {"targets":targets, "scores":{}}
    #this gets the names of the clusters in the original marker output
    for cluster in adata.uns['rank_genes_groups']['names'].dtype.names:
        #prepare blitzgsea input
        df = pd.DataFrame({"0":adata.uns['rank_genes_groups']['names'][cluster],
                           "1":adata.uns['rank_genes_groups']['scores'][cluster]})
        #possibly sort on absolute value of scores
        if absolute:
            df["1"] = np.absolute(df["1"])
            df = df.sort_values("1", ascending=False)
        #compute GSEA and store results/scores in output
        enrichment[cluster] = blitzgsea.gsea(df, targets, **kwargs)
        plot_gsea_args["scores"][cluster] = df
    #provide output
    if plot_args:
        return enrichment, plot_gsea_args
    else:
        return enrichment

def hypergeometric(adata, targets=None, nested=False, categories=None, pvals_adj_thresh=0.05, direction="both", corr_method="benjamini-hochberg", sep=","):
    '''
    Perform a hypergeometric test to assess the overrepresentation of gene group members 
    among marker genes computed for the original object.
    
    Returns a dictionary with clusters for which the original object markers were computed 
    as the keys, and data frames of test results sorted on q-value as the items.
    
    Input
    -----
    adata : ``AnnData``
        With marker genes computed via ``sc.tl.rank_genes_groups()`` in the original 
        expression space.
    targets : ``dict`` of lists of ``str``, optional (default: ``None``)
        The gene groups to evaluate. Can be targets of known drugs, GO terms, pathway 
        memberships, anything you can assign genes to. If ``None``, will use 
        ``d2c.score()`` output if present, and if not present load the ChEMBL-derived 
        drug target sets distributed with the package.
        
        Accepts two forms:
        
            - A dictionary with the names of the groups as keys, and the entries being the 
            corresponding gene lists.
            - A dictionary of dictionaries defined like above, with names of gene group 
            categories as keys. If passing one of those, please specify ``nested=True``.
        
    nested : ``bool``, optional (default: ``False``)
        Whether ``targets`` is a dictionary of dictionaries with group categories as keys.
    categories : ``str`` or list of ``str``, optional (default: ``None``)
        If ``targets=None`` or ``nested=True``, this argument  can be used to subset the 
        gene groups to one or more categories (keys of the original dictionary). In case 
        of the ChEMBL drug targets, these are ATC level 1/level 2 category codes.
    pvals_adj_thresh : ``float``, optional (default: ``0.05``)
        The ``pvals_adj`` cutoff to use on the ``sc.tl.rank_genes_groups()`` output to 
        identify markers.
    direction : ``str``, optional (default: ``"both"``)
        Whether to seek out up/down-regulated genes for the groups, based on the values 
        from ``scores``. Can be ``"up"``, ``"down"``, or ``"both"`` (for no selection).
    corr_method : ``str``, optional (default: ``"benjamini-hochberg"``)
        Which FDR correction to apply to the p-values of the hypergeometric test. Can be 
        ``"benjamini-hochberg"`` or ``"bonferroni"``.
    sep : ``str``, optional (default: ``","``)
        The delimiter that was used with ``d2c.score()`` for gene group storage.
    '''
    #get the universe of available genes, as a set for easy intersecting
    if adata.uns['rank_genes_groups']['params']['use_raw']:
        universe = set(adata.raw.var_names)
    else:
        universe = set(adata.var_names)
    #get {group:[targets]} form of gene groups to evaluate based on arguments
    #allow for target reconstruction for when this is ran after scoring
    targets = util.prepare_targets(
        adata,
        targets=targets,
        nested=nested,
        categories=categories,
        sep=sep,
        reconstruct=True
    )
    #intersect each group membership with the universe after turning it to a set
    for group in targets:
        targets[group] = set(targets[group]).intersection(universe)
    #kick out any empty keys, using dictionary comprehension
    targets = {k:v for k,v in targets.items() if v}
    #store the hypergeometric results in a dictionary, with groups as the keys
    overrepresentation = {}
    #this gets the names of the clusters in the original marker output
    for cluster in adata.uns['rank_genes_groups']['names'].dtype.names:
        #prepare overrepresentation output data frame
        results = pd.DataFrame(
            1, 
            index=list(targets.keys()),
            columns=['intersection','gene_group','markers','universe','pvals', 'pvals_adj']
        )
        #pre-type pvals as floats as otherwise pandas complains about deprecation
        results = results.astype({"pvals":"float64", "pvals_adj":"float64"})
        #identify the markers for this group from the original output
        #construct mask for significance
        mask = adata.uns['rank_genes_groups']['pvals_adj'][cluster] < pvals_adj_thresh
        #the sign of the score indicates up/down-regulation
        if direction == "up":
            mask = mask & (adata.uns['rank_genes_groups']['scores'][cluster] > 0)
        elif direction == "down":
            mask = mask & (adata.uns['rank_genes_groups']['scores'][cluster] < 0)
        #do this as a set because it will be easier to intersect later
        markers = set(adata.uns['rank_genes_groups']['names'][cluster][mask])
        #at this point we know how many genes we have as markers and the universe
        results['markers'] = len(markers)
        results['universe'] = len(universe)
        #perform hypergeometric test to assess overrepresentation
        #loop over the scored gene groups
        for ind in results.index:
            #retrieve gene membership and the intersection
            gene_group = targets[ind]
            common = gene_group.intersection(markers)
            results.loc[ind, 'intersection'] = len(common)
            results.loc[ind, 'gene_group'] = len(gene_group)
            #need to subtract 1 from the intersection length
            #https://alexlenail.medium.com/understanding-and-implementing-the-hypergeometric-test-in-python-a7db688a7458
            pval = hypergeom.sf(
                len(common)-1,
                len(universe),
                len(markers),
                len(gene_group)
            )
            results.loc[ind, 'pvals'] = pval
        #multiple testing correction
        #mirror sc.tl.rank_genes_groups() logic for consistency
        #just in case any NaNs popped up somehow, fill them to 1 so FDR works
        results = results.fillna(1)
        if corr_method == "benjamini-hochberg":
            results['pvals_adj'] = multipletests(results['pvals'], method="fdr_bh")[1]
        elif corr_method == "bonferroni":
            results['pvals_adj'] = np.minimum(results['pvals']*results.shape[0], 1.0)
        #sort on q-value and store output
        overrepresentation[cluster] = results.sort_values("pvals_adj")
    #that's it. return the dictionary
    return overrepresentation
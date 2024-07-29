import blitzgsea
import numpy as np

from . import data

from collections import ChainMap

def reconstruct_targets(adata, sep=","):
    '''
    Reconstruct the targets dictionary from an AnnData with ``d2c.score()`` output stored.
    '''
    targets = {}
    #loop over the group names, i.e. this object's feature space
    for group in adata.uns['drug2cell'].var_names:
        #split the stored gene list
        targets[group] = adata.uns['drug2cell'].var.loc[group,'all_genes'].split(sep)
    return targets

def prepare_targets(adata, targets=None, nested=False, categories=None, sep=",", reconstruct=False):
    '''
    A helper function that ensures the function calling it gets a {group:[targets]} 
    dictionary back.
    '''
    #no provided explicit target set
    if targets is None:
        #if there's .uns['drug2cell'] and we're instructed to reconstruct, do so
        if reconstruct and ('drug2cell' in adata.uns):
            targets = reconstruct_targets(adata, sep=sep)
        else:
            #load ChEMBL
            targets = data.chembl()
            #this is nested, regardless of what the input may say
            nested = True
    else:
        #copy so the original is unaltered
        targets = targets.copy()
    #subset to provided categories. this should only be provided for nested versions
    if categories is not None:
        #turn to list if need be
        if type(categories) is not list:
            categories = list(categories)
        #so that they can be iterated over here for subsetting
        targets = {k:targets[k] for k in categories}
    #do we have a dict of dicts, with the categories as keys?
    if nested:
        #turn the dict of dicts structure to just a basic dictionary
        #get all the drugs for all the remaining categories, making a list of dicts
        #turn that list to a single dict via ChainMap magic
        targets = dict(ChainMap(*[targets[cat] for cat in targets]))
    #at this point we have the desired form of {group:[targets]}, return it
    return targets

def prepare_plot_args(adata, targets=None, categories=None):
    '''
    Prepare the ``var_names``, ``var_group_positions`` and ``var_group_labels`` arguments 
    for scanpy plotting functions to display scored gene groups and group them nicely. 
    Returns ``plot_args``, a dictionary of the values that can be used with scanpy 
    plotting as ``**plot_args``.
    
    Input:
    ------
    adata : ``AnnData``
        Point the function to the ``.uns['drug2cell']`` slot computed by the ``score()`` 
        function earlier. It's required to remove gene groups that were not represented 
        in the data.
    targets : ``dict`` of lists of ``str``, optional (default: ``None``)
        The gene groups to evaluate. Can be targets of known drugs, GO terms, pathway 
        memberships, anything you can assign genes to. If ``None``, will load the 
        ChEMBL-derived drug target sets distributed with the package. Must be the 
        ``nested=True`` version of the input as described in the score function.
    categories : ``str`` or list of ``str``, optional (default: ``None``)
        If ``targets=None`` or ``nested=True``, this argument  can be used to subset the 
        gene groups to one or more categories (keys of the original dictionary). In case 
        of the ChEMBL drug targets, these are ATC level 1/level 2 category codes.
    '''
    #load ChEMBL dictionary if not specified
    if targets is None:
        targets = data.chembl()
    else:
        #copy so the original argument is unaltered
        targets = targets.copy()
    #subset to provided categories
    if categories is not None:
        #turn to list if need be
        if type(categories) is not list:
            categories = list(categories)
        #so that they can be iterated over here for subsetting
        targets = {k:targets[k] for k in categories}
    #this time we only care about the group names. kick out their targets
    for group in targets:
        targets[group] = list(targets[group].keys())
    #can commence constructing the plotting arguments now
    var_names = []
    var_group_positions = []
    var_group_labels = []
    #we'll begin at the first feature. this is zero indexed
    start = 0
    for group in targets:
        #intersect the gene group names with what's available in the object
        #i.e. what was actually scored and is available for plotting
        targets[group] = list(adata.var_names[np.isin(adata.var_names, targets[group])])
        #skip if empty
        if len(targets[group]) == 0:
            continue
        #if not empty, store!
        #append group names to list
        var_names = var_names + targets[group]
        #the newest group starts at start and ends at the end of the var_names
        var_group_positions = var_group_positions + [(start, len(var_names)-1)]
        var_group_labels = var_group_labels + [group]
        #the new start will be the next feature that goes in the list
        start = len(var_names)
    #that's it, return the things as a dict
    plot_args = {'var_names':var_names,
                 'var_group_positions':var_group_positions,
                 'var_group_labels':var_group_labels}
    return plot_args

def plot_gsea(enrichment, targets, scores, n=10, interactive_plot=True, **kwargs):
    '''
    Display the output of ``d2c.gsea()`` with blitzgsea's ``top_table()`` plot.
    
    The first ``d2c.gsea()`` output variable is ``enrichment``, and passing the second 
    ``d2c.gsea()`` output variable with a ``**`` in front of it provides ``targets`` and 
    ``scores``.
    
    Input
    -----
    enrichment : ``dict`` of ``pd.DataFrame``
        Cluster names as keys, blitzgsea's ``gsea()`` output as values
    targets : ``dict`` of list of ``str``
        The gene group memberships that were used to compute GSEA
    scores : ``dict`` of ``pd.DataFrame``
        Cluster names as keys, the input to blitzgsea
    n : ``int``, optional (default: ``10``)
        How many top scores to show for each group
    interactive_plot : ``bool``, optional (default: ``True``)
        If ``True``, will display the plots within a Jupyter Notebook. If ``False``, 
        will collect the figures into a list and return it at the end.
    kwargs
        Any additional arguments to pass to ``blitzgsea.plot.top_table()``.
    '''
    #optionally save output
    if not interactive_plot:
        figs = []
    #make a top_table() plot for each cluster
    for cluster in enrichment:
        #get a figure, passing the various arguments
        #interactive_plot prepares the figure so that it .show()s in a notebook
        fig = blitzgsea.plot.top_table(scores[cluster], 
                                       targets, 
                                       enrichment[cluster], 
                                       n=n, 
                                       interactive_plot=interactive_plot,
                                       **kwargs
                                      )
        #retitle accordingly
        fig.suptitle(cluster)
        #either display the plot or stash it in output
        if interactive_plot:
            fig.show()
        else:
            figs.append(fig)
    if not interactive_plot:
        return figs
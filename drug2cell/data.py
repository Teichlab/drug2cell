import pkg_resources
import pandas as pd

def chembl():
    '''
    Load the default ChEMBL drug target dictionary distributed with the package.
    
    Returns the drug target dictionary - ATC categories as keys, with each ATC category 
    a dictionary with corresponding drugs as keys.
    '''
    #this picks up the pickle shipped with the package
    stream = pkg_resources.resource_stream(__name__, 'drug-target_dicts.pkl')
    targets = pd.read_pickle(stream)
    return targets

def consensuspathdb():
    '''
    Load the ConsensusPathDB pathway gene memberships distributed with the package.
    
    Returns a dictionary with pathway names as keys and memberships as items.
    '''
    #this picks up the pickle shipped with the package
    stream = pkg_resources.resource_stream(__name__, 'cpdb_dict.pkl')
    targets = pd.read_pickle(stream)
    return targets
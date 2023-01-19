import numpy as np
import pandas as pd

def filter_activities(dataframe,
                      drug_max_phase=None,
                      add_drug_mechanism=True,
                      assay_type=None,
                      remove_inactive=True,
                      include_active=True,
                      pchembl_target_column=None,
                      pchembl_threshold=None
                        ):
    '''
    Perform a sequential set of filtering operations on the provided ChEMBL data frame. 
    The order of the filters matches the order of the arguments in the input description. 
    Returns a data frame with the rows fulfilling the resulting criteria.
    
    Input
    -----
    dataframe : ``pd.DataFrame``
        The ChEMBL data frame to perform filtering operations on.
    drug_max_phase : ``int`` or list of ``int``, optional (default: ``None``)
        Subset the data frame to drugs in the provided clinical stages:
        
            - Phase 1: Testing of drug on healthy volunteers for dose-ranging
            - Phase 2: Initial testing of drug on patients to assess efficacy and safety
            - Phase 3: Testing of drug on patients to assess efficacy, effectiveness and safety (larger test group)
            - Phase 4: approved
    
    add_drug_mechanism : ``bool``, optional (default: ``True``)
        Grant subsequent filtering immunity to rows with drug mechanism information 
        present.
    assay_type : ``str`` or list of ``str``, optional (default: ``None``)
        Subset the data frame based on assay type information:
        
            - Binding (B) - Data measuring binding of compound to a molecular target, e.g. Ki, IC50, Kd.
            - Functional (F) - Data measuring the biological effect of a compound, e.g. %cell death in a cell line, rat weight.
            - ADMET (A) - ADME data e.g. t1/2, oral bioavailability.
            - Toxicity (T) - Data measuring toxicity of a compound, e.g., cytotoxicity.
            - Physicochemical (P) - Assays measuring physicochemical properties of the compounds in the absence of biological material e.g., chemical stability, solubility.
            - Unclassified (U) - A small proportion of assays cannot be classified into one of the above categories e.g., ratio of binding vs efficacy.
    
    remove_inactive : ``bool``, optional (default: ``True``)
        Subset the data frame to remove inactive drug-target interactions.
    include_active : ``bool``, optional (default: ``True``)
        Grant subsequent filtering immunity to active drug-target interactions.
    pchembl_target_column : ``str``, optional (default: ``None``)
        Use the selected column in the data frame to dictate custom pChEMBL thresholds 
        for each unique value in the column.
    pchembl_threshold : ``float`` or ``dict`` of ``float``, optional (default: ``None``)
        Subset the data frame to this pChEMBL minimum. If a single ``float``, use that 
        value. If a ``dict`` provided in conjunction with a ``pchembl_target_column``, 
        have the unique values of the specified column as keys of the dictionary, with 
        entries being the desired threshold for that category.
    '''
    
    #create a helper column that will be set to True as needed to protect rows from deletion
    dataframe['keep']=False
    
    # check boolean inputs
    if type(add_drug_mechanism) != bool:
        raise TypeError("argument 'add_drug_mechanism' should be 'bool' type")
    if type(remove_inactive) != bool:
        raise TypeError("argument 'remove_inactive' should be 'bool' type")
    if type(include_active) != bool:
        raise TypeError("argument 'include_active' should be 'bool' type")
    
    # filter based on drug max phase
    if drug_max_phase!=None:
        if type(drug_max_phase)==list:
            dataframe = dataframe[dataframe['molecule_dictionary|max_phase'].isin(drug_max_phase)]
        elif type(drug_max_phase)==int:
            dataframe = dataframe[dataframe['molecule_dictionary|max_phase']==drug_max_phase]
        else:
            raise TypeError("argument 'drug_max_phase' should be 'int' or 'list' type")
    
    # grant subsequent filtering immunity in the event of drug mechanism presence
    if add_drug_mechanism:
        drugmech_index = dataframe[dataframe['drug_mechanism|molregno'].notnull()].index
        dataframe.loc[drugmech_index,'keep']=True  # to keep 'drug_mechanism' activity
    
    # filter based on assay type
    if assay_type!=None:
        if type(assay_type)==list:
            dataframe = dataframe[(dataframe['assays|assay_type'].isin(assay_type))| \
                                        (dataframe['keep'])]
        elif type(assay_type)==str:
            dataframe = dataframe[(dataframe['assays|assay_type']==assay_type)| \
                                        (dataframe['keep'])]
        else:
            raise TypeError("argument 'assay_type' should be 'str' or 'list' type")
    
    # filter based on activity
    ## remove inactive
    if remove_inactive:
        dataframe = dataframe[(dataframe['activities|activity_comment'].isin(['inactive',
                                                                                     'Not Active',
                                                                                     'Not Active (inhibition < 50% @ 10 uM and thus dose-reponse curve not measured)',
                                                                                     'Inactive'])==False)| \
                                   (dataframe['keep'])] # keep 'drug_mechanism' activity if 'add_drug_mechanism' is True
    
    ## grant subsequent filtering immunity to active
    if include_active:
        active_index = dataframe[dataframe['activities|activity_comment'].isin(['active','Active'])].index
        dataframe.loc[active_index,'keep']=True  # to keep 'active' activity
        
    ## pChEMBL thresholding
    if pchembl_threshold!=None:
        if type(pchembl_threshold) == dict:
            if set(dataframe[pchembl_target_column])==set(pchembl_threshold.keys()):
                dataframe['pchembl_active']=False
                for k in pchembl_threshold.keys():
                    eachclass_df = dataframe[dataframe[pchembl_target_column]==k]
                    dataframe.loc[eachclass_df.index,'pchembl_active']=eachclass_df['activities|pchembl_value']>=pchembl_threshold[k]
                    del eachclass_df
            else:
                raise KeyError("argument 'pchembl_threshold' should have all the classes in key")
        else:
            #single value, do a single thresholding
            dataframe['pchembl_active'] = dataframe['activities|pchembl_value']>pchembl_threshold
        dataframe = dataframe[dataframe['pchembl_active'] | dataframe['keep']] # keep 'active' and 'drug_mechanism' activity if including those are True
    
    return dataframe
    

def create_dict(dataframe):
    
    dic = {}
    for c in set(dataframe['molecule_dictionary|chembl_id']):
        table_c = dataframe[dataframe['molecule_dictionary|chembl_id']==c]
        c_name = list(set(table_c['molecule_dictionary|pref_name']))

        if len(c_name)>1:
            raise ValueEror(f'{c} has multimle name: {c_name}')
        else:
            c_name = c_name[0]

        l=[]
        for x in table_c['component_synonyms|component_synonym']:
            l=l+x.split('|')
        dic[f'{c}|{c_name}']=list(set(l))
        del l, table_c, c_name

    return dic


def create_drug_dictionary(dataframe,
                           drug_grouping=None,
                           atc_level=['level1','level2']
                            ):
    '''
    This function creates drug:target-genenames dictionary
    * dataframe: drug-target activity dataframe 
    * drug_grouping: 'ATC_level','drug_max_phase',or None (default: None)
    * atc_level: specify level to categorise drugs and make drug-target dictionary. list. (default: ['level1','level2'])
    '''

    drug_target_dict={}
 
    if drug_grouping=='ATC_level':
        for l in atc_level:
            dataframe.fillna(value={f'atc_classification|{l}':'No-category'},inplace=True)
            for group in dataframe[f'atc_classification|{l}'].unique():
                # filtering category and creating dictionary
                group_df = dataframe[dataframe[f'atc_classification|{l}']==group]
                drug_target_dict[group] = create_dict(group_df)
                del group_df
    elif drug_grouping=='drug_max_phase':
        dataframe.fillna(value={'molecule_dictionary|max_phase':'No-phase'},inplace=True)
        for group in set(dataframe['molecule_dictionary|max_phase']):
            group_df = dataframe[dataframe['molecule_dictionary|max_phase']==group]
            drug_target_dict[f'maxphase_{str(group)}'] = create_dict(group_df)
            del group_df
    elif drug_grouping!=None:
        raise KeyError("invalid key <drug_grouping>")
    else:
        #no grouping, just make a drug:[targets] dict
        drug_target_dict = create_dict(dataframe)
    return drug_target_dict
    
    
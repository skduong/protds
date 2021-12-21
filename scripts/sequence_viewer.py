'''
A simplified version of protds_v3 that focuses only on structure visualizations
Faster than protds_v3 since it does not deal with binding sites or calculations
'''
#pypdb: access Protein Data Bank
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators import text_operators
#general utilities
import os, icn3dpy 
import pandas as pd
from IPython.display import display

proteins = {}

#Save Proteins dictionary
def saveProteins(fname = 'proteins_list.pkl'):
    with open(os.path.join("saves",fname), 'wb') as handle:
        pickle.dump(proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
def loadProteins(fname = 'proteins_list.pkl'):
    global proteins
    with open(os.path.join("saves",fname), 'rb') as handle: 
        savedProteins = pickle.load(handle)
    proteins.update(savedProteins)

def searchPDB(id): #add a new id to the dictionary 
    global proteins
    if id not in proteins: 
        #search PDB 
        search_operator = text_operators.ExactMatchOperator(
            value= id, 
            attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession")
        return_type = ReturnType.ENTRY
        try: 
            results = perform_search(search_operator, return_type)
            proteins[id] = results
        except: 
            print('There were no results for', id, '\n')    
            return False

def selectPDB(protein): #ask user to pick a structure if > 1 exist
    pdbs = proteins[protein]
    if len(pdbs) > 1: 
        print("\nProtein ", protein, " has ", len(pdbs), " structures:\n", ''.join(i+', ' for i in pdbs)[:-2], '\n', sep='')
        select = input("Choose one either by name or number (ex: type 2 to get the 2nd structure): ")
        try:
            val = int(select)
            if val <= len(pdbs) and val > 0:
                return int(val-1)
            else:
                print("Invalid selection, the first structure will be selected:\n")
                return 0   
        except ValueError:
            if select.upper() in pdbs:
                return pdbs.index(select.upper())
            else:
                print("Invalid PDB.", pdbs[0], "will be selected by default:\n")
                return 0 
    else:
        return 0
            
def getPepView(protid, data, table=True, colors = ['red', 'yellow', 'blue', 'green', 'cyan'], first=False): 
    proteinID = protid.split('-')[0] #disregard character after the uniprotID
    if proteinID not in proteins:
        searchPDB(proteinID)
        
    df = data[data['ProteinID'].str.contains(protid)].copy()
    readings = df.filter(regex='[0-9]min').apply(pd.to_numeric, errors='coerce').fillna(-1).values
    df['Color'] = [colors[-1] if max(i)<0 else colors[:-1][list(i).index(max(i))%(len(colors)-1)] for i in readings]
    if table: display(df)
        
    try:
        if first:
            pdb = proteins[proteinID][0]
        else:
            pdb = proteins[proteinID][selectPDB(proteinID)] 
       
        if len(df)>0:
            cmd = 'color A9A9A9;'
            for i in df[['PeptideSequence', 'Color']].values: cmd += 'select :'+ i[0] + '; color '+ i[1]+';'
            return icn3dpy.view(q='mmdbid='+pdb, command = cmd+';toggle highlight; view annotations; set view detailed view; set background white;')
        else:
            print("The input data did not contain", proteinID)
    except:
        print("Protein Data Bank did not have results for", proteinID)
from IPython.display import HTML
from pypdb import *
from pypdb.clients.search.search_client import perform_search
from pypdb.clients.search.search_client import ReturnType
from pypdb.clients.search.operators import text_operators
from Bio.SeqFeature import FeatureLocation
from IPython.display import display
from Bio import SwissProt
import math
import urllib
import icn3dpy
import warnings
warnings.filterwarnings(action='once')

class Protein: #each protein has a unique UniProtId, a list of PDB structures, and a UniProt record 
    def __init__(self, name, pdbs):
        self.UniProtId = name
        self.PDBids = pdbs
        self.record = None
    
    def getFeatures(self):
        return self.record.features
        
    def getModView(self, mod, selection): #gives PDB name, selected 3D structure with site highlights and recolored ModifiedLocationNumber
        if mod:
            mod = ''.join(str(i)+',' for i in mod)[:-1] 
            highlight = 'select :' + mod+';color white; select :'
        else:
            highlight = 'select :'
        for i in self.record.features:
            if 'BIND' in i.type or 'ACT' in i.type or 'METAL' in i.type or ('bind' in str(i.qualifiers) and 'CHAIN' not in i.type):  #ex: 5B3Z
                loc = FeatureLocation(i.location.start+1, i.location.end)
                highlight = highlight + str(loc)[1:-1].replace(':','-') + ',' #to fix the extra selection range
        return self.PDBids[selection], icn3dpy.view(q='mmdbid='+self.PDBids[selection], command = highlight[:-1]+';view annotations; set view detailed view')

def df_mod(data): # keep entries with modified sequence location:
    df = data[data['ModifiedLocationNum'].notna()].reset_index()
    df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int)
    df = df.groupby('ProteinID').agg(lambda x: list(x)) 
    return df
        
#searches
def get_uniprot (query='',query_type='ACC'): #for querying UniProtDB; returns an opened url to be used as swiss-prot handle
    url = 'https://www.uniprot.org/uploadlists/'
    params = {
    'from':query_type,
    'to':'ACC',
    'format': 'txt',
    'query':query
    }
    data = urllib.parse.urlencode(params)
    data = data.encode('ascii')
    request = urllib.request.Request(url, data)
    return urllib.request.urlopen(request, timeout=20) 

def searchPDB(id, proteins): #add a new id to the dictionary
    if id not in proteins: 
        #search PDB 
        search_operator = text_operators.ExactMatchOperator(value= id, attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession")
        return_type = ReturnType.POLYMER_ENTITY
        try: 
            results = perform_search(search_operator, return_type)
            proteins[id] = Protein(id, list(map(lambda x: x.split('_')[0], results)))
        except: 
            print('There were no results for ', id)    
            return id 
        #get sites from UniProt
        handle = get_uniprot(query = id, query_type = 'ACC')
        try:
            proteins[id].record = SwissProt.read(handle)
        except:
            print("UniProt did not have results for ", id)
    
#visualizations
def checkDB(protid, proteins, data, df): #check to see if the ProteinID exists
    if protid in proteins:
        return proteins[protid].getModView(selectLoc(protid, data, df), select3D(protid, proteins))
    else:
        print("No results\n")
        return '', icn3dpy.view()

def select3D(protid, proteins): #ask user to choose from multiple structures if applicable
    pdbs = proteins[protid].PDBids
    if len(pdbs) > 1:
        print("Protein ", protid, " has ", len(pdbs), " structures:\n", ''.join(i+', ' for i in pdbs)[:-2], '\n')
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

def selectLoc(protid, df, mods): #ask user to choose a row for ProteinIDs showing in multiple rows
    data = df[df['ProteinID']==protid]#.reset_index()
    print("Here are the rows containing", protid)
    display(data)
    choice = input("Type a row number to view, or type \"all\": ")
    try: 
        a = int(choice)
        if a in data.index:
            locnum = data.loc[a]['ModifiedLocationNum']
            if not math.isnan(locnum):
                print('ModifiedLocationNum',locnum, 'will be selected.\n')
                return [int(locnum)]
            else:
                print("No modified location present")
                return []
    except:
        if choice.lower() == 'all':
            return mods.loc[protid]['ModifiedLocationNum']
        else:
            print("Invalid response. No ModifiedLocationNum selected.\n")
            return []
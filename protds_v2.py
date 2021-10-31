#pypdb: access Protein Data Bank
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators import text_operators
#biotite: tool for working with structures and their coordinates 
import biotite.sequence as seq
import biotite.structure as struc
import biotite.database.rcsb as rcsb
import biotite.structure.io.mmtf as mmtf
#Bio: access UniProt's site information
from Bio.SeqFeature import FeatureLocation
from Bio import SwissProt
#general utilities
import numpy as np
import pandas as pd
from IPython.display import display
import os, math, urllib, icn3dpy, warnings
from tempfile import gettempdir, NamedTemporaryFile, TemporaryDirectory
warnings.filterwarnings(action='once')

proteins = {}

class PDB: #each ProteinDataBank entry that is searched is stored as its own object
    def __init__(self, id, struct):
        self.PDBid = id #four character name
        self.structure = struct #structure info & coordinates

class Protein: #each unique protein in the dataset is represented by a Protein object
    def __init__(self, name, pdbs):
        self.UniProtId = name #UniProtID
        self.PDBids = pdbs #list of PDB names [Redundant, will remove in future update]
        self.record = None #Uniprot features and sites
        self.structures = None #list of PDB objects
    
    def getPDBs(self): #give full list of PDBs associated with this Protein
        return [i.PDBid for i in self.structures]
    
    def getSites(self): #list the position(s) and label of sites 
        sites = []; labels = []
        for i in self.record.features:
            if 'BIND' in i.type or 'ACT' in i.type or 'METAL' in i.type or ('bind' in str(i.qualifiers) and 'CHAIN' not in i.type): 
                loc = FeatureLocation(i.location.start+1, i.location.end)
                sites.append(loc)
                labels.append(i.type)
        return list(map(lambda x: list(range(x.start, x.end+1)), sites)), labels
        
    def getStrucs(self): #list the structure(s) coordinates
        return [i.structure for i in self.structures]

#Search
def get_uniprot (query='',query_type='ACC'): #for querying UniProtDB; returns an opened url to be used as swiss-prot handle
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from':query_type, 'to':'ACC', 'format': 'txt', 'query':query}
    data = urllib.parse.urlencode(params).encode('ascii')
    request = urllib.request.Request(url, data)
    return urllib.request.urlopen(request, timeout=20)
    
def get_mmtf(pdbid): #get atom coordinates of a pdb structure from first mmtf model; returns a biotite AtomArray 
    with TemporaryDirectory() as tempdir:
        mmtf_file_path = rcsb.fetch(pdbid, "mmtf", os.path.join(gettempdir(), tempdir))
        mmtf_file = mmtf.MMTFFile.read(mmtf_file_path)
        structure = mmtf.get_structure(mmtf_file) 
    return structure[0]  

def searchPDB(id): #add a new id to the dictionary 
    global proteins
    if id not in proteins:
        #search PDB 
        search_operator = text_operators.ExactMatchOperator(value= id, attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession")
        return_type = ReturnType.POLYMER_ENTITY
        try: 
            results = perform_search(search_operator, return_type)
            pdbids = list(map(lambda x: x.split('_')[0], results))
            proteins[id] = Protein(id, pdbids)
            #get sites from UniProt
            handle = get_uniprot(query = id, query_type = 'ACC')
            try:
                proteins[id].record = SwissProt.read(handle)
                if proteins[id].getSites()[0]: #to save time, structure coordinates are only taken if binding sites exist for the Protein 
                    #store structure coordinates
                    proteins[id].structures = list(map(lambda x: PDB(x, get_mmtf(x)), pdbids)) 
            except Exception as e:
                print("UniProt did not have results for ", id)
                print(e)
        except: 
            print('There were no results for', id, '\n')    
            return False

#Chain Info
def chainSeq(chain, struc): #chain: a string indicating which chain to look at; struc: a biotite structure for some PDB; returns the sequence of a structure's chain
    aaList = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', #AminoAcids
              'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    return seq.ProteinSequence(np.take(struc[(struc.chain_id == chain) & (np.isin(struc.res_name, aaList))].res_name, 
            np.unique(struc[(struc.chain_id == chain) & (np.isin(struc.res_name, aaList))].res_id, return_index=True)[1]))

def checkChains(pdb, modstr): #modstr: ModifiedSequence of a row; returns a list of valid chains for a given pdb
    chainIDs = np.unique(pdb.chain_id)
    if len(chainIDs) == 1: #single chain
        return [pdb.chain_id[0]] #assuming that single chains will always contain the subsequence
    else:
        chainList = []
        for j in chainIDs:
            chain = str(chainSeq(j, pdb))
            if modstr.split('[')[0] in chain and modstr.split(']')[-1] in chain: #ModifiedSequence is a subsequence of this chain
                chainList.append(j)
    return chainList
    
def entryChains(entry): #returns a list of valid chains for all PDBs associated with an entry's protein
    return map(lambda x: checkChains(x, entry[2]), proteins[entry[0]].getStrucs())
    
#Distance
def calcDist(structure, chain, entry): #returns the distances between a PDB structure's binding sites and the entry's ModifiedLocationNum for a given chain
    arr, lab = proteins[entry[0]].getSites()
    return map(lambda x,y: (x,y), lab, map(lambda x: struc.distance(struc.mass_center(structure[(structure.chain_id==chain) & (structure.res_id==entry[1])]), struc.mass_center(x)), [structure[(structure.chain_id==chain) & (np.isin(structure.res_id, j))] for j in arr]))
    
def entryDist(entry): #get distances for all PDBs of single entry row
    return map(lambda x: [(j, list(calcDist(proteins[entry[0]].getStrucs()[x[0]], j, entry))) for j in x[1]], list(enumerate(entryChains(entry))))
    
def printDists(entry): #present the results of entryDist() nicely in tables
    for n, pdb in enumerate(entryDist(entry)):
        print(proteins[entry[0]].PDBids[n])
        if len(pdb) > 0:
            testdf = pd.DataFrame(pdb[0][1], columns = ['Type', 'ChainA'])
            if len(pdb)==1:
                display(testdf)
            else:
                for i in pdb[1:]:
                    testdf = pd.concat([testdf, pd.Series(data = map(lambda x: x[1], i[1]), name='Chain'+i[0])], axis=1)
                testdf['Average'] = testdf.mean(numeric_only=True, axis=1) #last column contains average of all chains' distances
                display(testdf)
        else:
            display(pd.DataFrame())
    
#Visualization
def selectPDB(protein): #ask user to pick a structure if > 1 exist
    pdbs = protein.PDBids
    if len(pdbs) > 1: 
        print("\nProtein ", protein.UniProtId, " has ", len(pdbs), " structures:\n", ''.join(i+', ' for i in pdbs)[:-2], '\n', sep='')
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
        
def entryView(entry): #get highlighted 3D view for a selected row (entry); returns an iCn3D view
    if entry[0] in proteins.keys(): #if protein exists 
        pdb = proteins[entry[0]].PDBids[selectPDB(proteins[entry[0]])] 
        if proteins[entry[0]].structures: #if structures are available
            cmd = ''
            sites = ''
            for i in proteins[entry[0]].record.features: #format binding sites
                        if 'BIND' in i.type or 'ACT' in i.type or 'METAL' in i.type or ('bind' in str(i.qualifiers) and 'CHAIN' not in i.type):  #ex: 5B3Z
                            loc = FeatureLocation(i.location.start+1, i.location.end)
                            sites = sites + str(loc)[1:-1].replace(':','-') + ',' #to fix the extra selection range
            for i in checkChains(pdb.structure, entry[2]): #highlight only on applicable chains
                cmd += 'select .{chain}:{mod}; color white; select .{chain}:'.format(chain = i, mod = str(entry[1]))+sites[:-1]+'; color FF0; '
            return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd+';toggle highlight; view annotations; set view detailed view')
        else:
            return icn3dpy.view(q='mmdbid='+pdb, command = 'select :'+str(entry[1])+'; color white;'+';toggle highlight; view annotations; set view detailed view')   
    else:
        print("No results for", entry[0])
        return icn3dpy.view()
        
def selectLoc(protid, df, mods): #ask user to choose a row for ProteinIDs showing in multiple rows [unused; pending removal]
    data = df[df['ProteinID']==protid]
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
            
# notebook use
def modData(data): #look at only modified rows/ convert location floats to integer 
    df = data[data['ModifiedLocationNum'].notna()] #select entries with ModifiedLocations
    df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int) #convert ModifiedLocationNum to integers
    return df
    
def processRow(rowIndex, df): #simplify user experience by having them only query using row indices
    entry = getEntry(rowIndex, df)
    if entry[0] not in proteins.keys():
        try:
            searchPDB(entry[0])
            return entry
        except: 
            print("Row couldn't be processed")
            return []
    else:
        return entry

def getEntry(rowNum, moddata): 
    try:
        rowList = moddata.loc[rowNum][['ProteinID', 'ModifiedLocationNum', 'ModifiedSequence']].values.tolist()
        rowList.append(moddata.loc[rowNum].name)
        return rowList
    except:
        print("Invalid Row Entry")
        return ['NA', 'NA', 'NA', 'NA']
    
def getProteins(data): #automate processing row-by-row and get modified entries; input: the dataset; output: list of modified rows
    global proteins #populate the dictionary of Proteins
    df = data[data['ModifiedLocationNum'].notna()].reset_index()
    df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int)
    for i in data['ProteinID'].unique()[:4]: #[TESTING] only first 4 proteins; [Have not tried iterating through all of them]
        searchPDB(i) 
    return list(filter(lambda x: x[0] in list(proteins.keys()) and proteins[x[0]].getSites()[0], df[['ProteinID', 'ModifiedLocationNum', 'ModifiedSequence', 'index']].values.tolist()))

def getDistances(entry, data):
    try:
        display(data[data['ModifiedLocationNum'].notna()].loc[entry[-1]].to_frame().T) #the only purpose of data is this display...
        printDists(entry)
    except Exception as e:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        elif not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
        else:
            print(e, "Invalid row entry") 

def getView(entry, data):
    try:
        display(data[data['ModifiedLocationNum'].notna()].loc[entry[-1]].to_frame().T) 
        return entryView(entry)
    except:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        elif not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
            return entryView(entry) 
        else:
            print("Invalid row entry") 
            
def getView_old(protid): #search by protid [OLD; pending removal]
    if protid in proteins.keys(): 
        pdbs = proteins[protid].PDBids
        cmd = ''
        sites = ''
        
        select = 0
        if len(pdbs) > 1: 
            print("\nProtein ", proteins[protid].UniProtId, " has ", len(pdbs), " structures:\n", ''.join(i+', ' for i in pdbs)[:-2], '\n', sep='')
            select = input("Choose one either by name or number (ex: type 2 to get the 2nd structure): ")
            try:
                val = int(select)
                if val <= len(pdbs) and val > 0:
                    select = int(val-1)
                else:
                    print("Invalid selection, the first structure will be selected:\n")
            except ValueError:
                if select.upper() in pdbs:
                    select =  pdbs.index(select.upper())
                else:
                    print("Invalid PDB.", pdbs[0], "will be selected by default:\n")
        
        for i in proteins[protid].record.features:
                    if 'BIND' in i.type or 'ACT' in i.type or 'METAL' in i.type or ('bind' in str(i.qualifiers) and 'CHAIN' not in i.type):  #ex: 5B3Z
                        loc = FeatureLocation(i.location.start+1, i.location.end)
                        sites = sites + str(loc)[1:-1].replace(':','-') + ',' #to fix the extra selection range
        cmd += 'select :' +sites[:-1]+ '; color FF0; '

        return icn3dpy.view(q='mmdbid='+pdbs[select], command = cmd+'; view annotations; set view detailed view')
        
    else:
        print("No results for", protid)
        return icn3dpy.view()
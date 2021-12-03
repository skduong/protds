#pypdb: access Protein Data Bank
from pypdb.clients.search.search_client import QueryGroup, LogicalOperator
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.search_client import perform_search_with_graph
from pypdb.clients.search.operators import text_operators
#biotite: tool for working with structure coordinates and sequences
import biotite.sequence as seq
import biotite.structure as struc
import biotite.database.rcsb as rcsb
import biotite.sequence.align as align
import biotite.structure.io.mmtf as mmtf
#Bio: access UniProt's site information
from Bio.SeqFeature import FeatureLocation
from Bio import SwissProt
#general utilities
import pickle
import numpy as np
import pandas as pd
from IPython.display import display
import os, math, urllib, icn3dpy, warnings
from tempfile import gettempdir, NamedTemporaryFile, TemporaryDirectory
warnings.filterwarnings(action='once')

proteins = {}

class PDB:
    def __init__(self, id, struct):
        self.PDBid = id
        self.structure = struct
        self.bestChain = None #[chain with highest score, alignment score for said chain]

class Protein: #each unique protein in the dataset is represented by a Protein object
    def __init__(self, name):
        self.UniProtId = name 
        self.record = None #Uniprot features and sites
        self.structures = None
    
    def getPDBs(self):
        return [i.PDBid for i in self.structures] 
    
    def getSequence(self):
        return seq.ProteinSequence(self.record.sequence.replace('U', 'X').replace('O', 'X'))
    
    def getSites(self): #lists the position(s) and label of sites 
        sites = []; labels = []
        for i in self.record.features:
            if 'BIND' in i.type or 'METAL' in i.type or 'SITE' in i.type or ('bind' in str(i.qualifiers) and 'REGION' in i.type): 
                sites.append(FeatureLocation(i.location.start+1, i.location.end))
                labels.append(i.type)
        return list(map(lambda x: list(range(x.start, x.end+1)), sites)), labels
        
    def getStrucs(self):
        return [i.structure for i in self.structures]

#Save Proteins dictionary
def saveProteins(fname = 'proteins.pkl'):
    with open(os.path.join(os.getcwd(),"saves",fname), 'wb') as handle:
        pickle.dump(proteins, handle, protocol=pickle.HIGHEST_PROTOCOL)
        
def loadProteins(fname = 'proteins.pkl'):
    global proteins
    with open(os.path.join(os.getcwd(),"saves",fname), 'rb') as handle: 
        savedProteins = pickle.load(handle)
    proteins.update(savedProteins)
        
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
        structure = mmtf.get_structure(mmtf_file, model=1) 
    return structure 
    
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
            if len(results)>15: 
                results = []
                for i in np.arange(1,5,0.5):
                    try: 
                        resolution_operator = text_operators.ComparisonOperator(
                            value=i, 
                            attribute="rcsb_entry_info.resolution_combined", 
                            comparison_type=text_operators.ComparisonType.LESS)
                        both = QueryGroup(queries = [search_operator, resolution_operator], logical_operator = LogicalOperator.AND)
                        results += perform_search_with_graph(query_object=both, return_type=return_type)
                        if len(results)>5: break 
                    except:
                        continue            
            proteins[id] = Protein(id)
            #get sites from UniProt
            handle = get_uniprot(query = id, query_type = 'ACC')
            try:
                proteins[id].record = SwissProt.read(handle)
                proteins[id].structures = list(map(lambda x: PDB(x, get_mmtf(x)), map(lambda x: x.split('_')[0], results))) #need list?
            except Exception as e:
                print("UniProt did not have results for ", id)
                print(e)
        except: 
            print('There were no results for', id, '\n')    
            return False

#Chain Info
def chainSeq(chain, pdbStruc): #chain: a string indicating which chain to look at; pdbStruc: a biotite structure for some PDB
    aaList = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN', 'CYS', 'SEC', 
              'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE', 'TYR', 'TRP']
    try:
        structure = struc.array(np.take(pdbStruc[(pdbStruc.chain_id == chain) & (np.isin(pdbStruc.res_name, aaList))], 
                np.unique(pdbStruc[(pdbStruc.chain_id == chain) & (np.isin(pdbStruc.res_name, aaList))].res_id, return_index=True)[1]))
        return seq.ProteinSequence(structure.res_name), structure.res_id
    except:
        return [],[]
        
def checkChains(pdb, protid): #assign the chain with highest alignment score
    if not pdb.bestChain:
        chainIDs = pd.unique(pdb.structure.chain_id)
        if len(chainIDs) == 1: #single chain
            pdb.bestChain = chainIDs[0], align.align_optimal(proteins[protid].getSequence(), chainSeq(chainIDs[0], pdb.structure)[0], align.SubstitutionMatrix.std_protein_matrix())[0].score
        else:
            scores = []
            chainIDs = [i for i in chainIDs if (len(pdb.structure[pdb.structure.chain_id == i].res_name[0]) == 3)] 
            for i in chainIDs:
                seq = chainSeq(i, pdb.structure)[0] 
                if len(seq)>0: 
                    ali = align.align_optimal(proteins[protid].getSequence(), seq, align.SubstitutionMatrix.std_protein_matrix())[0]
                    if align.get_sequence_identity(ali) > 0.5: #filter out poorly aligned chains
                        scores.append(ali.score) 
                    else:
                        scores.append(-99999)
            pdb.bestChain = chainIDs[np.argmax(scores)], max(scores)
    return pdb.bestChain

def bestPDB(protid): #return the structure with the best alignment to the row's Protein sequence
    return proteins[protid].structures[np.argmax([i[1] for i in [checkChains(x, protid) for x in proteins[protid].structures]])]
    
#Alignment
def alignLoc(locs, protid, pdb): 
    seq = chainSeq(checkChains(pdb, protid)[0], pdb.structure)
    ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix())[0]
    try:
        mainSeq = [i[0] for i in ali.trace]
        alignedLocs = [[list(map(lambda x: ali.trace[mainSeq.index(x-1)][1], location)) for location in loc] for loc in locs]
        if len(locs)>1: 
            missing = [[np.array(i[0])[np.where(np.array(i[1]) == -1)[0]] for i in zip(locs[0], alignedLocs[0])],
                       [i[0] for i in zip(locs[1], alignedLocs[1]) if i[1]==[-1]]]
        else: 
            missing = []
        return ([list(map(lambda x: seq[1][x], [list(map(lambda x: ali.trace[mainSeq.index(x-1)][1], 
                filter(lambda x: ali.trace[mainSeq.index(x-1)][1] > -1, location))) for location in loc])) for loc in locs], missing)
    except Exception as e:
        print(e, "Something went wrong")
        return []

#Distance
def calcDist(pdb, chain, entry): #returns the distances between a PDB structure's binding sites and the entry's ModifiedLocationNum for a given chain
    structure = pdb.structure
    uniProtSites = proteins[entry[0]].getSites()
    alignedSites = alignLoc([uniProtSites[0], [[entry[1]]]], entry[0], pdb)
    if not alignedSites[1][1]: 
        return map(lambda x,y,z: (x,y,z), uniProtSites[1], map(lambda x: struc.distance(struc.mass_center(structure[(structure.chain_id==chain) & (structure.res_id==alignedSites[0][1][0])]), struc.mass_center(x)), [structure[(structure.chain_id==chain) & (np.isin(structure.res_id, j))] for j in alignedSites[0][0]]), alignedSites[1][0])
    else:
        return []
        
def printDists(entry, best): #display the results of entryDists with tables
    bestStruc = bestPDB(entry[0])
    for pdbDist in [[(i.bestChain[0], calcDist(i, i.bestChain[0], entry)), i.PDBid] for i in ([bestStruc] if best else proteins[entry[0]].structures)]:
        print(pdbDist[1])
        if isinstance(pdbDist[0][1], map) and len(pdbDist)>0:
            df = pd.DataFrame(pdbDist[0][1], columns = ['Type', 'Chain'+pdbDist[0][0], 'Missing'])
            df['Location'] = locToStr(proteins[entry[0]].getSites()[0])
            df = df[['Type', 'Location', 'Chain'+pdbDist[0][0], 'Missing']]
            display(df)
        else:
            print("***ModifiedLocationNum", entry[1], "is missing from the structure***")
            display(pd.DataFrame())

def dfDists(data): #calculate and append sites and distance information to the dataset [IN PROGRESS]
    df = data.copy()[['ProteinID', 'ModifiedLocation', 'ModifiedLocationNum']]
    noResults = [i for i in pd.unique(df['ProteinID']) if searchPDB(i)==False]
    df['HasResults'] = ~df['ProteinID'].isin(noResults)
    entries = df[['ProteinID', 'HasResults']].values.tolist() 
    #sites
    df['Sites'] = [[locToStr([i])[0] for i in proteins[e[0]].getSites()[0]] if e[1] else [] for e in entries]
    df['Type'] = [[i for i in proteins[e[0]].getSites()[1]] if e[1] else [] for e in entries]
    df['NumOfSites'] = [len(i) for i in df['Sites']]
    entries = df[['ProteinID', 'ModifiedLocationNum', 'ModifiedLocation', 'NumOfSites']].values.tolist()
    #distances
    missing=[]; dists=[]
    for entry in entries:
        if entry[3]>0: 
            struc = bestPDB(entry[0])
            dist = list(calcDist(struc, struc.bestChain[0], entry))
            if len(dist)>0:
                dists.append([round(i[1],3) for i in dist])
                missing.append([i[2] for i in dist])
            else: 
                dists.append([np.NaN for i in range(entry[3])])
                missing.append([[] for i in range(entry[3])]) 
                entry[2] += '*'
        else:
            dists.append([])
            missing.append([])
    df['Distances'] = dists
    df['ModifiedLocation'] = [i[2] for i in entries]
    df['MissingSites'] = missing
    return df.drop(columns=['ModifiedLocationNum'])[['ProteinID', 'ModifiedLocation', 'HasResults', 'NumOfSites',
                                                     'Type', 'Sites', 'Distances', 'MissingSites']]
            
#Visualization
def locToStr(locList): #locList = [[L1,L2], [L3], [etc]] 
    return ['-'.join(loc) for loc in [pd.unique([str(i[0]), str(i[-1])]).tolist() for i in locList if len(i)>0] if '[]' not in loc]

def pdbView(protid, loc1, loc2, chainSelect, full): #highlight iCn3D structure with locations from loc1(red) and loc2(black)
    if protid in proteins.keys(): 
        pdb = proteins[protid].structures[selectPDB(proteins[protid])]
        chainSelect = '.'+checkChains(pdb, protid)[0] if chainSelect=='best' else ''
        viewSelect = '' if full else 'select {chain}; show selection;'.format(chain = chainSelect) 
        settings = ';toggle highlight; view annotations; set view detailed view; set background white;'
        if len(loc1)>0 and loc2[0][0]: #both sites exist
            sites1, sites2 = [locToStr(i) for i in alignLoc([loc1, loc2], protid, pdb)[0]]
            if len(sites1)>0 and len(sites2)>0: #they both have sites after alignment
                cmd = viewSelect + 'select {chain}:{s1}; color red; select {chain}:'.format(chain = chainSelect, s1 = ','.join(
                [i for i in sites1]))+ ','.join([i for i in sites2])+'; color 000000;'
            elif len(sites1)>0:
                cmd = viewSelect + 'select {chain}:{s1}; color red;'.format(chain = chainSelect, s1 = ','.join([i for i in sites1]))
            elif len(sites2)>0:
                cmd = viewSelect + 'select {chain}:{s2}; color 000000;'.format(chain = chainSelect, s2 = ','.join([i for i in sites2]))
            else:
                cmd = ''
            return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd + settings)
            
        elif len(loc1)>0 or loc2[0][0]: #one of them exists
            sites = locToStr(alignLoc([loc1], protid, pdb)[0][0]) if len(loc1)>0 else locToStr(alignLoc([loc2], protid, pdb)[0][0])
            cmd = viewSelect + 'select {chain}:{s1}; color 000000;'.format(chain = chainSelect, s1 = ','.join([i for i in sites])) if len(sites)>0 else ''
            return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd + settings)
            
        else: #none exist, just load the structure with no highlights
            print("No sites detected")
            return icn3dpy.view(q='mmdbid='+pdb.PDBid)
    else:
        print("No results for", protid)
        return icn3dpy.view(command = 'set background white')
        
def selectPDB(protein): #ask user to pick a structure if > 1 exist
    pdbs = protein.getPDBs()
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

# notebook use
def modData(data): #look at only modified rows/ convert location floats to integer 
    df = data[data['ModifiedLocationNum'].notna()].copy()
    df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int)
    return df
    
def getProteins(data, load_pickle=False, save_pickle=False, fname = 'proteins.pkl'):
    if load_pickle: loadProteins() 
    for i in data['ProteinID'].unique(): 
        searchPDB(i)
    if save_pickle: saveProteins() 
    #df = data[data['ModifiedLocationNum'].notna()].reset_index()
    #df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int)
    #return list(filter(lambda x: x[0] in proteins.keys(), df[['ProteinID', 'ModifiedLocationNum', 'index']].values.tolist()))

def printSites(prot): #display active sites for a given protein
    display(pd.DataFrame(list(zip(prot.getSites()[1], prot.getSites()[0])), columns=['Type', 'Location']))

def checkRow(row): #verify rows and return [ProteinID, ModifiedLocationNum, Index] entries for easy processing
    try:
        if not math.isnan(row['ModifiedLocationNum']):
            entry = [row['ProteinID'], row['ModifiedLocationNum'].astype(int), row.name]
        else:
            entry = [row['ProteinID'], None, row.name]
        if entry[0] not in proteins.keys():
            try:
                searchPDB(entry[0])
                return entry
            except:
                print("Search for", entry[0], 'failed.')
                return ['NA', 'NA', 'NA']
        else: 
            return entry
    except Exception as e: 
        print(e, "Invalid row entry")
        return ['NA', 'NA', 'NA']    

def getPepView(protid, data, table=True, colors = ['red', 'yellow', 'blue', 'green', 'cyan'], first=False): 
    proteinID = protid.replace(';','-').split('-')[0] #disregard character after the uniprotID
    if proteinID not in proteins:
        searchPDB(proteinID)
        
    df = data[data['ProteinID'].str.contains(protid)].copy()
    readings = df.filter(regex='[0-9]min').apply(pd.to_numeric, errors='coerce').fillna(-1).values
    df['Color'] = [colors[-1] if max(i)<0 else colors[:-1][list(i).index(max(i))%(len(colors)-1)] for i in readings]
    if table: display(df)
        
    try:
        if first:
            pdb = proteins[proteinID].structures[0]
        else:
            pdb = proteins[proteinID].structures[selectPDB(proteins[proteinID])] 
       
        if len(df)>0:
            cmd = 'color A9A9A9;'
            for i in df[['PeptideSequence', 'Color']].values: cmd += 'select :'+ i[0] + '; color '+ i[1]+';'
            return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd+';toggle highlight; view annotations; set view detailed view; set background white;')
        else:
            print("The input data did not contain", proteinID)
    except:
        print("Protein Data Bank did not have results for", proteinID)
        
def getDistances(row, best=True): #error handling for user input before printing
    entry = checkRow(row)
    try:
        display(row.to_frame().T)
        if entry[1]:
            if not proteins[entry[0]].getSites()[0]:
                print(entry[0], "did not have binding sites listed in the UniProt databases")
            else:
                printDists(entry, best)
        else:
            print("No ModifiedLocationNum found")
    except Exception as e:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        else:
            print(e, "Invalid row entry")  
 
def getEntryView(row, chainSelect='best', full=True): 
    entry = checkRow(row)
    try:
        display(row.to_frame().T)
        if not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
        return pdbView(entry[0], proteins[entry[0]].getSites()[0], [[entry[1]]], chainSelect, full)
    except Exception as e:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        else:
            print(e,"Invalid row entry")

def getOutput(data, filename): #will most-likely fail for really large datasets [TESTING]
    moddata = data[data['ModifiedLocationNum'].notna()].copy() #include only modified entries
    moddata['ModifiedLocationNum'] = moddata['ModifiedLocationNum'].astype(int)
    try:
        outTable = dfDists(moddata)
        #long table
        outTable.explode(['Type', 'Sites', 'Distances', 'MissingSites']).to_csv(os.path.join(os.getcwd(),"output",filename+'_long.csv')) 
        #summary table
        outTable['MissingSites'] = [[x for x in i if len(x)>0] for i in outTable['MissingSites']]
        outTable.to_csv(os.path.join(os.getcwd(),"output",filename+'_summary.csv')) 
    except Exception as e:
        print(e)
            
#Email Requests
def emailRow(row):
    return [row['ProteinID'], row['ModifiedLocationNum'].astype(int), row.name]
def emailView(protid, loc1, loc2, best=True): #highlights all chains, assumes 2 sets of locations provided by datasets
    entry = [protid, 1, 2]
    if protid in proteins.keys(): 
        pdb = proteins[entry[0]].structures[selectPDB(proteins[entry[0]])]
        chainSelect = '.'+checkChains(pdb, protid)[0] if best else ''
        sites1, sites2 = [locToStr(i) for i in alignLoc([loc1, loc2], entry[0], pdb)[0]]
        if len(sites1)>0 and len(sites2)>0: #they both have sites after alignment
             cmd = 'select {chain}:{s1}; color 000000; select {chain}:'.format(chain = chainSelect, s1 = ','.join(
                [i for i in sites1]))+ ','.join([i for i in sites2])+'; color FFD700;'
        else:
            cmd = ''
        return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd+';toggle highlight; view annotations; set view detailed view; set background white')
    else:
        print("No results for", entry[0])
        return icn3dpy.view()
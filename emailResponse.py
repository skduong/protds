from protds_v2 import *

#New requests from email
#(Mostly minor modifications to take into consideration the lack of ModifiedSequence)
def emailView(protid, modlist): #accommodating the new task from the email; entry = [ProteinID, [list of locs1], [list of locs2]]
    if protid in proteins.keys(): #if protein exists 
        pdb = proteins[protid].PDBids[selectPDB(proteins[protid])]
        cmd = ''
        for i in modlist[0]: #modified locations listed in 1st dataset
            cmd += 'select :{modloc}; color white;'.format(modloc=i)
        for j in modlist[1]: #locations to highlight from 2nd dataset
            cmd += 'select :{highlightloc}; color yellow;'.format(highlightloc=j)
        return icn3dpy.view(q='mmdbid='+pdb, command = cmd +';toggle highlight; view annotations; set view detailed view')
    else:
        print("No results for", protid)
        return icn3dpy.view()
        
def checkChains_loc(pdb, modloc): #include chain if it's within range of modified location
    chainIDs = np.unique(pdb.chain_id)
    if len(chainIDs) == 1: #single chain
        return [pdb.chain_id[0]] #assuming that single chains will always contain the subsequence
    else:
        chainList = []
        for j in chainIDs:
            if modloc in pdb[pdb.chain_id == j].res_id:
                chainList.append(j)
    return chainList
    
def emailEntryChains(entry): #returns a list of valid chains for all PDBs associated with an entry's protein
    return map(lambda x: checkChains_loc(x, entry[1]), proteins[entry[0]].getStrucs())
        
def emailCalcDist(structure, chain, entry): #entry = [ProtID, modloc]
    arr, lab = proteins[entry[0]].getSites()
    return map(lambda x,y: (x,y), lab, map(lambda x: struc.distance(struc.mass_center(structure[(structure.chain_id==chain) & (structure.res_id==entry[1])]), struc.mass_center(x)), [structure[(structure.chain_id==chain) & (np.isin(structure.res_id, j))] for j in arr]))        

def emailEntryDist(entry): #get distances for all PDBs of single entry row
    return map(lambda x: [(j, list(emailCalcDist(proteins[entry[0]].getStrucs()[x[0]], j, entry))) for j in x[1]], list(enumerate(emailEntryChains(entry))))
 
def emailPrintDists(entry): #present the results of entryDist() nicely in tables
    for n, pdb in enumerate(emailEntryDist(entry)):
        print(proteins[entry[0]].PDBids[n])
        if len(pdb) > 0:
            df = pd.DataFrame(pdb[0][1], columns = ['Type', 'ChainA'])
            df['Location'] = ['-'.join(loc) for loc in [np.unique([str(i[0]), str(i[-1])]).tolist() for i in proteins[entry[0]].getSites()[0]]]
            df = df[['Type', 'Location', 'ChainA']]
            if len(pdb)==1:
                display(df)
            else:
                for i in pdb[1:]:
                    df = pd.concat([df, pd.Series(data = map(lambda x: x[1], i[1]), name='Chain'+i[0])], axis=1)
                df['Average'] = df.mean(numeric_only=True, axis=1) #last column contains average of all chains' distances
                display(df)
        else:
            display(pd.DataFrame())
            
def emailCheckRow(row):
    try:
        entry = [row['ProteinID'], row['ModifiedLocationNum'].astype(int), row.name]
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
    
def emailGetDistances(row): #error handling for user input before printing
    entry = emailCheckRow(row)
    try:
        display(row.to_frame().T)
        emailPrintDists(entry)
    except:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        elif not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
        else:
            print("Invalid row entry")     
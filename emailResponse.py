from protds_v2 import *

#New requests from email
#(Mostly minor modifications to take into consideration the lack of ModifiedSequence)
def emailGetEntry(rowNum, moddata): #returns [ProtID, modloc, index]
    try:
        rowList = moddata.loc[rowNum][['ProteinID', 'ModifiedLocationNum']].values.tolist()
        rowList.append(moddata.loc[rowNum].name)
        return rowList
    except:
        print("Invalid Row Entry")
        return ['NA', 'NA', 'NA']

def emailProcessRow(rowIndex, df): #simplify user experience by having them only query using row indices
    entry = emailGetEntry(rowIndex, df)
    if entry[0] not in proteins.keys():
        try:
            searchPDB(entry[0])
            return entry
        except: 
            print("Row couldn't be processed")
            return []
    else:
        return entry
 
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
        
def emailEntryChains(entry): #no chain info; assuming all 
    return map(lambda x: np.unique(x.chain_id).tolist(), proteins[entry[0]].getStrucs())
        
def emailCalcDist(structure, chain, entry): #entry = [ProtID, modloc]
    arr, lab = proteins[entry[0]].getSites()
    return map(lambda x,y: (x,y), lab, map(lambda x: struc.distance(struc.mass_center(structure[(structure.chain_id==chain) & (structure.res_id==entry[1])]), struc.mass_center(x)), [structure[(structure.chain_id==chain) & (np.isin(structure.res_id, j))] for j in arr]))        

def emailEntryDist(entry): #get distances for all PDBs of single entry row
    return map(lambda x: [(j, list(emailCalcDist(proteins[entry[0]].getStrucs()[x[0]], j, entry))) for j in x[1]], list(enumerate(emailEntryChains(entry))))
 
def emailPrintDists(entry): #present the results of entryDist() nicely in tables
    for n, pdb in enumerate(emailEntryDist(entry)):
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

def emailGetDistances(entry, data): #called by user to get distances
    try:
        display(data[data['ModifiedLocationNum'].notna()].loc[entry[-1]].to_frame().T) 
        emailPrintDists(entry)
    except:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        elif not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
        else:
            print("Invalid row entry") 
'''
This module provides get-functions and error handling to simplify running the program on Jupyter Notebook. 
'''
def modData(data): #look at only modified rows/ convert location floats to integer 
    df = data[data['ModifiedLocationNum'].notna()].copy()
    df['ModifiedLocationNum'] = df['ModifiedLocationNum'].astype(int)
    return df
    
def getProteins(data, load_pickle=False, save_pickle=False, fname = 'proteins.pkl'): 
    if load_pickle: loadProteins() 
    for i in data['ProteinID'].unique(): 
        searchPDB(i)
    if save_pickle: saveProteins() 

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

def getEntryView(row, showChain='best', full=True): #set showChain='' to select locations on every chain, or specify a chain letter to choose one other than the 'best' aligned one
    entry = checkRow(row)
    try:
        display(row.to_frame().T)
        if not proteins[entry[0]].getSites()[0]:
            print(entry[0], "did not have binding sites listed in the UniProt databases")
        return pdbView(entry[0], proteins[entry[0]].getSites()[0], [[entry[1]]], showChain, full)
    except Exception as e:
        if entry[0] not in proteins.keys():
            print("There were no results for", entry[0])
        else:
            print(e,"Invalid row entry")

def getOutput(data, filename, filepath=None):
    moddata = data[data['ModifiedLocationNum'].notna()].copy() #include only modified entries
    moddata['ModifiedLocationNum'] = moddata['ModifiedLocationNum'].astype(int)
    if not filepath: filepath = os.path.join(os.getcwd(),"output")
    try:
        outTable = dfDists(moddata)
        #long table
        outTable.explode(['Type', 'Sites', 'Distances', 'MissingSites']).to_csv(os.path.join(filepath,filename+'_long.csv'))
        #summary table
        outTable['MissingSites'] = [[x.tolist() for x in i if len(x)>0] for i in outTable['MissingSites']]
        outTable.to_csv(os.path.join(filepath,filename+'_summary.csv'))
    except Exception as e:
        print(e)
            
'''
Run as a Python program:
If users only want the output files, they can also run this program on terminal instead of opening a Notebook.
Creates and saves a "summary" and "long" csv file for the input data, in the same location
'''
def inputCheck(instr): #handle quotations in input responses
    checkstr = instr.replace('"', "'").split("'", 1)
    if len(checkstr)==1: 
        return checkstr[0]
    else:
        return checkstr[1][:-1]

def processData(filepath): 
    #input handling
    if filepath==None:
        filepath = inputCheck(input("Enter the path of the csv file's folder location: "))
        
    csvfiles = [i for i in os.listdir(filepath) if '.csv' in i]
    if len(csvfiles)==1:
        csvname = csvfiles[0]
    elif len(csvfiles)>1:
        csvname = inputCheck(input("More than one csv file found in that directory. Please enter the file's name: "))
        if '.csv' not in csvname: csvname += '.csv'
    else:
        print("Error. No csv files detected in", filepath)
        return 0
    
    try:
        data = pd.read_csv(os.path.join(filepath,csvname))
        getOutput(data, csvname[:-4], filepath)
        print("Done. Output is saved in", filepath)
    except Exception as e:
        print("Failed to process that data file.", e)
    
#(testing miscellaneous code and requests)
def emailRow(row):
    return [row['ProteinID'], row['ModifiedLocationNum'].astype(int), row.name]
def emailView(protid, loc1, loc2, bestChain=True, chooseStruc=False): #assumes 2 sets of locations provided by datasets
    if protid in proteins.keys(): 
        pdb = proteins[protid].structures[selectPDB(proteins[protid])] if chooseStruc else bestPDB(protid) 
        chainSelect = '.'+checkChains(pdb, protid)[0] if bestChain else ''
        sites1, sites2 = [locToStr(i) for i in alignLoc([loc1, loc2], protid, pdb)[0]]
        if len(sites1)>0 and len(sites2)>0: #they both have sites after alignment
             cmd = 'select {chain}:{s1}; color 000000; select {chain}:'.format(chain = chainSelect, s1 = ','.join(
                [i for i in sites1]))+ ','.join([i for i in sites2])+'; color FFD700;'
        else:
            cmd = ''
        return icn3dpy.view(q='mmdbid='+pdb.PDBid, command = cmd+';toggle highlight; view annotations; set view detailed view; set background white')
    else:
        print("No results for", protid)
        return icn3dpy.view()

if __name__ == "__main__":
    from protds import *
    path = None #can hard-code a path to avoid being asked every run
    #loadProteins()
    processData(path)
    
else: #imported by a Jupyter Notebook
    from protds.protds import *
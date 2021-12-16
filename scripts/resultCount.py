'''
If you encounter: ModuleNotFoundError: No module named 'pypdb', install pypdb by entering this command in the terminal:
pip install git+git://github.com/williamgilpin/pypdb
'''
#pypdb: access Protein Data Bank
from pypdb.clients.search.search_client import perform_search, ReturnType
from pypdb.clients.search.operators import text_operators
import pandas as pd
import os

def inputCheck(instr): #handle quotations in input responses
    checkstr = instr.replace("'", '"').split('"')
    if len(checkstr)==1: 
        return checkstr[0]
    else:
        return checkstr[1]
        
def getResultCount():
    '''
    Input: a FASTA file
    Output: a csv with each ProteinID's number of results found in ProteinDataBank
            results are saved in the same directory as the FASTA file
    '''
    #input handling
    pathname = inputCheck(input("Enter the path of the FASTA file's folder location: "))
    fastafiles = [i for i in os.listdir(pathname) if '.fasta' in i]

    if len(fastafiles)==1:
        fastaname = fastafiles[0]
    elif len(fastafiles)>1:
        fastaname = inputCheck(input("More than one FASTA file found in that directory. Please enter the file's name: "))
        if '.fasta' not in fastaname: fastaname += '.fasta'
    else:
        print("Error. No FASTA files detected in", pathname)
        return 0

    try:
        data = open(os.path.join(pathname,fastaname), "r").read()
    except Exception as e:
        print("Failed to read that file.", e)
        return 0
            
    #getting result count
    pids = [i[:6] for i in data.split(">sp|")][1:]
    numResults = {}
    for i in pids:
        try:
            search_operator = text_operators.ExactMatchOperator(value= i, attribute="rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession")
            return_type = ReturnType.ENTRY
            results = perform_search(search_operator, return_type)
            numResults[i] = len(results)
        except:
            numResults[i] = 0

    #save as csv
    dfResults = pd.DataFrame(numResults.items(), columns=['ProteinID', 'NumOfResults'])
    dfResults.to_csv(os.path.join(pathname,fastaname[:-6]+'_count.csv'), index=False)

    #print results
    total = len(dfResults)
    print("Summary:")
    print("no results:", len(dfResults[dfResults['NumOfResults']==0]),'/',total,'=', len(dfResults[dfResults['NumOfResults']==0])/total)
    counts = [0,10,15,20,50,100]
    for i in range(len(counts)-1):
        match = len(dfResults[(dfResults['NumOfResults']>counts[i]) & (dfResults['NumOfResults']<=counts[i+1])])
        print(counts[i]+1, 'to', counts[i+1], "results:", match,'/', total, '=', match/total)
    print("over 100 results:", len(dfResults[dfResults['NumOfResults']>100]),'/',total, '=', len(dfResults[dfResults['NumOfResults']>100])/total)

    print("\nCumulative:")
    for i in [10,15,20,50,100]:
        print("<", i, "results:", len(dfResults[dfResults['NumOfResults']<i]),'/', total,'=', len(dfResults[dfResults['NumOfResults']<i])/total)
    #print(">=100 results:", len(dfResults[dfResults['NumOfResults']>=100]),'/',total,'=', len(dfResults[dfResults['NumOfResults']>=100])/total)
    
if __name__ == "__main__":
    getResultCount()
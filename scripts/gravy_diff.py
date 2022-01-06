import os, urllib
import pandas as pd
from Bio import SeqIO, SwissProt
import biotite.sequence as seq

def get_uniprot (query='',query_type='ACC'): #for querying UniProtDB; returns an opened url to be used as swiss-prot handle
    url = 'https://www.uniprot.org/uploadlists/'
    params = {'from':query_type, 'to':'ACC', 'format': 'txt', 'query':query}
    data = urllib.parse.urlencode(params).encode('ascii')
    request = urllib.request.Request(url, data)
    return urllib.request.urlopen(request, timeout=20)

def fetchSequence(upid, seqlist): #get sequences from UniProt
    if upid not in seqlist.keys():
        try:
            seqlist[upid] = SwissProt.read(get_uniprot(query = upid, query_type = 'ACC')).sequence
            return seqlist[upid]
        except:
            return None
    else:
        return seqlist[upid]
    
def getSequence(data, fasta_file_path=''): #look for FASTA file to pull sequences from first, if it doesn't exists, then get them from UniProt records instead
    df = data.copy()
    ids = [i.split(';')[0] for i in df['ProteinID'].values]
    seqs = {}
    try:
        with open(fasta_file_path) as fasta_file: 
            record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        record_dict = dict(zip([i.split('|')[1] for i in record_dict.keys()], [str(j.seq) for j in record_dict.values()]))
        df['ProteinSequence'] = [record_dict[i] if i in record_dict else fetchSequence(i, seqs) for i in ids]
    except Exception as e: #not recommended to rely on fetching sequences from UniProt; this is a fail-safe
        print("Unable to use FASTA, so sequences will be taken from UniProt records. This process will take a while.", e)
        df['ProteinSequence'] = [seqs[id] if id in seqs.keys() else fetchSequence(id, seqs) for id in ids]
    return df

def calculate_gravy(sequence): #GRAVY score is sum of hydropathy divided by length of the sequence
    hydropathy_index = {'I':4.5, 'V':4.2, 'L':3.8,'F':2.8, 'C':2.5,'M': 1.9, 'A': 1.8, 'G': -.4, 'T': -.7, 'W': -.9, 'S': -.8,
                        'Y':-1.3, 'P':-1.6, 'H':-3.2, 'E':-3.5, 'Q':-3.5, 'D': -3.5, 'N':-3.5, 'K': -3.9, 'R': -4.5}
    try:
        hydropathy = sum(hydropathy_index[letter] for letter in sequence if letter in hydropathy_index.keys())
        gravy_score = hydropathy/len(sequence) 
        return gravy_score
    except:
        return None
    
def process(data, fastaPath, protIDs = "PG.UniProtIds", pepSeqs = "PEP.StrippedSequence"):
    '''
    The input data should have columns named "ProteinID" and "PeptideSequence" 
    By default, it will be assumed that "PG.UniProtIds" and "PEP.StrippedSequence" are possible column names (change these to accommodate your dataset)
    this function gives the data columns a consistant naming convention and adds a column of ProteinSequence if it doesn't exist
    Output: dataset sorted by ProteinIDs with the following additional columns:
        - PeptideGRAVY: the GRAVY index of the peptide subsequence
        - SequenceGRAVY: the GRAVY index of the ProteinSequence
    '''
    #preprocessing 
    if "ProteinID" not in data.columns or "PeptideSequence" not in data.columns: #handle the lack of "ProteinID" or "PeptideSequence"
        df = data.rename(columns={protIDs: "ProteinID", pepSeqs: "PeptideSequence"})
    else: df = data.copy()
    if 'ProteinSequence' not in df.columns: #add column of "ProteinSequence" if not available
        df = getSequence(df, fastaPath)
    #GRAVY 
    df['SequenceGRAVY'] = [calculate_gravy(seq) for seq in df['ProteinSequence'].values]
    df['PeptideGRAVY'] = [calculate_gravy(seq) for seq in df['PeptideSequence'].values]
    
    return df.sort_values(by='ProteinID')
    
def gravyDiff(sortedData): #GRAVY of proteinSequence-peptideSubsequences
    gravDiff = [] 
    for p in sortedData.groupby('ProteinID'):
        gravDiff += [p[1]['SequenceGRAVY'].values[0] - sum(p[1]['PeptideGRAVY']) for i in range(len(p[1]))]
    sortedData['GRAVYdifference'] = gravDiff
    return sortedData
    
def gravyDiff2(table1, table2): #GRAVY of PeptideSequences from Table1 minus the ones also found in Table2
    #assuming table1 has gone through process(table1), where its columns have been renamed, its GRAVY has been found, and it has been sorted by ProteinID
    #assuming table2's columns will have the following names: (edit as needed, but it's recommended to just have a "ProteinID" and "PeptideSequence" column)
    table2Proteins = "PG.UniProtIds"
    table2PeptideSeq = "Stripped.Sequence"
    if "ProteinID" in table2.columns: table2Proteins = "ProteinID"
    if "PeptideSequence" in table2.columns: table2PeptideSeq = "PeptideSequence"

    #handling isomers
    table1['UPID'] = [i.split(';')[0] for i in table1["ProteinID"].values]
    table2['UPID'] = [i.split(';')[0] for i in table2[table2Proteins].values]
    
    d2 = dict(list(table2.groupby('UPID')))
    diff2 = []
    for i in table1.groupby("UPID"):
        if i[0] in d2.keys():
            minustable2 = i[1][~i[1]['PeptideSequence'].isin(d2[i[0]][table2PeptideSeq].values)]
            diff2 += [sum(minustable2['PeptideGRAVY']) for j in range(len(i[1]))]
        else:
            diff2 += [sum(i[1]['PeptideGRAVY']) for j in range(len(i[1]))]
    table1['GRAVYdifference2'] = diff2
    return table1.drop('UPID', axis=1)
    
def pepPositions(df): #adds a column of PepMid: the position number of the Peptide's middle letter (closer to the left) on the full sequence
    seqs = df[['PeptideSequence', 'ProteinSequence']].values
    start=[]; end=[]
    for i in seqs: 
        begin = seq.find_subsequence(seq.ProteinSequence(i[1].replace('U','X')), seq.ProteinSequence(i[0].replace('U','X'))) if i[1] else pd.array([])
        if begin.size>0:
            stop = begin+len(i[0])
            start.append(begin+1)
            end.append(stop)
        else:
            start.append(None)
            end.append(None)
    df['PepStart'] = start
    df['PepEnd'] = end
    #If a PeptideSequence appears on multiple locations, they will be looked at individually (one start|end per row)
    df = df[df['PepStart'].notna()]
    df = df.explode(['PepStart', 'PepEnd']) 
    df[['PepStart','PepEnd']] = df[['PepStart','PepEnd']].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
    df['PepMid'] = [i[0]+int((i[1]-i[0])/2) for i in df[['PepStart', 'PepEnd']].values]
    
    return df.drop(['PepStart', 'PepEnd'], axis=1)
    
def inputCheck(instr): #handle quotations in input responses
    checkstr = instr.replace('"', "'").split("'", 1)
    if len(checkstr)==1: 
        return checkstr[0]
    else:
        return checkstr[1][:-1]
        
def getGRAVYdiffs(fastaPath, table1Path, table2Path, save = True):
    if ".fasta" not in fastaPath: fastaPath+=".fasta"
    if ".csv" not in table1Path: table1Path+=".csv"
    if ".csv" not in table2Path: table2Path+=".csv"
    
    try: #get ProteinSequences, rename columns, and get GRAVY indexes for Table 1
        data = pd.read_csv(table1Path)
        df1 = process(data, fastaPath)
    except Exception as e:
        print("Failed processing Table 1.", e)
        return 0
        
    try: #get GRAVYdifference2 using Table1 minus Table2's peptide sequences 
        data2 = pd.read_csv(table2Path)
        df2 = gravyDiff2(df1, data2) 
    except Exception as e:
        print("Failed calculating GRAVYdifference2.", e)
        df2 = df1
        
    try: #get the middle positions of Table1's peptide sequences (to determine any repeats) and proceed to getting GRAVYdifference = full sequence - peptide sequences
        df3 = gravyDiff(pepPositions(df2)) 
    except Exception as e:
        print("Failed calculating GRAVYdifference.", e)
        df3 = df2
    
    if save:
        df3.to_csv(table1Path[:-4]+"_GRAVY.csv")
        print("GRAVY calculations successful. New data table saved to", os.path.dirname(table1Path))
        
    return df3
    
if __name__ == "__main__":
    fastaPath = inputCheck(input("Please enter the full path of the FASTA file:"))
    table1Path = inputCheck(input("Please enter the full path of Table 1: "))
    table2Path = inputCheck(input("Please enter the full path of Table 2: "))
    
    table1diff = getGRAVYdiffs(fastaPath, table1Path, table2Path)   
    
    #append results to table2:
    data2 = pd.read_csv(table2Path)
    if "ProteinID" not in data2.columns: data2 = data2.rename(columns={"PG.UniProtIds": "ProteinID", "Stripped.Sequence": "PeptideSequence"})
    #handling isomers:
    table1diff["UPID"] = [i.split(';')[0] for i in table1diff["ProteinID"].values]
    data2["UPID"] = [i.split(';')[0] for i in data2["ProteinID"].values]

    subset = table1diff.drop_duplicates("UPID")[["UPID", "ProteinSequence", "SequenceGRAVY", "GRAVYdifference", "GRAVYdifference2"]]
    table2diff = pd.merge(data2,subset, on=["UPID"], how='left')
    table2diff.drop('UPID', axis=1).to_csv(table2Path[:-4]+'_GRAVY.csv', index=False)
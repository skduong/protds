'''
Get distances between a set of PeptideSequences and their overall center of mass
Each PeptideSequence is represented by its middle letter position (closest to the left)
'''
from protds_17 import *

def preProcess(data):
    '''
    The starting dataset has a column of "ProteinID" and "PeptideSequence"
    If it does not provide "ProteinSequence", then one will be made
    If a PeptideSequence appears on multiple locations, they will be looked at individually (one location per row)
    Output: filtered and sorted dataset of ProteinIDs that have results on ProteinDataBank with peptide start, end, and middle positions
    '''
    #search for the ProteinSequence and structures (takes a while)
    ids = [i.replace(';','-').split('-')[0] for i in pd.unique(data['ProteinID'])]
    noResults = [i for i in ids if searchPDB(i)==False]
    df = data[~data['ProteinID'].isin(noResults)].copy()
    if 'ProteinSequence' not in data.columns:
        df['ProteinSequence'] = [proteins[i.replace(';','-').split('-')[0]].record.sequence for i in df['ProteinID'].values]
    #start, end, and middle positions
    seqs = df[['PeptideSequence', 'ProteinSequence']].values
    start=[]; end=[]
    for i in seqs: 
        begin = seq.find_subsequence(seq.ProteinSequence(i[1].replace('U','X')), seq.ProteinSequence(i[0].replace('U','X')))
        if begin.size>0:
            stop = begin+len(i[0])
            start.append(begin+1)
            end.append(stop)
        else:
            begin = seq.find_subsequence(seq.ProteinSequence(i[1].replace('U','X')), 
                                         seq.ProteinSequence(i[0][:int(len(i[0])/2)].replace('U','X')))
            stop = begin+len(i[0])
            start.append(begin+1)
            end.append(stop)
    df['PepStart'] = start
    df['PepEnd'] = end
    df = df.explode(['PepStart', 'PepEnd'])
    #dfDiscard = data[data['ProteinID'].isin(noResults)] concat with df[df['PepStart'].isna()]
    df = df[df['PepStart'].notna()]
    df[['PepStart','PepEnd']] = df[['PepStart','PepEnd']].apply(pd.to_numeric, errors='coerce', downcast='integer', axis=1)
    entry = df[['PepStart', 'PepEnd']].values
    df['PepMid'] = [i[0]+int((i[1]-i[0])/2) for i in entry]
    return df.sort_values(by='ProteinID')
    
def pepDistances(sortedData):
    protGroup = df.groupby('ProteinID')
    dists = []
    for p in protGroup:
        entries = p[1][['ProteinID', 'PepMid']].values
        protid= entries[0][0].replace(';', '-').split('-')[0]
        pdb = bestPDB(protid)
        aligned = alignLoc([[[i[1]] for i in entries]], protid, pdb)
        if sum([len(i) for i in aligned[0][0]])>0:
            structure=pdb.structure 
            com = [struc.mass_center(structure[(structure.chain_id==pdb.bestChain[0]) & 
                                    (structure.res_id==i[0])]) if len(i)>0 else 'Missing' for i in aligned[0][0]]
            fullset = [structure[(structure.chain_id==pdb.bestChain[0]) & (structure.res_id==i[0])] for i in aligned[0][0] if len(i)>0]
            master = struc.mass_center(struc.array([i for sub in fullset for i in sub]))
            dists += [struc.distance(k, master) if k!='Missing' else k for k in com]
        else: 
            dists+=['NA' for i in range(len(p[1]))]
    sortedData['DistanceToCenter'] = dists
    return sortedData
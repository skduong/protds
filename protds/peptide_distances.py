from protds.protds_v3 import *
from Bio import SeqIO

def calculate_gravy(sequence):
    hydropathy_index = {'I':4.5, 'V':4.2, 'L':3.8,'F':2.8, 'C':2.5,'M': 1.9, 'A': 1.8, 'G': -.4, 'T': -.7, 'W': -.9, 'S': -.8,
                        'Y':-1.3, 'P':-1.6, 'H':-3.2, 'E':-3.5, 'Q':-3.5, 'D': -3.5, 'N':-3.5, 'K': -3.9, 'R': -4.5}
    try:
        hydropathy = sum(hydropathy_index[letter] for letter in sequence if letter in hydropathy_index.keys())
        gravy_score = hydropathy/len(sequence) #GRAVY score is sum of hydropathy divided by length of the sequence
        return gravy_score
    except:
        return None

def preProcess(data, fasta_file_path=''):
    '''
    The starting dataset has a column of "ProteinID" and "PeptideSequence"
    If it does not provide "ProteinSequence", then one will be made
    If a PeptideSequence appears on multiple locations, they will be looked at individually (one location per row)
    Output: filtered and sorted dataset of ProteinIDs that have results on ProteinDataBank with peptide middle positions
            along with GRAVY index of full sequence and GRAVY of peptide subsequence
    '''
    #search for the ProteinSequence and structures (takes a while)
    noResults = [i for i in pd.unique(data['ProteinID']) if searchPDB(i.replace(';','-').split('-')[0])==False]
    df = data[~data['ProteinID'].isin(noResults)].copy()
    if 'ProteinSequence' not in data.columns: #add column of "ProteinSequence" if not available
        try: #will look for FASTA file to pull sequences from first, if it doesn't exists, then get them from UniProt records instead
            with open(fasta_file_path) as fasta_file: 
                record_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
            record_dict = dict(zip([i.split('|')[1] for i in record_dict.keys()], [str(j.seq) for j in record_dict.values()]))
            df['ProteinSequence'] = [record_dict[i.split(';')[0]] if i.split(';')[0] in record_dict else None for i in df['ProteinID'].values]
        except Exception as e:
            df['ProteinSequence'] = [proteins[i.replace(';','-').split('-')[0]].record.sequence for i in df['ProteinID'].values]
            print("Unable to use FASTA. Sequences will be taken from UniProt records. ", e)
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
    #GRAVY 
    df['SequenceGRAVY'] = [calculate_gravy(seq) for seq in df['ProteinSequence'].values]
    df['PeptideGRAVY'] = [calculate_gravy(seq) for seq in df['PeptideSequence'].values]
    return df.drop(['PepStart', 'PepEnd'], axis=1).sort_values(by='ProteinID')
    
def pepDistances(sortedData):
    '''
    Get distances between a set of PeptideSequences and their overall center of mass
    Each PeptideSequence is represented by its middle letter position (closest to the left)
    '''
    protGroup = sortedData.groupby('ProteinID')
    dists = []; means = []; stds = []; gravDiff = []
    for p in protGroup:
        protid= p[0].replace(';', '-').split('-')[0]
        pdb = bestPDB(protid)
        aligned = alignLoc([[[i] for i in p[1]['PepMid'].values]], protid, pdb)
        if sum([len(i) for i in aligned[0][0]])>0:
            structure=pdb.structure 
            com = [struc.mass_center(structure[(structure.chain_id==pdb.bestChain[0]) & 
                                    (structure.res_id==i[0])]) if len(i)>0 else 'Missing' for i in aligned[0][0]]
            fullset = [structure[(structure.chain_id==pdb.bestChain[0]) & (structure.res_id==i[0])] for i in aligned[0][0] if len(i)>0]
            combined = struc.mass_center(struc.array([i for sub in fullset for i in sub]))
            distances = [i for i in [struc.distance(k, combined) if k!='Missing' else k for k in com] if not isinstance(i, str)]
            dists += [struc.distance(k, combined) if k!='Missing' else k for k in com]
            means += [np.array(distances).mean() if i!='Missing' else i for i in com]
            stds += [np.array(distances).std() if i!='Missing' else i for i in com]
        else: 
            dists+= ['NA' for i in range(len(p[1]))]
            means += ['NA' for i in range(len(p[1]))]
            stds += ['NA' for i in range(len(p[1]))]
        #GRAVY of proteinSequence-peptideSubsequences
        gravDiff += [p[1]['SequenceGRAVY'].values[0] - sum(p[1]['PeptideGRAVY']) for i in range(len(p[1]))]
    sortedData['GRAVYdifference'] = gravDiff
    sortedData['DistanceToCenter'] = dists
    sortedData['MeanDistances'] = means
    sortedData['StdDistances'] = stds
    return sortedData
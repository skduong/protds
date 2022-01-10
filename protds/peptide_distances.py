def alignPep(proteinGroup, getCenters = True): #getCenters: return only the center of mass for each Peptide
    #getting best aligned PDB structure and its sequence
    protid = proteinGroup[0].split(';')[0]
    pdb = bestPDB(protid)
    seq = chainSeq(checkChains(pdb, protid)[0], pdb.structure)
    
    #aligning "PepMid" positions to the structure
    protLen = len(proteins[protid].getSequence())
    if len(seq[0]) in range(int(protLen-protLen*.2), int(protLen+protLen*.2)):
        ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix())[0]
    else:
        ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix(), local=True)[0]
    mainseq = [i[0] for i in ali.trace]
    alignedseq = [seq[1][ali.trace[mainseq.index(i)][1]] if i in mainseq else -1 for i in proteinGroup[1]["PepMid"].values]
    alignedseq = map(lambda x: seq[1][x] if x!=-1 else x, [ali.trace[mainseq.index(i)][1] if i in mainseq else -1 for i in proteinGroup[1]["PepMid"].values])
    structure = pdb.structure
    coords = [structure[(structure.chain_id==pdb.bestChain[0]) & (structure.res_id==i)] for i in alignedseq]
    
    if getCenters:
        return [struc.mass_center(i) if len(i)!=0 else i for i in coords]
    else:
        return coords

def pepCenterDist(sortedData, getNoResults=False):
    '''
    Get distances between a set of PeptideSequences and their overall center of mass
    Each PeptideSequence is represented by its middle letter position (closest to the left)
    '''
    #search for Protein structures (takes a while)
    uniqueIDs = sortedData.drop_duplicates(["ProteinID"])[["ProteinID", "ProteinSequence"]].values
    noResults = [i[0] for i in uniqueIDs if searchPDB(i[0].split(';')[0], False, i[1])==False]
    sortedData = sortedData[~sortedData['ProteinID'].isin(noResults)]
    protGroup = sortedData.groupby('ProteinID')
    dists = []; means = []; stds = []
    for p in protGroup:
        pepCoords = alignPep(p, False)
        if sum([len(i) for i in pepCoords]) > 0:
            combined = struc.mass_center(struc.array([i for sub in pepCoords for i in sub]))
            pepCom = [struc.mass_center(i) if len(i)!=0 else i for i in pepCoords]
            distances = [struc.distance(i, combined) for i in pepCom]
            dists += [i if isinstance(i,np.float32) else 'Missing' for i in distances]
            complete = np.array([i for i in distances if isinstance(i,np.float32)])
            means += [complete.mean() if len(i)!=0 else 'Missing' for i in pepCoords]
            stds += [complete.std() if len(i)!=0 else 'Missing' for i in pepCoords]
        else: 
            dists+= ['NA' for i in range(len(p[1]))]
            means += ['NA' for i in range(len(p[1]))]
            stds += ['NA' for i in range(len(p[1]))]
    sortedData['DistanceToCenter'] = dists
    sortedData['MeanDistances'] = means
    sortedData['StdDistances'] = stds
   
    if getNoResults: return sortedData, noResults
    else: return sortedData
        
def getDistances():
    return 0 #working on it
    
if __name__ == "__main__":
    from protds_v3 import *
    from ..scripts.gravy_diff.py import *
    getDistances()
    
else: #imported by a Jupyter Notebook
    from protds.protds_v3 import *
    from scripts.gravy_diff import *
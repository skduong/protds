def pepDistances(sortedData):
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
        protid= p[0].split(';')[0]
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
    sortedData['DistanceToCenter'] = dists
    sortedData['MeanDistances'] = means
    sortedData['StdDistances'] = stds
    return sortedData
        
def getDistances():
    return 0
    
if __name__ == "__main__":
    from protds_v3 import *
    from ..scripts.gravy_diff.py import *
    getDistances()
    
else: #imported by a Jupyter Notebook
    from protds.protds_v3 import *
    from scripts.gravy_diff import *
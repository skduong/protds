def alignPep(proteinGroup, getCenters = True): #getCenters: return only the center of mass for each Peptide
    protid = proteinGroup[0].split(';')[0]
    if proteins[protid].structures: #do alignment on PDB structures if they exist
        pdb = bestPDB(protid)
        seq = chainSeq(checkChains(pdb, protid)[0], pdb.structure)
        #aligning "PepMid" positions to the structure
        protLen = len(proteins[protid].getSequence())
        if len(seq[0]) in range(int(protLen-protLen*.2), int(protLen+protLen*.2)):
            ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix())[0]
        else:
            ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix(), local=True)[0]
            
        mainseq = [i[0] for i in ali.trace]
        try:
            alignedseq = map(lambda x: seq[1][x] if x!=-1 else -111, [ali.trace[mainseq.index(i-1)][1] if i in mainseq else -1 for i in proteinGroup[1]["PepMid"].values])
        except ValueError:
            ali = align.align_optimal(proteins[protid].getSequence(), seq[0], align.SubstitutionMatrix.std_protein_matrix())[0]
            mainseq = [i[0] for i in ali.trace]
            alignedseq = map(lambda x: seq[1][x] if x!=-1 else -111, [ali.trace[mainseq.index(i-1)][1] if i in mainseq else -1 for i in proteinGroup[1]["PepMid"].values]) 
        structure = pdb.structure
        coords = [structure[(structure.chain_id==pdb.bestChain[0]) & (structure.res_id==i)] for i in alignedseq]
        structName = pdb.PDBid
        if any([i.shape[0]==0 for i in coords]): #missing coordinates
            pred = proteins[protid].getPredictedStrucs()
            if pred: 
                coords = [pred[pred.res_id==i] for i in proteinGroup[1]["PepMid"].values]
                structName = "AlphaFold"
    else: #look for AlphaFold predictions if PDB structures aren't available
        pred = proteins[protid].getPredictedStrucs()
        coords = [pred[pred.res_id==i] for i in proteinGroup[1]["PepMid"].values] if pred else [[]]
        structName = "AlphaFold"
    
    if getCenters:
        return [struc.mass_center(i) if len(i)!=0 else i for i in coords], structName
    else:
        return coords, structName

def pepCenterDist(sortedData, getNoResults=False):
    '''
    Get distances between a set of PeptideSequences and their overall center of mass
    Each PeptideSequence is represented by the coordinates of its middle letter position (closest to the left)
    '''
    #search for Protein structures (takes a while)
    uniqueIDs = sortedData.drop_duplicates(["ProteinID"])[["ProteinID", "ProteinSequence"]].values
    noResults = [i[0] for i in uniqueIDs if searchPDB(i[0].split(';')[0], False, i[1])==False]
    sortedData = sortedData[~sortedData['ProteinID'].isin(noResults)]
    protGroup = sortedData.groupby('ProteinID')
    dists = []; means = []; stds = []; strucs=[]
    for p in protGroup:
        pepCoords, pdb = alignPep(p, False)
        strucs += [pdb]*len(p[1])
        if sum([len(i) for i in pepCoords]) > 0:
            combined = struc.mass_center(struc.array([i for sub in pepCoords for i in sub]))
            pepCom = [struc.mass_center(i) if len(i)!=0 else i for i in pepCoords]
            distances = [struc.distance(i, combined) for i in pepCom]
            dists += [i if isinstance(i,np.float32) else 'Missing' for i in distances]
            complete = np.array([i for i in distances if isinstance(i,np.float32)])
            means += [complete.mean() if len(i)!=0 else 'Missing' for i in pepCoords]
            stds += [complete.std() if len(i)!=0 else 'Missing' for i in pepCoords]
        else: 
            dists+= ['NA']*len(p[1])
            means += ['NA']*len(p[1])
            stds += ['NA']*len(p[1])
    sortedData['DistanceToCenter'] = dists
    sortedData['MeanDistances'] = means
    sortedData['StdDistances'] = stds
    sortedData['Structure'] = strucs
   
    if getNoResults: return sortedData, noResults
    else: return sortedData
    
def perpDist(params, xyz, abs=True):
    a, b, c, d = params
    x, y, z = xyz
    if abs:
        return np.abs(a * x + b * y + c * z + d) / np.sqrt(a**2 + b**2 + c**2)
    else:
        return (a * x + b * y + c * z + d) / np.sqrt(a**2 + b**2 + c**2)

def angle(n1, n2, complement=False): #input: normal of plane1 (a1,b1,c1) and normal of plane2 (a2,b2,c2)
    dot = np.abs(np.dot(n1,n2)) if complement else np.dot(n1,n2)
    angle = np.arccos(dot / (np.linalg.norm(n1) * np.linalg.norm(n2)))
    return angle*180/np.pi #return angle in degrees
    
def planeSVD(X, p=[]): #X is a set of (X,Y,Z) points; p is a point the plane passes through
    """
    Singular value decomposition method.
    Source: https://gist.github.com/lambdalisue/7201028
    """
    # Find the average of points (centroid) along the columns
    C = np.average(X, axis=0)
    # Create CX vector (centroid to point) matrix
    CX = X - C
    # Singular value decomposition
    U, S, V = np.linalg.svd(CX)
    # The last row of V matrix indicate the eigenvectors of smallest eigenvalues (singular values)
    N = V[-1]

    # Extract a, b, c, d coefficients.
    x0, y0, z0 = p if len(p) !=0 else C
    a, b, c = N
    d = -(a * x0 + b * y0 + c * z0) 
    
    return a, b, c, d

def pepPlane(proteinGroup, customPoint=False, returnP=False):
    xyz, pdb = alignPep(proteinGroup)
    points = [i for i in xyz if len(i)>0 and not np.isnan(sum(i))] 
    if len(points) == 0:
        return None, xyz, pdb
    if len(points) >3 and customPoint: #(if there are <3 points, the optimal plane goes through every point)
        try:
            pindex = np.nanargmax(np.delete(proteinGroup[1]['1st 15min'].tolist(), [i[0] for i in enumerate(xyz) if len(i[1])==0]))
            p = points[pindex]
            points = [points[i] for i in range(len(points)) if i!=pindex]
        except:
            print(proteinGroup[0],"had no '1st 15min' readings. Failed to find point of highest intensity.")
            p=[]
    else: p = []
    
    #get optimized plane with SingularValueDecomposition
    x,y,z = [np.array(i) for i in zip(*points)]
    if returnP: 
        return planeSVD(np.array([x, y, z]).T, p), xyz, (points,p), pdb
    else:
        return planeSVD(np.array([x, y, z]).T, p), xyz, pdb
 
def pepPlaneDist(sortedData, customPoint=False):
    '''
    Get distances between a set of PeptideSequences and a Plane fitted through them
    Each PeptideSequence is represented by the coordinates of its middle letter position (closest to the left)
    The Plane is optimized through SingularValueDecomposition, which minimizes the perpendicular distances between it and the Peptides
    If customPoint=True, then the Plane will pass through the PeptideSequence with the highest 1st 15min intensity
    '''
    #searching for Protein results/coordinates
    uniqueIDs = sortedData.drop_duplicates(["ProteinID"])[["ProteinID", "ProteinSequence"]].values
    noResults = [i[0] for i in uniqueIDs if searchPDB(i[0].split(';')[0], False, i[1])==False]
    sortedData = sortedData[~sortedData['ProteinID'].isin(noResults)]
    protGroup = sortedData.groupby('ProteinID')

    planeDists=[]; strucs=[]
    for g in protGroup:
        plane = pepPlane(g, customPoint)
        if not plane[0]:
            planeDists += ['NA' for i in plane[1]]
            strucs += [plane[-1] for i in plane[1]]
            continue
        else:
            allx, ally, allz = [np.array(i) for i in zip(*[i if len(i)>0 else np.array([np.nan, np.nan, np.nan]) for i in plane[1]])]
            planeDists += [i if not np.isnan(i) else "Missing" for i in perpDist(plane[0], (allx,ally,allz))]
            strucs += [plane[-1] for i in perpDist(plane[0], (allx,ally,allz))]
   
    sortedData["DistanceToPlane"] = planeDists
    sortedData["Structure"] = strucs
    return sortedData 

def plotPlane(proteinGroup, customPoint=False):
    from mpl_toolkits.mplot3d import Axes3D
    import matplotlib.pyplot as plt
    #from ipywidgets import *
    #%matplotlib notebook
    if customPoint:
        plane = pepPlane(proteinGroup, True, True)
        x,y,z = [np.array(i) for i in zip(*[i for i in plane[2][0] if len(i)>0])]
        p = plane[2][1]
    else:
        plane = pepPlane(proteinGroup)
        x,y,z = [np.array(i) for i in zip(*[i for i in plane[1] if len(i)>0])]
        p = []
    
    plt.figure()
    ax = plt.subplot(111, projection='3d')
    ax.scatter(x, y, z, color='b')
    if len(p)>0: ax.scatter(p[0], p[1], p[2], color='r')

    # plot plane
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    X,Y = np.meshgrid(np.arange(xlim[0], xlim[1]),
                      np.arange(ylim[0], ylim[1]))
                      
    normal = plane[0][:-1]
    Z = (-normal[0] * X - normal[1] * Y - plane[0][-1]) * 1. /normal[2]

    ax.plot_wireframe(X,Y,Z, color='k')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.show()
        
def getCenterDistByTime(filename, df, times):
    #All Times:
    alltimes = pepCenterDist(peptidePercentage(gravyDiff(df)), True) 
    alltimes[0].drop("PepMid", axis=1).to_csv(os.path.join("output","AllTimes"+filename))
    summary = pd.DataFrame(pd.unique(alltimes[0]["ProteinID"]), columns=["ProteinID"])

    #Separate Times:
    df = df[~df['ProteinID'].isin(alltimes[1])]
    for i in times:
        subset = pepCenterDist(peptidePercentage(gravyDiff(df[df[i].notna()])))
        subset.drop([t for t in times if t!=i]+["PepMid"], axis=1).to_csv(os.path.join("output",i+filename))
        subset.drop_duplicates(subset=['ProteinID'])
        #summary table:
        for col in ['GRAVYdifference', 'PeptidePercentage', 'MeanDistances', 'StdDistances']:
            summary[i+col] = summary.ProteinID.map(dict(zip(subset["ProteinID"], subset[col].values)))
    for col in ['GRAVYdifference', 'PeptidePercentage', 'MeanDistances', 'StdDistances']:
        summary["AllTimes"+col] = summary.ProteinID.map(dict(zip(alltimes[0]["ProteinID"], alltimes[0][col].values)))
    summary.to_csv(os.path.join("output","Summary"+filename), index=False)
    print("Center of Mass distances for all times have been successfully calculated.")
        
def getAngles(data1, data2):
    '''
    data1 corresponds to the dataset used for Plane1, which crossed through the point of highest 1st15min intensity
    data2 corresponds to the dataset used for Plane2, which was not required to pass through any particular point
    '''
    #handling isomers:
    if "UPID" not in data1.columns or "UPID" not in data2.columns:
        data1["UPID"] = [i.split(';')[0] for i in data1["ProteinID"].values]
        data2["UPID"] = [i.split(';')[0] for i in data2["ProteinID"].values]
    both = data1[data1["UPID"].isin(data2["UPID"])]
    
    angles = {}
    for protein in pd.unique(both["UPID"]):
        if searchPDB(protein, False, data1[data1["UPID"]==protein].ProteinSequence.values[0]) != False:
            p1 = pepPlane((protein, data1[data1["UPID"]==protein]), True)[0]
            p2 = pepPlane((protein, data2[data2["UPID"]==protein]), False)[0]
            try:
                if sum(np.cross(p1[:3], p2[:3])==[0,0,0]) == 3: #parallel planes
                    angles[protein] = "NAN"
                else: 
                    angles[protein] = angle(p1[:3], p2[:3])
            except:
                if not p1 and not p2: angles[protein] = "Planes1&2 NA"
                elif not p1: angles[protein] = "Plane1 NA"
                elif not p2: angles[protein] = "Plane2 NA"
                else: angles[protein] = "NA"    
        else:
            angles[protein] = 'No PDB result'
            
    data1['Angle'] = data1.UPID.map(angles)
    return data1.drop("UPID", axis=1)
        
def getDistances():
    return 0 #working on it
    
if __name__ == "__main__":
    from protds import *
    from ..scripts.gravy_diff.py import *
    getDistances()
    
else: #imported by a Jupyter Notebook
    from protds.protds import *
    from scripts.gravy_diff import *
import numpy as np
import pandas as pd
import urllib, requests, io, gzip, os
from Bio import SwissProt
#biotite: tool for working with structure coordinates and sequences
import biotite.sequence as seq
import biotite.structure as struc
import biotite.database.rcsb as rcsb
import biotite.sequence.align as align
import biotite.structure.io.mmtf as mmtf
import biotite.structure.io.pdbx as pdbx

ALPHAFOLD_DOWNLOADS = "path/to/downloaded/alphafolds"

def alphafold_download_folder(path):
    global ALPHAFOLD_DOWNLOADS
    ALPHAFOLD_DOWNLOADS = path

def get_alphafold(upid, struc_format='cif'):
    #get alphafold structures from downloaded files
    try:
        fname = np.array(os.listdir(ALPHAFOLD_DOWNLOADS))[[(upid in i) and (struc_format in i) for i in os.listdir(ALPHAFOLD_DOWNLOADS)]][0]
        path = os.path.join(ALPHAFOLD_DOWNLOADS, fname)
    except IndexError:
        print(upid, "not found in", ALPHAFOLD_DOWNLOADS)
        return False
    
    with gzip.open(path, mode="rt") as f:
        file_content = f.read()
        file = io.StringIO(file_content)
        cif_file = pdbx.PDBxFile.read(file)
        structure = pdbx.get_structure(cif_file)
    return structure 

def get_swissprot(pid):
    file_link = "https://rest.uniprot.org/uniprotkb/"+pid+".txt"
    with urllib.request.urlopen(file_link) as response:
        swissprot = SwissProt.read(response)
    return swissprot

def get_distance_to_features(upid, mod_num):
    features = get_swissprot(upid).features
    structure = get_alphafold(upid)
    
    if not structure: 
        return ['no alphafold'],[0]
    
    modlocs_pts = structure[0][structure.res_id == mod_num]
    feature_types=[]; feature_locs=[]; feature_dists=[] 
    
    for f in features:
        featlocs = [i for i in range(f.location.start, f.location.end+1)]
        feature_pts = structure[0][[i in featlocs for i in structure.res_id]]
        feature_types.append(f.type)
        feature_locs.append(str(f.location).replace(":", "-")[1:-1])
        feature_dists.append(struc.distance(struc.mass_center(feature_pts), struc.mass_center(modlocs_pts)))
                             
    return feature_types, feature_locs, feature_dists

# def get_distance_to_features(upid, mod_num):
#     features = get_swissprot(upid).features
#     structure = get_alphafold(upid)
    
#     if not structure: 
#         return ['no alphafold'],[0]
    
#     modlocs_pts = structure[0][structure.res_id == mod_num]
#     feature_types=[]; feature_dists=[]
    
#     for f in features:
#         featlocs = [i for i in range(f.location.start, f.location.end+1)]
#         feature_pts = structure[0][[i in featlocs for i in structure.res_id]]
#         feature_types.append(f.type)
#         feature_dists.append(struc.distance(struc.mass_center(feature_pts), struc.mass_center(modlocs_pts)))
                             
#     return feature_types, feature_dists

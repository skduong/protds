import numpy as np
import pandas as pd
import icn3dpy
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

def get_alphafold(pid, struc_format='cif'):
    #get alphafold structures from downloaded files
    try:
        fname = np.array(os.listdir(ALPHAFOLD_DOWNLOADS))[[(pid in i) and (struc_format in i) for i in os.listdir(ALPHAFOLD_DOWNLOADS)]][0]
        path = os.path.join(ALPHAFOLD_DOWNLOADS, fname)
    except IndexError:
        print(pid, "not found in", ALPHAFOLD_DOWNLOADS)
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
    upid = upid.split('-')[0]
    features = get_swissprot(upid).features
    structure = get_alphafold(upid)
    
    if not structure: 
        return ['no alphafold'],['no alphafold'],['no alphafold'],[-1]
    
    modlocs_pts = structure[0][structure.res_id == mod_num]
    feature_types=[]; feature_locs=[]; feature_dists=[]; feature_notes=[] 
    
    for f in features:
        feature_types.append(f.type)
        
        if f.location.start+1 == f.location.end:
            feature_locs.append(str(f.location.end))
        else:
            feature_locs.append(str(f.location.start+1) +'-'+ str(f.location.end))
            
        try:
            description = list(f.qualifiers.items())[0]
        except IndexError:
            description = ('note', '')
        feature_notes.append(description[0]+' = '+description[1])
        
        try:
            featlocs = range(f.location.start+1, f.location.end+1)
            feature_pts = structure[0][[i in featlocs for i in structure.res_id]]
            feature_dists.append(struc.distance(struc.mass_center(feature_pts), struc.mass_center(modlocs_pts)))
        except TypeError: #UnknownPosition()
            feature_dists.append(-1)
         
    return feature_types, feature_notes, feature_locs, feature_dists

def highlight_alphafold(prot_id, loc1, loc2, col1='FF0', col2='F00'):
    highlight1 = "; select " + " or ".join([".A:"+l for l in loc1]) + "; color "+col1
    highlight2 = "; select " + " or ".join([".A:"+l for l in loc2]) + "; color "+col2
    settings = "; select .A; toggle highlight; view annotations; set view detailed view; set background white;"
    
    return icn3dpy.view(command="load af "+prot_id+ highlight1 + highlight2 + settings)
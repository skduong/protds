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

def get_alphafold(pid, alphafold_path, struc_format='cif'): #get alphafold structures from downloaded files
    try:
        fname = np.array(os.listdir(alphafold_path))[[(pid in i) and (struc_format in i) for i in os.listdir(alphafold_path)]][0]
        path = os.path.join(alphafold_path, fname)
    except IndexError:
        print(pid, "not found in", alphafold_path)
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

def highlight_alphafold(prot_id, loc1, loc2, col1='FF0', col2='F00'):
    highlight1 = "; select " + " or ".join([".A:"+l for l in loc1]) + "; color "+col1
    highlight2 = "; select " + " or ".join([".A:"+l for l in loc2]) + "; color "+col2
    settings = "; select .A; toggle highlight; view annotations; set view detailed view; set background white;"
    
    return icn3dpy.view(command="load af "+prot_id+ highlight1 + highlight2 + settings)


def uniprot_record(pid, record_dict):
    if pid not in record_dict.keys():
        record_dict[pid] = get_swissprot(pid)
    return record_dict[pid]

def center_of_mass_dist(mod_loc, feat_locations, structure):
    mod_coords = structure[0][structure.res_id == mod_loc]
    res_range = range(feat_locations.start+1, feat_locations.end+1)
    feature_coords = structure[0][[i in res_range for i in structure.res_id]]
    return struc.distance(struc.mass_center(feature_coords), struc.mass_center(mod_coords))

def min_max_dists(mod_loc, feat_locations, structure):
    mod_coords = structure[0][structure.res_id == mod_loc]
    res_range = np.arange(feat_locations.start+1, feat_locations.end+1)
    
    dists_to_feature = [struc.distance(struc.mass_center(mod_coords), 
                                       struc.mass_center(structure[0][structure.res_id == res]))
                        for res in res_range]
    
    min_i, max_i = np.argmin(dists_to_feature), np.argmax(dists_to_feature)
    return  res_range[min_i], dists_to_feature[min_i], res_range[max_i], dists_to_feature[max_i]

def get_distance_to_features(mod_num, features, structure):
    if not structure: 
        return ['no alphafold'],['no alphafold'],['no alphafold'],[-1],['no alphafold'],[-1],['no alphafold'],[-1]
    
    feature_types=[]; feature_locs=[]; feature_notes=[]
    feature_dists=[]; min_locs=[]; min_dists=[]; max_locs=[]; max_dists=[]
    
    for f in features:
        feature_types.append(f.type)

        try:
            description = list(f.qualifiers.items())[0]
        except IndexError:
            description = ('note', '')
        feature_notes.append(description[0]+' = '+description[1])
        
        if f.location.start+1 == f.location.end:
            feature_locs.append(str(f.location.end))
            feature_dists.append(center_of_mass_dist(mod_num, f.location, structure))
            min_locs.append(''); min_dists.append(''); max_locs.append(''); max_dists.append('')
        else:
            feature_locs.append(str(f.location.start+1) +'-'+ str(f.location.end))
            try:
                if mod_num in range(f.location.start+1, f.location.end+1):
                    feature_dists.append(0)
                    min_locs.append(''); min_dists.append(''); max_locs.append(''); max_dists.append('')
                else:
                    feature_dists.append(center_of_mass_dist(mod_num, f.location, structure))
                    min_max = min_max_dists(mod_num, f.location, structure)
                    min_locs.append(min_max[0]); min_dists.append(min_max[1])
                    max_locs.append(min_max[2]); max_dists.append(min_max[3])
            except TypeError: #UnknownPosition()
                feature_dists.append(-1)
                min_locs.append(''); min_dists.append(''); max_locs.append(''); max_dists.append('')
         
    return feature_types, feature_notes, feature_locs, feature_dists, min_locs, min_dists, max_locs, max_dists

def process_data(df:pd.DataFrame, alphafold_downloads=''):
    uniprot_records = {}
    dist_list=[]; feat_list=[]; feat_loc_list=[]; feat_note_list=[]
    min_loc_list=[]; min_dist_list=[]; max_loc_list=[]; max_dist_list=[]

    for i in df.index:
        entry = df.loc[i]
        upid = entry.ProteinID.split('-')[0]
        features = uniprot_record(upid, uniprot_records).features
        structure = get_alphafold(upid, alphafold_downloads)
        
        feats, feat_notes, feat_locs, dists, min_locs, min_dists, max_locs, max_dists = get_distance_to_features(entry.ModifiedLocationNum, features, structure)
        dist_list.append(dists)
        feat_list.append(feats)
        feat_loc_list.append(feat_locs)
        feat_note_list.append(feat_notes)
        min_loc_list.append(min_locs)
        min_dist_list.append(min_dists)
        max_loc_list.append(max_locs)
        max_dist_list.append(max_dists)

    df["Feature"] = feat_list
    df["Feature_Description"] = feat_note_list
    df["Feature_Location"] = feat_loc_list
    df["Distance_to_ModifiedLocationNum"] = dist_list
    df["Min_Feature_Location"] = min_loc_list
    df["Minimal_Distance"] = min_dist_list
    df["Max_Feature_Location"] = max_loc_list
    df["Maximal_Distance"] = max_dist_list
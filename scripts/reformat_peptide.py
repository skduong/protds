from tkinter import Tk 
from tkinter.filedialog import askopenfilename
import pandas as pd
import os

def process(data):
    if all(label in data.columns for label in ["Assigned Modifications", "Start", "Protein ID", "Gene", "Protein Description"]):
        return process1(data)
    elif all(label in data.columns for label in ["EG.ProteinPTMLocations", "EG.ModifiedSequence", "PEP.PeptidePosition", "PG.UniProtIds", "PG.Genes"]):
        return processDIA(data)
    else:
        return None

def process1(data):
    """
    Please make sure the data has the following column names:
    “Assigned Modifications”, "Start", "Protein ID", "Gene", "Protein Description"
    """
    #drop peptides with no modifications:
    df = data.copy(deep=False) #avoid SettingWithCopyWarning
    df = df[df["Assigned Modifications"].notna()]
    #separate same modifications as different rows:
    df["Assigned Modifications"] = df["Assigned Modifications"].str.replace(' ', '').str.split(',')
    df = df.explode("Assigned Modifications")
    #delete entries with Assigned Modifications of 15.9949, 57.0215, or 42.0106:
    df = df[~df["Assigned Modifications"].str.contains('15.9949|15.9949|42.0106')]
    #“Modification Type”, “Modification Amino Acid”, and “Modification Location”:
    df["Modification Type"] = df["Assigned Modifications"].map(lambda x: x.split('(')[1][:-1])
    def getMods(seq, start, mod):
        if "N-term" in mod:
            return seq[0], start
        else:
            letter = next(filter(str.isalpha, mod))
            return letter, start-1+int(mod[:mod.find(letter)])
    df[["Modification Amino Acid", "Modification Location"]] = list(map(getMods, df["Peptide Sequence"], df["Start"], df["Assigned Modifications"]))
    #aggregate intensities
    cols = [i for i in df.columns if "Intensity" in i and "MaxLFQ" not in i]
    group = df.groupby(["Protein ID", "Gene", "Protein Description", "Modification Type", "Modification Location", "Modification Amino Acid"], as_index=False)
    return group[cols].agg('sum') 
    
def processDIA(data):
    """
    Please make sure the data has the following column names:
    "EG.ProteinPTMLocations", "EG.ModifiedSequence", "PEP.PeptidePosition", "PG.UniProtIds", "PG.Genes"
    """
    intensity_col_keywords = ["PEP.Quantity"] #[All columns containing any intensity_col_keywords will be considered]
    
    #drop peptides with no modifications:
    df = data.copy(deep=False) #avoid SettingWithCopyWarning
    df = df[df["EG.ProteinPTMLocations"].notna()]
    #include N-term modifications, without "Acetyl":
    nterms = (df["EG.ModifiedSequence"].str.contains("N-term")) & ~(df["EG.ModifiedSequence"].str.contains("Acetyl"))
    sub = df[["EG.ProteinPTMLocations", "EG.ModifiedSequence", "PEP.PeptidePosition"]].loc[nterms]
    df["EG.ProteinPTMLocations"].loc[nterms] = list(map(lambda x: x[0][:-1]+ f",{x[1][1]}{x[2]})", sub.values))
    #separate same modifications as different rows:
    df["EG.ProteinPTMLocations"] = df["EG.ProteinPTMLocations"].apply(lambda x: x[1:-1].replace(' ','').split(','))
    df = df.explode("EG.ProteinPTMLocations")
    #split “EG.ProteinPTMLocations” column to “Modification Amino Acid” and “Modification Location”:
    df[["Modification Amino Acid", "Modification Location"]] = list(map(lambda x: (x[0], x[1:]), df["EG.ProteinPTMLocations"]))
    #modification type:
    df["Modification Type"] = list(map(lambda x: x[0][x[0].find(x[1]+'[')+2:].split(']')[0], df[["EG.ModifiedSequence", "Modification Amino Acid"]].values))
    #delete rows containing “Carbamidomethyl” or “Oxidation” in the “Modification type” column:
    df = df[~df["Modification Type"].str.contains("Carbamidomethyl|Oxidation")]
    #aggregate intensities:
    cols = [i for i in df.columns if any(word in i for word in intensity_col_keywords)]
    group = df.groupby(["PG.UniProtIds", "PG.Genes", "Modification Type", "Modification Location", "Modification Amino Acid"], as_index=False)
    return group[cols].agg('sum') 
    
def getOutput(filename):
    filepath, extension = os.path.splitext(filename)
    try:
        if extension.lower() == ".tsv":
            process(pd.read_csv(filename, sep='\t')).to_csv(filepath+"_reformatted.tsv", sep='\t', index=False)
        elif extension.lower() == ".csv":
            process(pd.read_csv(filename)).to_csv(filepath+"_reformatted.csv", index=False)
        else: 
            process(pd.read_excel(filename)).to_excel(filepath+'_reformatted'+extension, index=False)
        print(filepath.split('/')[-1]+extension, "has been successfully reformatted!")
    except AttributeError:
        print("Error: Invalid column names.\n")
    except ValueError:
        print("Error: Unrecognized file type. Try csv, tsv, or excel files.\n")

if __name__ == "__main__":
    Tk().withdraw() 
    print("Choose a file for reformating\n")
    filename = askopenfilename()
    getOutput(filename)
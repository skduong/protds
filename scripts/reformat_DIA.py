from tkinter import Tk 
from tkinter.filedialog import askopenfilename
import pandas as pd
import os

"""
Please make sure the data has the following column names:
"EG.ProteinPTMLocations", "EG.ModifiedSequence", "PEP.PeptidePosition", "PG.UniProtIds", "PG.Genes"
"""
#all columns containing any intensity_col_keywords will be considered:
intensity_col_keywords = ["PEP.Quantity"]

def process(data):
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
    except:
        print("Unrecognized file type. Try csv, tsv, or excel files.")
        return 0

if __name__ == "__main__":
    Tk().withdraw() 
    print("Choose a file for reformating\n")
    filename = askopenfilename()
    getOutput(filename)
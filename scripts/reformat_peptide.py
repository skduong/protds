from tkinter import Tk 
from tkinter.filedialog import askopenfilename
import pandas as pd
import os

def process(data):
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
    
def getOutput(filename):
    filepath, extension = os.path.splitext(filename)
    try:
        if extension.lower() == ".tsv":
            process(pd.read_csv(filename, sep='\t')).to_csv(filepath+"_reformated.tsv", sep='\t', index=False)
        elif extension.lower() == "csv":
            process(pd.read_csv(filename)).to_csv(filepath+"_reformated.csv", index=False)
        else: 
            process(pd.read_excel(filename)).to_excel(filepath+'_reformated'+extension, index=False)
        print(filepath.split('/')[-1]+extension, "has been sucessfully reformated!")
    except:
        print("Unrecognized file type. Try csv, tsv, or excel files.")
        return 0

if __name__ == "__main__":
    Tk().withdraw() 
    print("Choose a file for reformating\n")
    filename = askopenfilename()
    getOutput(filename)
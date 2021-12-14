import sys 
sys.path.append('..')
from protds.peptide_distances import *

def Process_File():
    fasta_file_path = None
    file_directory = None
    
    if fasta_file_path == None: fasta_file_path = input("Please enter path to fasta file without quotations: ")
    if not ".fasta" in fasta_file_path: fasta_file_path += ".fasta"
    if file_directory == None: file_directory = input("Please enter path to file folder without quotations: ")
    
    try:
        project_files = os.listdir(file_directory)
        for file in project_files:
            if file[-4:] == '.csv':
                try:
                    df = preProcess(pd.read_csv(os.path.join(file_directory, file)), fasta_file_path)
                    result = pepDistances(df)
                    result.to_csv(os.path.join(file_directory, file[:-4]+'_distances.csv'))
                    #break
                except Exception as e:
                    print(e)
                    continue
        print('\nFile processed successfully!')
        
    except Exception as e:
        print('\nCould not complete operation.', e)

if __name__ == "__main__":
    Process_File()
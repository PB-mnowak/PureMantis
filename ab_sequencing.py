# Stop / Amber mutation detection
# Connection with database and identification of unique binders
# Report template

from os import getcwd, system
from Bio import BiopythonWarning
from Bio.SeqIO.FastaIO import SimpleFastaParser
from openpyxl import load_workbook
from os import listdir, makedirs
from os.path import isfile, join
from src.ab_data import *
from src.ab_classes import *
import warnings


def get_fasta(path):
    fasta_files = [f for f in listdir(path) if isfile(join(path, f)) and f.lower().endswith("fasta")]
    return fasta_files

def get_scf(path):
    # try:
        path = join(path, "scf_files")
        scf_files = [f for f in listdir(path) if isfile(join(path, f)) and f.lower().endswith("scf")]
        return scf_files
    # except:
    #     print('No scf_files folder found')

def mode_select():
    print("Select analysis mode based on type of sequences:")
    print("1. Nucleotide sequences")
    print("2. Aminoacid sequences")
    while True:
        m = input("Mode number: ")
        if m == "1":
            return "nt"
        elif m == "2":
            return "aa"
        else:
            print(" ---< Mode not recognized. Try again >---")

def scf_file_sort(ub_list, path, file_dir):
    # ub_main_cont = [(rfs, lib), [len(seqs), seqs]
    scf_files = get_scf(path=path)
    n = len(ub_list)
    if scf_files:
        for i in range(n):
            makedirs(f'{file_dir}\\Variant {i+1}', exist_ok=True)
            for seq_id in ub_list[i][1][1]:
                if "_(reversed)" in seq_id:
                    seq_id = seq_id.replace("_(reversed)", "")
                if '.scf' not in seq_id:
                    seq_id = seq_id + '.scf'
                if seq_id in scf_files:
                    system(f'move ".\\scf_files\\{seq_id}" "{file_dir}\\Variant {i+1}\\" >nul')
                    scf_files.remove(seq_id)
        for file in scf_files:
            makedirs(f'{file_dir}\\Not classified', exist_ok=True)
            system(f'move ".\\scf_files\\{file}" "{file_dir}\\Not classified\\" >nul')
    else:
        
        print("No scf files detected")

def print_frame(func):
    def inner(*args, **kwargs):
        system("cls")
        print('/ PureMantis \\'.center(80, '_'))
        func(*args, **kwargs)
        print(''.center(80, '‾'))
    return inner

@print_frame
def main():

    warnings.filterwarnings('ignore', category=BiopythonWarning, module='Bio')

    # m = input("Mode - nt / aa: ")
    # m = "nt"
    # m = mode_select()
    m = "nt"
    mypath = getcwd()
            
    fasta_files = get_fasta(mypath)

    # Main data container of the analysis
    rf_main_cont = []
    ub_main_cont = {}
    binder = Binder()

    # Excel template file path and name
    template_dir = "P:\\_research group folders\\AB Antibodies and Phage Display\\_CDR analysis script\\Templates"
    wb = load_workbook(join(template_dir, "rf_template.xlsx"))
    ws_rf = wb["RF analysis"]
    ws_ub = wb["Unique binders"]

    for fasta in fasta_files:
        with open(join(mypath, fasta), "r") as handle:
            for values in SimpleFastaParser(handle):

                # nt sequence extraction
                seq_id = values[0]
                seq_nt = values[1].replace("-", "")
                binder.load_seq_nt(seq=seq_nt)
                
                # !!! Library recognition -> start/end seq or universal !!!
                
                # No sequence recognized
                if binder.seq:
                    binder.extract_rfs(mode=m)
                    if m == "nt":
                        binder.translate_rfs()
                    rf_extracted = binder.rfs
                    
                    # Library detection and mutation search
                    mutations = []
                    
                    if binder.lib:
                        frames = lib_frame_dict[binder.lib]
                        for i in range(8):
                            if frames[i] not in seq_nt:
                                mutations.append(fr_dict[i])
                            else:
                                mutations.append("")
                    else:
                        binder.lib = "Not recognized"
                        mutations = ["-" for _ in range(8)]

                    seq_aa = binder.seq
                    if m == "nt":
                        seq_aa = seq_aa.translate()
                    
                    # Stop codon search
                    codons = [str(binder.seq[i:i+3]) for i in range(0, len(binder.seq), 3)]
                    amber = codons.count("TAG")
                    ochre = codons.count("TAA")
                    opal = codons.count("TGA")

                    stop_mut = [n if n else "" for n in [amber, ochre, opal]]

                    # Appending data
                    rf_main_cont.append([seq_id, rf_extracted, binder.lib, mutations, stop_mut, seq_aa])
                    
                    if all(rf_extracted):
                        ub_key = ("_".join(map(str, rf_extracted)), binder.lib)
                        seq_ids = ub_main_cont.get(ub_key, False)
                        if seq_ids:
                            ub_main_cont[ub_key].append(seq_id)
                        else:
                            ub_main_cont[ub_key] = [seq_id]
                
                else:
                    # 0 - id, 1[0:6] - rf, 2 - binder.lib, 3[0:8] - mut, 4 - stop, 5 - seq_aa
                    rf_main_cont.append([seq_id, ["" for _ in range(6)], "Not recognized", ["-" for _ in range(8)], [""] * 3, "-"])

    # Sort based on RF1H
    rf_main_cont.sort(key=lambda x: "-".join(map(str, x[1:7])))

    for i, row in enumerate(rf_main_cont, 3):
        
        # Add sequence ID
        ws_rf.cell(column=2, row=i).value = row[0]
        # Add RFs
        for j in range(6):
            ws_rf.cell(column=3+j, row=i).value = str(row[1][j])
        ws_rf.cell(column=11, row=i).value = str(row[2])
        # Add mutations of NRFs
        for j in range(8):
            ws_rf.cell(column=12+j, row=i).value = row[3][j]
        # Add stop codon detection
        for j in range(3):
            ws_rf.cell(column=20+j, row=i).value = row[4][j]
        # Add full AA sequence
        ws_rf.cell(column=23, row=i).value = str(row[5])
        

    # Change dict into list and add number of IDs per RF
    ub_main_cont = [[ub_key, [len(seqs), seqs]] for ub_key, seqs in ub_main_cont.items()]
    
    # Sort descending by number of IDs sharing RF
    ub_main_cont.sort(key=lambda x: x[1][0], reverse=True)

    # Write RFs of given variant and IDs sharing it
    for i, ub in enumerate(ub_main_cont, start=3):
        ws_ub.cell(column=i, row=2).value = ub[0][1]    # Library
        
        rfs_split = ub[0][0].split("_")     # List of RFs
        for j in range(6):
            ws_ub.cell(column=i, row=j+3).value = rfs_split[j]
        ws_ub.cell(column=i, row=9).value = ub[1][0]    # Number of seqs
        for j, seq in enumerate(ub[1][1]):              # Seq IDs
            ws_ub.cell(column=i, row=11+j).value = seq

    # name = input("Save as: ")
    # file_dir = join(mypath, name)
    # makedirs(f'{file_dir}', exist_ok=True)

    # Sort scf files
    # scf_file_sort(ub_main_cont, mypath, file_dir)

    # Save results  
    # filename = join(file_dir, name + ".xlsx")
    # wb.save(filename=filename)
    
    # Move fasta
    # for fasta in fasta_files:
    #     system(f'move "{fasta}" "{file_dir}" >nul')


if __name__ == "__main__":
    main()
    system('pause')
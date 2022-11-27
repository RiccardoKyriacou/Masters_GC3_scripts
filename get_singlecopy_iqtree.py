import os
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import call as unix

'''
Script to replace Gene ID
headers with species names
for single copy orthogroups
'''

parse = argparse.ArgumentParser()

parse.add_argument("-s", "--singlecopy",type=str, help="Path 1 to 1 orthogroup data to build tree from (use copied file)",required=True)
parse.add_argument("-p", "--proteins",type=str, help="Path to primary transcrips proteins file containing species within OG",required=True)#

args = parse.parse_args()

if args.singlecopy[-1] != '/':
    args.singlecopy = args.singlecopy + '/'

#Matching species names to protein sequences 

sp_id_dict = {}
for proteins in glob.glob(args.proteins + "*.fa"):
    sp_name = proteins.split("/")[-1].split("-")[0]
    if sp_name.startswith("Drosophila_melanogaster"): 
        try:
            sp_name = sp_name.split(".")[0]
        except:
            sp_name = sp_name
    with open(proteins, "r") as f:
        for record in SeqIO.parse(f, 'fasta'):            
            ID_pep = record.id

            sp_id_dict[ID_pep] = sp_name #Dictionary of each gene ID and corresponding species 

for orthogroup in glob.glob(args.singlecopy+ "*.fa"):
    with open(orthogroup, "r") as f, open(orthogroup + ".tree", "w") as OutF:
        for record in SeqIO.parse(f, 'fasta'):            
            ID = record.id
            pep_seq = str(record.seq)

            species_header = sp_id_dict[ID]

            OutF.write(">" + species_header + "\n" + pep_seq + "\n")

# Aligning single copy OGs individually 
print("Running mafft...")
for orthogroup in glob.glob(args.singlecopy + "*.tree"):
    OGname = orthogroup.split("/")[-1].split(".")[0]
    unix("mafft " + orthogroup + " > " + args.singlecopy + OGname + ".fa.tree.mafft", shell=True)

#Creating supermatrix
print("Creating supermatrix...")

os.chdir(args.singlecopy)
unix("ls -1 *.mafft >> fasta_list.txt", shell=True)

unix("phykit create_concat -a fasta_list.txt -p phykit_concat", shell=True)

#IQTREE
print("Running IQTREE...")
unix("iqtree -s phykit_concat.fa ", shell=True)



            

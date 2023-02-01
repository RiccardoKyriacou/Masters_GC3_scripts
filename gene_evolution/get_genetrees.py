#!/usr/bin/env python3
import glob
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import call as unix

'''
This script first moves only OGs which contain 
a high GC3 gene (i.e. genes in heatmap) into 
a new directiory. Script will then perform 
trimal alignment before running IQ tree on
each OG file, creating gene trees. 
'''

# Usage: python3 get_genetrees.py --path ~/GC_analysis/GC_Lep/GC_cds_Lep/results/GC3/1to1_orthologues/outlier_GC3-comparison.tsv --out ../data/1to1_highGC3_OGs/ --species_tree ../data/Species_Tree/SpeciesTree_rooted_node_labels.txt

parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to outlier_GC3-comparison.tsv file to pull relevant OGS", required = True)
parse.add_argument("--out",type=str, help="path to output for OGs", required = True)
parse.add_argument("--species_tree",type=str, help="path to orthofinder species tree", required = True)

args = parse.parse_args()

#Parse tsv of all OGs which containg gene at GC3 > 85% 
def parse_tsv(tsv):
    OG_lst = []
    with open(tsv) as f: 
        for line in f:
            if line.startswith("Orthogroup"): #Skip header of tsv file
                continue
            else: 
                OG_name = line.split("\t")[0]
                OG_lst.append(OG_name)
    #Create a set of each OG appearing once in the list
    return set(OG_lst)

#Copy each file specified in OG-lst to new dir 
def move_files(OG_lst, output_dir):
    for orthogroup in OG_lst:
        filename = args.path.split("/")[-1]
        OG_path = args.path.split(filename)[0]
        aligned_OG = OG_path + orthogroup + ".fa.mafft"
        shutil.copy(aligned_OG, output_dir, follow_symlinks=True)

#Run trimal withut backtranslation to trim OGs and prduce timmed amin-acid alignment
def run_trimal(fn):
    for alignment in glob.glob(fn + "*.mafft"):
        unix(f"trimal -in {alignment} -out {alignment}.trimal -gappyout", shell=True)

#Create dictionary of gene_id : sp_name
def get_OG_sp_dict(tsv):
    ID_sp_dict = {}
    with open(tsv) as f: 
        for line in f:
            if line.startswith("Orthogroup"): #Skip header of tsv file
                continue
            else: 
                sp_name = line.split("\t")[1]
                gene_id = line.split("\t")[2]

                ID_sp_dict[gene_id] = sp_name

    return ID_sp_dict

#Replace header such that:        
#>GENE_ID 
#Is re-written to 
#>sp_name
def rename_trimal_files(trimal_files, ID_sp_dict):
    for alignment in glob.glob(trimal_files + "*.trimal"): #Open alignemnt files
        with open(alignment) as f, open(f"{alignment}.renamed", "w") as outf: #Create .renamed folders 
            #Parse trimal folders to get gene ID and seq
            for record in SeqIO.parse(f, 'fasta'): 
                Gene_ID = record.id
                seq = str(record.seq)
                sp_name = ID_sp_dict[Gene_ID] #Match gene_ID to seq

                outf.write(f">{sp_name}\n{seq}\n") #Re-write.rename file

#Run IQtree to get gene trees for each .trimal file
def run_gene_iqtree(trimal_files, species_tree):
    for alignment in glob.glob(trimal_files + "*.renamed"):
        OG_name = alignment.split("/")[-1].split(".fa")[0] #Getting OG name
        print(f"Running gene IQtree for {OG_name}")
        #Running IQtree on command line
        unix(f"iqtree -s {alignment} -te {species_tree} --prefix {OG_name}", shell=True) 

def main():
    #Parse tsv and copy files
    OG_lst = parse_tsv(args.path)
    move_files(OG_lst, args.out)
    print(f"Copied {len(OG_lst)} OGs to {args.out}")  
    #Running trimal
    print(f"Running trimal (gappyout) for {len(OG_lst)} OGs...")
    run_trimal(args.out) 
    #Currently, trimal files have gene ID as headers but species tree has species names as headers 
    #However, in order to run gene tree, need to match gene tree headers to species tree headers
    #Therefore next two fucntions get dict of gene_ID : sp_name and use it to re-write alignemts as .rename
    ID_sp_dict = get_OG_sp_dict(args.path)
    rename_trimal_files(args.out, ID_sp_dict)
    #Once files are appropriatley names, able to runIQ tree for each gene from the command line 
    run_gene_iqtree(args.out, args.species_tree)

if __name__ == "__main__":
    main()
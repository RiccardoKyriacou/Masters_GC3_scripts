#!/usr/bin/env python3
import os
import re
import glob
import shutil
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from subprocess import call as unix

'''
This script:
1) first moves only OGs which contain a high GC3 gene (i.e. genes in heatmap) into a new directiory. 
2) Performs trimal alignment 
3) Runs IQ tree on each OG file, creating gene trees.
4) Finally, script will create binary tree files needed for RER converge
'''

# Usage: python3 get_genetrees.py --path ~/GC_analysis/GC_Lep/GC_cds_Lep/results/GC3/1to1_orthologues/outlier_GC3-comparison.tsv --out ../data/1to1_highGC3_OGs/ --species_tree ../data/Species_Tree/SpeciesTree_rooted_node_labels.txt

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to outlier_GC3-comparison.tsv file to pull relevant OGS", required = True)
parse.add_argument("-o", "--out",type=str, help="path to output for OGs", required = True)
parse.add_argument("-s", "--species_tree",type=str, help="path to orthofinder species tree", required = True)

args = parse.parse_args()

#Parse tsv of all OGs which containg gene at GC3 > 95% (Lep) or 85% (other clades) 
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
        #Running IQtree on command line (-keep-indent used to keep identical sequences)
        unix(f"iqtree -s {alignment} -te {species_tree} --prefix {OG_name} -keep-ident", shell=True) 

#This fucntion simply concatonates all the gene trees into an appropraite RER files
def create_RER_input(gene_tree_files):
    with open("combined_genetrees.txt", "w") as outf:
        for tree in glob.glob(gene_tree_files + "*.treefile"):
            with open(tree) as f:
                OG_name = tree.split("/")[-1].split(".treefile")[0]
                for line in f:
                    tree_data = line.rstrip()

                    outf.write(f"{OG_name}\t{tree_data}\n")

#From here, script is looking to create a binary tree file for each OG 
#These are needed for binary trait analysis in RER converge 
#Each Binary tree file contains a binary tree with 1s for species with high GC3 
#And 0s for species with lower GC3

#First fucntion uses RegEx to replace each tree branch length value with a 0 
def get_binary_tree_skeleton(combined_genetree_path):
    with open(f"{combined_genetree_path}") as f:
        #Get first tree from cmbined gene-tree output
        first_line = f.readline()
        tree = first_line.split("\t")[1]
        #Match to every decimal place and replace with 0 
        new_tree = re.sub("\d+\.\d+", "0", tree)
        #Replace residual numbers 
        new_tree = re.sub("\d", "0", new_tree)

    return str(new_tree) 

#Dictionary of each OG with list of species it is high GC3 >95% in 
def get_highGC3_OG_dict(tsv, OG_lst):
    highGC3_sp_OG_dict = defaultdict(list)
    with open(tsv) as f: 
        for line in f:
            matched_OG = line.split("\t")[0]
            GC3 = line.split("\t")[4]
            sp_name = line.split("\t")[1]
            for OG in OG_lst:
                if (OG == matched_OG) and (float(GC3) > 95):
                    highGC3_sp_OG_dict[OG].append(sp_name) 
    
    return highGC3_sp_OG_dict

#Getting all binary trees
def get_all_binary_trees(tree_skeleton, highGC3_sp_OG_dict):  
    for OG, species in highGC3_sp_OG_dict.items():
        with open(f"{OG}.binarytree", "w") as outf: #Open new binary tree file for each OG
            binary_tree = tree_skeleton #Reset binary tree to skeleton
            for sp in species: #For each high GC3 species
                
                binary_tree = binary_tree.replace(f"{sp}:0", f"{sp}:1")  #Replace 0 with 1 for high GC3 sp

            outf.write(binary_tree)

def main():
    #Parse tsv to get list of relevant OGS
    OG_lst = parse_tsv(args.path)
    move_files(OG_lst, args.out) #copy files to new directory
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
    #Tree is run with constraint to match toplology of the species tree for RER analysis 
    run_gene_iqtree(args.out, args.species_tree)
    #Then, must combine all these gene trees in tsv file for RER analysis 
    #File in format
    #OG_name \t tree_data
    print("Combining trees into file combined_genetrees.txt")
    current_path = f"{os.getcwd()}/" #IQTREE outputs to current dir so must get current path for create_RER_input
    create_RER_input(current_path)
   
    # RER requires "binary trees" for binary trait calculations 
    # Therefore the rest of this scrip will now look to etting all 
    # binary trees for RER analysis
    print(f"Creating binary trees for {len(OG_lst)} OGs with GC3 outliers...")
    #Firstly, get a tree "skeleton" with same topology as gene tree 
    tree_skeleton = get_binary_tree_skeleton(f"combined_genetrees.txt")
    #Make dictionary of each high GC3 OG : list of species in which gene is high GC
    highGC3_OG_dict = get_highGC3_OG_dict(args.path, OG_lst)
    #Main function to get trees
    get_all_binary_trees(tree_skeleton, highGC3_OG_dict)


if __name__ == "__main__":
    main()

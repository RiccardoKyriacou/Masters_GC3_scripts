#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
from Bio import SeqIO
from collections import defaultdict
from subprocess import call as unix
from joblib import Parallel, delayed

"""
Script to translate aligned protein files for high GC3 containing BUSCOs to 
codon aware nucleotide alignments taht can be used in hyphy selection analyisis.

Requires: HyPhy and HyPhy-analysis to be built in the same directiory running script 

Usage: python3 get_hyphy_input.py --tsv ../../BUSCO_Lep/results/outlier_GC3-comparison.tsv --busco ../../BUSCO_Lep/data/fasta/ --out ../data/aligned_protein_sequences/ --species ../data/lep_tree.nwk

"""
#Parse tsv of all busco_ids which containg gene at GC3 > 95% (Lep) or 85% (other clades) 
def parse_tsv(tsv):
    BUSCO_lst = []
    with open(tsv) as f: 
        for line in f:
            if line.startswith("BUSCO_id"): #Skip header of tsv file
                continue
            else: 
                BUSCO_name = line.split("\t")[0]
                BUSCO_lst.append(BUSCO_name)
    #Create a set of each busco_id appearing once in the list
    return set(BUSCO_lst)

#Copy each file specified in busco_id-lst to new dir 
def copy_files(BUSCO_lst, busco_path, output_dir):
    for busco in BUSCO_lst:
        #Copy aligned proteins files
        protein_aligned_busco = f"{busco_path}/proteins/{busco}.faa.mafft"
        shutil.copy(protein_aligned_busco, output_dir, follow_symlinks=True)
        #Copy nucelotide fasta files
        nucleotide_fasta = f"{busco_path}/nucleotides/{busco}.fna"
        shutil.copy(nucleotide_fasta, output_dir, follow_symlinks=True)

#Run pal2nal for each gene in the buco fasta file 
def run_pal2nal(output_dir, BUSCO_lst):
    #Loop through each aligned protein file and run through pal2nal providing correct nucelotide sequence 
    for busco in BUSCO_lst: 
        busco_aligned_pep = f"{output_dir}{busco}.faa.mafft"
        busco_nuc = f"{output_dir}{busco}.fna"
        unix(f"pal2nal.pl {busco_aligned_pep} {busco_nuc} -output fasta > {output_dir}{busco}.fasta", shell=True)
    #Change headers to just sp_name
    for codon_alignment in glob.glob(f"{output_dir}*.fasta"):
        busco_id = codon_alignment.split("/")[-1].split("_")[0]
        with open(codon_alignment) as f, open(f"{output_dir}{busco_id}.aligned_codon", "w") as outf:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description #Header in form LimLuna.OU830598.1:6261747-6269071
                sp_name = header.split(".")[0] 
                seq = str(record.seq)
                
                outf.write(f">{sp_name}\n{seq}\n")


#Getting labelled forground species (high GC3) for RELAX HyPhy analysis 
def get_labeled_tree(tsv, species_tree, output_dir, src_path):
    #Dict of Busco ID: [list of high GC3 species]
    high_GC3_species_dict = defaultdict(list)
    #Parse tsv to get a list of the species which are high GC3 for each BUSCO
    with open(tsv) as f: 
        for line in f:
            if line.startswith("BUSCO_id"): #Skip header of tsv file
                continue
            else: 
                BUSCOid = line.split("\t")[0]
                sp_name = line.split("\t")[1]
                GC3 = float(line.split("\t")[2])
                #Threshold to consider gene "high GC3"
                if GC3 > 85:
                    high_GC3_species_dict[BUSCOid].append(sp_name)

    #Iterate through dictionary and write list to label trees 
    for busco_id, species_list in high_GC3_species_dict.items():
        #Need to make output list.txt file to feed into LabelTree
        with open(f"{output_dir}{busco_id}_species_list.txt", "w") as outf:
            for sp in species_list:
                outf.write(f"{sp}\n")

    #Can now pass each species list to HyPhy LabelTrees 
    #NOTE: hyphy and hyphy-analysis must be built in the same directiory that you are running the script 
    for species_list in glob.glob(f"{output_dir}*.txt"): #Open alignment files 
        busco_id = species_list.split("/")[-1].split("_")[0]
        unix(f"{src_path}/hyphy/hyphy LIBPATH={src_path}/hyphy/res {src_path}/hyphy-analyses/LabelTrees/label-tree.bf --tree {species_tree} --list {species_list} --output {output_dir}{busco_id}_labeled_tree.nwk --internal-nodes None", shell=True)


#Function to run HyPhy Relax for one gene family / BUSCO set 
def run_hyphy_relax(codon_alignment, labelled_tree):
    #Command line argument to run hyphy (may need to conda install hyphy)
    unix(f"hyphy relax --alignment {codon_alignment} --tree {labelled_tree} --test Foreground", shell=True)
      
def parallel_relax_run(output_dir, n):
    job_dict = {} 
    for codon_alignment in glob.glob(f"{output_dir}*.aligned_codon"):
        for labelled_tree in glob.glob(f"{output_dir}*.nwk"):
            #Make dictionary of jobs to run 
            job_dict[codon_alignment] =labelled_tree

    #Iterate though dict to run jobs
    Parallel(n_jobs=n)(delayed(run_hyphy_relax)(alignment, tree) for alignment, tree in job_dict.items())   

def main():
    #Move files into a new dir W
    BUSCO_lst = parse_tsv(args.tsv) #Get list of relevant busco_ids
    copy_files(BUSCO_lst, args.busco, args.out) #Copy files to new directory
    print(f"Copied {len(BUSCO_lst)} BUSCOs to {args.out}") 
    #Run pal2nal 
    run_pal2nal(args.out, BUSCO_lst)
    #Label trees
    current_path = f"{os.getcwd()}" #IQTREE outputs to current dir so must get current path for create_RER_input
    get_labeled_tree(args.tsv, args.species, args.out, current_path)
    #Run HyPhy relax in parallel 
    print("Running HyPhy RELAX...")
    #Get number of runs fo HyPhy (number of unique BUSCO genes of goven cutoff)
    n = int(len(BUSCO_lst))
    parallel_relax_run(args.out, n)

if __name__ == "__main__":
    parse = argparse.ArgumentParser()

    parse.add_argument("-t", "--tsv",type=str, help="path to outlier_GC3-comparison.tsv file", required = True)
    parse.add_argument("-b", "--busco",type=str, help="path to BUSCO fastas", required = True)
    parse.add_argument("-o", "--out",type=str, help="path to output for busco_ids", required = True)
    parse.add_argument("-s", "--species",type=str, help="species tree", required = True)

    args = parse.parse_args()

    if args.busco[-1] != '/':
        args.busco += '/'

    if args.out[-1] != '/':
        args.out += '/'

    main()

#Citation: 
#Wertheim, JO et al. "RELAX: detecting relaxed selection in a phylogenetic framework." Mol. Biol. Evol. 32, 820–832 (2015).clear

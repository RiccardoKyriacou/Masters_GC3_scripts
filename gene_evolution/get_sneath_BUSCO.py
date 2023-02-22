#!/usr/bin/env python3
import os
import glob
import shutil
import argparse
from Bio import SeqIO
from collections import defaultdict
from subprocess import call as unix

"""
This script is stand alone

Method:
1) Copies all busco_ids for the order of intrest and renames protein alignments 
2) Requires .txt file to get species in the form outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3 
   with each speices separated by a tab in the format:
   outgroup_species_1 \t outgroup_species_2 \t outgroup_species_3 \t ingroup1 \t ingroup2 \t ingroup3
3) Caculates sneath value. This is done on an gene by gene basis for all busco_ids. 
3a) First, a fucntion checks if each amino-acid residue in the gene is the same for all three outgroup species
3b) If the residue is the same, representing the ancestral state, then the equivalent aa is found for each ingroup speices
3c) Adjusted Sneath value is calculated by comaping ancestral residues to the residues ofr each ingroup species and dividng 
    by alignemnt length 

Usage: python3 get_sneath_BUSCO.py --trimmed ../results/trimmedGC3_BUSCO.tsv --busco ../data/fasta/proteins/ --species ../results/sneath_distribution_run1/sp_list_run1 --out ../data/sneath_distribution_run1/

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
        aligned_busco = busco_path + busco + ".faa.mafft"
        shutil.copy(aligned_busco, output_dir, follow_symlinks=True)

#Run trimal withut backtranslation to trim busco_ids and prduce timmed amin-acid alignment
def run_trimal(fn):
    for alignment in glob.glob(fn + "*.mafft"):
        unix(f"trimal -in {alignment} -out {alignment}.trimal -gappyout", shell=True)

#Replace header such that:        
#>PamFasc.OU452285
#Is re-written to 
#>PamFasc
def rename_trimal_files(trimal_files):
    for alignment in glob.glob(trimal_files + "*.trimal"): #Open alignemnt files
        with open(alignment) as f, open(f"{alignment}.renamed", "w") as outf: #Create .renamed folders 
            #Parse trimal folders to get gene ID and seq
            for record in SeqIO.parse(f, 'fasta'): 
                header = record.description
                sp_name = header.split(".")[0]
                seq = str(record.seq)

                outf.write(f">{sp_name}\n{seq}\n") #Re-write.rename file

#Calculate sneath distribution...

#This fucntion reads the species.txt file of species to use in analysis
def get_species(species_file):
    #Species in file should be separated by \t
    with open(species_file) as f:
        for line in f:
            outgroup1 = line.split("\t")[0].strip()
            outgroup2 = line.split("\t")[1].strip()
            outgroup3 = line.split("\t")[2].strip()
            ingroup1 = line.split("\t")[3].strip()
            ingroup2 = line.split("\t")[4].strip()
            ingroup3 = line.split("\t")[5].strip()
    
    return outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3

#Function to parse the .renamed alignemtn files, creating new .sneath files containing only the species of intrest 
def parse_alignemts(output_dir, outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3):
    for alignment in glob.glob(f"{output_dir}*.renamed"): #Open alignemnt files
        species_of_intrest_lst = [outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3]
        with open(alignment) as f, open(f"{alignment}.sneath", "w") as outf: #Create 6 species file for sneath analysis 
            #Parse trimal folders to get gene ID and seq
            for record in SeqIO.parse(f, 'fasta'): 
                sp_name = record.description
                seq = str(record.seq)
                #For each species of intrest 
                if sp_name in species_of_intrest_lst:
                    outf.write(f">{sp_name}\n{seq}\n") #Re-write alignemnt to only include species of intrest 

# Main function to comapre outgroup residues and calculate the Sneath index
# Uses following get_ingroup_conserevd_residues() and get_sneath_value() functions defined below
# 1) Fucntion essentially opens each busco_id file and makes a list of outgroup residues and dicitonaty of species : ingroup residues
# 2) Then for each ingroup species residue in the dictionary, it is passed through get_ingroup_conserevd_residues() to get conserved 
#    outgroup residues and corresponding ingroup residues
# 3) Finally, residues for each ingroup species is comapred to outgroup residues to get adjusted Sneath value 
def comapare_sneath_residues(output_dir, outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3):
    with open("sneath_values.tsv", "w") as outf:  
        for sneath_alignment in glob.glob(f"{output_dir}*.sneath"):
            busco_id_name = sneath_alignment.split("/")[-1].split(".fa")[0] #Getting name of busco_id for output
            with open(sneath_alignment) as f:
                #List of sequences for genes in the outgroups
                outgroup_seq_lst = []   
                #Dictionary of species name : sequences for genes in the ingroups 
                ingroup_dict = {}
                for record in SeqIO.parse(f, 'fasta'): 
                    sp_name = record.description            
                    
                    #Appending outgroup genes to a list (species name not needed)
                    if (sp_name == outgroup1) or (sp_name == outgroup2) or (sp_name == outgroup3):
                        #Get the gene sequnce
                        seq = str(record.seq)
                        #Store these sequences in a list        
                        outgroup_seq_lst.append(seq)   
                    
                    #For the ingroup species, append species name and genes to dictionary 
                    if (sp_name == ingroup1) or (sp_name == ingroup2) or (sp_name == ingroup3): 
                        seq = str(record.seq)      
                        ingroup_dict[sp_name] = seq
                
                #Iterate through each ingroup residue for each species
                #This allows us to retrieve the correct species name and corresponding Sneath in the output 
                for ingroup_species, ingroup_residues in ingroup_dict.items():
                    
                    #Use get_ingroup_conserevd_residues() to get conserved outgroup and ingroup residues 
                    conserved_aa_outgroups, aa_ingroup = get_ingroup_conserevd_residues(outgroup_seq_lst, ingroup_residues)
        
                    #Get sneath value for each ingroup comapred to conserved outgroup residues 
                    adjusted_sneath = get_sneath_value(conserved_aa_outgroups, aa_ingroup)
                    
                    outf.write(f"{busco_id_name}\t{ingroup_species}\t{adjusted_sneath}\n")

#This function takes a list of the outgroup sequences and outgroup sequences and retunrs a list of only the "consevred aa"
#Conserved aas are residues that are the same in all three ingroup species
def get_ingroup_conserevd_residues(outgroup_seq_lst, ingroup_gene):          
    #Unpack each gene outgroup individually from outgroup sequence list 
    outgroup_gene1 = [*outgroup_seq_lst[0]]
    outgroup_gene2 = [*outgroup_seq_lst[1]]
    outgroup_gene3 = [*outgroup_seq_lst[2]]
   
    #Unpack aa for the ingroup species 
    ingroup_gene = [*ingroup_gene]

    conserved_aa_outgroups = [] #List of ancestral aa residues 
    corresponding_aa_ingroup = [] #List of the corresponding ingroup residues 

    #Compare each amino-acid in parallel using zip()
    for aa in zip(outgroup_gene1, outgroup_gene2, outgroup_gene3, ingroup_gene):  
        #aa[0] = aa for outgroup 1 and so on...
        if (aa[0] == aa[1]) and (aa[0] == aa[2]):
            #If its the same aa in all three outgroups,  aa represents ancestral state 
            #Make a list of the aa in this ancestral state for the gene 
            conserved_aa_outgroups.append(aa[0])
            #Make list of these amino-acids for each of the ingroup species 
            corresponding_aa_ingroup.append(aa[3])
            #Returns 2 lists 
            #list of the conserved amino acids between the outgroup species
            #list of the same amino-acids for ingroup species
    
    return conserved_aa_outgroups, corresponding_aa_ingroup 

#Can comapre the conserevd residues of the outgroups to the ingroups to get Sneath value for  ancestral sites
#Therefore If the amino acid changes in the ingroup, we can calcualte Sneath index 
def get_sneath_value(gene_lst_outgroup, gene_lst_ingroup):
    
    sneath_matrix = {
    "L": {"L": 0, "I": 5, "V": 9, "G": 24, "A": 15, "P": 23, "Q": 22, "N": 20, "M": 20, "T": 23, "S": 23, "C": 24, "E": 30, "D": 25, "K": 23, "R": 33, "Y": 30, "F": 19, "W": 30, "H": 25},
    "I": {"L": 5, "I": 0, "V": 7, "G": 25, "A": 17, "P": 24, "Q": 24, "N": 23, "M": 22, "T": 21, "S": 25, "C": 26, "E": 31, "D": 28, "K": 24, "R": 34, "Y": 34, "F": 22, "W": 34, "H": 28},
    "V": {"L": 9, "I": 7, "V": 0, "G": 19, "A": 12, "P": 20, "Q": 25, "N": 23, "M": 23, "T": 17, "S": 20, "C": 21, "E": 31, "D": 28, "K": 26, "R": 36, "Y": 36, "F": 26, "W": 37, "H": 31},
    "G": {"L": 24, "I": 25, "V": 19, "G": 0, "A": 9, "P": 17, "Q": 32, "N": 26, "M": 34, "T": 20, "S": 19, "C": 21, "E": 37, "D": 33, "K": 31, "R": 43, "Y": 36, "F": 29, "W": 39, "H": 34},
    "A": {"L": 15, "I": 17, "V": 12, "G": 9, "A": 0, "P": 16, "Q": 26, "N": 25, "M": 25, "T": 20, "S": 16, "C": 13, "E": 34, "D": 30, "K": 26, "R": 37, "Y": 34, "F": 26, "W": 36, "H": 29},
    "P": {"L": 23, "I": 24, "V": 20, "G": 17, "A": 16, "P": 0, "Q": 33, "N": 31, "M": 31, "T": 25, "S": 24, "C": 25, "E": 43, "D": 40, "K": 31, "R": 43, "Y": 37, "F": 27, "W": 37, "H": 36},
    "Q": {"L": 22, "I": 24, "V": 25, "G": 32, "A": 26, "P": 33, "Q": 0, "N": 10, "M": 13, "T": 24, "S": 21, "C": 22, "E":14, "D": 22, "K": 21, "R": 23, "Y": 29, "F": 24, "W": 31, "H": 27},
    "N": {"L": 20, "I": 23, "V": 23, "G": 26, "A": 25, "P": 31, "Q": 10, "N": 00, "M": 21, "T": 19, "S": 15, "C": 19, "E": 19, "D": 14, "K": 27, "R": 31, "Y": 28, "F":24, "W": 32, "H": 24},
    "M": {"L": 20, "I": 22, "V": 23, "G": 34, "A": 25, "P": 31, "Q": 13, "N": 21, "M": 0, "T": 25, "S": 22, "C": 17, "E": 26, "D": 31, "K": 24, "R": 28, "Y": 32, "F": 24, "W": 31, "H" : 30},
    "T": {"L": 23, "I": 21, "V": 17, "G": 20, "A": 20, "P": 25, "Q": 24, "N": 19, "M": 25, "T": 0, "S": 12, "C": 19, "E": 34, "D": 29, "K": 34, "R": 38, "Y": 32, "F": 28, "W": 38, "H": 34},
    "S": {"L": 23, "I": 25, "V": 20, "G": 19, "A": 16, "P": 24, "Q": 21, "N": 15, "M": 22, "T": 12, "S": 0, "C": 13, "E": 29, "D": 25, "K": 31, "R": 37, "Y": 29, "F": 25, "W": 35, "H": 28},
    "C": {"L": 24, "I": 26, "V": 21, "G": 21, "A": 13, "P": 25, "Q": 22, "N": 19, "M": 17, "T": 19, "S": 13, "C": 0, "E": 33, "D": 28, "K": 32, "R": 36, "Y": 34, "F": 29, "W": 37, "H": 31},
    "E": {"L": 30, "I": 31, "V": 31, "G": 37, "A": 34, "P": 43, "Q": 14, "N": 19, "M": 26, "T": 34, "S": 29, "C": 33, "E": 0, "D": 7, "K": 26, "R": 31, "Y": 34, "F": 35, "W": 43, "H": 27},
    "D": {"L": 25, "I": 28, "V": 28, "G": 33, "A": 30, "P": 40, "Q": 22, "N": 14, "M": 31, "T": 29, "S": 25, "C": 28, "E": 7, "D": 0, "K": 34, "R": 39, "Y": 34, "F": 35, "W": 45, "H": 35},
    "K": {"L": 23, "I": 24, "V": 26, "G": 31, "A": 26, "P": 31, "Q": 21, "N": 27, "M": 24, "T": 34, "S": 31, "C": 32, "E": 26, "D": 34, "K": 0, "R": 14, "Y": 34, "F": 28, "W": 34, "H": 27},
    "R": {"L": 33, "I": 34, "V": 36, "G": 43, "A": 37, "P": 43, "Q": 23, "N": 31, "M": 28, "T": 38, "S": 37, "C": 36, "E": 31, "D": 39, "K": 14, "R": 0, "Y": 36, "F": 34, "W": 36, "H": 31},
    "Y": {"L": 30, "I": 34, "V": 36, "G": 36, "A": 34, "P": 37, "Q": 29, "N": 28, "M": 32, "T": 32, "S": 29, "C": 34, "E": 34, "D": 34, "K": 34, "R": 36, "Y": 0, "F": 13, "W": 21, "H": 23},
    "F": {"L": 19, "I": 22, "V": 26, "G": 29, "A": 26, "P": 27, "Q": 24, "N": 24, "M": 24, "T": 28, "S": 25, "C": 29, "E": 35, "D": 35, "K": 28, "R": 34, "Y": 13, "F": 0, "W": 13, "H": 18},
    "W": {"L": 30, "I": 34, "V": 37, "G": 39, "A": 36, "P": 37, "Q": 31, "N": 32, "M": 31, "T": 38, "S": 35, "C": 37, "E": 43, "D": 45, "K": 34, "R": 36, "Y": 21, "F": 13, "W": 0, "H": 25},
    "H": {"L": 25, "I": 28, "V": 31, "G": 34, "A": 29, "P": 36, "Q": 27, "N": 24, "M": 30, "T": 34, "S": 28, "C": 31, "E": 27, "D": 35, "K": 27, "R": 31, "Y": 23, "F": 18, "W": 25, "H": 0},
    }
    
    S = 0
    for i, j in zip(gene_lst_outgroup, gene_lst_ingroup): #Compare each amino-acid in parallel 
        if (i == "-") or (j == "-"): 
            #If the aa is missing (indicated by -) do not consider for Sneath calculating 
            continue
        if i != j: #If residue is different between ancestral state and ingroup
            original_amino_acid = sneath_matrix[i] #Look for the values in the original amino-acid 
            dissimilarity = original_amino_acid[j] #Find the value of the dfference
            S += dissimilarity #Add this to the sneath value
    
    if len(gene_lst_outgroup) == len(gene_lst_ingroup):
        #Get "adjusted sneath value" by dividing by protein alignemnt length 
        adjusted_sneath_value = S / len(gene_lst_outgroup)
    else:
        print("Genes not of same length")
        #Divide by average of the two genes 
        adjusted_sneath_value = S / ( ( len(gene_lst_outgroup) + len(gene_lst_ingroup) ) / 2 ) 
        
    return adjusted_sneath_value

#These next two fucntions: make_busco_id_dict() and get_R_output simply re-format the sneath_values.tsv file 
#to a format that can be plotted in R whereby one species on x and anotehr on y, plotting the coorrds of each point
#as the sneath for (species_x, species_y)

def make_busco_id_dict(sneath_output):
    #Create a dictionary to store the values for each busco_id and species
    busco_id_dict = defaultdict(dict)
    with open(sneath_output) as f:
        for row in f:
            #Loop through the input table and add the values to the dictionary
            busco_id = row.split("\t")[0]
            sp = row.split("\t")[1]
            sneath = row.split("\t")[2].strip()
            
            busco_id_dict[busco_id][sp] = sneath
    
    return busco_id_dict    


# To get GC3 alongside the snetah values for each species, need to make a dictionary for every species with GC3 value
def get_ingroup_id_GC3_dict(tsv, sp_name):
    id_GC3_dict = {}
    with open(tsv) as f:        
        for line in f:
            if line.startswith("BUSCO_id"): #Skip header of tsv file
                continue
            else:
                busco_id = line.split("\t")[0]
                species = line.split("\t")[1]
                GC3 = line.split("\t")[2]

                if species == sp_name:
                    id_GC3_dict[busco_id] = GC3
    
    return id_GC3_dict

# Output new file as a table in format
# busco_id   SP1 SP2 SP3
# busco_id1  s   s   s
# busco_id2  s   s   s
def create_tsv(busco_id_dict, ingroup1, ingroup2, ingroup3):
    with open("sneath_R_output.tsv", "w") as outf:
        #Open and write headers for table
        outf.write("{:<9}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\n".format("BUSCO_id", ingroup1, f"{ingroup1}_GC3", ingroup2, f"{ingroup2}_GC3", ingroup3, f"{ingroup3}_GC3"))
        
        #MAke three dictionaries, each of the busco_id : GC3 for each ingroup species
        ingroup1_dict = get_ingroup_id_GC3_dict(args.trimmed, ingroup1)
        ingroup2_dict = get_ingroup_id_GC3_dict(args.trimmed, ingroup2)
        ingroup3_dict = get_ingroup_id_GC3_dict(args.trimmed, ingroup3)

        for busco_id, species_sneath in busco_id_dict.items():
            
            #Retieve sneath value ffor species 1
            species1 = species_sneath.get(ingroup1, "")
            #Get the gC3 for ingroup species 1 
            ingroup1_GC3 = ingroup1_dict[busco_id]       
            
            species2 = species_sneath.get(ingroup2, "")
            ingroup2_GC3 = ingroup2_dict[busco_id]

            species3 = species_sneath.get(ingroup3, "")
            ingroup3_GC3 = ingroup3_dict[busco_id]

            outf.write("{:<9}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\t{:<23}\n".format(busco_id, species1, ingroup1_GC3, species2, ingroup2_GC3, species3, ingroup3_GC3))


def main():
    #Using fucntions from get_genetree.py to move files into a new dir 
    BUSCO_lst = parse_tsv(args.trimmed) #Get list of relevant busco_ids
    copy_files(BUSCO_lst, args.busco, args.out) #Copy files to new directory
    print(f"Copied {len(BUSCO_lst)} BUSCOs to {args.out}") 
    #Running trimal to get amino-acid alignemnt 
    print(f"Running trimal (gappyout) for {len(BUSCO_lst)} busco_ids...")
    run_trimal(args.out) 
    #Get species name to make visual comparisons easier
    rename_trimal_files(args.out)

    #Calculating Sneath value: 
    #Parse species.txt file to get names of species in analysis 
    outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3 = get_species(args.species)
    #Create alignment .sneath file containing only the 6 species of intrest
    parse_alignemts(args.out, outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3)
    #Main sneath function - get sneath outputput using the comapare_sneath_residues() function 
    comapare_sneath_residues(args.out, outgroup1, outgroup2, outgroup3, ingroup1, ingroup2, ingroup3)

    #Get appropriate output for R analysis:
    busco_id_dict = make_busco_id_dict("sneath_values.tsv")
    create_tsv(busco_id_dict, ingroup1, ingroup2, ingroup3)

if __name__ == "__main__":
    parse = argparse.ArgumentParser()

    parse.add_argument("-t", "--trimmed",type=str, help="path to trimmedGC3_ortbusco_idroups.tsvfile to pull info", required = True)
    parse.add_argument("-b", "--busco",type=str, help="path to BUSCO fastas", required = True)
    parse.add_argument("-s", "--species",type=str, help="list of species for analysis", required = True)
    parse.add_argument("-o", "--out",type=str, help="path to output for busco_ids", required = True)

    args = parse.parse_args()

    if args.busco[-1] != '/':
        args.busco += '/'

    main()
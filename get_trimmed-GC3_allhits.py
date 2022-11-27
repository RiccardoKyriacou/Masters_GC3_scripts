import os
import re
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq

'''
Script to calculate GC3 for
aligned cds fasta files
'''

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="Path to directiory containing .trimal alignments",required=True)
parse.add_argument("-t", "--tsv",type=str, help="Path to combined_outlier_genes.tsv from blast",required=True)

args = parse.parse_args()

os.makedirs(args.input, exist_ok=True)
if args.input[-1] != '/':
    args.input = args.input + '/'

#Creating dictionary of gene names
gene_dict = {}
with open(args.tsv, "r") as f:
    for line in f:
        gene_id = line.split("\t")[0].split("_")[0]
        gene_drome = re.findall(".+description:(.+?)\t", line) #Extracting genes with headers containing "description: [gene name]"
        if len(gene_drome) > 0:
            gene_hits = "\t".join(gene_drome)
            gene_dict[gene_id] = gene_hits 
        gene_cds = re.findall(".+cds (.+?)\t", line) #Extracting genes with headers containing "cds [gene name]"
        if len(gene_cds) > 0:
            gene_cds_hits = "\t".join(gene_cds)
            gene_dict[gene_id] = gene_cds_hits
        gene_protein = re.findall("protein (.+?)\t", line) #Extracting genes with headers containing "cds [gene name]"
        if len(gene_protein) > 0:
            gene_protein_hits = "\t".join(gene_protein)
            gene_dict[gene_id] = gene_protein_hits
 
#Calculating GC3
outF1 = open("allspecies_trimmed_GC3.tsv", "w") 
outF1.write("Species_Name_ID" + "\t" + "Blast_Hit" + "\t"  + "GC3" + "\t" + "GC3_Length" + "\n") 
for trimmed in glob.glob(args.input + "*.trimal"): 
    sp_name = trimmed.split("/")[-1].split("_blast")[0] 
    with open(trimmed) as f:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.description.split("_")[0]
            seq = str(record.seq)
            GC3 = seq[2::3]
            GC3_count = GC3.count('G') + GC3.count('g') + GC3.count('C') + GC3.count('c')
            GC3_percent = ( GC3_count / len(GC3) ) * 100
            if not str(ID).startswith( ("Bombyx", "Drosophila", "Heliconius") ):
                gene_name = gene_dict[str(ID)]
                outF1.write(sp_name + "\t" + gene_name + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\n")
            elif str(ID).startswith("Bombyx"):
                outF1.write("Bombyx_mori" + "\t" + "Gene_hit" + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\n")
            elif str(ID).startswith("Drosophila"):
                outF1.write("Drosophila_melanogaster" + "\t" + "Gene_hit" + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\n")
            elif str(ID).startswith("Heliconius"):
                outF1.write("Heliconius_melpomene" + "\t" + "Gene_hit" + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\n")
                
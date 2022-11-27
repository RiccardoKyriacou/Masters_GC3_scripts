import os
import re
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from subprocess import call as unix

'''
Script to calculate aligned GC3 
for orthofinder orthogroups 
given cds primary transcripts
'''

parse = argparse.ArgumentParser()

parse.add_argument("-s", "--singlecopy",type=str, help="path to 1to1_orthologues directory (copied from Single_Copy_Orthologue_Sequences from orthofinder run)",required=True)
parse.add_argument("-c", "--cds",type=str, help="Path to cds primary transcripts file",required=True)
#Note, primary transcripts folder must be within the cds folder downloaded from ensemble 
parse.add_argument("-g", "--genome",type=str, help="Path to genomes file",required=True)

args = parse.parse_args()

if args.singlecopy[-1] != '/':
    args.singlecopy = args.singlecopy + '/'

if args.cds[-1] != '/':
    args.cds = args.cds + '/'

if args.genome[-1] != '/':
    args.genome = args.genome + '/'

#Dictionary of cds IDS : cds seq
cds_dict = {}
for cds in glob.glob(args.cds + "*.fa"):
    with open(cds) as f:
        for record in SeqIO.parse(f, 'fasta'):            
            ID_cds = record.id
            cds_seq = str(record.seq)
            cds_dict[ID_cds] = cds_seq

#retun cds given orthogroup ID
for orthogroup in glob.glob(args.singlecopy + "*.fa"):
    OGname = orthogroup.split("/")[-1].split(".")[0]
    with open(orthogroup) as f1, open(args.singlecopy + OGname + ".fa.cds", "w") as OutF1:
        for record in SeqIO.parse(f1, 'fasta'):            
            ID_OG = record.id
            cds_OG = cds_dict[ID_OG]
            OutF1.write(">" + ID_OG + "\n" + cds_OG+ "\n")

print("Running mafft...")
for orthogroup in glob.glob(args.singlecopy + "*.fa"):
    OGname = orthogroup.split("/")[-1].split(".")[0]
    unix("mafft " + orthogroup + " > " + args.singlecopy + OGname + ".fa.mafft", shell=True)

print("Running trimal (gappyout) with backtranslation...")
for alignemnt in glob.glob(args.singlecopy + "*.mafft"):
    OGname = alignemnt.split("/")[-1].split(".")[0]
    unix("trimal -in " + alignemnt + " -out " + alignemnt + ".trimal" + " -gappyout -backtrans " + args.singlecopy + OGname + ".fa.cds"+ " -ignorestopcodon", shell=True)

print("Calculating aligned GC3...")

#Make dict of species name and ID
name_dict = {}
for cds in glob.glob(args.cds + "*.fa"):
    sp_name = cds.split("/")[-1].split("-")[0]
    with open(cds) as f:
        for record in SeqIO.parse(f, 'fasta'):
                ID_cds = record.id
                name_dict[ID_cds] = sp_name

#Dict of dicts which has for each species, the length of each chromosome 
species_genome_dict = defaultdict(dict)
for files in glob.glob(args.genome + "*.fasta"):
    genome_dict = {}
    if "Drosophila_melanogaster.BDGP6" in files: #/home/zoo/quee4075/GC_Chrm_runs/GC_genome_run-1/data/genomes/Drosophila_melanogaster.BDGP6.32.dna.toplevel.fasta
        with open(files) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.id
                sp_name = "Drosophila_melanogaster"
                seq = str(record.seq)
                len_chrm = len(seq)
                genome_dict[header] = len_chrm
    elif "BomTerr1" in files:
        with open(files) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                sp_name = "Bombus_terrestris"
                if "chromosome" in header:
                    chrm_no = header.split(" ")[4]
                    seq = str(record.seq)
                    len_chrm = len(seq)
                    genome_dict[chrm_no] = len_chrm
                elif "unplaced" in header:
                    unplaced_id = record.id
                    seq = str(record.seq)
                    len_chrm = len(seq)
                    genome_dict[unplaced_id] = len_chrm
    else:
        with open(files) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                sp_name = header.split(" ")[1] + "_" + header.split(" ")[2]
                if "chromosome:" in header:
                    chrm_no = header.split(": ")[1]
                    seq = str(record.seq)
                    len_chrm = len(seq)
                    genome_dict[chrm_no] = len_chrm
                elif "contig:" in header:
                    contig_id = record.id 
                    seq = str(record.seq)
                    len_chrm = len(seq)
                    genome_dict[contig_id] = len_chrm

    species_genome_dict[sp_name] = genome_dict

#Two dicts = cds:start_pos to get position and cds:chrm_no to get chromsome number 
pos_dict = {}
chrm_no_dict = {}
cds_full = args.cds.split("primary")[0] #Get one dir back to full cds 
for cds in glob.glob(cds_full + "*.fa"):
    if "Drosophila_melanogaster" in cds:
        with open(cds) as f:
            for record in SeqIO.parse(f, 'fasta'): #primary_assembly:BDGP6.32:3R:16375284:16403690:-1
                header = record.description
                ID = record.description.split("gene:")[1].split(" ")[0] #Gene ID not transcript ID
                chrm_no = header.split(" ")[2].split(":")[2]
                start_pos = header.split(" ")[2].split(":")[3]
                seq = str(record.seq)
                pos_dict[ID] = start_pos
                chrm_no_dict[ID] = chrm_no
    else:
        with open(cds) as f:
            for record in SeqIO.parse(f, 'fasta'):
                header = record.description
                ID = record.description.split("gene:")[1].split(" ")[0]
                chrm_no = header.split(" ")[2].split(":")[1]
                start_pos = header.split(" ")[2].split(":")[2]
                seq = str(record.seq)
                pos_dict[ID] = start_pos
                chrm_no_dict[ID] = chrm_no

with open(args.singlecopy + "trimmedGC3_ortogroups.tsv", "w") as OutF1, open(args.singlecopy + "outliers_trimmedGC3.tsv", "w") as OutF2:
    OutF1.write("Orthogroup\tspecies\tgene_ID\tGC3_percent\tGC3_length\tchrm_no\tstart_position\tchrm_size\n")
    OutF2.write("Orthogroup\tspecies\tgene_ID\tGC3_percent\tGC3_length\tchrm_no\tstart_position\tchrm_size\n")
    for cds in glob.glob(args.singlecopy + "*.trimal"):
        OGname = cds.split("/")[-1].split(".")[0] 
        with open(cds) as f:
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id             
                chrm_name = chrm_no_dict[ID]
                start_pos = pos_dict[ID]

                species = name_dict[ID]
                species_chrm = species_genome_dict[species]
                chrm_size = species_chrm[chrm_name]

                seq = str(record.seq)
                GC3 = seq[2::3]
                GC3_count = GC3.count('G') + GC3.count('g') + GC3.count('C') + GC3.count('c')
                GC3_percent = ( GC3_count / len(GC3) ) * 100
                OutF1.write(OGname + "\t" + species + "\t" + ID  + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\t" + chrm_name + "\t" + start_pos + "\t" + str(chrm_size) +"\n")

                if GC3_percent > 95:
                    OutF2.write(OGname + "\t" + species + "\t" + ID  + "\t" + str(GC3_percent) + "\t" + str(len(GC3)) + "\t" + chrm_name + "\t" + start_pos + "\t" + str(chrm_size) +"\n")
                else:
                    continue 

#TODO Add comparison code under here 

orthogroup_lst = []
with open(args.singlecopy + "outliers_trimmedGC3.tsv") as f:
    for line in f:
        if line.startswith("Orthogroup"):
            continue
        else:
            gene_id = line.split("\t")[2]
            orthogroup = line.split("\t")[0]
            orthogroup_lst.append(orthogroup)

with open(args.singlecopy + "trimmedGC3_ortogroups.tsv") as f1, open("outlier_GC3-comparison.tsv", "w") as OutF:
    OutF.write("Orthogroup" + "\t" + "Species" + "\t" + "gene_ID" + "\t" + "chrm_no" + "\t" + "GC3" + "\t" + "Outlier" + "\t" + "Position" + "\n")   
    for line in f1:
        if line.startswith("Orthogroup"):
            continue
        else: 
            orthogroup = line.split("\t")[0]
            if orthogroup in set(orthogroup_lst):
                gene_ID = line.split("\t")[2]
                sp = line.split("\t")[1]
                GC3 = line.split("\t")[3]
                chrm_no = line.split("\t")[5]
                start_position = int(line.split("\t")[6])
                chrm_size = int(line.split("\t")[7])

                window = int(1_000_000) 
                near_telomere = ""
                if start_position < window:
                    near_telomere += "Near_Telomere"
                elif start_position > (chrm_size - window):
                    near_telomere += "Near_Telomere"
                else:
                    near_telomere = near_telomere

                if float(GC3) > 95:
                    OutF.write(orthogroup + "\t" + sp + "\t" + gene_ID + "\t" + chrm_no + "\t" + GC3 + "\t" + "Outlier_GC3" + "\t" + near_telomere + "\n")
                elif float(GC3) > 85: 
                    OutF.write(orthogroup + "\t" + sp + "\t" + gene_ID + "\t" + chrm_no + "\t" + GC3 + "\t" + "High_GC3" + "\t" + near_telomere + "\n")
                else:
                    OutF.write(orthogroup + "\t" + sp + "\t" + gene_ID + "\t" + chrm_no + "\t" + GC3 + "\t" + "" + "\t" + near_telomere + "\n")




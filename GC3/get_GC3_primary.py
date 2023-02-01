#!/usr/bin/env python3
import glob
import argparse
from Bio import SeqIO
from collections import defaultdict


'''
Script to calculate GC3 content per
chromosome/scaffold in a genome with
a given cutoff size. Able to use cds 
and primary tenplate data
'''

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to cds fasta files",required=True)
parse.add_argument("-c", "--cutoff",type=int, help="cutoff for GC3 to be included as outliers ",required=True)
parse.add_argument("-o", "--outfile",type=str, help="name of output file",required=True)
parse.add_argument("-i", "--info",type=str, help="tsv file with species phylum information",required=False)
parse.add_argument("-n", "--nucleotides",type=int, help="amount of nucleotides allowed for outlier genes", required=False)

args = parse.parse_args()

#If species-clade info file provided, store data
def store_info():
    sp_group = {}
    with open(args.info) as f:
        for line in f:
            lines = line.split('\t')
            species = lines[0]
            group = lines[1].strip()
            sp_group[species] = group
    return sp_group

#Prints out the number of genomes calcualting GC for 
def genome_list():
    genome_list = []
    for fasta in glob.glob(args.path + '*'):
        if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                genome_list.append(fasta)
    print('Calculating GC content for', len(genome_list), 'genomes...\n')

#Fucntion to create a file of outleir genes  - not really useful for GC but nice to have to get an idea 
def outlier_cutoff(GC_content, sp_name, header, seq, seq_len, chrom_info, outf1, chrm):
    if (seq_len > args.nucleotides) and (GC_content > args.cutoff): #Looking cutooff of 950 3rd codon positions, not nucelotides
        if seq.startswith("ATG"):
            outf1.write(f"{sp_name}\t{chrm}\t{header}_{chrom_info}\n")
        else:
            outf1.write(f"{sp_name}\t{chrm}\t{header}_{chrom_info}\tNo_ATG_StartCodon\n")
    
def count_GC3(sp_group):
    with open(args.outfile, 'w') as outf, open(f"outlierGC3_cds_c{args.cutoff}_n{args.nucleotides}.tsv", 'w') as outf1:
        for fasta in glob.glob(args.path + '*'):
            if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                genome = fasta.split('/')[-1]
                sp_name = genome.split('.')[0].split('-')[0] # Lineus_longissimus-GCA_910592395.2-2022_03-cds.fa
                if args.info:
                    group = sp_group[sp_name]

                chr_cds_seq = defaultdict(list) #Dictonary of list of cds
                print(sp_name)
                with open(fasta) as f:
                    first_line = f.readline()
                    if ' ' in first_line:
                        continue
                    else:
                        cds_path = args.path.split('primary')[0]
                        cds_file = glob.glob(cds_path + genome)
                        cds_file = cds_file[0]
                        geneID_header = {}
                        with open(cds_file) as f1:
                            for record in SeqIO.parse(f1, 'fasta'):
                                header = record.description
                                geneID = header.split(' ')[3].split(':')[1]
                                geneID_header[geneID] = header

                    for record in SeqIO.parse(f, 'fasta'): #Parsing cds files to get gene info
                        header = record.description
                        try:
                            chrom = header.split(':')[1]
                            chrom_info = header.split(" ")[2].split(":")[1:4]
                            chrom_info = ":".join(chrom_info)
                        except:
                            chrom_full = geneID_header[header]
                            chrom = chrom_full.split(':')[1]
                            chrom_info = chrom_full.split(" ")[2].split(":")[1:4]
                            chrom_info = ":".join(chrom_info)

                        seq = str(record.seq)
                        seq_len = len(seq)
                    
                        GC3 = seq[2::3] #Caculating GC3 

                        GC3_count = GC3.count('G') + GC3.count('g') + GC3.count('C') + GC3.count('c')#Counting GC content to find cutoff
                        GC3_len = len(GC3) 
                        GC_content = (GC3_count/GC3_len)*100
                        
                        outlier_cutoff(GC_content, sp_name, header, seq, seq_len, chrom_info, outf1, chrom)
            
                        chr_cds_seq[chrom].append(GC3)#Append all cds sequences to each chromosome in a dictionary

                for chrm, cds_regions in chr_cds_seq.items():#For each chromosome and all cds sequences in the dictionary
                    cds_joined = ''.join(cds_regions)#Join all cds sequences for a given chromosome into one big sequence
                    len_cds_region = len(cds_joined)
                    number_genes = len(cds_regions)
                    GC3_count = cds_joined.count('G') + cds_joined.count('g') + cds_joined.count('C') + cds_joined.count('c')
                    GC_content = (GC3_count/len_cds_region)*100

                    #Output file: sp    group   chrm    Gc_content  number of genes
                    if args.info:
                        outf.write(f"{sp_name}\t{group}\t{chrm}\t{str(GC_content)}\t{str(number_genes)}\n")
                    else:
                        outf.write(f"{sp_name}\t{chrm}\t{str(GC_content)}\t{str(number_genes)}\n")

def main():
    if args.info:
        sp_group = store_info()
        genome_list()
        count_GC3(sp_group)
    else:
        sp_group = None
        genome_list()
        count_GC3(sp_group)

if __name__ == "__main__":
    main()

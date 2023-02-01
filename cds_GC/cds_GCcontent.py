#!/usr/bin/env python3
import glob
import argparse
import statistics
from Bio import SeqIO
from collections import defaultdict

'''
Script to calculate cds GC content 
per chromosome/scaffold in a genome
with a given cutoff size.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to cds fasta files",required=True)
parse.add_argument("-c", "--cutoff",type=int, help="GC cutoff",required=True)
parse.add_argument("-o", "--outfile",type=str, help="name of output file",required=True)
parse.add_argument("-i", "--info",type=str, help="tsv file with species phylum information",required=False)

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
def outlier_cutoff(GC_content, cutoff, sp_name, header):
    with open('outlier_cds.tsv', 'w') as outf1:
        if GC_content > cutoff: #Find outlier genes with GC content above a certain threshold
            outf1.write(f"{sp_name}\t{header}\n")

def count_GC3(sp_group):
    with open(args.outfile, 'w') as outf:
        for fasta in glob.glob(args.path + '*'):
            if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                genome = fasta.split('/')[-1]
                sp_name = genome.split('.')[0].split('-')[0] # Split file name to get sp name  - Lineus_longissimus-GCA_910592395.2-2022_03-cds.fa
                if args.info:
                    group = sp_group[sp_name]

                chr_cds_seq = defaultdict(list) #Dictonary of list of cds
                print(sp_name)
                with open(fasta) as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        header = record.description
                        chrom = header.split(':')[1]
                        seq = str(record.seq)

                        GC_count = seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')#Counting GC content to find cutoff
                        seq_len = len(seq)
                        GC_content = (GC_count/seq_len)*100

                        outlier_cutoff(GC_content, args.cutoff, sp_name, header)        

                        chr_cds_seq[chrom].append(seq)#Append all cds sequences to each chromosome in a dictionary

                for chrm, cds_regions in chr_cds_seq.items():#For each chromosome and all cds sequences in the dictionary
                    cds_joined = ''.join(cds_regions)#Join all cds sequences for a given chromosome into one big sequence
                    len_cds_region = len(cds_joined)
                    number_genes = len(cds_regions)
                    GC_count = cds_joined.count('G') + cds_joined.count('g') + cds_joined.count('C') + cds_joined.count('c')
                    GC_content = (GC_count/len_cds_region)*100

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

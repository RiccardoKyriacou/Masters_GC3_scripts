import glob
import argparse
import statistics
from Bio import SeqIO
from collections import defaultdict

'''
Script to calculate GC content per
chromosome/scaffold in a genome
with a given cutoff size.
'''

parse = argparse.ArgumentParser()

parse.add_argument("-p", "--path",type=str, help="path to genome fasta files",required=True)
parse.add_argument("-c", "--cutoff",type=int, help="scaffold size cutoff",required=True)
parse.add_argument("-o", "--outfile",type=str, help="name of output file",required=True)
parse.add_argument("-i", "--info",type=str, help="tsv file with species phylum information",required=False)

args = parse.parse_args()

#If species-clade info file provided, store data
def store_info(args_info):
    sp_names = {}
    sp_group = {}
    with open(args_info) as f:
        for line in f:
            lines = line.split('\t')
            species = lines[0]
            group = lines[2].strip()
            sp_short= lines[1]
            sp_names[sp_short] = species
            sp_group[sp_short] = group
    return sp_names, sp_group

def genome_list(args_path):
    genome_list = []            
    for fasta in glob.glob(args_path + '*'):
        if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                genome_list.append(fasta)
    print('Calculating GC content for', len(genome_list), 'genomes...\n')
    return genome_list
  
#Fucntion to calcualte GC
def calculate_GC(args_outfile, args_path, sp_names, sp_group):
    with open(args_outfile, 'w') as outF:
        for fasta in glob.glob(args_path + '*'):
            if fasta.endswith(('.fa', '.fasta', '.fas', '.fna')):
                genome = fasta.split('/')[-1]
                sp_name = genome.split('.')[1].split('_')[1]
                if args.info:
                    species = sp_names[sp_name]
                    group = sp_group[sp_name]
                print(sp_name)
                with open(fasta) as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        seq = str(record.seq)
                        seq_len = len(seq)
                        GC_count = seq.count('G') + seq.count('g') + seq.count('C') + seq.count('c')
                        GC_content = (GC_count/seq_len)*100 #Measure GC content by dividing number of Gs and Cs by the scaffold length
                        if seq_len >= args.cutoff:
                            if args.info:
                                outF.write(f"{species}\t{group}\t{str(seq_len)}\t{str(GC_content)}\n")
                            else:
                                outF.write(f"{sp_name}\t{str(seq_len)}\t{str(GC_content)}\n")
                            
def main():
    if args.info:
        sp_names, sp_group = store_info(args.info)
        genome_list(args.path)
        calculate_GC(args.outfile, args.path, sp_names, sp_group)
    else:
        sp_names = None
        sp_group = None
        genome_list(args.path)
        calculate_GC(args.outfile, args.path, sp_names, sp_group)

if __name__ == "__main__":
    main()

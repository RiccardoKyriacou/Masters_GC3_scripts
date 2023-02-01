import os
import re
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from subprocess import call as unix
from joblib import Parallel, delayed

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="name of input tsv file to get OGs",required=True)
parse.add_argument("-f", "--OGfasta",type=str, help="path to orthogroup fasta files",required=True)
parse.add_argument("-o", "--output",type=str, help="name of directory to output files",required=True)
parse.add_argument("-b", "--protdb",type=str, help="path to directory with protdb files",required=True)

args = parse.parse_args()

def check_args():
    os.makedirs(args.output, exist_ok=True)
    if args.output[-1] != '/':
        args.output = args.output + '/'

    if args.OGfasta[-1] != '/':
        args.OGfasta = args.OGfasta + '/'

    args.protdb = os.path.abspath(args.protdb)

    if args.protdb[-1] != '/':
        args.protdb = args.protdb + '/'

def parse_tsv(fn): #Parse input tsv file and get OGs
    sp_OG_dict = defaultdict(list)
    with open(fn) as f:
        for line in f:
            if line.startswith("Orthogroup"):
                continue
            else: 
                lines = line.split('\t') 
                OG = lines[0]
                species = lines[1].strip()
                sp_OG_dict[OG].append(species) #Create dict of OGs : (list of species)
    return sp_OG_dict

def blast_OGs(sp_OG_dict, args_OGfasta, args_protdb, args_output): #Blast each outlier OG
    with open(args_output + "outliers_BLAST.tsv", "w") as outf:
        for OG, sp in sp_OG_dict.items(): #
            OG_fasta = glob.glob(args_OGfasta + OG + '.fa') 
            OG_fasta = OG_fasta[0] 
            print(f"BLASTing {OG} for {sp}")
            unix('blastp -query ' + OG_fasta + ' -db ' + args_protdb + '/insect_pep -out ' + args_output + OG + '_blastoutput.tsv -outfmt "6 qseqid sseqid stitle evalue pident bitscore qstart qend qlen sstart send slen" -evalue 1e-30 -max_target_seqs 5', shell=True)
            for file in glob.glob(args_output + OG + "*"): #Open blast output files
                with open(file) as f:
                    gene_lst = [] #Make list of the set of genes 
                    for line in f:
                        if "description:" in line: #Parse gene name from blast file for each gene
                            gene_name = line.split("description:")[1].split("/t")[0].split("\t")[0]
                            gene_lst.append(gene_name)

                    gene_lst = set(gene_lst) #Make set such that only unique gene hits for orthogroups

                    print(f"{', '.join(gene_lst)}\t{', '.join(sp)}\t{OG}\n")
                    outf.write(f"{', '.join(gene_lst)}\t{', '.join(sp)}\t{OG}\n")

def main():
    OG_lst = parse_tsv(args.input)
    blast_OGs(OG_lst, args.OGfasta, args.protdb, args.output)

if __name__ == "__main__":
    check_args()
    main()


#python3 orthofinder_outlier_blast.py --input ../results/GC3/1to1_orthologues/outliers_trimmedGC3.tsv --OGfasta ../results/GC3/1to1_orthologues/



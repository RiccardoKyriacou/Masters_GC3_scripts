import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from subprocess import call as unix

'''
Script to calculate aligned GC3 
for BUSCO orthogroups 

Usage: python3 get_BUSCO_GC3.py --busco ../data/fasta/ -o 95
'''

# Check arguments 
def check_args():
    if args.busco[-1] != '/':
        args.busco += '/'

# First step, protein level alignemnt, therefore need to run mafft
# Run mafft to align BUSCOs
def run_mafft(busco_dir):
    print("Running mafft...")
    proteins = f"{busco_dir}proteins/" #Path to protein dir
    for BUSCO in glob.glob(proteins + "*.faa"):
        BUSCO_name = BUSCO.split("/")[-1].split(".")[0]
        unix(f"mafft {BUSCO} > {proteins}{BUSCO_name}.faa.mafft", shell=True)

# Run trimal to trim alignemnts and backtrans to get DNA seq to coint GC3
def run_trimal(busco_dir):
    print("Running trimal (gappyout) with backtranslation...")
    proteins = f"{busco_dir}proteins/" #Path to protein dir
    nucleotides = f"{busco_dir}nucleotides/" #path to nuceotide dir 
    for alignemnt in glob.glob(proteins + "*.mafft"):
        BUSCO_name = alignemnt.split("/")[-1].split(".")[0]
        unix(f"trimal -in {alignemnt} -out {alignemnt}.trimal -backtrans {nucleotides}{BUSCO_name}.fna -ignorestopcodon -gappyout", shell=True)

# Get three outputs based on GC3 counts:
# 1) trimmedGC3_BUSCO.tsv - tsv file of all BUSCOS with GC3 counts
# 2) outlierGC3_95.tsv - tsv file of only outliers in same format 
def get_trimmed_GC3_output(busco_dir):
    with open(f"{busco_dir}trimmedGC3_BUSCO.tsv", "w") as outf1, open(f"{busco_dir}outlierGC3_{args.outlier}.tsv", "w") as outf2:
        #Write heading for output files 
        outf1.write("BUSCO_id\tspecies\tGC3_percent\tGC3_length\tchrm_name\tstart_position\n")
        outf2.write("BUSCO_id\tspecies\tGC3_percent\tGC3_length\tchrm_name\tstart_position\n")
        #Path to protein dir
        proteins = f"{busco_dir}proteins/" 
        for cds in glob.glob(proteins + "*.trimal"):
            BUSCO_name = cds.split("/")[-1].split(".")[0] 
            print(f"Calcuating GC3 for {BUSCO_name} ...")
            with open(cds) as f:
                for record in SeqIO.parse(f, 'fasta'):
                    header = record.id
                    #Headers for busco in form species name + . + chromosome name + : + gene start + - + gene end
                    species_name = header.split(".")[0]      
                    
                    # TODO: fpor some reason trimal removes headers
                    # chrm_name = header.split(":")[0].split(".")[1]  
                    # start_pos = header.split("-")[0].split(":")[1]

                    seq = str(record.seq)
                    GC3 = seq[2::3]
                    GC3_count = GC3.count('G') + GC3.count('g') + GC3.count('C') + GC3.count('c')
                    GC3_percent = ( GC3_count / len(GC3) ) * 100
                    
                    outf1.write(f"{BUSCO_name}\t{species_name}\t{GC3_percent}\t{str(len(GC3))}\n")

                    if GC3_percent > args.outlier:
                        outf2.write(f"{BUSCO_name}\t{species_name}\t{GC3_percent}\t{str(len(GC3))}\n")

# Also want to get comparison outlook to plot heatmap
# Therefore first get a list of outlier BUSCO 
def get_outlier_BUSCO_list(busco_dir):
    BUSCO_list = []
    with open(f"{busco_dir}outlierGC3_{args.outlier}.tsv") as f:
        for line in f:
            if line.startswith("BUSCO_id"):
                continue
            else:
                busco_id = line.split("\t")[0]
                BUSCO_list.append(busco_id)

    return set(BUSCO_list) #Get set of BUSCO outliers

# For each outlier BUSCO, find that BUSCO for each species and extract with GC3 
def get_outlier_GC3_comparison(busco_dir, BUSCO_list):
    with open(f"{busco_dir}trimmedGC3_BUSCO.tsv") as f, open(f"{busco_dir}outlier_GC3-comparison.tsv", "w") as outf:
        outf.write(f"BUSCO_id\tspecies\tGC3_percent\tOutlier\n")
        for line in f:
            if line.startswith("Orthogroup"): #get all trimmed BUSCO GC3 
                continue
            else: 
                busco_id = line.split("\t")[0]
                #For every outlier BUSCO, bull out for every species 
                if busco_id in BUSCO_list:
                    sp = line.split("\t")[1]
                    GC3 = line.split("\t")[2]

        
                    if float(GC3) > 95:
                        outf.write(f"{busco_id}\t{sp}\t{GC3}\tOutlier_GC3\n")
                    elif float(GC3) > 85: 
                        outf.write(f"{busco_id}\t{sp}\t{GC3}\tHigh_GC3\n")
                    else:
                        outf.write(f"{busco_id}\t{sp}\t{GC3}\t\t\n")

def main():
    # run_mafft(args.busco)
    # run_trimal(args.busco)
    get_trimmed_GC3_output(args.busco)
    BUSCO_list = get_outlier_BUSCO_list(args.busco)
    get_outlier_GC3_comparison(args.busco, BUSCO_list)

if __name__ == "__main__":

    parse = argparse.ArgumentParser()

    parse.add_argument("-b", "--busco",type=str, help="path to BUSCO directory",required=True)
    parse.add_argument("-o", "--outlier",type=int, help="Outlier for GC/AT %",required=True)

    args = parse.parse_args()

    main()

#TODO add way to get start positions and maybe telomere analysis 
#!/usr/bin/env python3
import glob
import argparse
from joblib import Parallel, delayed
from subprocess import call as unix

# Author: Peter Mulhair
# Date: 07/06/2021
# Usage python3 busco_run.py --path ../genome/ --lineage /home/zoo/quee4075/GC_Insects/GC3/data/genome/busco_downloads/lineages/insecta_odb10

parse = argparse.ArgumentParser()

parse.add_argument("--path",type=str, help="path to genomes for BUSCO",required=True)
parse.add_argument("--lineage",type=str, help="path to /busco_downloads/lineages/[lineage] ",required=True)

args = parse.parse_args()

def busco_run(genome):
    sp_name = genome.split('/')[-1].split('.')[1]
    spName = sp_name.replace('1', '').replace('2', '').replace('3', '').replace('4', '').replace('5', '')
    spName = spName[2:]
    print(spName)
    unix(f"busco -i {genome} -l {args.lineage} -o {spName} -m genome -c 10 --offline --out_path busco_outdirs", shell=True)

def create_genome_list(genomes):
    genome_list = []
    for genome in glob.glob(genomes + "*.fasta"):
        genome_list.append(genome)
    return genome_list

def main():
    genome_list = create_genome_list(args.path)
    Parallel(n_jobs=1)(delayed(busco_run)(fasta) for fasta in genome_list)

if __name__ == "__main__":
    main()

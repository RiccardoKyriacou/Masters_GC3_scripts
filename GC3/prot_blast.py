import os
import re
import glob
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from subprocess import call as unix
from joblib import Parallel, delayed

'''
Script to blast all GC3 outlier genes 
specfied in the get_GC3_primary.py 
output
'''

parse = argparse.ArgumentParser()

parse.add_argument("-i", "--input",type=str, help="name of input tsv file to get gene IDs (oultier gene file)",required=True)
parse.add_argument("-c", "--cds",type=str, help="path to directory with cds fasta files",required=True)
parse.add_argument("-o", "--output",type=str, help="name of directory to output files",required=True)
parse.add_argument("-b", "--blastdb",type=str, help="name of directory with blastdb files",required=True)

args = parse.parse_args()

def sense_check_args(): #Check to make sure command line arguments are correctly input
    os.makedirs(args.output, exist_ok=True)
    if args.output[-1] != '/':
        args.output += '/'

    if args.cds[-1] != '/':
        args.cds += '/'

    args.blastdb = os.path.abspath(args.blastdb)

    if args.blastdb[-1] != '/':
        args.blastdb += '/'

#Parse input tsv file
def parse_input(fn):
    sp_genes_dict = defaultdict(list)
    with open(fn) as f:
        for line in f:
            lines = line.split('\t') 
            species = lines[0]
            geneID = lines[2].split("_")[0]
            sp_genes_dict[species].append(geneID) #Create dict of species : (list of gene_ids)
    
    return sp_genes_dict

#Create fasta files of query genes
def create_fasta(sp_genes_dict):
    for sp, genes in sp_genes_dict.items():
        geneID_list = {}
        for gene in genes: #create dictionary of gene IDs from sp_genes_dict
            geneID = gene.split('_')[0]
            geneID_list[geneID] = gene #dictionary of gene name: Gene_ID
        cds_fasta = glob.glob(args.cds + sp + '*')
        cds_fasta = cds_fasta[0]
        with open(args.output + sp + '_query_prot.fasta', 'w') as outf, open(cds_fasta) as f: #write output fatsas 
            for record in SeqIO.parse(f, 'fasta'):
                ID = record.id
                if ID in geneID_list.keys():
                    gene_info = geneID_list[ID]
                    seq = str(record.seq)
                    nuc_seq = Seq(seq)
                    prot_seq = nuc_seq.translate()
                    outf.write(f">{gene_info}\n{str(prot_seq)[:-1]}\n")    

#Create prot_fasta_list
def get_fasta_list():
    prot_fasta_list = []
    for fasta in glob.glob(args.output + '*fasta'):
        if 'query_prot' in fasta:
            prot_fasta_list.append(fasta)
    
    return prot_fasta_list

#Blast command line functionality 
def prot_blastnr(fasta):
    sp_fas = fasta.split('.fas')[0]
    sp_name = fasta.split("/")[-1].split("_query")[0]
    print(f"Performing BLAST for {sp_name}") #pblast_GC3_95outliers/Bombus_pratorum_query_prot.fasta
    unix('blastp -query ' + fasta + ' -db ' + args.blastdb + 'prot_db/insect_pep -out ' + sp_fas + '_blastoutput.tsv -outfmt "6 qseqid sseqid stitle evalue pident bitscore qstart qend qlen sstart send slen" -evalue 1e-30 -max_target_seqs 5', shell=True)

#Get dictionary of protein ID : sequence 
def get_prot_id():
    prot_ID_seq_dict = {} 
    with open(args.blastdb + 'prot_db/insect_pep.fasta') as f:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.id
            seq = str(record.seq)
            prot_ID_seq_dict[ID] = seq

    return prot_ID_seq_dict

#get dictionary of species : dictionary 
def get_prot_sp():
    sp_dict_prot = {} 
    for fasta in glob.glob(args.blastdb + "prot_db/*fa"):
        species_name = fasta.split("/")[-1].split(".")[0]
        with open(fasta) as Fprot:
            for line in Fprot:
                if line.startswith(">"):
                    gene_id = line.split(">")[1].split(" ")[0] #Get hit_id
                    sp_dict_prot[gene_id] = species_name
                else: 
                    continue
    
    return sp_dict_prot

#Output tsv of all protein hits 
def get_output_tsv(prot_ID_seq_dict, sp_dict_prot):
    for tsv in glob.glob(args.output + '*blastoutput.tsv'):
        sp = tsv.split('/')[-1].split('_query')[0]
        searchtype = tsv.split('/')[-1].split('_')[3]
        sp_fasta = tsv.split('_blastout')[0]
        if searchtype == 'prot':
            with open(tsv) as f, open(sp_fasta + '.fasta') as f1:
                sp_seq_info = {}
                for record in SeqIO.parse(f1, 'fasta'):
                    header = record.description
                    sp_seq = str(record.seq)
                    sp_seq_info[header] = sp_seq
                    
                prot_hits = defaultdict(list)
                for line in f:
                    lines = line.split('\t')
                    sp_query = lines[0]
                    subject_header = lines[2]
                    prot_hits[sp_query].append(subject_header) 
                    
                for prot_query, hits_info in prot_hits.items():
                    with open(args.output + sp + '_' + prot_query.split('.')[0] + '_blast_prot.fasta', 'w') as outf:
                        sp_prot_seq = sp_seq_info[prot_query]
                        outf.write(f">{prot_query}\n{sp_prot_seq}\n")
                        for hit_info in hits_info:
                            hit_ID = hit_info.split(' ')[0]
                            transcript_id = hit_info.split("transcript:")[1].split(" ")[0]
                            sp_name = sp_dict_prot[hit_ID]
                            hit_ID_seq = prot_ID_seq_dict[hit_ID]
                            print(f">{sp_name}_{transcript_id}\n{hit_ID_seq}\n")
                            outf.write(f">{sp_name}_{transcript_id}\n{hit_ID_seq}\n")
        else:
            continue 

#At this point we have done a protein BLAST to retrieve protein hits 
#However we now want to output the nucleotide sequence for these genes in order to count GC3

#Get a dictionary of cds_IDS and cds sequences from --cds
def get_cds_dict():
    cds_dict = {}
    for cds in glob.glob(args.cds + "*fa"):
        with open(cds) as f:
            for record in SeqIO.parse(f, "fasta"):
                cds_ID = record.id
                cds_seq = str(record.seq)
                cds_dict[cds_ID] = cds_seq 

    return  cds_dict

#Dictionary of gene id : sp name for cds balstdb
def get_cds_sp_dict():
    sp_dict_cds = {} 
    for fasta in glob.glob(args.blastdb + "cds_db/*fa"):
        species_name = fasta.split("/")[-1].split(".")[0]
        with open(fasta) as Fcds:
            for line in Fcds:
                if line.startswith(">"):
                    gene_id = line.split(">")[1].split(" ")[0] #Get hit_id
                    sp_dict_cds[gene_id] = species_name
                else: 
                    continue    
    return sp_dict_cds

#Dictionary of ID : seq
def get_nuc_ID():
    nuc_ID_seq_dict = {}
    with open(args.blastdb + 'cds_db/insect_nuc.fasta') as f:
        for record in SeqIO.parse(f, 'fasta'):
            ID = record.id
            seq = str(record.seq)
            nuc_ID_seq_dict[ID] = seq 

    return nuc_ID_seq_dict 
     
#Main function to output all aligned cds files 
def get_aligned_cds(cds_dict, sp_dict_cds, nuc_ID_seq_dict ):        
    for tsv in glob.glob(args.output + "*tsv"):
        sp_name = tsv.split("/")[-1].split("_query")[0]
        with open(tsv) as f:
            cds_hits = defaultdict(list)
            for line in f:
                lines = line.split('\t')
                gene_query = lines[0]
                transcript_header = lines[2].split("transcript:")[1].split(" ")[0]
                cds_hits[gene_query].append(transcript_header)     

            for gene_query, transcript_header in cds_hits.items():
                with open(args.output + sp_name + '_' + gene_query.split('.')[0] + '_blast_cds.fasta', 'w') as outf:
                    cds_seq = cds_dict[gene_query.split("_")[0]]
                    outf.write('>' + gene_query + '\n' + cds_seq + '\n')
                    for hit in transcript_header:
                        hit_sp_name = sp_dict_cds[hit]
                        transcript_id = hit
                        hit_cds_seq = nuc_ID_seq_dict[hit]
                        outf.write(f">{hit_sp_name}_{transcript_id}\n{hit_cds_seq}\n")

#Creating combined TSV file   
def create_combined_output():
    unix("cat " + args.output + " *blastoutput.tsv >> " + args.output + "combined_outlier_genes.tsv", shell=True)

#Extracting only gene_id and Gene name from combined tsv 
def parse_combined_output():
    with open("combined_outlier_genes.tsv", "r") as f, open(args.output + "all_gene_hits.tsv", "w") as outf:
        outf.write("Gene_ID" + "\t" + "Putative_hits" + "\t" + "E-value" + "\t" + "Identities" + "\n")
        for line in f:
            gene_id = line.split("\t")[0].split("_")[0]
            E_value = line.split("\t")[3]
            Identities = line.split("\t")[4]
            gene_drome = re.findall(".+description:(.+?)\t", line) #Extracting genes with headers containing "description: [gene name]"
            if len(gene_drome) > 0:
                gene_hits = "\t".join(gene_drome)
                outf.write(f"{gene_id}\t{gene_hits}\t{E_value}\t{Identities}\n")     
            gene_cds = re.findall(".+cds (.+?)\t", line) #Extracting genes with headers containing "cds [gene name]"
            if len(gene_cds) > 0:
                gene_cds_hits = "\t".join(gene_cds)
                outf.write(f"{gene_id}\t{gene_cds_hits}\t{E_value}\t{Identities}\n")
            gene_protein = re.findall("protein (.+?)\t", line) #Extracting genes with headers containing "cds [gene name]"
            if len(gene_protein) > 0:
                gene_protein_hits = "\t".join(gene_protein)
                outf.write(f"{gene_id}\t{gene_protein_hits}\t{E_value}\t{Identities}\n")

def main():
    #First step is to create the alignment fasta files
    print(f"Creating alignment fasta files")
    sp_genes_dict = parse_input(args.input)   
    create_fasta(sp_genes_dict)  
    
    #Blast these alignements 
    print(f"Performing BLAST")
    prot_fasta_list = get_fasta_list()  
    Parallel(n_jobs=60)(delayed(prot_blastnr)(fas) for fas in prot_fasta_list)     
    
    #Outputting protein blast results as a series of TSV files
    print(f"Outputting BLAST TSV")
    prot_ID_seq_dict = get_prot_id()
    sp_dict_prot = get_prot_sp()
    get_output_tsv(prot_ID_seq_dict, sp_dict_prot)
    
    #Get nucleotide sequences from BLAST outputs and output to TSV as aligned fastas
    print(f"Aligning nucleotide")
    cds_dict = get_cds_dict()
    sp_dict_cds = get_cds_sp_dict
    nuc_ID_seq_dict = get_nuc_ID
    get_aligned_cds(cds_dict, sp_dict_cds, nuc_ID_seq_dict )
    
    #Create combined output file of blast hits 
    print(f"Creating combined outlier file")
    create_combined_output()
    
    #Extract only the gene IDs and names of hits from this file to make more readable output
    parse_combined_output() 

if __name__ == "__main__":
    sense_check_args()
    main()

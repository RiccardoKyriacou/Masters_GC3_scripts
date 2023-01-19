import argparse

'''
Plot the GC3 comparison 
file according to order 
of the phylogenetic tree
'''

parse = argparse.ArgumentParser()

parse.add_argument("-s", "--species",type=str, help="Path to ordered list of species",required=True)
parse.add_argument("-t", "--tsv",type=str, help="Path to outlier_GC3-comparison.tsv file from orthofinder script",required=True)

args = parse.parse_args()

def parse_list(file):
    sp_lst = []
    with open(file) as f:
        for line in f:
            sp = line.rstrip("\n")
            sp_lst.append(sp)
    return sp_lst

def order_tsv(tsv, sp_lst):
    with open("ordered_GC3-comparison.tsv", "w") as outF:
        outF.write(f"Orthogroup\tSpecies\tgene_ID\tchrm_no\tGC3\tOutlier\tPosition\n")
        for sp in sp_lst:
            with open(tsv) as f:
                for line in f:
                    if line.startswith("Orthogroup"):
                        continue
                    else:
                        Orthogroup = line.split("\t")[0]	
                        Species = line.split("\t")[1]	
                        gene_ID = line.split("\t")[2]	
                        chrm_no = line.split("\t")[3]	
                        GC3 = line.split("\t")[4]	
                        Outlier = line.split("\t")[5]
                        Position = line.split("\t")[6]

                        if Species == sp:
                            outF.write(f"{Orthogroup}\t{Species}\t{gene_ID}\t{chrm_no}\t{GC3}\t{Outlier}\t{Position}")

def order_by_species():
    sp_list = parse_list(args.species)
    order_tsv(args.tsv, sp_list)

if __name__ == "__main__":
    order_by_species()

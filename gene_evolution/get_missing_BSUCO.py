import glob
import argparse

'''
Script to get output of missing genes for each species
and output of which genes are missing for each species 

Usage: python3 get_missing_BSUCO.py --busco ../results/busco_outdirs/ --links ../results/buscoIDs_to_genes.txt
'''

# Get dictionary of BUSCO IDs : gene name
def get_buscoID_name_dict(fn):
    buscoID_name_dict ={}
    with open(fn) as f:
        for line in f:
            buscoID = line.split("\t")[0]
            gene_name = line.split("\t")[1]
            
            buscoID_name_dict[buscoID] = gene_name
    
    return buscoID_name_dict

# Get a list of the BUSCO missing genes for each species
def get_missing_buscos(tsv, buscoID_name_dict):
    gene_names = [] 
    for line in tsv:
        if line.startswith("#"):
            continue
        else:
            busco_id = line.rstrip()
            # Search dictionary to output gene name when given BUSCO ID
            gene_name = buscoID_name_dict[busco_id]
            # Can also output len(gene_name) in output to get a count for missing genes 
            gene_names.append(gene_name)
  
    return gene_names

# Main function 
def path_to_missing_busco_list(busco_outdirs, links):
    # Parse links file to get dictionary
    buscoID_name_dict = get_buscoID_name_dict(links)
    with open("missing_buscos.tsv", "w") as outf:
        for species_name in glob.glob(busco_outdirs + "*"):
            # Get sp_name from file path
            sp_name = species_name.split("/")[-1]
            # Path to the output busco file 
            missing_busco_tsv = f"{species_name}/run_lepidoptera_odb10/missing_busco_list.tsv"
            with open(missing_busco_tsv) as tsv:
                gene_names = get_missing_buscos(tsv, buscoID_name_dict)
                # Output file in form of sp_name \t missing_gene_count \t missing_gene_name
                outf.write(f"{sp_name}\t{len(gene_names)}\t{', '.join(gene_names)}\n")

if __name__ == "__main__":
    parse = argparse.ArgumentParser()

    parse.add_argument("-b", "--busco",type=str, help="path to busco_outdirs", required = True)
    parse.add_argument("-l", "--links",type=str, help="path to file for busco IDs to gene name", required = True)

    args = parse.parse_args()

    if args.busco[-1] != '/':
        args.busco += '/'

    path_to_missing_busco_list(args.busco, args.links)
import argparse
import json
import glob
from collections import defaultdict 

"""
Script to parse JSON output from HyPhy relax and output to 
more readible file

Usage: 

"""

# def parse_species_kvalues(relax_dir):
#     OG_sp_k = {}
#     for output in glob.glob(f"{relax_dir}*.RELAX.json"):
#         with open(output) as f:
#             relax_output = json.load(f)

#             branch_attributes = relax_output["branch attributes"]
#             species_tests = branch_attributes["0"]

#             for species, relax_stats in species_tests.items():
#                 print(species)
#                 for stat_name, value in relax_stats.items():
#                     if ("k (general descriptive)") in stat_name:
#                         k = stat_name
#                         k_value = value

def get_highGC3_species(relax_dir):
    buscoID_sp_dict = defaultdict(list)
    for sp_lst in glob.glob(f"{relax_dir}*_species_list.txt"):
        buscoID = sp_lst.split("/")[-1].split("_")[0]
        with open(sp_lst) as f:
            for line in f:
                sp_name = line.strip()
                buscoID_sp_dict[buscoID].append(sp_name)

    return buscoID_sp_dict
        

def parse_test_results(relax_dir, buscoID_sp_dict):
    with open("HyPhy_RELAX_result.tsv", "w") as outf:
        outf.write(f"BUSCO_name\tk_value\tp_value\tresult\tforeground_species\n")
        for output in glob.glob(f"{relax_dir}*.RELAX.json"):
            BUSCO_name = output.split("/")[-1].split(".")[0]
            with open(output) as f:
                #Parse and read .json output
                relax_output = json.load(f)
                test_attributes = relax_output["test results"]

                #Search busco dict to get species 
                high_GC3_species = buscoID_sp_dict[BUSCO_name]

                #Loop through test results
                for stat_name, value in test_attributes.items():
                    #Get p-value
                    if ("p-value") in stat_name:
                        p_value = value
                    #Get k
                    if ("relaxation or intensification parameter") in stat_name:
                        k_value = value
                    #Compare P-value to k to see if selection intensified or relaxed 
                        result = ""        
                        if (k_value > 1) and (p_value < 0.05):
                            result += "INTENSIFIED"
                        if (k_value < 1) and (p_value < 0.05):
                            result += "RELAXED"
                        else:
                            result += ""
                        
                        outf.write(f"{BUSCO_name}\t{k_value}\t{p_value}\t{len(high_GC3_species)}\t{result}\t{high_GC3_species}\n")

                    
if __name__ == "__main__":
    parse = argparse.ArgumentParser()

    parse.add_argument("-r", "--relax",type=str, help="path to RELAX output directory", required = True)

    args = parse.parse_args()

    if args.relax[-1] != '/':
        args.relax += '/'

    buscoID_sp_dict = get_highGC3_species(args.relax)
    parse_test_results(args.relax, buscoID_sp_dict)



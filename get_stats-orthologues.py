import argparse
from scipy import stats
from collections import defaultdict
from scipy.stats import fisher_exact

'''
Script to calcualte Fisher'sexact 
test and Student's T for genes 
within 1MB window of telomeres
'''

parse = argparse.ArgumentParser()

parse.add_argument("-t", "--tsv",type=str, help="Path to trimmedGC3_orthogroups file from pblast script",required=True)

args = parse.parse_args()

sp_list = []
with open(args.tsv) as f: #Create list of all species in tsv file
    for line in f:
        if line.startswith("Orthogroup"): #Skip past header line 
            continue
        else:
            species = line.split("\t")[1]
            sp_list.append(species)

sp_list = set(sp_list) #Make list of each species appearing only once 

sp_chrms = defaultdict(list)
for sp in sp_list: #Create dictionary of speceis and all chromosomes 
    with open(args.tsv) as f:
        for line in f:
            if line.startswith("Orthogroup"): 
                continue
            else:
                species = line.split("\t")[1]
                chrm_no = line.split("\t")[5]
                sp_chrms[species].append(chrm_no)

def fisher_exact_sliding_window(telomere_highGC3, telomere_lowGC3, nontelomere_highGC3, nontelomere_lowGC3):
    output = ""
    oddsratio, pvalue = fisher_exact([[telomere_highGC3, telomere_lowGC3], [nontelomere_highGC3, nontelomere_lowGC3]]) 
    Pvalue = str(pvalue)
    if pvalue < 0.05:
        output =  Pvalue + "\t" + "Result is significant" 
        return output
    else:
        output = Pvalue + "\t" + "" 
        return output

def t_test(GC3_telomere, GC3_nontelomere):
    statistic, pvalue = stats.ttest_ind(GC3_telomere, GC3_nontelomere)
    mean_GC3_telomere = sum(GC3_telomere) / len(GC3_telomere) 
    mean_GC3_nontelomere = sum(GC3_nontelomere) / len(GC3_nontelomere) 
    if pvalue < 0.05:
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t" + "Significant"
    else:
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t"                        

with open("Fischers-test_allspecies.txt", "w") as OutF, open("T-test_allspecies.txt", "w") as OutF1:    
    OutF.write("Species" + "\t" + "Chromosome" + "\t" + "Chromosome_length" + "\t" + "P-value" + "\t" + "Significance" + "\t" + "High_GC3_genes" + "\n")
    OutF1.write("Species" + "\t" + "Chromosome" + "\t" + "Chromosome_length" + "\t" + "Mean_GC3_telomere" + "\t" + "Mean_GC3_non-telomere" + "\t" + "P-value" + "\t" + "Significance" + "\n")
    for sp, chrms in sp_chrms.items():
        chrm_size_dict = {} #Make a normal dict of chromosome number to size
        chrm_position_GC3 = defaultdict(dict)#Create a nested dictionary where chromosome : position_GC3 dictionary
        for chrm in set(chrms):
            position_GC3 = {}
            with open(args.tsv) as f:
                for line in f:
                    if line.startswith("Orthogroup"): 
                        continue
                    else:
                        species = line.split("\t")[1]
                        gene_id = line.split("\t")[2]
                        GC3 = float(line.split("\t")[3])
                        chrm_no = line.split("\t")[5]
                        position = float(line.split("\t")[6])
                        
                        if (chrm == chrm_no) and (sp == species):
                            chrm_size = float(line.split("\t")[7]) #Chromosome size needs to be defined when your defined chrm matches the chromosome number in the tsv file, not before it
                            chrm_size_dict[chrm] = chrm_size #Add chromosome number and size to new dictionary defined above
                            position_GC3[position] = GC3   

            chrm_position_GC3[chrm] = position_GC3 #Add chromosome number as key and position_GC3 as value
        
        for chromos, position_GC3_dict in chrm_position_GC3.items():#Loop through the nested dictionary i.e. for chromosome and all of its gene position and GC3s 
            telomere_highGC3 = 0
            telomere_lowGC3 = 0
            nontelomere_highGC3 = 0
            nontelomere_lowGC3 = 0
            GC3_telomere = []
            GC3_nontelomere = []
            chrm_len = chrm_size_dict[chromos]#get the length of the chromosome from new chrm_size_dict dictionary
            for position, GC3 in position_GC3_dict.items():
                window = 1_000_000 
                if position < window: #if gene is in first window then count genes above 85% GC3
                    GC3_telomere.append(GC3)
                    if GC3 >=85:
                        telomere_highGC3 += 1
                    else:
                        telomere_lowGC3 += 1
                if position > (chrm_len - window): # if gene is in last window then count genes above 85% GC3
                    GC3_telomere.append(GC3)
                    if GC3 >= 85:
                        telomere_highGC3 += 1 
                    else:
                        telomere_lowGC3 += 1
                if (position > window) and (position < (chrm_len - window)): #Genes between the first and last window
                    GC3_nontelomere.append(GC3)
                    if GC3 >= 85:
                        nontelomere_highGC3 += 1
                    else:
                        nontelomere_lowGC3 += 1 

            print("Calculating Fisher's exact and T test for " + sp + " " + "Chromosome: " + chromos + "...")
            if telomere_highGC3 > 0: 
                result = fisher_exact_sliding_window(telomere_highGC3, telomere_lowGC3, nontelomere_highGC3, nontelomere_lowGC3)
                OutF.write(sp + "\t" + chromos + "\t" + str(chrm_len) + "\t" + result  + "\t" + str(nontelomere_highGC3 + telomere_highGC3) + "\n" )
            else:
                OutF.write(sp + "\t" + chromos + "\t" + str(chrm_len) + "\t" + "" +  "\t" + "" + "\t" + str(nontelomere_highGC3 + telomere_highGC3) + "\n")

            if (len(GC3_telomere) > 0) and (len(GC3_nontelomere) > 0 ):
                ttest = t_test(GC3_telomere, GC3_nontelomere) 
                OutF1.write(sp + "\t" + chromos + "\t" + str(chrm_len) + "\t" + ttest + "\n" )
            else:
                OutF1.write(sp + "\t" + chromos + "\t" + str(chrm_len) + "\t" + "\n" )
               

            
                            



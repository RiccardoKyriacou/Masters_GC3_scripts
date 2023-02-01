import argparse
from scipy import stats
from collections import defaultdict
from scipy.stats import fisher_exact

'''
Script to calcualte Fisher's exact 
test and Student's T for genes 
within 1MB window of telomeres
'''

parse = argparse.ArgumentParser()

parse.add_argument("-t", "--tsv",type=str, help="Path to trimmedGC3_orthogroups file from orthofinder script",required=True)

args = parse.parse_args()

#Create list of all species in tsv file
sp_list = []
with open(args.tsv) as f: 
    for line in f:
        if line.startswith("Orthogroup"): #Skip past header line 
            continue
        else:
            species = line.split("\t")[1]
            sp_list.append(species)

sp_list = set(sp_list) #Make list of each species appearing only once 

#Create dictionary of speceis and all chromosomes
sp_chrms = defaultdict(list)
for sp in sp_list: 
    with open(args.tsv) as f:
        for line in f:
            if line.startswith("Orthogroup"): 
                continue
            else:
                species = line.split("\t")[1]
                chrm_no = line.split("\t")[5]
                sp_chrms[species].append(chrm_no)

def fisher_exact_sliding_window(telomere_highGC3, telomere_lowGC3, nontelomere_highGC3, nontelomere_lowGC3):
    telomereput = ""
    oddsratio, pvalue = fisher_exact([[telomere_highGC3, telomere_lowGC3], [nontelomere_highGC3, nontelomere_lowGC3]]) 
    Pvalue = str(pvalue)
    if pvalue < 0.05:
        telomereput =  Pvalue + "\t" + "Result is significant" 
        return telomereput
    else:
        telomereput = Pvalue + "\t" + "" 
        return telomereput

def t_test(GC3_telomere, GC3_nontelomere):
    statistic, pvalue = stats.ttest_ind(GC3_telomere, GC3_nontelomere)
    mean_GC3_telomere = sum(GC3_telomere) / len(GC3_telomere) 
    mean_GC3_nontelomere = sum(GC3_nontelomere) / len(GC3_nontelomere) 
    if (pvalue < 0.05) and (pvalue > 0.01) :
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t" + "*"
    elif (pvalue < 0.01) and (pvalue > 0.001):
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t" + "**"
    elif pvalue < 0.001:
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t" + "***"
    else:
        return str(mean_GC3_telomere) + "\t" +str(mean_GC3_nontelomere) + "\t" + str(pvalue) + "\t" 
    # If a p-value is less than 0.05, it is flagged with one star (*) 
    # If a p-value is less than 0.01, it is flagged with 2 stars (**) 
    # If a p-value is less than 0.001, it is flagged with three stars (***).

def sliding_window(position, chrm_size ):
    telomere = ""
    window = 1_000_000 
    if position < window: #if gene is in first window
        telomere += "Near_Telomere"
    if position > (chrm_size - window): # if gene is in last window then count genes above 85% GC3
        telomere += "Near_Telomere"  
    if (position > window) and (position < (chrm_size - window)): #Genes between the first and last window
       telomere += "Non_telomere"   
    return telomere 

#main
with open("Fischers-test_allspecies.tsv", "w") as telomeref, open("T-test_genome.tsv", "w") as telomeref1:    
    telomeref.write(f"Species\tChromosome\tChromosome_length\tP-value\tSignificance\tHigh_GC3_genes\n")
    telomeref1.write(f"Species\tMean_GC3_telomere\tMean_GC3_non-telomere\tP-value\tSignificance\n")
    for sp, chrms in sp_chrms.items():
        chrm_size_dict = {} #Make a normal dict of chromosome number to size
        chrm_position_GC3 = defaultdict(dict) #Create a nested dictionary where {chromosome : {position: GC3} } 
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
        genome_GC3_total = []
        GC3_total_telomere = []
        GC3_total_nontelomere = []      

        for chromos, position_GC3_dict in chrm_position_GC3.items():#Loop through the nested dictionary i.e. for chromosome and all of its gene position and GC3s 
            telomere_highGC3 = 0
            telomere_lowGC3 = 0
            nontelomere_highGC3 = 0
            nontelomere_lowGC3 = 0
            GC3_telomere = []
            GC3_nontelomere = []
            chrm_len = chrm_size_dict[chromos] #get the length of the chromosome from new chrm_size_dict dictionary
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
        
            #Perform Fisher's exact test for each chromosome for each species 
            if telomere_highGC3 > 0: 
                result = fisher_exact_sliding_window(telomere_highGC3, telomere_lowGC3, nontelomere_highGC3, nontelomere_lowGC3)
                telomeref.write(f"{sp}\t{chromos}\t{str(chrm_len)}\t{result}\t{str(nontelomere_highGC3 + telomere_highGC3)}\n")
            else:
                telomeref.write(f"{sp}\t{chromos}\t{str(chrm_len)}\t\t\t{str(nontelomere_highGC3 + telomere_highGC3)}\n")
       
            #Calculate total GC3 for each chromosome
            GC3_total = GC3_telomere + GC3_nontelomere 

            #Create list of total, telomere and non-telomere GC3 for all chromsomes for each pecies 
            genome_GC3_total.append(GC3_total) 
            GC3_total_telomere.append(GC3_telomere)
            GC3_total_nontelomere.append(GC3_nontelomere)

        #Sum total GC3 for each species 
        genome_GC3_total = sum(genome_GC3_total, [])  
        GC3_total_telomere = sum(GC3_total_telomere, [])
        GC3_total_nontelomere = sum(GC3_total_nontelomere, [])
        
        #Checking for normally distributed data to perform T-test
        shapiro_test = stats.shapiro(genome_GC3_total)
        pvalue = shapiro_test.pvalue
        if pvalue <  0.05: #If genome is normally distributed...
            print(f"{sp} genome GC3 normally distributed. Performing T-test...")
            if (len(GC3_total_telomere) > 0) and (len(GC3_total_nontelomere) > 0 ): #Perform T test...
                telomeref1.write(f"{sp}\t{t_test(GC3_total_telomere, GC3_total_nontelomere)}\n")
            else:     
                telomeref1.write(f"{sp}\n")                              
        else:
            print(f"{sp} Genome GC3 NOT normally distributed")

#Getting boxplot output:

#SLiding window function
def sliding_window(position, chrm_size):
    telomere = ""
    window = 1_000_000 
    if position < window: #if gene is in first window
        telomere += "Telomere"
    elif position > (chrm_size - window): # if gene is in last window then count genes above 85% GC3
        telomere += "Telomere"  
    elif (position > window) and (position < (chrm_size - window)): #Genes between the first and last window
       telomere += "Non_telomere"   
    return telomere 

#Parsing tsv  
def create_boxplot_data(file): 
    with open(file) as f, open("boxplot_data.tsv", "w") as outf2:    
        outf2.write(f"Orthogroup\tspecies\tgene_ID\tGC3_percent\tchrm_size\ttelomere\n")
        for line in f:
            if line.startswith("Orthogroup"): 
                continue
            else:
                orthogroup = line.split("\t")[0]
                species = line.split("\t")[1]
                gene_id = line.split("\t")[2]
                GC3 = float(line.split("\t")[3])
                position = float(line.split("\t")[6])
                chrm_size = int(line.split("\t")[7])
            
                telomere = sliding_window(position, chrm_size)
                outf2.write(f"{orthogroup}\t{species}\t{gene_id}\t{GC3}\t{chrm_size}\t{telomere}\n")

create_boxplot_data(args.tsv)


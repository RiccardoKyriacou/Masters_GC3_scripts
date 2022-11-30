# Useful_scripts
Scripts to count GC3

Current script pipeline:
* get_GC3_priamry.py 
  * (calculate GC3 content per chromosome)
* prot_blast.py 
  * (BLAST all outliers to find candidate GOIs)
* get_trimmed_allhits.py 
  * (align and trim GOIs)
* orthofinder_GC3.py 
  * (retrieve GC3 and chromosome information for orthofinder orthogroup)
* get_stats-orthologues.py 
  * (run F and T test on data) 
* get_singlecopy_iqtree.py
  * (Create an species tree for all speices used) 
* order_by_species-tree.py
  * (Order comparison tsv file according to the produced species tree)

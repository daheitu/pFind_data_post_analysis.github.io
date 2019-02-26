# coding = utf8
import re

f = open("./DSSO_CID_MS2_CID_MS3_T1_CIDFT_rep_pair.txt", 'r')
tab = f.readlines()

for line in tab:
    line_list = line.strip().split("\t")
    if len(line_list) == 5:
        prec_mz = line_list[2]
        chrg = line_list[1]
        pairs = line_list[-1].split(";")
        if len(pairs) == 1:
            sites = line_list[-1][2:-2].split("), (")
            site1_list = sites[0].split(",")
            print(site1_list)
            

        
    
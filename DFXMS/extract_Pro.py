# -*- coding: utf-8 -*-
"""
Created on Sun Nov  4 21:47:36 2018

@author: Yong
"""

import os
import re
from bs4 import BeautifulSoup
import urllib3
import time
import random


def get_SubLocation(Uniprot_ID):
    url = "https://www.uniprot.org/uniprot/" + Uniprot_ID + ".txt"
    urllib3.disable_warnings()
    # url = "https://www.uniprot.org/uniprot/P0A953.txt"

    http = urllib3.PoolManager()
    respon = http.request('GET', url)
    soup = BeautifulSoup(respon.data, "html5lib")
    soup_txt = str(soup).split("\n")
    sub_loc = ""
    for line in soup_txt:
        # print(line)
        if line[:29] == "CC   -!- SUBCELLULAR LOCATION":
            sub_loc = line[30:-1] 
            print(sub_loc)
            break
        else:
            continue
    return sub_loc


sleep_time = random.uniform(1, 3.5)
f = open("./selected_protein_result.txt", 'r').readlines()
print(f[1].split("\t"))
b = open("extrated_loc.txt", 'w')
for line in f[1:]:
    line_list = line.split("\t")
    psm_count = line_list[4]
    uni_id = line_list[1]
    if int(psm_count) > 100:
        time.sleep(sleep_time)
        sub_loc = get_SubLocation(uni_id)
        b.write("\t".join([uni_id, sub_loc]) + "\n")
    else:
        break

b.close()

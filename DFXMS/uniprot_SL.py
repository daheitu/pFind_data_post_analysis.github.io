# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 11:15:28 2018

@author: Yong
"""


from bs4 import BeautifulSoup
import urllib3
import time

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

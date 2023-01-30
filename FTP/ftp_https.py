# coding = utf-8
import requests
from bs4 import BeautifulSoup
import time
import re

url = "https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/04/PXD030983/"

b = open("http_19483.txt", 'w')

def getLinkfromHTML(url):
    r = requests.get(url)
    print(r.status_code)
    demo = r.text
    soup = BeautifulSoup(demo, 'html.parser')
    # print(soup.prettify())
    text = ""
    for link in soup.find_all('a'):
        raw_name = link.string
        # print(raw_name)
        if ".raw" in raw_name:
            print(url + raw_name)
            b.write(url + raw_name + "\n")
    # print(text)
    
    # return "None sentence","None sentence", "None sentence", "None sentence"

print(getLinkfromHTML(url))
b.close()
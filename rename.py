import os

# /nibs/home/ycao/Collaboration/TC/liuzhu_5diub/A/BS3/quant
# /nibs/home/ycao/Collaboration/TC/YJ/quant/YJ_2QUANT/BEST
"""
for i in range(4,7):
    path = "/nibs/home/ycao/Collaboration/TC/YJ/quant/YJ_" + str(i) + "QUANT/"
    os.chdir(path)
    os.system("rename pQuant.proteins pQuant.proteins_YJ" + str(i) + " pQuant.proteins.list")
"""

os.system("cd /nibs/home/ycao/Collaboration/TC/liuzhu_5diub")
os.system("mkdir Quant")
os.system("cd Quant")
os.chdir("/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/Quant")
for s in ["B","C","D","E"]:
    os.system("mkdir "+s)
    os.system("cd "+s)
    os.chdir("/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/Quant/"+s)    
    for linker in ["BS3", "BS2G", "DST", "EGS"]:
        os.system("mkdir "+ linker)
        os.system("cd "+linker)
        os.chdir("/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/Quant/"+s+"/"+ linker)
        path = "/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/" + s + "/" + linker + "/quant"
        os.system("cp -r "+ path +" ./")
        os.system("cd ..")
        os.chdir("/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/Quant/"+s)
    os.system("cd ..")
    os.chdir("/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/Quant/")    
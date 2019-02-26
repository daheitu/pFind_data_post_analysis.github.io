import os
os.chdir(r"E:\MASS SEPTRA DATA\jianhua\DOPA3_16DOT1")
f = open("pFind.protein").readlines()
b=open("filter.txt",'w')
c=open("freq.txt",'w')
mass_list=[]
mass_dic={}
for i in range(3,len(f)):
    #print(f[i])
    modi= f[i].strip().split("\t")[6]
    freq = int(f[i].strip().split("\t")[-1])
    if "PFIND_DELTA_" in modi:
        b.write(f[i])
        modi_mass = modi.split(",")[1]
        mass = round(float(modi_mass[modi_mass.find("PFIND_DELTA_")+12:modi_mass.find(";")]),2)
        modi_residue = f[i].strip().split("\t")[1][int(modi.split(",")[0])-1]
        if mass > 300 and modi_residue == "K":
            if mass not in mass_dic:
                mass_dic[mass] = freq
            else:
                mass_dic[mass] = mass_dic[mass] + freq
"""          
for i in mass_list:
    if mass_list.count(i)>0:
        mass_dic[i]=mass_list.count(i)
"""
print(mass_dic)
for key in mass_dic:
    if mass_dic[key] >1:
        c.write("\t".join([str(key),str(mass_dic[key])]))
        c.write("\n")
b.close()
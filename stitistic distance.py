import os

os.chdir(r"C:\Users\Yong\Desktop\DISTANCE")
f = open("list.txt").readlines()
b = open("dis_sditri.txt", "w")
dis_dic = {
    "<4": 0,
    "4-8": 0,
    "8-12": 0,
    "12-16": 0,
    "16-20": 0,
    "20-24": 0,
    "24-28": 0,
    "28-32": 0,
    "32-36": 0,
    "36-40": 0,
    ">40": 0
}
for line in f:
    num = float(line.strip())
    if num < 4:
        dis_dic["<4"] += 1
    elif num < 8:
        dis_dic["4-8"] += 1
    elif num < 12:
        dis_dic["8-12"] += 1
    elif num < 16:
        dis_dic["12-16"] += 1
    elif num < 20:
        dis_dic["16-20"] += 1
    elif num < 24:
        dis_dic["20-24"] += 1
    elif num < 28:
        dis_dic["24-28"] += 1
    elif num < 32:
        dis_dic["28-32"] += 1
    elif num < 36:
        dis_dic["32-36"] += 1
    elif num < 40:
        dis_dic["36-40"] += 1
    else:
        dis_dic[">40"] += 1

print(dis_dic)
for key in dis_dic:
    b.write("\t".join([key, str(dis_dic[key])]))
    b.write("\n")

b.close()
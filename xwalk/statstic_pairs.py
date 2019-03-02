import os
os.chdir(r'E:\Script\xwalk_output\KR')
cut_off = 32.2
# f = open('3zrj_RR.txt', 'r').readlines()
b = open('KR_stac.csv', 'w')
num_dic = {}
total_num = 0
for file in os.listdir(os.getcwd()):
    if file[-6:] == "RK.txt":
        f = open(file, 'r').readlines()
        num = 0
        for line in f:
            line_list = line.strip().split("\t")
            sdad_dis = line_list[6]
            if float(sdad_dis) < cut_off:
                num += 1
                total_num += 1
        num_dic[file] = num

for key in num_dic:
    b.write(",".join([key, str(num_dic[key])])+"\n")

b.close()


print(total_num)
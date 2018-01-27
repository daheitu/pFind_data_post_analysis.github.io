import os
import statistics
os.chdir(r"D:\E\Collabaration\TC\lIUZHU\A")
file_list = os.listdir(os.getcwd())
for file in file_list:
    if file[-5:] == ".list":
        new_file = file[:-5] + ".txt"
        f = open(file, 'r').readlines()
        n15_ratio_list = []
        r = open(new_file, 'w')
        for line in f:
            line_list = line.strip().split("\t")
            if ";" in line_list[-2]:
                n15 = line_list[-2][:-1].split(';')
                for ratio_sigma in n15:
                    n15_ratio_list.append(float(ratio_sigma.split(',')[0]))
        print(file)
        mdi = statistics.median(n15_ratio_list)
        for line in f:
            line_list = line.strip().split("\t")
            if line_list[3][0].isdigit():
                ratio = str(round(float(line_list[3]) / mdi, 2)) + "; " + str(
                    round(float(line_list[10]) / mdi, 2)) + "; " + str(
                        round(float(line_list[17]) / mdi, 2))
                sigma = line_list[4] + "; " + line_list[11] + "; " + line_list[18]
                line_simp = [
                    line_list[0], line_list[2], line_list[8], ratio, sigma
                ]
                r.write("\t".join(line_simp))
                r.write("\n")
        r.close()
        
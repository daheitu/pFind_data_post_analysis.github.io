import os



def merger_dis_fl(fl_name):
    b = open(fl_name + "merge.csv", 'w')
    for root, dirs, fls in os.walk(r"G:\msData\20200419"):
        # print(root)
        for fl in fls:
            if fl == fl_name:
                print(os.path.join(root, fl))
                f = open(os.path.join(root, fl)).readlines()
                for line in f:
                    lineList = line.rstrip().split("\t")
                    site, dist = lineList[:2]
                    if dist != "no stru":
                        b.write(",".join([site, dist])+"\n")
    b.close()
# print(os.getcwd())

for name in ["BSMEGoutput.txt", "DSSoutput.txt", "DSSOoutput_score.txt"]:
    merger_dis_fl(name)
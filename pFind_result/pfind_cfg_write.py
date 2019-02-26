import os

path = os.popen("pwd").readlines()[0].rstrip("\n") + "/"
f = os.popen("ls *.pf2").readlines()
total_file = len(f)
para = open("/nibs/home/ycao/scripts/pFind_Win.cfg", 'r').readlines()
new_para = open("pFind_Win.cfg", 'w')
i = 0 
while i < len(para):
    if para[i].strip() == "[datalist]":
        new_para.write(para[i])
        new_para.write("msmsnum=" + str(total_file)+"\n")
        n = 1
        for line in f:
            file_name = line
            write_list = ["msmspath", str(n), "=", path, file_name]
            new_para.write("".join(write_list))
            n += 1
    else:
        new_para.write(para[i])

    i += 1
new_para.close()



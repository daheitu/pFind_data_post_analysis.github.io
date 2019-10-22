import os
os.chdir(r"C:\Users\Yong\Downloads")


def remove_blank_digtal(str):
    seq = ""
    for word in list(str):
        if word in ["a", "g", "u", "c"]:
            seq += word
    return seq


f = open("miRNA.dat", 'r').readlines()
begin_list = []
end_list = []
for i in range(len(f)):
    if f[i][:2] == "ID":
        begin_list.append(i)
    elif f[i][:2] == "//":
        end_list.append(i)
    else:
        pass

print(len(begin_list), len(end_list))
b = open("miRNA.txt", 'w')
for i in range(len(begin_list)):
    seq =""
    for j in range(begin_list[i], end_list[i]):
        if f[j][:2] == "ID":
            ID_name = f[j].rstrip("\n").split("   ")[1]
            # print(ID_name)
        elif f[j][:13] == "SQ   Sequence":
            for n in range(j+1,end_list[i]):
                seq += f[n].strip()
            seq = remove_blank_digtal(seq)

    position_list = []
    for m in range(begin_list[i], end_list[i]):
        if f[m][:10] == "FT   miRNA":
            position_list.append(f[m].rstrip("\n").split("           ")[1])

    wtite_list = [ID_name, seq]
    for i in range(len(position_list)):
        up_site = int(position_list[i].split("..")[0])
        down_site = int(position_list[i].split("..")[1])
        wtite_list.append(
            "[" + position_list[i] + "] " + seq[up_site - 1:down_site - 1])

    b.write("\t".join(wtite_list))
    b.write("\n")

b.close()

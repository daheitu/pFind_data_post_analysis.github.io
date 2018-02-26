import os
os.chdir(r"C:\Users\Yong\Desktop\sgc\saga")


def get_linked_site_inform(linked_site):
    m = linked_site.find("(")
    n = linked_site.find(")")
    p = linked_site.find("-")
    x = linked_site.find("(", p)
    y = linked_site.find(")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:x].strip()
    position1 = linked_site[m + 1:n]
    position2 = linked_site[x + 1:y]
    return [protein1, protein2, position1, position2]


f = open("list.txt", 'r').readlines()
b = open("report.txt", 'w')
for line in f:
    a = get_linked_site_inform(line.strip())

    b.write("\t".join(a))
    b.write("\n")

b.close()

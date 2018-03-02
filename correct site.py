import os
os.chdir(
    r"C:\Users\Yong\Documents\pLink\search_task_2018.01.23.18.20.37_SAGA_DSS\reports"
)


def site_correct(link_pair):
    m = link_pair.find("(")
    n = link_pair.find(")")
    p = link_pair.find("-")
    x = link_pair.find("(", p)
    y = link_pair.find(")", p)
    protein1 = link_pair[:m].strip()
    protein2 = link_pair[p + 1:x].strip()
    position1 = int(link_pair[m + 1:n])
    position2 = int(link_pair[x + 1:y])
    if position1 > position2:
        site = link_pair[p + 1:] + "-" + link_pair[:p]
    else:
        site = link_pair
    return site


f = open("linksite.txt", 'r').readlines()

tmp = open("pLink2 site.txt", 'w')
for line in f:
    tmp.write(site_correct(line.strip()))
    tmp.write("\n")

tmp.close()

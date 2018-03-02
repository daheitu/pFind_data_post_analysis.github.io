import re

s = "ATG12(3)-AT256(123)"

print(re.match("\((\d*)\)-", s))
print(re.match(".*-.*\((\d*)\)$", s))

print(re.findall("\((\d*)\)", s)[0])

# pattern = \(\d+\)

def get_linked_site_inform(linked_site):
    pos_list = re.findall("\((\d*)\)", linked_site)
    position1 = pos_list[0]
    position2 = pos_list[1]
    p = linked_site.find("-")
    m = linked_site.find("(" + position1 + ")-")
    n = linked_site.find("(" + position2 + ")", p)
    protein1 = linked_site[:m].strip()
    protein2 = linked_site[p + 1:n].strip()
    return protein1, protein2, position1, position2


print(get_linked_site_inform(s))

f = open(r"D:\FastaDatabase\mg1655small.fasta").readlines()
b = open(r"D:\FastaDatabase\mg1655small_2.fasta", 'w')
for line in f:
    if not line.startswith("+") and not line.startswith("-"):
        b.write(line)
b.close()
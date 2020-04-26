import os

# os.chdir(r"D:\fasta database")
InputFasta = "SAGA.fasta"


def formatFast(fastaPath):
    f = open(fastaPath, 'r').readlines()
    Formated_fast = os.path.basename(fastaPath).replace(".fasta", "_format.fasta")
    b = open(Formated_fast, "w")
    for line in f:
        if line[0] == ">":
            b.write(line)
        else:
            pro_length = len(line.strip())
            n, p = pro_length.__divmod__(60)
            if p == 0:
                for i in range(n):
                    b.write(line[60 * i:60 * (i + 1)])
                    b.write("\n")
            else:
                for i in range(n):
                    b.write(line[60 * i:60 * (i + 1)])
                    b.write("\n")
                b.write(line[-p - 1:])
    b.close()
                

#if __name__ == "__main__":
    #formatFast(InputFasta)

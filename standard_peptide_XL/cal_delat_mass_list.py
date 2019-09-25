# coding = utf-8

input_com = "[I|P|Y]-[A|W]	[T|Q|G]-[N|F]	5.00543"


def loadAAmassTable():
    mpMassTable={'A':71.037114,'R':156.101111,'N':114.042927,'D':115.026943,'C':103.009185, \
'E':129.042593,'Q':128.058578,'G':57.021464,'H':137.058912,'I':113.084064, \
'L':113.084064,'K':128.094963,'M':131.040485,'F':147.068414,'P':97.052764, \
'S':87.032028,'T':101.047679,'U':150.95363,'W':186.079313,'Y':163.06332,'V':99.068414, \
'H2O':18.01056,'Proton':1.0072766}
    return mpMassTable


def loadHydrophoicData():
    aaHydroDic = {}
    with open(
            r"F:\OneDrive\github\pFind_data_post_analysis.github.io\standard_peptide_XL\AA_hyfrophbic.ini"
    ) as file:
        f = file.readlines()
    for line in f:
        aa = line[-2]
        hydroVaule = float(line[:-3].split(":")[1].strip())
        aaHydroDic[aa] = hydroVaule
    return aaHydroDic


def seqPreTreat(seq):
    return seq[1:-1].split("|")


def getSeq(unit):
    seq1, seq2 = unit.split("-")
    return seqPreTreat(seq1), seqPreTreat(seq2)


def calhydro(seq):
    aaHydroDic = loadHydrophoicData()
    hydroVaule = 0
    for aa in seq:
        hydroVaule += aaHydroDic[aa]
    return hydroVaule


mpMassTable = loadAAmassTable()
inList = input_com.split("\t")
list1, list2 = getSeq(inList[0])
list3, list4 = getSeq(inList[1])

seqMassDic = {}
for aa1 in list1:
    for aa2 in list2:
        for aa3 in list3:
            for aa4 in list4:
                seq = aa1 + aa2 + aa3 + aa4
                seqMass = 0
                for aa in seq:
                    seqMass += mpMassTable[aa]
                seqMassDic[seq] = seqMass
                #print(seq, seqMass)
b = sorted(seqMassDic.items(), key=lambda x: x[1])
print(b)
for i in range(len(b) - 1):
    dMass = b[i + 1][1] - b[i][1]
    if dMass < 7:
        print(b[i + 1], b[i], (calhydro(b[i + 1][0]) - calhydro(b[i][0])))

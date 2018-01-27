import os
sample_list = ['A', 'B', 'C', 'D', 'E']

crossliner = ['BS3', 'BS2G', 'DST', 'EGS']

for sample in sample_list:
    for linker in crossliner:
        path = "/nibs/home/ycao/Collaboration/TC/liuzhu_5diub/" + sample + '/' + linker + '/quant'
        os.system('cd ' + path)
        os.chdir(path)
        os.system("mkdir Mono")
        os.system("mv *_"+linker+"* ./Mono")
        template = open("pQuant.cfg", 'r').readlines()
        tmp = open("pQuant_best.cfg", "w")
        for line in template:
            if line[:18] == "TYPE_PEPTIDE_RATIO":
                tmp.write(line[:19] + "1;")
                tmp.write("\n")
            elif line[:22] == 'NUMBER_SCANS_HALF_CMTG':
                tmp.write(line[:23] + "180;")
                tmp.write("\n")
            else:
                tmp.write(line)
        tmp.close()

        os.system(
            "java -jar /nibs/home/yhding/program/pQuant20160120/pQuant_64bit.jar ./pQuant_best.cfg"
        )
        os.system('rename pQuant.proteins pQuant.proteins_' + linker +
                  ' pQuant.proteins*')
        os.system('rename pQuant.spectra pQuant.spectra_' + linker +
                  ' pQuant.spectra*')
        os.system("mkdir Best")
        os.system("mv *_"+linker+"* ./Best")

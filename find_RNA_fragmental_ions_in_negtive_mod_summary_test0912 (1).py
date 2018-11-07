import os

os.chdir(r"C:\Users\Yong\Desktop\ZMQ")

precosur_mass_tol = 500
considered_ion_type = []
fregment_ion_tol = 50
match_rate_cutoff = 0.01
ion_types = ['A-B', 'A', 'B', 'C', 'D', 'D-H2O', 'W', 'X', 'Y', 'Z']

residues_mass_dic = {
    "A": 329.0525182,
    "G": 345.0474332,
    "C": 305.0412852,
    "U": 306.0253016,
    "P": 409.115115,
    "Q": 425.11003,
    "R": 385.103882,
    "T": 386.0878984
}

nutral_loss_dic = {
    "A": 135.05452,
    "G": 151.04942,
    "C": 111.04232,
    "U": 112.02732,
    "P": 409.115115,
    "Q": 425.11003,
    "R": 385.103882,
    "T": 386.0878984
}

# P:alkynyl_A,   Q:alkynyl_G     R:alkynyl_C,   T:alkynyl_U


def generate_ion_mass_range(num, tol):
    deta = num * tol / 1000000
    return [num - deta, num + deta]


def seq_mass_cal(seq):
    mass = 18.0105642 - 79.96633
    for resi in residues_mass_dic:
        mass += list(seq).count(resi) * residues_mass_dic[resi]

    return round(mass, 4)


def seq_reader(seq_fasta_path):
    seq_dict = dict()  #value: description key: sequence
    with open(seq_fasta_path, 'r') as seqreader:
        contents = seqreader.readlines()
        now_seq = ''
        description = ''
        for line in contents:
            if '>' in line:
                if now_seq != '':
                    seq_dict[now_seq] = description
                now_seq = ''
                description = line.strip()[1:]
            else:
                now_seq = now_seq + line.strip()

        if now_seq != '':
            seq_dict[now_seq] = description
    print(seq_dict)
    return seq_dict


# print(seq_reader("ten_RNA.fasta"))


def find_condi_pep(precu_mass, peptide_dic):
    condi_pep_list = []
    for pep in peptide_dic:
        pep_mass = seq_mass_cal(pep)
        [minm, maxm] = generate_ion_mass_range(pep_mass, precosur_mass_tol)
        if precu_mass > minm and precu_mass < maxm:
            condi_pep_list.append(pep)
        else:
            continue

    return condi_pep_list


def generate_series_ion_mass(seq):
    ion_dic = {}
    seq_mass = seq_mass_cal(seq)
    for ion_type in ion_types:
        ion_dic[ion_type] = []
    i = 0
    ion_dic['A'].append(residues_mass_dic[seq[i]] - 78.95850 - 1.00782)
    ion_dic['A-B'].append(ion_dic['A'][i] - nutral_loss_dic[seq[i]])
    ion_dic['W'].append(seq_mass - ion_dic['A'][i])
    ion_dic['B'].append(ion_dic['A'][i] + 15.99491)
    ion_dic['X'].append(seq_mass - ion_dic['B'][i])
    ion_dic['C'].append(ion_dic['B'][i] + 62.96359 + 1.00782)
    ion_dic['Y'].append(seq_mass - ion_dic['C'][i])
    ion_dic['D'].append(ion_dic['C'][i] + 18.01056)
    ion_dic['D-H2O'].append(ion_dic['D'][i] - 18.01056)
    ion_dic['Z'].append(seq_mass - ion_dic['D'][i])
    i += 1
    while i < len(seq) - 1:
        ion_dic['A'].append(ion_dic['A'][i - 1] + residues_mass_dic[seq[i]])
        ion_dic['A-B'].append(ion_dic['A'][i] - nutral_loss_dic[seq[i]])
        ion_dic['B'].append(ion_dic['A'][i] + 15.99491)
        ion_dic['C'].append(ion_dic['B'][i] + 62.96359 + 1.00782)
        ion_dic['D'].append(ion_dic['C'][i] + 18.01056)
        ion_dic['D-H2O'].append(ion_dic['D'][i] - 18.01056)
        ion_dic['W'].insert(0, seq_mass - ion_dic['A'][i])
        ion_dic['X'].insert(0, seq_mass - ion_dic['B'][i])
        ion_dic['Y'].insert(0, seq_mass - ion_dic['C'][i])
        ion_dic['Z'].insert(0, seq_mass - ion_dic['D'][i])
        i += 1
    return ion_dic


# print(seq_mass_cal("ACUGAGGUC"))
# print(generate_series_ion_mass("ACUG"))


def compare_therION_MS2(ms2_list, therION_list):
    count_list = [0] * len(therION_list)
    max_ion_list = []
    ion_range_dic = {}
    for ion in therION_list:
        ion_range_dic[ion] = generate_ion_mass_range(ion, fregment_ion_tol)
        max_ion_list.append(ion_range_dic[ion][1])

    ion_idx = 0
    ms2_idx = 0
    while ion_idx < len(therION_list) and ms2_idx < len(ms2_list):
        now_ion_max = max_ion_list[ion_idx]
        if now_ion_max < float(ms2_list[ms2_idx]):
            ion_idx += 1
        elif ion_range_dic[therION_list[ion_idx]][0] <= float(
                ms2_list[ms2_idx]):
            therION = therION_list[ion_idx]
            realION = ms2_list[ms2_idx]
            count_list[ion_idx] += 1
            ion_idx += 1
        else:
            ms2_idx += 1
    return count_list


def matchANDjudge_pep_ms2(pep, ms2_list, charge):
    ion_dic = generate_series_ion_mass(pep)
    match_type_dic = {}
    for ion_type in ion_dic:
        match_charge_dic = {}
        for i in range(1, charge + 1):
            ther_mass_list = []
            for ion in ion_dic[ion_type]:
                ther_mass_list.append((ion + 1.00782 * i) / i)
            match_list = compare_therION_MS2(ms2_list, ther_mass_list)
            match_charge_dic[i] = match_list
        final_match_list = []
        for m in range(len(ion_dic[ion_type])):
            type_total_match_list = []
            for i in range(1, charge + 1):
                type_total_match_list.append(match_charge_dic[i][m])
            final_match_list.append(sum(type_total_match_list))
        match_type_dic[ion_type] = final_match_list
    total_unmatch_num = 0
    total_ion_num = 0
    for ion_type in match_type_dic:
        total_unmatch_num += match_type_dic[ion_type].count(0)
        total_ion_num += len(match_type_dic[ion_type])

    match_rate = 1 - float(total_unmatch_num) / total_ion_num

    if match_rate > match_rate_cutoff:
        return match_rate, match_type_dic
    else:
        return match_rate, match_type_dic


file_list = os.listdir(os.getcwd())
# ms2_total_ion = {}  #=dict() #loaded ms2 mz&intesity
b = open("report.txt", 'w')
#loading process
for name in file_list:
    if name[-5:] == 'fasta':
        seq_path = name
        seq_dict = seq_reader(seq_path)
    elif name[-8:] != "opic.mgf":  #other files
        all = open(name, 'r').readlines()
        i = 0
        while i < len(all):
            if all[i].strip() == "BEGIN IONS":
                title = all[i + 1][5:all[i + 1].find(" File:")]
                # pepmz = float(all[i + 3].strip().split("=")[1])
                precusor_charge = int(all[i + 4].strip().split("=")[1][:-1])
                if " " in all[i + 3]:
                    pepmz = float(
                        all[i + 3].strip().split(" ")[0].split("=")[1])
                else:
                    pepmz = float(all[i + 3].strip().split("=")[1])
                pepMass = (pepmz + 1.00782) * precusor_charge
            else:
                print("wrong")

            condi_pep_list = find_condi_pep(pepMass, seq_dict)

            p = i + 5
            if condi_pep_list == []:
                while p < len(all):
                    if all[p].strip() == "END IONS":
                        break
                    else:
                        pass
                    p += 1
            else:
                # print(p)
                print(condi_pep_list)
                ms2_mass_list = []
                while p < len(all):
                    if all[p].strip() == "END IONS":
                        break
                    else:
                        mass = all[p].split(" ")[0]
                        itsty = all[p].split(" ")[1]
                        ms2_mass_list.append(mass)
                    p += 1

                for pep in condi_pep_list:
                    match_rate, match_type_dic = matchANDjudge_pep_ms2(
                        pep, ms2_mass_list, precusor_charge)
                    if match_rate > match_rate_cutoff:
                        output = "\t".join([title, pep, str(match_rate)])
                        b.write(output+ "\n")

                        for ion_type in match_type_dic:
                            write_line = [ion_type]
                            for num in match_type_dic[ion_type]:
                                write_line.append(str(num))
                            b.write("\t".join(write_line) + "\n")
                    else:
                        continue
            print(p)
            i = p + 1
b.close()


"""
        
        for i in range(len(all)):
            if all[i].strip() == "BEGIN IONS":
                ion_be.append(i)
            elif all[i].strip() == "END IONS":
                ion_end.append(i)
        # print(len(ion_be))
        if len(ion_be) != len(ion_end):
            print("erro")
        else:  #target files
            report_file_name = name[:name.find(".mgf")] + "_mathced_ions.mgf"
            report_file = open(report_file_name, 'w')
            for k in ion_be:
                title = all[k + 1].strip().split(' ')[0].split('=')[1]
                ms2_ion_list = []
                ms2_ion_list_intensity = []
                ms2_ion = {}
                #print k+5,ion_end[ion_be.index(k)]
                for m in range(k + 5, ion_end[ion_be.index(k)]):
                    ms2_ion[all[m].strip().split(" ")[0]] = all[
                        m].strip().split(" ")[1]
                #ms2_ion = dict(zip(ms2_ion_list, ms2_ion_list_intensity))

                ms2_total_ion[title] = ms2_ion
                #print(len(ms2_ion), len(ms2_total_ion))
    elif name[-5:] == 'fasta':
        seq_path = name
        seq_dict = seq_reader(seq_path)

seq_thero_mass_dic = {}
for seq in seq_dict:
    seq_thero_mass_dic[seq] = seq_mass_cal(seq)


for sequence in seq_dict:
    #generate theoretical ions
    ##sequence = 'CCAUGGGAG'    
    ion_types = ['A-B', 'A', 'B', 'C', 'D', 'D-H2O', 'W', 'X', 'Y', 'Z']
    the_ions = {}  # theoretical ions records by ion types

    #initialzation
    for ion_type in ion_types:
        the_ions[ion_type] = []

    sequence_mass_list = [residues_mass_dic[k] for k in sequence]

    charge_list = [-1, -2, -3, -4, -5, -6, -7]
    mz_charge = {}  # key:mz value:[charge, position]

    for charge in charge_list:  # different charge
        n = 1
        while n < len(sequence_mass_list):
            five_sequence_mass_list = sequence_mass_list[0:n]
            sum_five = sum(five_sequence_mass_list)
            sum_five_1 = sum(five_sequence_mass_list[0:n - 1])
            three_sequence_mass_list = sequence_mass_list[
                -(len(sequence_mass_list) - n):]
            sum_three = sum(three_sequence_mass_list)
            # 5'OH, 3'OH theoretical mass, no modification
            w = (sum_three + 17.0027396 + 1.0078246 + 1.0078246 * charge) / (
                -charge)  # w = sum_3+ OH+ H
            y = (sum_three + 17.0027396 + 1.0078246 - 79.9663326 +
                 1.0078246 * charge) / (-charge)  # y = w -HPO3
            x = (sum_three + 17.0027396 + 1.0078246 - 15.994915 +
                 1.0078246 * charge) / (-charge)  # x = w - O
            z = (sum_three + 17.0027396 + 1.0078246 - 95.9612476 +
                 1.0078246 * charge) / (-charge)  # z = w -HPO4
            d = (sum_five + 17.0027396 + 1.0078246 + 1.0078246 * charge) / (
                -charge)  # w = sum_five +OH+H
            d_h2O = (sum_five + 17.0027396 + 1.0078246 - 18.0105642 +
                     1.0078246 * charge) / (-charge)
            c = (
                sum_five + 17.0027396 + 1.0078246 - 15.994915 - 1.0078246 * 2 +
                1.0078246 * charge
            ) / (
                -charge
            )  #c = d-O-H-H Don't know why minus 2H, just for be aggred with RNAModMapper
            b = (sum_five + 17.0027396 + 1.0078246 - 79.9663326 +
                 1.0078246 * charge) / (-charge)  #b = d - HPO3
            a = (sum_five + 17.0027396 + 1.0078246 - 95.9612476 +
                 1.0078246 * charge) / (-charge)  #a = d - HPO4
            a_B = (sum_five_1 + 17.0027396 + 1.0078246 + 96.0204358 +
                   1.0078246 * charge) / (-charge)

            the_ions['A'].append(a)
            the_ions['A-B'].append(a_B)
            the_ions['B'].append(b)
            the_ions['C'].append(c)
            the_ions['D'].append(d)
            the_ions['D-H2O'].append(d_h2O)
            the_ions['W'].append(w)
            the_ions['X'].append(x)
            the_ions['Y'].append(y)
            the_ions['Z'].append(z)

            mz_charge[a] = [charge, n]
            mz_charge[a_B] = [charge, n]
            mz_charge[b] = [charge, n]
            mz_charge[c] = [charge, n]
            mz_charge[d] = [charge, n]
            mz_charge[d_h2O] = [charge, n]
            mz_charge[w] = [charge, n]
            mz_charge[x] = [charge, n]
            mz_charge[y] = [charge, n]
            mz_charge[z] = [charge, n]

            n = n + 1
    #print (the_ions)

    #sort to compare in ascending order
    for ion_type in ion_types:
        the_ions[ion_type].sort()

    # ion match
    for title in ms2_total_ion:

        #load ions of one experimental spectrum
        ms2_ion = ms2_total_ion[title]  #ions of one spectrum

        mz_convert = {}  #key: float number mz value: string number mz
        ms2_mz = []  #list of float number mz
        for mz in ms2_ion:
            ms2_mz.append(float(mz))
            mz_convert[float(mz)] = mz

        ms2_mz.sort()

        #matching process
        for ion_type in ion_types:

            the_ion = the_ions[ion_type]  #ions of one ion type
            the_ion_range = {}  #key: maximum mz of one ion value: [maximum, original, minimum] mz of one ion value
            the_ion_max = []  # list of maximum mz of ion values

            for ion in the_ion:
                [minmz, maxmz] = generate_ion_mass_range(ion)
                the_ion_range[maxmz] = [minmz, ion, maxmz]
                the_ion_max.append(maxmz)

            the_ion_idx = 0
            exp_ion_idx = 0

            #matching by using the merge sort technique
            while the_ion_idx < len(the_ion_max):

                if exp_ion_idx >= len(ms2_mz): break

                the_now_ion_max = the_ion_max[the_ion_idx]

                # exp ion > max the ion
                if the_now_ion_max < ms2_mz[exp_ion_idx]:
                    the_ion_idx += 1

                # exp ion in range of  the ion
                elif the_ion_range[the_now_ion_max][0] <= ms2_mz[exp_ion_idx]:
                    now_the_ion_mz = the_ion_range[the_now_ion_max][1]
                    now_exp_ion_mz = mz_convert[ms2_mz[exp_ion_idx]]
                    output='\t'.join([title, str(now_exp_ion_mz), ms2_ion[now_exp_ion_mz], \
                          str(mz_charge[now_the_ion_mz][0]), ion_type+ "_" +str(mz_charge[now_the_ion_mz][1]), \
                            seq_dict[sequence]])
                    report_file.write(output + '\n')
                    the_ion_idx += 1

                # exp ion < min the ion
                else:
                    exp_ion_idx += 1

report_file.close()
f.close()
"""
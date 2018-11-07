import os

os.chdir(r"G:\DSSO_1008\CID_BASED_Methods\pparse")


def generate_ion_mass_range(num, tol):
    deta = num * tol / 1000000
    return [num - deta, num + deta]


def judge_mass_diff(mz1, mz2, chrg, delta_ms, tol):
    m1 = mz1 * chrg
    m2 = mz2 * chrg
    if m1 < m2:
        pass
    else:
        m1, m2 = m2, m1
    min_m, max_m = generate_ion_mass_range(m1 + delta_ms, tol)
    if m2 > min_m and m2 < max_m:
        return True
    else:
        return False


file_list = os.listdir(os.getcwd())
for fl in file_list:
    if fl[-4:] == ".ms3":
        fl_name = fl[:-3]
        rep_name = fl_name + "extract.ms3"
        B = open(rep_name, 'w')
        f = open(fl, 'r').readlines()
        i = 0
        while i < len(f):
            if f[i][0].isdigit():
                i += 1
            else:
                if f[i][0] == 'S':
                    B.write(f[i])
                    i += 1
                elif f[i][:8] == "I	Filter":
                    B.write(f[i])
                    i += 1
                elif f[i][:15] == "I	PrecursorScan":
                    B.write(f[i])
                    i += 1
                elif f[i][:3] == "Z	0":
                    B.write(f[i])
                    i += 4
                elif f[i][0] == "Z":
                    B.write(f[i])
                    i += 1
                else:
                    i += 1
        B.close()
        p = open(rep_name, 'r').readlines()
        b = open("report_pair_sacn.txt", 'w')
        total_row = len(p)
        total_ms3_scan = total_row / 4
        i = 0
        while i < total_ms3_scan - 1:
            print(i)
            scan = p[4 * i + 0].split("\t")[1]
            precus_scan = p[4 * i + 2].strip().split("\t")[-1]
            precus_mz = p[4 * i + 0].strip().split("\t")[-1]
            chrg = p[4 * i + 3].split("\t")[1]

            scan_next = p[4 * (i + 1) + 0].split("\t")[1]
            precus_scan_next = p[4 * (i + 1) + 2].strip().split("\t")[-1]
            precus_mz_next = p[4 * (i + 1) + 0].strip().split("\t")[-1]
            chrg_next = p[4 * (i + 1) + 3].split("\t")[1]
            print(scan, scan_next)
            scan_bool = int(scan_next) == (int(scan) + 1)
            precus_scan_bool = precus_scan == precus_scan_next
            # if "NA" in [chrg, ]
            chrg_bool = chrg == chrg_next
            ms_diff_bool = judge_mass_diff(
                float(precus_mz), float(precus_mz_next), int(chrg), 31.9, 30)
            none_pair_num = 0
            if scan_bool and precus_scan_bool and ms_diff_bool:

                wt_list = [
                    scan, precus_scan, precus_mz, chrg, scan_next,
                    precus_scan_next, precus_mz_next, chrg_next
                ]
                b.write("\t".join(wt_list) + "\n")
                i += 2
            else:
                none_pair_num += 1
                i += 1
    else:
        continue

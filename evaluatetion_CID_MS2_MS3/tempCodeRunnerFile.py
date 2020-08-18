    delta_32_2plus = 0
            p = i + 1
            while p < len(f):
                if f[p].startswith("#masspairs"):
                    break
                else:
                    sub_list = f[p].strip().split("\t")
                    if sub_list[0] == "all_doublet":
                        charge = int(sub_list[4])
                        if charge > 1:
                            delta_32_2plus += 1
                    p += 1
            
            title_num_dic[title] = delta_32_2plus
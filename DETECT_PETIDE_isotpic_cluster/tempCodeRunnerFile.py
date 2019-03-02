    while k < len(ms2_ms_list):
            if float(ms2_ms_list[k]) > up_ms:
                return False, -1
            elif float(ms2_ms_list[k]) >= low_ms:
                return True, k
            else:
                k += 1
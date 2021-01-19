import numpy as np

mat = np.zeros([9,9], int)
mat[0,2] = 2
mat[0,3] = 5
mat[0,5] = 1
mat[0,8] = 6

mat[1,2] = 7
mat[1,7] = 4
mat[1,8] = 9

mat[2,0] = 3
mat[2,1] = 9
mat[2,3] = 2
mat[2,6] = 5


mat[3,1] = 6
mat[3,4] = 8
mat[3,6] = 9

mat[4,1] = 8
mat[4,5] = 9
mat[4,8] = 4

mat[5,2] = 9
mat[5,3] = 7
mat[5,4] = 1
mat[5,7] = 2

mat[6,0] = 6
mat[6,1] = 7
mat[6,6] = 8

mat[7,0] = 9
mat[7,4] = 5
mat[7,7] = 6

mat[8,0] = 1
mat[8,3] = 6
mat[8,4] = 7
mat[8,5] = 3
mat[8,6] = 4









def find_sub_max(m, n, mat):
    m0 = 3 * (m // 3) 
    n0 = 3 * (n // 3)
    return np.reshape(mat[m0: m0+3, n0: n0+3], (1,9))[0]


print(mat)
print(find_sub_max(2,5,mat))


def the_left_nums(raw_list, column_list, sub_list):
    raw_set = set(list(raw_list))
    column_set = set(list(column_list))
    sub_set = set(list(sub_list))
    nums = set(list(range(1,10)))
    return nums - (raw_set & column_set & sub_set)


def interpreter_mat(mat):
    for m in range(9):
        for n in range(9):
            if mat[m, n] == 0:
                raw_list = mat[m, ...]
                column_list = mat[..., n]
                sub_list = find_sub_max(m, n, mat)
                num_left_set = the_left_nums(raw_list, column_list, sub_list)
                if len(num_left_set) == 1:
                    mat[m, n] = list(num_left_set)[0]
                else:
                    
    return mat


def main(mat):
    while np.count_nonzero(mat) != 81:
        mat_new = interpreter_mat(mat)
        if (mat_new == mat).all():
            break
        else:
            mat = mat_new
        print(mat)


main(mat)
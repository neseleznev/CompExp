import copy

N = 4
M = 20 - N
A = 0.1 * M + 0.01 * N
B = 0.2 * M + 0.02 * N + 0.001*M*N + 0.001*N*N

x_star = [
[1],
[2],
[0.1*M + 0.01*N]]

set_of_equations = [
[1.2345, 3.1415, 1,           7.5175 + A],
[2.3456, 5.9690, 0,           14.2836],
[3.4567, 2.1828, 2 + 0.1 * N, 7.8223 + B]
]

ACCURACY_AFTER_DOT = 0

def round(value, accuracy):
    return int((value * (10 ** accuracy))) / (10 ** accuracy)

def get_indexes_of_no_zero_items(list):
    indexes = []
    for i in range(len(list)):
        if abs(round(list[i], ACCURACY_AFTER_DOT - 1)) > 0:
            indexes.append(i)
    return indexes
    
    
def print_matrix(matrix):
    s = [[str(e) for e in row] for row in matrix]
    lens = [max(map(len, col)) for col in zip(*s)]
    fmt = '\t'.join('{{:{}}}'.format(x) for x in lens)
    table = [fmt.format(*row) for row in s]
    print('\n'.join(table))
        
def matrixmult (A, B):
    rows_A = len(A)
    cols_A = len(A[0])
    rows_B = len(B)
    cols_B = len(B[0])

    if cols_A != rows_B:
      print("Cannot multiply the two matrices. Incorrect dimensions.")
      return

    # Create the result matrix
    # Dimensions would be rows_A x cols_B
    C = [[0 for row in range(cols_B)] for col in range(rows_A)]

    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += round(A[i][k] * B[k][j], ACCURACY_AFTER_DOT)
    return C

def calculate_cell_by_gauss(i, j, matrix_a, matrix_b, matrix_c):
    sum = 0
		
    for k in range(0, i):
        sum += round(matrix_b[i][k] * matrix_c[k][j], ACCURACY_AFTER_DOT)
    
    for k in range(i + 1, len(matrix_b[1])):
        sum += round(matrix_b[i][k] * matrix_c[k][j], ACCURACY_AFTER_DOT)
    
    if (j > i): # this element is in C
        matrix_c[i][j] = round((matrix_a[i][j] - sum)/matrix_b[i][i], ACCURACY_AFTER_DOT)
	
    else: # element in B			
        matrix_b[i][j] = round((matrix_a[i][j] - sum)/matrix_c[i][i], ACCURACY_AFTER_DOT)


def calculate_matrix_B_nad_C_by_gauss(matrix_a):
    matrix_b = [
    [0, 0, 0],
    [0, 0, 0],
    [0, 0, 0]
    ]

    matrix_c = [
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0]
    ]
    
    height = len(matrix_a)
    weedth = len(matrix_a[0])

    try:
        for k in range(0, weedth):
            for i in range(k, height):
                calculate_cell_by_gauss(i, k, matrix_a, matrix_b, matrix_c)
            for j in range(k + 1, weedth):
                calculate_cell_by_gauss(k, j, matrix_a, matrix_b, matrix_c)
    except ZeroDivisionError:
        raise Exception("Method can not find solution")
        
    return matrix_b, matrix_c
    
def solve_up_triangle_system(matrix, b_column):
    x_column = []
    
    for i in range(len(b_column)):
        x_column.append([0])
    
    for i in range(len(b_column))[::-1]:
        sum = 0
        for k in range(i):
            sum += round(matrix[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)
        
        for k in range(i+1, len(b_column)):
            sum += round(matrix[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)
        
        x_column[i][0] = round((b_column[i][0] - sum)/matrix[i][i], ACCURACY_AFTER_DOT)
        
    return x_column
        
    
def solve_any_transformed_system(set_of_equations):
    x_column = []
    for i in range(len(set_of_equations[0]) - 1):
        x_column.append([0])
        
    
    computed_elements_of_x = []
    num_of_computed_x_on_last_iteartion = 0
    no_zero_indexes_on_last_iteration = []
    
    while (len(computed_elements_of_x) != len(x_column)):
        
        for i in range(len(set_of_equations)):
            indexes_of_no_zero_items = get_indexes_of_no_zero_items(set_of_equations[i][:-1])
            
            new_no_zero_elements = list(set(indexes_of_no_zero_items) - set(no_zero_indexes_on_last_iteration))
            
            if (len(new_no_zero_elements) != 1):
                continue
            
            j = new_no_zero_elements[0]

            sum = 0
            for k in range(0, j):
                sum += round(set_of_equations[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)
            
            for k in range(j + 1, len(set_of_equations[0]) - 1):
                sum += round(set_of_equations[i][k] * x_column[k][0], ACCURACY_AFTER_DOT)
                
            x_column[j][0] = round((set_of_equations[i][-1] - sum)/set_of_equations[i][j], ACCURACY_AFTER_DOT)
            
            computed_elements_of_x.append(j)
            num_of_computed_x_on_last_iteartion += 1
            no_zero_indexes_on_last_iteration = new_no_zero_elements

            break
        
        if (num_of_computed_x_on_last_iteartion == len(computed_elements_of_x) - 1):
            raise Exception("Method can not find solution!!!!")
    return x_column
    
def compact_gauss_scheme(set_of_equations):
    matrix_b, matrix_c = calculate_matrix_B_nad_C_by_gauss(set_of_equations)
	
    y_column = []
    for i in range(0, len(matrix_c)):
        y_column.append([matrix_c[i][-1]])
        

    matrix_c = [line[:-1] for line in matrix_c]
    
    x_by_method = solve_up_triangle_system(matrix_c, y_column)
	
    print("Matrix B:")
    print_matrix(matrix_b)
    print()
    
    print("Matrix C:")
    print_matrix(matrix_c)
    print()
    
    print("column y:")
    print_matrix(y_column)
    print()
    
    return x_by_method

def find_max_cell(matrix, start_string_num):
    max_val = matrix[start_string_num][0]
    res_i = start_string_num
    res_j = 0
    for i in range(start_string_num, len(matrix)):
        for j in range(0, len(matrix[i])):
            if abs(matrix[i][j]) > max_val:
                
                max_val = matrix[i][j]
                res_i = i
                res_j = j
                
    return res_i, res_j
   
def swap_matrix_lines(matrix, string_num_1, string_num_2):
    temp = matrix[string_num_1]
    matrix[string_num_1] = matrix[string_num_2]
    matrix[string_num_2] = temp
   
   
def do_forward_step_of_gauss_by_main_element(matrix):
    set_of_equations = copy.deepcopy(matrix)
    
    for start_string_index in range(len(set_of_equations)):
            matrix_a = [line[:-1] for line in set_of_equations]
            
            main_i, main_j = find_max_cell(matrix_a, start_string_index)
            swap_matrix_lines(set_of_equations, start_string_index, main_i)
            
            for i in range(start_string_index + 1, len(set_of_equations)):
                multiplicator = round(-set_of_equations[i][main_j] / set_of_equations[start_string_index][main_j], ACCURACY_AFTER_DOT)

                for j in range(len(set_of_equations[i])):
                    set_of_equations[i][j] = round(set_of_equations[i][j] + set_of_equations[start_string_index][j] * multiplicator, ACCURACY_AFTER_DOT)
 
    return set_of_equations
   
def main_item_gauss(matrix):
        set_of_equations = do_forward_step_of_gauss_by_main_element(matrix)
        
        print("Set of equations after forward step of method:")
        print_matrix(set_of_equations)
        print()
        
        x_column = solve_any_transformed_system(set_of_equations)
        return x_column
        
        
 
    

    
def do_method_investigation(method):
    for accuracy in [2, 4, 6]:
        global ACCURACY_AFTER_DOT
        ACCURACY_AFTER_DOT = accuracy
        
        copy_of_system = []
        
        for i in range(len(set_of_equations)):
            copy_of_system.append([])
            for j in range(len(set_of_equations[i])):
                copy_of_system[i].append(round(set_of_equations[i][j], ACCURACY_AFTER_DOT))
        
        print("\naccuracy: ", accuracy)
        
        try:
            x_by_method = method(copy_of_system)
        except Exception as e:
            print(e)
            continue
        
        print("answer:")
        print_matrix(x_by_method)
        print()
        
	

    
def main():
    print("real answer(X*):")
    print_matrix(x_star)
    print("===============================================================")
    print()
    
    print("Compact Gauss Scheme: ")
    do_method_investigation(compact_gauss_scheme)	
    print("===============================================================")
    
    print("Gauss Scheme with main element: ")
    do_method_investigation(main_item_gauss)	

main()
#main_item_gauss(set_of_equations)







	

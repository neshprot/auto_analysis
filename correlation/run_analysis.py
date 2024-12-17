import os
from run_correlation import *
from top_mutations_seq import *
from make_xl import make_exel_file

file_name = 'log_within_6'
merged_dict = read_and_parse(f'../logouts/{file_name}')
write_to_file(merged_dict, f'{file_name}')


# top muts part
top_number = 40

top_num_first_param, max_first_param, min_first_param = get_top_num(merged_dict, 0, top_number)
top_num_second_param, max_second_param, min_second_param = get_top_num(merged_dict, 1, top_number)

max_desc = np.array([max_first_param, max_second_param])
min_desc = np.array([min_first_param, min_second_param])
top_fit = get_top_fitness_num(merged_dict, top_number, max_desc, min_desc, weight_one=0.5, weight_two=0.5)


top_first_dict = dict(top_num_first_param)
top_second_dict = dict(top_num_second_param)
top_fit_dict = dict(top_fit)

merged_dict = {**top_first_dict, **top_second_dict, **top_fit_dict}

top_fit_25_75 = get_top_fitness_num(merged_dict, top_number, max_desc, min_desc, weight_one=0.25, weight_two=0.75)
top_fit_75_25 = get_top_fitness_num(merged_dict, top_number, max_desc, min_desc, weight_one=0.75, weight_two=0.25)


histogram_multi_muts('his_1', merged_dict, 1)
histogram_multi_muts('his_2', merged_dict, 2)
histogram_multi_muts('his_3', merged_dict, 3)

# end of top param part

correlation_first_sites, correlation_second_sites = correlation_between_sites(f'out_{file_name}')
correlation_first_sites = sorted(correlation_first_sites, key=lambda x: x[2])
correlation_second_sites = sorted(correlation_second_sites, key=lambda x: x[2])

correlation_first, correlation_second, values, columns_per_mutation, unique_sites = correlation_between_muts(f'out_{file_name}')
correlation_first = sorted(correlation_first, key=lambda x: x[2])
correlation_second = sorted(correlation_second, key=lambda x: x[2])

first_desc = [x[0] for x in values]
second_desc = [x[1] for x in values]
model_coef_first = lasso(columns_per_mutation, first_desc, 0.1, unique_sites)
model_coef_second = lasso(columns_per_mutation, second_desc, 0.1, unique_sites)

common_elements_first, mut_to_cor_first = compare_correl_lasso(correlation_first, model_coef_first)
common_elements_second, mut_to_cor_second = compare_correl_lasso(correlation_second, model_coef_second)


increase_value_first, decrease_value_first = cut_to_dec_inc(common_elements_first, mut_to_cor_first)
increase_value_second, decrease_value_second = cut_to_dec_inc(common_elements_second, mut_to_cor_second)

if increase_value_first:
    print('Mutation increasing delta pKa')
    for mut in increase_value_first:
        print(mut[0])
if decrease_value_first:
    print('Mutation decreasing delta pKa')
    for mut in decrease_value_second:
        print(mut[0])
if increase_value_second:
    print('Mutation increasing dihedral angle')
    for mut in increase_value_second:
        print(mut[0])
if decrease_value_second:
    print('Mutation decreasing dihedral angle')
    for mut in decrease_value_second:
        print(mut[0])



make_exel_file([top_num_first_param, top_num_second_param, top_fit],
               file_name, [increase_value_first, decrease_value_first, increase_value_second, decrease_value_second])
os.remove('his_1.png')
os.remove('his_2.png')
os.remove('his_3.png')
os.remove(f'out_{file_name}')
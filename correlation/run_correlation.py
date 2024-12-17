import numpy as np
from scipy.stats import pearsonr
from scipy.sparse import csr_matrix
from sklearn.linear_model import Lasso
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.feature_selection import SelectFromModel
from top_mutations_seq import *

def lasso(columns_per_mutation, y, alpha, unique_sites):
    X = np.column_stack(columns_per_mutation)

    # Разбиение на обучающую и тестовую выборки
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    best_lasso_model = Lasso(alpha=alpha)
    best_lasso_model.fit(X_train, y_train)

    # Обучение модели LASSO с помощью SelectFromModel для оценки важности признаков
    lasso = Lasso(alpha=alpha) # alpha - настраиваемый параметр 0.1
    sfm = SelectFromModel(lasso)
    sfm.fit(X_train, y_train)

    # Выбор признаков на основе важности
    X_train_selected = sfm.transform(X_train)
    X_test_selected = sfm.transform(X_test)

    lasso_model = Lasso(alpha=alpha)
    lasso_model.fit(X_train_selected, y_train)

    # Предсказания
    y_pred = lasso_model.predict(X_test_selected)

    # Метрики
    mse = mean_squared_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)

    # Оценка важности признаков (с учетом отбора)
    model_coef = []
    for i, importance in enumerate(sfm.estimator_.coef_):
        if abs(importance) != 0:
            model_coef.append([unique_sites[i], importance])

    return model_coef

def read_and_parse(input_file):
    value_desc = None
    merged_dict = {}
    not_to_skip = False
    with open(input_file, 'r') as inf:
        for line in inf.readlines():
            words = line.split()
            if not_to_skip:
                new_line = line.split('\n')[0]
                if new_line == '':
                    new_line = 'wt'
                if new_line not in merged_dict:
                    merged_dict[new_line] = value_desc
                not_to_skip = False
            if 'descriptors' in words:
                value_desc = [float(w) for w in words[2:-4]]
                not_to_skip = True
    del merged_dict['wt']
    return merged_dict


def read_file_only_sites(input_file):
    mutation_sets = []
    values = []
    actually_muts = {}
    muts_values = {}
    with open(input_file, 'r') as inf:
        for line in inf.readlines():
            mutation = [int(i) for i in line.split('/')[1::2]]
            mut_name = " ".join([str(i) for i in mutation])
            value = [float(j) for j in line.split('/')[-1].split()[1:]]
            mut = []
            for k in line.split('/')[::2]:
                mut.extend(k.split())
            muts_values[" ".join(mut[1:-2:2])] = value
            if mutation in mutation_sets:
                idx = mutation_sets.index(mutation)
                values[idx][0] += value[0]
                values[idx][1] += value[1]
                values[idx][2] += 1
                actually_muts[mut_name] += [mut[1:-2:2]]
            else:
                mutation_sets.append(mutation)
                values.append(value + [1])
                actually_muts[mut_name] = [mut[1:-2:2]]
    values = [[x[0]/x[2], x[1]/x[2]] for x in values]
    return mutation_sets, values, actually_muts, muts_values


def read_file_sites_and_muts(input_file):
    mutation_sets = []
    values = []
    actually_muts = {}
    muts_values = {}
    with open(input_file, 'r') as inf:
        for line in inf.readlines():
            mutation = [i for i in line.split()[:-2:]]
            mut_name = " ".join([str(i) for i in mutation])
            value = [float(j) for j in line.split('/')[-1].split()[1:]]
            mut = []
            for k in line.split('/')[::2]:
                mut.extend(k.split())
            muts_values[" ".join(mut[1:-2:2])] = value
            if mutation in mutation_sets:
                idx = mutation_sets.index(mutation)
                values[idx][0] += value[0]
                values[idx][1] += value[1]
                values[idx][2] += 1
                actually_muts[mut_name] += [mut[1:-2:2]]
            else:
                mutation_sets.append(mutation)
                values.append(value + [1])
                actually_muts[mut_name] = [mut[1:-2:2]]
    values = [[x[0] / x[2], x[1] / x[2]] for x in values]
    return mutation_sets, values, actually_muts, muts_values

def numerical_to_onehot(numerical_arrays):
    """
    Converts multiple numerical arrays to one-hot encoded sparse arrays.

    Args:
        numerical_arrays: A list or tuple of NumPy arrays containing numerical values.

    Returns:
        A list of csr_matrix objects, representing the one-hot encoded arrays.
    """

    def combine_onehot_max(arrays):
        """
        Combines multiple one-hot arrays into a single array using maximum reduction.

        Args:
          arrays: A list of NumPy arrays (must be 1-dimensional).

        Returns:
          A NumPy array representing the combined one-hot encoding. Returns None if input is invalid.
        """

        return np.maximum.reduce(arrays)

    unique_values = set()
    for arr in numerical_arrays:
        unique_values.update(arr)

    value_to_index = {val: idx for idx, val in enumerate(sorted(list(unique_values)))}
    onehot_arrays = []
    for arr in numerical_arrays:
        row_indices = np.arange(len(arr))
        col_indices = np.array([value_to_index[val] for val in arr])
        data = np.ones(len(arr))  # All values are 1 in one-hot
        num_rows = len(arr)
        num_cols = len(value_to_index)
        sparse_onehot = csr_matrix((data, (row_indices, col_indices)), shape=(num_rows, num_cols))
        onehot_arrays.append(combine_onehot_max(sparse_onehot.toarray()))

    return onehot_arrays, sorted(list(unique_values))

def create_column_arrays(arrays):
    """
    Creates new arrays from the columns of multiple input arrays.

    Args:
        arrays: A list of NumPy arrays of the same length.

    Returns:
        A list of NumPy arrays, each representing a column from the input arrays. Returns None if input is invalid.
    """
    stacked_array = np.column_stack(arrays)
    return [stacked_array[i, :] for i in range(stacked_array.shape[0])]

def find_correlation(columns_per_mutation, values, unique_sites, alpha):
    correlation_first = []
    correlation_second = []
    first_desc = [x[0] for x in values]
    second_desc = [x[1] for x in values]
    for idx, mutation in enumerate(columns_per_mutation):
        correlation_f, p_value_f = pearsonr(mutation, first_desc)
        correlation_s, p_value_s = pearsonr(mutation, second_desc)
        if p_value_f < alpha:
            correlation_first.append([unique_sites[idx], correlation_f, p_value_f])
        if p_value_s < alpha:
            correlation_second.append([unique_sites[idx], correlation_s, p_value_s])
    return correlation_first, correlation_second


def correlation_between_sites(inp_file):
    mut_descs, values, actually_muts, muts_values = read_file_only_sites(inp_file)
    onehot_encoded, unique_sites = numerical_to_onehot(mut_descs)
    columns_per_mutation = create_column_arrays(onehot_encoded)

    # correlation between mutants and values
    correlation_first, correlation_second = find_correlation(columns_per_mutation, values, unique_sites, alpha=0.05)
    return correlation_first, correlation_second


def correlation_between_muts(inp_file):
    mut_descs, values, actually_muts, muts_values = read_file_sites_and_muts(inp_file)
    onehot_encoded, unique_sites = numerical_to_onehot(mut_descs)
    columns_per_mutation = create_column_arrays(onehot_encoded)
    correlation_first, correlation_second = find_correlation(columns_per_mutation, values, unique_sites, alpha=0.05)
    return correlation_first, correlation_second, values, columns_per_mutation, unique_sites


def write_to_file(merged_dict, file_name):
    with open(f'out_{file_name}', 'w') as ouf:
        for muts, values in merged_dict.items():
            ouf.write(f'{muts} {values[0]} {values[1]}\n')


def compare_correl_lasso(cor, las):
    dict1 = {tuple1[0]: tuple1 for tuple1 in cor}
    dict2 = {tuple2[0]: tuple2 for tuple2 in las}
    common_keys = set(dict1.keys()) & set(dict2.keys())
    common_elements = list(common_keys)
    return common_elements, dict1


def cut_to_dec_inc(common_elements, mut_to_cor):
    increase_value = []
    decrease_value = []
    for i in common_elements:
        values = mut_to_cor[i]
        if values[1] > 0:
            increase_value.append((values[0], values[1]))
        if values[1] < 0:
            decrease_value.append((values[0], values[1]))
    return increase_value, decrease_value

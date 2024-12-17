import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from itertools import combinations


def read_and_parse_out_file(filename):
    merged_dict = {}
    with open(filename, 'r') as inf:
        for line in inf.readlines():
            words = line.split()
            merged_dict[" ".join(words[:-2:])] = [float(x) for x in words[-2::]]
    return merged_dict

def get_top_num(muts_dic, param, number):
    sorted_muts = sorted(muts_dic.items(), key=lambda item: item[1][param], reverse=True)
    max_param = sorted_muts[0][1][param]
    min_param = sorted_muts[-1][1][param]
    top_num = sorted_muts[:number]
    return top_num, max_param, min_param


def get_top_fitness_num(muts_dic, number, max_descr, min_descr, weight_one, weight_two):
    def fitness(vals):
        def norm(vals):
            descriptors = []
            for idx, descr in enumerate(vals):
                if max_descr[idx] != min_descr[idx]:
                    descriptors.append((descr - min_descr[idx]) / (max_descr[idx] - min_descr[idx]))
                else:
                    descriptors.append(0)
            return descriptors

        def linear(vals):
            return np.dot(np.array(vals), np.array([weight_one, weight_two]))
        fit = linear(norm(vals))
        return fit

    sorted_muts = sorted(muts_dic.items(), key=lambda item: fitness(item[1]), reverse=True)

    top_num = sorted_muts[:number]
    return top_num


def plot_two_params(muts_val, color):
    x_values = [value[1][0] for value in muts_val]
    y_values = [value[1][1] for value in muts_val]
    # Создание графика
    plt.scatter(x_values, y_values, label='label', marker='o', linestyle='-', color=color)

def histogram_multi_muts(fig_name, muts_desc, num):
    sep_muts = []
    pairs_list = []

    for mutation in muts_desc:
        sep_muts.append(mutation)

    for string in sep_muts:
        words = string.split()
        for combo in combinations(words, num):
            pairs_list.append(tuple(sorted(combo)))

    muts_counts = Counter(pairs_list)

    sorted_muts = sorted(muts_counts.items(), key=lambda item: item[1], reverse=True)
    total = sum(item[1] for item in sorted_muts)
    labels = []
    values = []
    for mut in sorted_muts[:40]:
        labels.append(" ".join(mut[0]))
        values.append(mut[1] / total * 100)

    plt.figure(figsize=(7.2, 3.6))
    # Создание гистограммы
    plt.bar(labels, values)

    # Настройка отображения
    plt.xlabel("Мутации", fontsize=12, fontname='Calibri')  # Подпись оси OX
    plt.ylabel("Частота встречаемости(%)", fontsize=12, fontname='Calibri')  # Подпись оси OY
    plt.title("Гистограмма", fontsize=12, fontname='Calibri')  # Заголовок
    plt.xticks(rotation=45, ha="right", fontsize=9.6, fontname='Calibri')  # Поворот подписей для лучшей читаемости, если они длинные

    # Отображение гистограммы
    plt.tight_layout()  # Подгонка компоновки, чтобы подписи не перекрывались
    plt.savefig(fig_name)

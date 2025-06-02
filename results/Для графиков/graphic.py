import matplotlib.pyplot as plt
import numpy as np
import csv
import ast
import os
import glob

l1 = 0.0
l2 = 0.0
d1 = 0.0
d2 = 0.0
k=5

folder_path = "./Данные для графиков/Функция расслоения"
csv_files = glob.glob(os.path.join(folder_path, "*.csv"))

def create_graphic(values, name_file, l, d):
    x = np.linspace(-1, 1, 20)
    plt.plot(x, [c.real for c in values], marker='o', color='b', label = "Real")
    plt.plot(x, [c.imag for c in values], marker='o', color='r', label = "Imaginary")
    plt.title(f'l = {l}, d = {d}, k = {k}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(name_file)
    plt.close()

for file1 in csv_files:
    with open(file1, "r") as file:
        reader = csv.reader(file, delimiter=';')
        l1 = float(next(reader)[0])
        d1 = float(next(reader)[0])
        l2 = float(next(reader)[0])
        d2 = float(next(reader)[0])
        k  = float(next(reader)[0])
        data = next(reader)
    complex_data = [complex(*ast.literal_eval(value)) for value in data]
    first_20 = complex_data[:20]
    last_20 = complex_data[-20:]
    create_graphic(first_20, f"Графики/Функция расслоения/k{k}l{l1}d{d1}.png", l1, d1)
    create_graphic(last_20, f"Графики/Функция расслоения/k{k}l{l2}d{d2}.png", l2, d2)
import matplotlib.pyplot as plt
import numpy as np
import csv
import ast

l1 = 0.0
l2 = 0.0
d1 = 0.0
d2 = 0.0
k = 0.0

# Читаем данные из файла
with open("../данные/функция расслоения/", "r") as file: #указать название файла
    reader = csv.reader(file, delimiter=';')
    l1 = float(next(reader)[0])
    d1 = float(next(reader)[0])
    l2 = float(next(reader)[0])
    d2 = float(next(reader)[0])
    k = float(next(reader)[0])
    data = next(reader)

complex_data = [complex(*ast.literal_eval(value)) for value in data]

first_20 = complex_data[:20]
last_20 = complex_data[-20:]

def create_graphic(values, name_file, l, d):
    x = np.linspace(-1, 1, 20)
    plt.plot(x, [c.real for c in values], marker='o', color='b', label = "Real")
    plt.plot(x, [c.imag for c in values], marker='o', color='r', label = "Imaginary")
    plt.title(f'l = {l}, d = {d}, k = {k}')
    plt.legend()
    plt.tight_layout()
    plt.savefig(name_file)
    plt.close()

create_graphic(first_20, f"../графики/функция расслоения/k{k}l{l1}d{d1}.png", l1, d1)
create_graphic(last_20, f"../графики/функция расслоения/k{k}{l2}d{d2}.png", l2, d2)
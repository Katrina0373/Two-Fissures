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
k = 0.0
data = []

folder_path = "../данные"
csv_files = glob.glob(os.path.join(folder_path, "*.csv"))


def create_graphic(values, name_file):
    x = np.linspace(0.0, 10.0, 1000)
    plt.plot(x, [c.real for c in values], color='b', label = "Real")
    plt.plot(x, [c.imag for c in values], color='r', label = "Imaginary")
    plt.legend()
    plt.grid(True)
    plt.title(f"l1 = {l1}, d1 = {d1}, l2 = {l2}, d2 = {d2}")
    plt.tight_layout()
    plt.savefig(name_file)
    plt.close()


for file1 in csv_files:
# Читаем данные из файла
    with open(file1, "r") as file:
        reader = csv.reader(file, delimiter=';')
        l1 = float(next(reader)[0])
        d1 = float(next(reader)[0])
        l2 = float(next(reader)[0])
        d2 = float(next(reader)[0])
        k =  float(next(reader)[0])
        data = next(reader)

    complex_data = [complex(*ast.literal_eval(value)) for value in data]
    create_graphic(complex_data, f"../графики/классика/k{k}l1{l1}d1{d1}l2{l2}d2{d2}.png")

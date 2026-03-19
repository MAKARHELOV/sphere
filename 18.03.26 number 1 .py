import numpy as np
import scipy.special as sp
import json
import toml
import matplotlib.pyplot as plt
from pathlib import Path


class SphereRCSCalculator:
    def __init__(self, diameter):
        self.diameter = diameter
        self.radius = diameter / 2.0
        self.c = 299792458.0  # Скорость света

    def calculate_an_bn(self, n, x):


        # Значения функций
        jn = sp.spherical_jn(n, x)
        yn = sp.spherical_yn(n, x)

        # Производные
        djn = sp.spherical_jn(n, x, derivative=True)
        dyn = sp.spherical_yn(n, x, derivative=True)

        # Функция Ханкеля 2-го рода (outgoing wave): h2 = j - i*y
        hn = jn - 1j * yn
        dhn = djn - 1j * dyn

        # Вычисление производных от (x * f(x))
        # d/dx (x * f) = f + x * f'
        d_xjn = jn + x * djn
        d_xhn = hn + x * dhn

        an = - d_xjn / d_xhn
        bn = - jn / hn

        return an, bn

    def calculate_rcs(self, frequency):
        wavelength = self.c / frequency
        k = 2 * np.pi / wavelength
        x = k * self.radius

        # Предел суммирования
        n_max = int(np.ceil(x + 4 * np.power(x, 1 / 3) + 1)) + 15
        if n_max < 10: n_max = 10

        series_sum = 0.0 + 0.0j

        for n in range(1, n_max + 1):
            an, bn = self.calculate_an_bn(n, x)
            # Формула (1): сумма (-1)^n * (n + 0.5) * (bn - an)
            term = ((-1) ** n) * (n + 0.5) * (bn - an)
            series_sum += term

        rcs = (wavelength ** 2 / np.pi) * (np.abs(series_sum) ** 2)
        return rcs

    def run_calculation(self, f_min, f_max, points=600):
        freqs = np.logspace(np.log10(f_min), np.log10(f_max), points)
        rcs_vals = []
        lambdas = []

        print(f"Расчет (D={self.diameter}м)...")
        for i, f in enumerate(freqs):
            rcs = self.calculate_rcs(f)
            rcs_vals.append(rcs)
            lambdas.append(self.c / f)
            if i % 100 == 0:
                print(f"  Точка {i}: f={f / 1e9:.2f} ГГц, RCS={rcs:.6f}")

        return freqs, np.array(lambdas), np.array(rcs_vals)


class ResultExporter:
    @staticmethod
    def save_json_type3(freqs, lambdas, rcs, filename):
        data = {
            "freq": freqs.tolist(),
            "lambda": lambdas.tolist(),
            "rcs": rcs.tolist()
        }
        with open(filename, 'w') as f:
            json.dump(data, f, indent=2)
        print(f"Файл {filename} сохранен.")

    @staticmethod
    def plot(freqs, rcs, diameter):
        plt.figure(figsize=(10, 6))
        x_axis = freqs / 1e9
        plt.plot(x_axis, rcs, label=f'D={diameter}м', linewidth=1.5)

        geo_cs = np.pi * (diameter / 2) ** 2
        plt.axhline(y=geo_cs, color='red', linestyle='--', alpha=0.7, label=f'Геом. сечение ({geo_cs:.4f})')

        plt.title(f'ЭПР сферы (Вариант 12)\nD={diameter}м')
        plt.xlabel('Частота, ГГц')
        plt.ylabel('ЭПР, м²')
        plt.grid(True, which='both', linestyle='--', alpha=0.6)
        plt.legend()
        plt.savefig('result_plot.png', dpi=300)
        plt.show()


def main():
    try:
        with open('task_rcs_02.toml', 'r') as f:
            config = toml.load(f)
        v12 = config['data']['variant_12']
        D = float(v12['D'].replace('"', ''))
        fmin = float(v12['fmin'].replace('"', ''))
        fmax = float(v12['fmax'].replace('"', ''))
    except Exception as e:
        print(f"Ошибка: {e}")
        return

    calc = SphereRCSCalculator(D)
    freqs, lambdas, rcs = calc.run_calculation(fmin, fmax)

    exporter = ResultExporter()
    exporter.save_json_type3(freqs, lambdas, rcs, 'result.json')
    exporter.plot(freqs, rcs, D)


if __name__ == '__main__':
    main()
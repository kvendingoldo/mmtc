# -*- coding: utf-8 -*-
# @Author: Alexander Sharov

import scipy.optimize
import scipy
import matplotlib.pyplot as plt

from mpmath import *

# constants
a = 1
b = 1
c = 1
k = 1
x_01 = 3
t_k = 2


def system(time, vars):
    """
    исходная система ДУ
    """
    x1, x2, Ψ1, Ψ2 = vars[0], vars[1], vars[2], vars[3]
    fun1 = x2
    fun2 = Ψ2 / (2 * a) + (-c * x1)
    fun3 = 2 * x1 * a + Ψ2 * c
    fun4 = 2 * x2 * b
    return fun1, fun2, fun3, fun4


def FU(vars):
    """
    граничные условия
    """
    Ψ1, Ψ2 = vars[0], vars[1]
    # решение задачи Коши
    f = odefun(system, 0, [x_01, 0, Ψ1, Ψ2])
    x1_k, x2_k, Ψ1k, Ψ2k = f(t_k)
    # учет условий 4.4 и 4.12
    return x1_k, x2_k, Ψ1k, Ψ2k

# выбор начальных условий, откуда найдем Ψ1_0, Ψ2_0
# вектор из едениц - начальное приближение
x = scipy.optimize.fsolve(FU, [0, 0, 0, 0])
print(x)

# t = tk
print(FU(x))

# решаем задачи Коши с начальными условиями, которые мы нашли
# x[0] = Ψ1_0, x[1] = Ψ2_0
f = odefun(system, 0, [x_01, 0, x[0], x[1]])

N = 10
x1 = []
x2 = []
Ψ1 = []
Ψ2 = []
time = []
dx2dt = []

for i in range(0, N):
    cur_time = i * t_k / N
    f_t = f(cur_time)
    x1.append(f_t[0])
    x2.append(f_t[1])
    Ψ1.append(f_t[2])
    Ψ2.append(f_t[3])
    time.append(cur_time)
    dx2dt.append(f_t[3] / 2 * a)

print("x1=", x1)
print("x2=", x2)
print("Ψ1=", Ψ1)
print("Ψ2=", Ψ2)
print("time=", time)

# зависимость x1(t)
plt.subplot(2, 2, 1)
plt.plot(time, x1)
plt.grid(True)
plt.xlabel('time')
plt.ylabel('x1')
ax = plt.gca()
ax.annotate('', xy=(0, 1), xytext=(-15, 2), ha='left', va='top', xycoords='axes fraction',
            textcoords='offset points', fontsize=20)
ax.annotate('', xy=(0.85, 0.05), ha='left', va='top', xycoords='axes fraction', fontsize=15)

# зависимость x2(t)
plt.subplot(2, 2, 2)
plt.plot(time, x2)
plt.grid(True)
plt.xlabel('time')
plt.ylabel('$x_2$')

# зависимость x1, x2
plt.subplot(2, 2, 3)
plt.plot(x1, x2)
plt.grid(True)
plt.xlabel('$x_1$')
plt.ylabel('$x_2$')

# зависимость
plt.subplot(2, 2, 4)
plt.plot(time, dx2dt)
plt.grid(True)
plt.xlabel('time')
plt.ylabel('$\\dfrac{dx_2}{dt}$')
plt.savefig('plot.png')
plt.show()

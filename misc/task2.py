import scipy.optimize
import scipy
import math
import matplotlib.pyplot as plt

from mpmath import *

# α = const > 0
α = 3
# конечный момент времени
tk = 10
# вектор скорости
v = 150
# ускорение свободного падения
g = 9.8
# масштабы длинны и времени
L = v ** 2 / g
T = v / g

xk1 = 0.7
xk2 = 0.5
xk3 = 1


def system(time, vars):
    x1, x2, x3, Ψ1, Ψ2, Ψ3 = vars[0], vars[1], vars[2], vars[3], vars[4], vars[5]
    fun1 = cos(x3)
    fun2 = sin(x3)
    fun3 = Ψ3 / (2 * α)
    fun4 = 0
    fun5 = 0
    fun6 = Ψ1 * sin(x3) - Ψ2 * cos(x3)
    return fun1, fun2, fun3, fun4, fun5, fun6


def FU(y):
    Ψ1, Ψ2, Ψ3 = y[0], y[1], y[2]
    f = odefun(system, 0, [1, 0, 1, Ψ1, Ψ2, Ψ3])
    x1_k, x2_k, x3_k, Ψ1k, Ψ2k, Ψ3k = f(tk)
    return (xk2 + x1_k * tan(xk3), x3_k, Ψ1k + Ψ2k * tan(xk3))


x = scipy.optimize.fsolve(FU, [0.5, 0.02, 1])
print("x=", x)
print("FU=", FU(x))

f = odefun(system, 0, [1, 0, 1, x[0], x[1], x[2]])
N = 10
j = []
x1 = []
x2 = []
x3 = []
phy1 = []
phy2 = []
phy3 = []
time = []
x2t = []
for i in range(0, N):
    t = i * tk / N
    f_t = f(t)
    j.append(f_t)
    x1.append(f_t[0])
    x2.append(f_t[1])
    x3.append(f_t[2])
    phy1.append(f_t[3])
    phy2.append(f_t[4])
    phy3.append(f_t[5])
    time.append(t)
    x2t.append(f_t[5] / 2 * α)

print("x1=", x1)
print("x2=", x2)
print("x3=", x3)
print("phy1=", phy1)
print("phy2=", phy2)
print("phy3=", phy3)
print("time=", time)

# plt.subplot(2,2,1)
# plt.plot(time, x1)
# plt.grid(True)
##plt.xlabel('x1')
##plt.ylabel('time')
# ax = plt.gca()
# ax.annotate('$x1$', xy=(0, 1), xytext=(-15,2), ha='left', va='top', xycoords='axes fraction', textcoords='offset points', fontsize=20)
# ax.annotate('$time$', xy=(0.85, 0.05), ha='left', va='top', xycoords='axes fraction', fontsize=15)
# plt.savefig("1.png")

# plt.subplot(2,2,2)
# plt.plot(time, x2)
# plt.grid(True)
# plt.xlabel('x2')
# plt.ylabel('time')
# plt.savefig("1.png")

# plt.subplot(2,2,3)
# plt.plot(x1, x2)
# plt.grid(True)
# plt.xlabel('x1')
# plt.ylabel('x2')
# plt.savefig("1.png")

# plt.subplot(2,2,4)
# plt.plot(time, x2t)
# plt.grid(True)
# plt.xlabel('time')
# plt.ylabel('x2t')
# plt.savefig("1.png")
# plt.show()

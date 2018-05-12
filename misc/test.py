from mpmath import *
import scipy
import scipy.optimize
import matplotlib.pyplot as plt
#x01 = Symbol('x01')
a = 4
b = 2
c = 3
d = 5
m = 6
tk = 1

def system(t, x):
    return [(x[3]*c)/m, -x[2], x[1], (0.5*x[3])/a + (d-c*x[0])/m]

def fun1(psi):
    psi10, psi20 = psi[0], psi[1]
    f = odefun(system, 0, [1, 0, psi10, psi20]) # решение задачи Коши
    #f0 = f(t0)
    fk = f(tk)
    return [a*fk[0] - b*fk[1], b*fk[2] + a*fk[3]] # условия 4.4 и 4.12

x = scipy.optimize.fsolve(fun1, [1, 1]) # выбор начальных условий
#найдем psi10, psi20
print(x)
#print(f(0))
print(fun1(x)) # при t = tk
f = odefun(system, 0, [1, 0, x[0], x[1]]) # решение задачи Коши x[0] = psi10, x[1] = psi20

ox3 = []
oy3 = []
ox5 = []
oy5 = []
#график phi = a*fk[0] - b*fk[1]
for i in arange(0,tk,0.01):
    fk = f(i)
    ox3.append(fk[0])
    oy3.append(fk[1])

for i in arange(-2,4,0.1):
    ox5.append(i)
    oy5.append(a * i/b)

plt.axis([-2.0, 2.0, -1.0, 4.0])
plt.plot(ox3,oy3)
plt.plot(ox5,oy5)
plt.grid(True)
plt.savefig("1.png")


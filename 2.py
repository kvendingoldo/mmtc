# -*- coding: utf-8 -*-
# @Author: Alexander Sharov

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as optimize

from scipy.integrate import odeint
from sympy import integrate, symbols
from math import radians, degrees

g = 9.8
velocity = 150

# масштаб длинны
length_scale = velocity ** 2 / g

# α = const > 0
α = 1

# inital values
y0 = 4000
x0 = 3000
θ_0 = radians(42)


# final values
yk0 = 4000
xk2 = yk0 / length_scale
θ_k = radians(50)
xk3 = θ_k

initial_approx = [1.9, -1.6, 1.25, 15.0]


def calculate_hamiltonian(vars):
    x1, x2, x3, Ψ1, Ψ2, Ψ3 = vars
    management = Ψ3 / (2 * α) # 5.18
    # 5.11
    return -1 - α * management ** 2 + Ψ1 * np.cos(x3) + Ψ2 * np.sin(x3) + Ψ3 * management


def system(vars, time):
    x1, x2, x3, Ψ1, Ψ2, Ψ3 = vars
    dx1dt = np.cos(x3)
    dx2dt = np.sin(x3)
    dx3dt = Ψ3 / 2. * α
    dΨ1dt = 0
    dΨ2dt = 0
    dΨ3dt = Ψ1 * np.sin(x3) - Ψ2 * np.cos(x3)
    array = [dx1dt, dx2dt, dx3dt, dΨ1dt, dΨ2dt, dΨ3dt]
    return array


def calculate_ivp(vars):
    Ψ1, Ψ2, Ψ3, Ψ4 = vars
    time = np.linspace(0, Ψ4, 1000)
    init_conditionals = [x0 / length_scale, y0 / length_scale, θ_0, Ψ1, Ψ2, Ψ3]
    return odeint(system, init_conditionals, time)


def calculate_final_condition(vars):
    x1, x2, x3, Ψ1, Ψ2, Ψ3 = vars
    return [xk2 + x1 * np.tan(xk3) - x2,  # 5.13
            xk3 - x3,  # 5.13
            Ψ1 + Ψ2 * np.tan(xk3),  # transversality condition (5.15)
            calculate_hamiltonian(vars)]


def calculate_bvp(vars):
    def subordinate_solution_final_condition(vars):
        ivp_solution = calculate_ivp(vars)
        # the last element of the array of solution values is [-1]
        return calculate_final_condition(ivp_solution[-1])

    return optimize.fsolve(subordinate_solution_final_condition, vars)


def calculate_functional(management, time):
    integrand_function = 1 + α * management ** 2 # 5.4
    return integrate(integrand_function, (symbols('t'), 0, time))


bvp_solution = calculate_bvp(initial_approx)
print('Solution of boundary value problem: %s \n' % bvp_solution)

ivp_solution = calculate_ivp(bvp_solution)
print('Solution of initial value problem: %s \n' % ivp_solution)

check = calculate_final_condition(ivp_solution[-1])
print('Check: %s \n' % check)

time = np.linspace(0, bvp_solution[-1], 1000)
x1 = ivp_solution[:, 0]
x2 = ivp_solution[:, 1]
Ψ3 = ivp_solution[:, 5]

management = []
for element in Ψ3:
    u = element / (2 * α) # 5.1.1
    management.append(degrees(np.arctan(u)))

x22 = []
for element in x1:
    x22.append(xk2 + element * np.tan(xk3))

functional = []
for ind in range(len(management)):
    value = calculate_functional(management[ind], bvp_solution[-1])
    functional.append(value)

min_value = min(functional)
print('Functional: %f' % min_value)

fig = plt.figure()

ax1 = plt.subplot2grid((2, 2), (0, 0))
ax1.plot(time, x1)
ax1.set_xlabel('$t$')
ax1.set_ylabel('$x_1(t)$')

ax2 = plt.subplot2grid((2, 2), (1, 0))
ax2.plot(time, x2)
ax2.set_xlabel('$t$')
ax2.set_ylabel('$x_2(t)$')

ax3 = plt.subplot2grid((2, 2), (0, 1))
ax3.set_title('optimal management')
ax3.plot(time, management)
ax3.set_xlabel('$t$')
ax3.set_ylabel('$y(t)$')

ax4 = plt.subplot2grid((2, 2), (1, 1))
ax4.set_title('mechanical trajectory of system')
ax4.plot(x1, x2)
ax4.plot(x1, x22, '-b', label='diversity')
ax4.set_xlabel('$x_1$')
ax4.set_ylabel('$x_2(x_1)$')

for ax in fig.axes:
    ax.grid(True)

plt.tight_layout()
plt.show()

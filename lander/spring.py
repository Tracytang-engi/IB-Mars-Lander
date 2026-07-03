

import numpy as np
import matplotlib.pyplot as plt

# mass, spring constant, initial position and velocity
m = 1
k = 1
x = 0
v = 1
xver = 0
vver = 1

# Save initial conditions for exact analytical solution
x0 = x
v0 = v
# simulation time, timestep and time
t_max = 1000
dt = 2
t_array = np.arange(0, t_max, dt)

# initialise empty lists to record trajectories
x_list = []
v_list = []

# Euler integration
for t in t_array:

    # append current state to trajectories
    x_list.append(x)
    v_list.append(v)

    # calculate new position and velocity
    a = -k * x / m
    x = x + dt * v
    v = v + dt * a

# convert trajectory lists into arrays, so they can be sliced (useful for Assignment 2)
x_array = np.array(x_list)
v_array = np.array(v_list)


# plot the position-time graph
plt.figure(1)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_array, label='Euler x (m)')
plt.plot(t_array, v_array, label='Euler v (m/s)')
plt.legend()
plt.show()

#Verlet Method：
xver_list = []
vver_list = []

# MODIFIED: Initial conditions for Verlet
x_prev = xver  # x(0)
a0 = -k * x_prev / m
x_curr = x_prev + dt * vver + 0.5 * (dt ** 2) * a0  # x(dt), via Taylor/Euler-like start

# Record t=0 state
xver_list.append(x_prev)
vver_list.append(vver)

# MODIFIED: Main Verlet loop starts from second time point
for i in range(1, len(t_array)):
    # position Verlet: x(t+dt) = 2x(t) - x(t-dt) + dt^2 * a(t)
    a_curr = -k * x_curr / m
    x_next = 2 * x_curr - x_prev + (dt ** 2) * a_curr

    # central-difference velocity at current point
    v_curr = (x_next - x_prev) / (2 * dt)

    xver_list.append(x_curr)
    vver_list.append(v_curr)

    x_prev = x_curr
    x_curr = x_next

# MODIFIED: Convert Verlet trajectories to arrays for plotting
xver_array = np.array(xver_list)
vver_array = np.array(vver_list)

plt.figure(2)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, xver_array, label='xver (m)')
plt.plot(t_array, vver_array, label='vver (m/s)')
plt.legend()
plt.show()


# Exact analytical solution for simple harmonic oscillator:
# x(t) = x0 cos(wt) + (v0/w) sin(wt)
# v(t) = -x0 w sin(wt) + v0 cos(wt), where w = sqrt(k/m)
omega = np.sqrt(k / m)
x_exact = x0 * np.cos(omega * t_array) + (v0 / omega) * np.sin(omega * t_array)
v_exact = -x0 * omega * np.sin(omega * t_array) + v0 * np.cos(omega * t_array)

# Figure 3: exact analytical solution only
plt.figure(3)
plt.clf()
plt.xlabel('time (s)')
plt.grid()
plt.plot(t_array, x_exact, label='Exact x (m)')
plt.plot(t_array, v_exact, label='Exact v (m/s)')
plt.legend()
plt.show()
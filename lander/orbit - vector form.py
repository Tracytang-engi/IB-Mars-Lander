import numpy as np
import matplotlib.pyplot as plt

# Mars gravitational constants
G = 6.673e-11
M = 6.42e23

#initial condition
t_max = 20000
dt = 1.0
t_array = np.arange(0, t_max, dt)

# Initial state (circular-like orbit in x-y plane)
r0 = np.array([3386000.0 + 200000.0, 0.0, 0.0], dtype=float)
v_circ = np.sqrt(G * M / np.linalg.norm(r0))
v0 = np.array([0.0, v_circ, 0.0], dtype=float)


def accel_gravity(r):
    r_norm = np.linalg.norm(r)
    return -(G * M / (r_norm ** 3)) * r


# Euler integration (vector form)
r_e = r0.copy()
v_e = v0.copy()
r_e_list = []
v_e_list = []
for _ in t_array:
    r_e_list.append(r_e.copy())
    v_e_list.append(v_e.copy())
    a_e = accel_gravity(r_e)
    r_e = r_e + dt * v_e
    v_e = v_e + dt * a_e

r_e_array = np.array(r_e_list)
v_e_array = np.array(v_e_list)


# Verlet integration (vector form)
r_prev = r0.copy()
a0 = accel_gravity(r_prev)
r_curr = r_prev + dt * v0 + 0.5 * (dt ** 2) * a0
r_v_list = [r_prev.copy()]
v_v_list = [v0.copy()]

for _ in range(1, len(t_array)):
    a_curr = accel_gravity(r_curr)
    r_next = 2.0 * r_curr - r_prev + (dt ** 2) * a_curr
    v_curr = (r_next - r_prev) / (2.0 * dt)
    r_v_list.append(r_curr.copy())
    v_v_list.append(v_curr.copy())
    r_prev = r_curr
    r_curr = r_next

r_v_array = np.array(r_v_list)
v_v_array = np.array(v_v_list)


# Plot orbital trajectories (x-y plane)
plt.figure(1)
plt.clf()
plt.xlabel("x (m)")
plt.ylabel("y (m)")
plt.axis("equal")
plt.grid()
plt.plot(r_e_array[:, 0], r_e_array[:, 1], label="Euler orbit")
plt.plot(r_v_array[:, 0], r_v_array[:, 1], label="Verlet orbit")
plt.legend()
plt.show()

# Plot radius vs time for stability comparison
plt.figure(2)
plt.clf()
plt.xlabel("time (s)")
plt.ylabel("|r| (m)")
plt.grid()
plt.plot(t_array, np.linalg.norm(r_e_array, axis=1), label="Euler |r|")
plt.plot(t_array, np.linalg.norm(r_v_array, axis=1), label="Verlet |r|")
plt.legend()
plt.show()
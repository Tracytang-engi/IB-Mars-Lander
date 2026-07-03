from typing import Any

import numpy as np
import matplotlib.pyplot as plt
from numpy import ndarray, dtype, float64

# Physical constants for Mars
G = 6.673e-11
M_MARS = 6.42e23
R_MARS = 3386000.0
MU = G * M_MARS


def gravity_acceleration(position):
    r = np.linalg.norm(position)#r is the length of the position vector
    return -(MU / (r ** 3)) * position


def integrate_euler(position0, velocity0, dt, t_max):
    t_array: ndarray[tuple[int], dtype[float64 | Any]] | ndarray[tuple[int], dtype[Any]] = np.arange(0.0, t_max, dt)
    position = position0.copy()#copy is a method for arrays
    velocity = velocity0.copy()

    pos_hist = []
    vel_hist = []

    for t in t_array:#dummy variable
        pos_hist.append(position.copy())
        vel_hist.append(velocity.copy())
        acc = gravity_acceleration(position)
        position = position + dt * velocity
        velocity = velocity + dt * acc

    return t_array, np.array(pos_hist), np.array(vel_hist)


def integrate_verlet(position0, velocity0, dt, t_max):
    t_array = np.arange(0.0, t_max, dt)

    r_prev = position0.copy()
    a0 = gravity_acceleration(r_prev)
    r_curr = r_prev + dt * velocity0 + 0.5 * (dt ** 2) * a0

    pos_hist = [r_prev.copy()]
    vel_hist = [velocity0.copy()]

    for _ in range(1, len(t_array)):
        a_curr = gravity_acceleration(r_curr)
        r_next = 2.0 * r_curr - r_prev + (dt ** 2) * a_curr
        v_curr = (r_next - r_prev) / (2.0 * dt)

        pos_hist.append(r_curr.copy())
        vel_hist.append(v_curr.copy())

        r_prev = r_curr
        r_curr = r_next

    return t_array, np.array(pos_hist), np.array(vel_hist)


def build_scenario(scenario_id):
    altitude0 = 200000.0
    r0 = np.array([R_MARS + altitude0, 0.0, 0.0], dtype=float)

    v_circ = np.sqrt(MU / np.linalg.norm(r0))
    v_esc = np.sqrt(2.0 * MU / np.linalg.norm(r0))

    if scenario_id == 1:
        # Straight down descent: zero initial velocity
        v0 = np.array([0.0, 0.0, 0.0], dtype=float)
        t_max = 1500.0
        dt = 0.1
        title = "Scenario 1: Straight down descent"
    elif scenario_id == 2:
        # Circular orbit
        v0 = np.array([0.0, v_circ, 0.0], dtype=float)
        t_max = 12000.0
        dt = 1.0
        title = "Scenario 2: Circular orbit"
    elif scenario_id == 3:
        # Elliptical orbit (below circular speed)
        v0 = np.array([0.0, 0.8 * v_circ, 0.0], dtype=float)
        t_max = 12000.0
        dt = 1.0
        title = "Scenario 3: Elliptical orbit"
    elif scenario_id == 4:
        # Hyperbolic escape (above escape speed)
        v0 = np.array([0.0, 1.1 * v_esc, 0.0], dtype=float)
        t_max = 12000.0
        dt = 1.0
        title = "Scenario 4: Hyperbolic escape"
    else:
        raise ValueError("scenario_id must be 1..4")

    return r0, v0, dt, t_max, title


def plot_scenario(scenario_id):
    r0, v0, dt, t_max, title = build_scenario(scenario_id)

    t_e, r_e, _ = integrate_euler(r0, v0, dt, t_max)
    t_v, r_v, _ = integrate_verlet(r0, v0, dt, t_max)

    if scenario_id == 1:
        alt_e = np.linalg.norm(r_e, axis=1) - R_MARS
        alt_v = np.linalg.norm(r_v, axis=1) - R_MARS

        plt.figure()
        plt.title(title)
        plt.xlabel("time (s)")
        plt.ylabel("altitude (m)")
        plt.grid()
        plt.plot(t_e, alt_e, label="Euler altitude")
        plt.plot(t_v, alt_v, label="Verlet altitude")
        plt.legend()
    else:
        plt.figure()
        plt.title(title)
        plt.xlabel("x (m)")
        plt.ylabel("y (m)")
        plt.axis("equal")
        plt.grid()
        plt.plot(r_e[:, 0], r_e[:, 1], label="Euler trajectory")
        plt.plot(r_v[:, 0], r_v[:, 1], label="Verlet trajectory")
        plt.scatter([0.0], [0.0], s=30, label="Mars center")
        plt.legend()


if __name__ == "__main__":
    for sid in (1, 2, 3, 4):
        plot_scenario(sid)
    plt.show()

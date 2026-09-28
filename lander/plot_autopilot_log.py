import pathlib
import numpy as np
import matplotlib.pyplot as plt


log_path = pathlib.Path("autopilot_descent_log_s1_0016,02,035.txt")
if not log_path.exists():
    raise FileNotFoundError(
        "autopilot_descent_log.txt not found. "
        "Run the simulator with autopilot enabled in Scenario 1 or 5 first."
    )

data = np.loadtxt(log_path, skiprows=1)
if data.ndim == 1:
    data = data.reshape(1, -1)

time_s = data[:, 0]
altitude_m = data[:, 1]
v_r_mps = data[:, 2]
v_target_mps = data[:, 3]
throttle = data[:, 4]

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Assignment 5 requested comparison: actual vs target descent rate against altitude
axes[0].plot(altitude_m, v_r_mps, label="Actual $v\\cdot e_r$")
axes[0].plot(altitude_m, v_target_mps, "--", label="Target $v\\cdot e_r$")
axes[0].set_xlabel("Altitude h (m)")
axes[0].set_ylabel("Radial velocity (m/s)")
axes[0].set_title("Autopilot: actual vs target")
axes[0].grid(True)
axes[0].legend()

# Helpful diagnostic plot
axes[1].plot(altitude_m, throttle, color="tab:red", label="Throttle")
axes[1].set_xlabel("Altitude h (m)")
axes[1].set_ylabel("Throttle")
axes[1].set_title("Throttle vs altitude")
axes[1].grid(True)
axes[1].legend()

fig.suptitle("Assignment 5 log analysis")
fig.tight_layout()
plt.show()

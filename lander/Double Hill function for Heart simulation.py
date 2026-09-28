import numpy as np
import matplotlib.pyplot as plt

# ---------- 生理参数设置（均在正常范围）----------
E_min = 0.05          # 舒张末期弹性 (mmHg/mL)
E_max = 1.0           # 收缩末期峰值弹性 (mmHg/mL)，按你的要求设为1

# 控制收缩上升沿的参数（时间常数和斜率）
tau1 = 0.12           # 收缩开始时间 (s)
m1   = 2.0            # 收缩陡峭度

# 控制舒张下降沿的参数
tau2 = 0.30           # 舒张开始时间 (s)
m2   = 2.5            # 舒张陡峭度

T = 0.8               # 心动周期 (s)，对应心率 75 bpm
dt = 0.001            # 时间步长 (s)
t = np.arange(0, T, dt)   # 时间数组

# ---------- 计算双 Hill 函数 ----------
# 第一部分：收缩上升支 (S形曲线)
hill1 = (t / tau1) ** m1 / (1 + (t / tau1) ** m1)

# 第二部分：舒张下降支 (反向S形曲线)
hill2 = 1 / (1 + (t / tau2) ** m2)

# 未归一化的乘积
product = hill1 * hill2

# 计算归一化因子 A（使峰值等于1）
A = np.max(product)

# 最终的弹性函数
E_t = E_min + (E_max - E_min) / A * product

# ---------- 可视化 ----------
plt.figure(figsize=(10, 5))
plt.plot(t, E_t, linewidth=2, color='darkred')
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Elastance E(t) (mmHg/mL)', fontsize=12)
plt.title('Ventricular Elastance over a Cardiac Cycle (Double Hill Function)', fontsize=14)
plt.grid(True, linestyle='--', alpha=0.6)
plt.xlim(0, T)

# 标注关键阶段
plt.axvline(x=0, color='gray', linestyle=':', alpha=0.5)
plt.axvline(x=tau1, color='blue', linestyle=':', alpha=0.5, label='τ₁ (收缩上升)')
plt.axvline(x=tau2, color='green', linestyle=':', alpha=0.5, label='τ₂ (舒张下降)')
plt.legend()
plt.tight_layout()
plt.show()

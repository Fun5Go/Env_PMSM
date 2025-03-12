import gymnasium as gym
import numpy as np
import matplotlib.pyplot as plt
from stable_baselines3 import DDPG
from stable_baselines3.common.monitor import Monitor
from stable_baselines3.common.evaluation import evaluate_policy

# 1️⃣ 选择环境
from environments import EnvLoad3RL  # 确保这个文件包含自定义环境

# 2️⃣ 设定系统参数
sys_params_dict = {
    "dt": 1 / 10e3,  # 采样时间 [s]
    "r": 1,          # 电阻 [Ohm]
    "l": 1e-2,       # 电感 [H]
    "vdc": 500,      # DC 母线电压 [V]
    "we_nom": 200 * 2 * np.pi,  # 额定角速度 [rad/s]
}

# 计算最大电流 [A]
def idq_max_norm(vdq_max, we, r, l):
    return vdq_max / np.sqrt(r**2 + (we * l)**2)

sys_params_dict["i_max"] = idq_max_norm(sys_params_dict["vdc"]/2, sys_params_dict["we_nom"],
                                        sys_params_dict["r"], sys_params_dict["l"])
sys_params_dict["reward"] = "quadratic"  # 设定奖励函数类型

# 3️⃣ 创建 Gym 环境
env = EnvLoad3RL(sys_params=sys_params_dict)
env = Monitor(env)  # 监视环境，记录训练数据

# 4️⃣ 创建 DDPG RL 代理
model = DDPG("MlpPolicy", env, verbose=1, learning_rate=1e-3, buffer_size=50000, batch_size=64)

# 5️⃣ 训练 RL 代理
TIMESTEPS = 1000  # 训练步数
model.learn(total_timesteps=TIMESTEPS, log_interval=10)

# 6️⃣ 评估训练效果
mean_reward, std_reward = evaluate_policy(model, env, n_eval_episodes=10)
print(f"评估结果：平均奖励 = {mean_reward:.2f} ± {std_reward:.2f}")

# 7️⃣ 可视化 RL 训练过程
obs, _ = env.reset()
rewards = []
states = []
actions = []

for _ in range(500):  # 运行 500 步
    action, _ = model.predict(obs, deterministic=True)  # 选择最佳动作
    obs, reward, done, truncated, _ = env.step(action)  # 执行动作
    states.append(obs[:2])  # 记录状态（id, iq）
    actions.append(action)  # 记录动作
    rewards.append(reward)  # 记录奖励
    if done or truncated:
        break

# 8️⃣ 绘制结果
fig, axs = plt.subplots(3, 1, figsize=(10, 8))

axs[0].plot([s[0] for s in states], label="i_d")
axs[0].plot([s[1] for s in states], label="i_q")
axs[0].set_ylabel("Current (A)")
axs[0].legend()
axs[0].set_title("电流轨迹")

axs[1].plot(actions, label=["V_d", "V_q"])
axs[1].set_ylabel("Voltage (V)")
axs[1].legend()
axs[1].set_title("电压控制")

axs[2].plot(rewards, label="Reward")
axs[2].set_ylabel("Reward")
axs[2].legend()
axs[2].set_title("奖励变化")

plt.xlabel("Time Steps")
plt.show()

import numpy as np
from pyssa import Simulation
import matplotlib.pyplot as plt

# zero order
V_r = np.array([[0]])
V_p = np.array([[1]])
X0 = np.array([100], dtype=np.int64)
k_det = np.array([1.1])
sim = Simulation(V_r, V_p, X0, k_det)
sim.simulate()
sim.plot(disp=False)
plt.savefig("images/ex_0.png")

# first order
V_r = np.array([[1], [0]])
V_p = np.array([[0], [1]])
X0 = np.array([100, 20], dtype=np.int64)
k_det = np.array([1.1])
sim = Simulation(V_r, V_p, X0, k_det)
sim.simulate()
sim.plot(disp=False)
plt.savefig("images/ex_1a.png")
sim.simulate(algorithm="tau_leaping", n_rep=20)
sim.plot(plot_indices=[1], names=["B"], disp=False)
plt.savefig("images/ex_1b.png")

# enzyme kinetics
V_r = np.array([[1, 0, 0], [1, 0, 0], [0, 1, 1], [0, 0, 0]])
V_p = np.array([[0, 1, 0], [0, 1, 1], [1, 0, 0], [0, 0, 1]])
X0 = np.array([200, 50, 0, 0], dtype=np.int64)
k_det = np.array([0.006, 0.005, 0.1])
sim = Simulation(V_r, V_p, X0, k_det)
sim.simulate(max_t=50, n_rep=10)
sim.plot(disp=False, names=["S", "E", "SE", "P"])
plt.savefig("images/ex_2a.png")
sim = Simulation(V_r, V_p, X0, k_det, volume=5.0)
sim.simulate(max_t=50, n_rep=10)
sim.plot(disp=False, names=["S", "E", "SE", "P"])
plt.savefig("images/ex_2b.png")
# sim = Simulation(V_r, V_p, X0, k_det, chem_flag=True)
# sim.simulate(max_t=50, n_rep=10)
# sim.plot(disp=False, names=["S", "E", "SE", "P"])
# plt.savefig("images/ex_2c.png")

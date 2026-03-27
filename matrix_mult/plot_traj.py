import numpy as np
import matplotlib.pyplot as plt

# ---- Load trajectory file ----
data = np.loadtxt("traj_1.dat")

t = data[:,0]
x = data[:,1]
y = data[:,2]

# ---- Plot trajectory in xy ----
plt.figure()
plt.plot(x, y, "-")
plt.gca().set_aspect("equal")

plt.xlabel("x")
plt.ylabel("y")
plt.title("Trajectory of particle 1")

plt.show()

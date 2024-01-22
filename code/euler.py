import matplotlib
from matplotlib import pyplot
import numpy as np

rho = 1
L = 1
u = 1
Gamma = 0.02
phi_0 = 0
phi_L = 1
num_nodes = 50
MAX_T = 1.00
method = 2 # 1: Central Difference; 2: Backward in Convection

Pe = rho * u * L / Gamma

# Initial Mesh
dx = L / num_nodes
x = np.linspace(0, L, num_nodes + 1)

# Theoretical Solution
phi_theory = phi_0 + (np.exp(Pe * x / L) - 1) / (np.exp(Pe) - 1) * (phi_L - phi_0)

# Calculate Matrix
t = 0.00
dt = 0.1
while t < MAX_T:
    d = Gamma * dt / rho / dx / dx
    c = u * dt / dx
    A_E = np.zeros(num_nodes + 1)
    A_W = np.zeros(num_nodes + 1)
    A_P = np.zeros(num_nodes + 1)
    for i in range(1, num_nodes):
        A_E[i] = c / 2 - d
        A_W[i] = -d - c / 2
        A_P[i] = 1 - (A_E[i] + A_W[i])
    # Boundary condition
    Q_0 = -A_W * phi_0
    Q_N = -A_E * phi_L
    phi = np.zeros(num_nodes + 1)
    phi[1] = phi_0
    phi[num_nodes] = phi_L
    # Solve
    A = np.diag(A_P * np.ones(num_nodes - 1)) + np.diag(A_W * np.ones(num_nodes - 2), k=-1) + np.diag(A_E * np.ones(num_nodes - 2), k=1)
    A = np.full_like(A, A)
    Q = np.zeros((num_nodes - 1, num_nodes - 1))
    Q[1, 1] = (Q_0)
    Q[-1, 1] = (Q_N)
    phi[1:-1] = np.linalg.solve(A, Q)
    # Plot
    pyplot.plot(x, phi)
        

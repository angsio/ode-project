from compute_components import *

omega = 2*np.pi*100
ang_vel = 30

R = 0.05
r = 0.0323
a = 0.02
b = 0.02

d = 0.003
sigma = 5.79e7

B_0 = 0.3

x = np.linspace(0, a, 300)
y = np.linspace(0, b, 300)

# infitesimal dimensions
dx = x[1] - x[0]
dy = y[1] - y[0]
dA = dx*dy

# Unchanging Quantities
X, Y = np.meshgrid(x[1:-1], y[1:-1], indexing='ij')
R = position_field(X, Y, a, b, r)

Ex = E_x(a, b, r, ang_vel, B_0, X, Y)
Ey = E_y(a, b, r, ang_vel, B_0, X, Y)

v = velocity_field(a, b, r, ang_vel, X, Y)

T = []
times = np.linspace(0, 2*np.pi/omega, 200)

for t in times:
    E = electric_field(Ex, Ey, omega, t)
    B = magnetic_field(X, B_0, omega, t)
    J = current_density_field(E, v, B, sigma, d)
    tau_d = torque_density(R, J, B)
    tau_dz = tau_d[..., 2]
    T.append(np.sum(tau_dz) * dA)

T = np.array(T)
T_rms = np.sqrt(np.mean(T**2))
print(T_rms)


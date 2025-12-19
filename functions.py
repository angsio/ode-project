import numpy as np

# Calculate x-component of DC Electric Field at every point in PPA
def E_x(a, b, r, ang_vel, B_0, X, Y):

    term_1 = -(b + 2*r)*np.atan((a - Y)/(b - X))
    term_2 = (b - 2*r)*np.atan((a - Y)/X)
    term_3 = -(b + 2*r)*np.atan(Y/(b - X))
    term_4 = (b - 2*r)*np.atan(Y/X)

    # Long Term
    numerator = (4*(X**2) + 4*(Y**2)) * (4*(a**2) + 4*(X**2) - 8*a*Y + 4*(Y**2))
    denominator = (4*(b**2) + 4*(X**2) - 8*b*X + 4*(Y**2)) * (4*(a**2) + 4*(b**2) - 8*b*X - 8*a*Y + 4*(X**2) + 4*(Y**2))
    long_term = 0.5*a*np.log(numerator/denominator)

    final_sum = ((B_0*ang_vel)/(4*np.pi)) * (term_1 + term_2 + term_3 + term_4 + long_term)

    return final_sum

# Calculate y-component of DC Electric Field at every point in PPA
def E_y(a, b, r, ang_vel, B_0, X, Y):

    term_1 = (2*a)*np.atan((b - X)/(a - Y))
    term_2 = (2*a)*np.atan(X/(a - Y))
    term_3 = -(2*a)*np.atan(X/Y)
    term_4 = (2*a)*np.atan((X - b)/Y)

    # Long Term 1
    numerator_1 = (X**2 + Y**2) * (a**2 + b**2 - 2*b*X + X**2 - 2*a*Y + Y**2)
    denominator_1 = (b**2 + X**2 - 2*b*X + Y**2) * (a**2 - 2*a*Y + Y**2 + X**2)
    long_term_1 = 2*r*np.log(numerator_1/denominator_1)

    # Long Term 2
    numerator_2 = (a**2 - 2*a*Y + Y**2 + X**2) * (a**2 + b**2 - 2*b*X + X**2 - 2*a*Y + Y**2)
    denominator_2 = (b**2 + X**2 - 2*b*X + Y**2) * (X**2 + Y**2)
    long_term_2 = b*np.log(numerator_2/denominator_2)

    final_sum = -((B_0*ang_vel)/(8*np.pi)) * (term_1 + term_2 + term_3 + term_4 + long_term_1 + long_term_2)

    return final_sum

# Calculate velocity at every given point in PPA
def velocity_field(a, b, r, ang_vel, X, Y):
    
    x_comp = -ang_vel*(Y - b/2)
    y_comp = ang_vel*(X + r - a/2)
    z_comp = np.zeros_like(X)

    v = np.stack((x_comp, y_comp, z_comp), axis=-1)

    return v # Shape (200, 200, 3)

# Calculate magnetic field at every given point in PPA
def magnetic_field(X, B_0, omega, t):

    x_comp = np.zeros_like(X)
    y_comp = np.zeros_like(X)
    z_comp = B_0*np.sin(omega*t) * np.ones_like(X)

    B = np.stack((x_comp, y_comp, z_comp), axis=-1)

    return B # Shape (200, 200, 3)

# Calculate electric field at every given point in PPA
def electric_field(Ex, Ey, omega, t):

    x_comp = Ex*np.cos(omega*t)
    y_comp = Ey*np.cos(omega*t)
    z_comp = np.zeros_like(Ex)

    E = np.stack((x_comp, y_comp, z_comp), axis=-1)

    return E # Shape (200, 200, 3)

# Calculate current density field at every given point in PPA
def current_density_field(E, v, B, sigma, d):

    J = sigma*d*(E + np.cross(v, B, axis=-1))

    return J

# Define position field for torque lever arm relative to disk center
# This formula translates the disk center to the previously defined PPA Origin
def position_field(X, Y, a, b, r):

    x_comp = X + r - a/2
    y_comp = Y - b/2
    z_comp = np.zeros_like(X)

    R = np.stack((x_comp, y_comp, z_comp), axis=-1)

    return R

# Calculate "torque density" at every point in PPA
# Intended for surface integration
def torque_density(R, J, B):
    
    f = np.cross(J, B, axis=-1)
    tau_density = np.cross(R, f, axis=-1)

    return tau_density








    


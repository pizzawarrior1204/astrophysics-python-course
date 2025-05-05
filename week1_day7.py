import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def equations_of_motion(t, state, mu_sun, mu_jup, a_jup):
    """
    Compute derivatives for Earth + Jupiter orbiting Sun (Sun fixed at origin).
    state = [x_e, y_e, vx_e, vy_e, x_j, y_j, vx_j, vy_j] (Earth, Jupiter)
    mu_sun: G * M_sun
    mu_jup: G * M_jupiter
    a_jup: Jupiter's semi-major axis (assumed circular orbit)
    """
    # Unpack state
    x_e, y_e, vx_e, vy_e, x_j, y_j, vx_j, vy_j = state
    
    # Earth's position and distance from Sun
    r_e = np.sqrt(x_e**2 + y_e**2)
    
    # Jupiter's position (circular orbit approximation)
    # omega_j = sqrt(mu_sun / a_jup^3) from Kepler's third law
    omega_j = np.sqrt(mu_sun / a_jup**3)
    x_j = a_jup * np.cos(omega_j * t)
    y_j = a_jup * np.sin(omega_j * t)
    # Jupiter's velocity (circular orbit)
    vx_j = -a_jup * omega_j * np.sin(omega_j * t)
    vy_j = a_jup * omega_j * np.cos(omega_j * t)
    
    # Distance between Earth and Jupiter
    r_ej = np.sqrt((x_e - x_j)**2 + (y_e - y_j)**2)
    
    # Earth's acceleration
    # Sun's gravity on Earth
    ax_e = -mu_sun * x_e / r_e**3
    ay_e = -mu_sun * y_e / r_e**3
    # Jupiter's perturbation on Earth
    ax_e -= mu_jup * (x_e - x_j) / r_ej**3
    ay_e -= mu_jup * (y_e - y_j) / r_ej**3
    
    # Jupiter's acceleration (circular, so we don't integrate it, but include for completeness)
    ax_j = -mu_sun * x_j / (a_jup**3)  # Simplified for circular orbit
    ay_j = -mu_sun * y_j / (a_jup**3)
    
    return [vx_e, vy_e, ax_e, ay_e, vx_j, vy_j, ax_j, ay_j]

def simulate_perturbed_orbit(a_earth=1.496e9, e_earth=0.0167, a_jup=7.785e9, n=2000):
    """
    Simulate Earth's orbit perturbed by Jupiter.
    a_earth: Earth's semi-major axis (m)
    e_earth: Earth's eccentricity
    a_jup: Jupiter's semi-major axis (m)
    n: number of time points
    """
    # Gravitational parameters
    mu_sun = 1.32712440018e20  # G * M_sun (m^3 s^-2)
    m_jup = 1.898e27  # Jupiter's mass (kg)
    m_jup *= 1000
    G = 6.67430e-11
    mu_jup = G * m_jup  # G * M_jupiter
    
    
    # Orbital period of Earth
    T_earth = 2 * np.pi * np.sqrt(a_earth**3 / mu_sun)
    
    # Initial conditions for Earth at perihelion
    r0_e = a_earth * (1 - e_earth)
    v0_e = np.sqrt(mu_sun * (2/r0_e - 1/a_earth))
    state0_e = [r0_e, 0, 0, v0_e]
    
    # Initial conditions for Jupiter (circular orbit, t=0)
    omega_j = np.sqrt(mu_sun / a_jup**3)
    x0_j = a_jup
    y0_j = 0
    vx0_j = 0
    vy0_j = a_jup * omega_j
    state0_j = [x0_j, y0_j, vx0_j, vy0_j]
    
    # Combined initial state (Earth + Jupiter)
    state0 = state0_e + state0_j
    
    # Time span (one Earth year)
    t_span = (0, T_earth*12)
    t_eval = np.linspace(0, T_earth*12, n)
    
    # Integrate
    sol = solve_ivp(equations_of_motion, t_span, state0, args=(mu_sun, mu_jup, a_jup),
                    t_eval=t_eval, method='RK45', rtol=1e-8, atol=1e-8)
    
    # Extract solution
    x_e, y_e = sol.y[0], sol.y[1]
    x_j, y_j = sol.y[4], sol.y[5]
    
    # Plot
    fig, ax = plt.subplots()
    ax.plot(x_e, y_e, 'b-', label="Earth's Perturbed Orbit")
    ax.plot(x_j, y_j, 'g-', label="Jupiter's Orbit (Circular)")
    ax.scatter([0], [0], color='red', label='Sun')
    ax.set_title("Earth's Orbit Perturbed by Jupiter")
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.axis('equal')
    ax.legend()
    ax.grid(True)
    plt.show()

if __name__ == "__main__":
    simulate_perturbed_orbit()
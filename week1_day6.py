import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

def equations_of_motion(t, state, mu):
    """
    Compute derivatives for a 2D orbit.
    state = [x, y, vx, vy]
    mu = G * M (gravitational parameter)
    """
    x, y, vx, vy = state
    r = np.sqrt(x**2 + y**2)
    ax = -mu * x / r**3
    ay = -mu * y / r**3
    return [vx, vy, ax, ay]

def simulate_orbit(a=1.496e9, e=0.5, n=2000):
    """
    Simulate orbit numerically.
    a: semi-major axis (meters)
    e: eccentricity
    n: number of time points
    """
    # Gravitational parameter (Sun: M = 1.989e30 kg, G = 6.67430e-11 m^3 kg^-1 s^-2)
    mu = 1.32712440018e20  # G * M_sun (m^3 s^-2)
    
    # Orbital period (Kepler's third law: T = 2pi * sqrt(a^3 / mu))
    T = 2 * np.pi * np.sqrt(a**3 / mu)
    print(f"Orbital period T = {T} seconds")
    
    # Initial conditions at perihelion (theta = 0)
    r0 = a * (1 - e)  # Distance at perihelion
    v0 = np.sqrt(mu * (2/r0 - 1/a))  # Velocity at perihelion
    state0 = [r0, 0, 0, v0]  # [x, y, vx, vy]
    
    # Time span
    t_span = (0, T)
    t_eval = np.linspace(0, T, n)
    
    # Integrate with tighter tolerances
    sol = solve_ivp(equations_of_motion, t_span, state0, args=(mu,), t_eval=t_eval, 
                    method='RK45', rtol=1e-8, atol=1e-8)
    
    # Extract solution
    x, y = sol.y[0], sol.y[1]
    
    # Plot
    fig, ax = plt.subplots()
    ax.plot(x, y, 'b-', label='Numerical Orbit')
    ax.scatter([0], [0], color='red', label='Star')
    ax.set_title('Numerically Integrated Orbit')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.axis('equal')
    ax.legend()
    ax.grid(True)
    plt.show()

if __name__ == "__main__":
    simulate_orbit(e=0.5)
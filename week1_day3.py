import numpy as np
import matplotlib.pyplot as plt

def plot_orbit(a=1.496e9,e=0.0167, n=100):
    """
    plot elliptical orbit.
    a: semi-major axis (meters)
    e: eccentricity
    n: number of points
    """
    # Generate theta values
    theta = np.linspace(0, 2 * np.pi, n)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))  # Polar equation
    x = r * np.cos(theta)  # Convert to Cartesian coordinates
    y = r * np.sin(theta)

    plt.plot(x, y,label='Orbit', color='blue')
    plt.scatter([0], [0], color='red', label='Star')
    plt.title('Elliptical Orbit')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    plot_orbit() # earth orbit a=1 au, e=0.0167
    
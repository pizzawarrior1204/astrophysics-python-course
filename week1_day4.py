import numpy as np
import matplotlib.pyplot as plt

def plot_orbit_with_planet(a=1.496e9,e=0.6167,theta_planet=0,n=100):
    """
    plot elliptical orbit with planet
    a: semi-major axis (meters)
    e: eccentricity
    theta_planet: angle of planet (randians)
    n: number of points to plot
    """
    theta = np.linspace(0, 2*np.pi, n)
    r = a * (1-e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    
    #planet position
    r_p = a * (1-e**2) / (1 + e * np.cos(theta_planet))
    x_p = r_p * np.cos(theta_planet)
    y_p = r_p * np.sin(theta_planet)

    plt.plot(x, y, label='Orbit',color='blue')
    plt.scatter([0],[0],color='red',label='Star')
    plt.scatter(x_p,y_p,color='green',label='Planet')
    plt.title('Elliptical Orbit with Planet')
    plt.xlabel('x (m)')
    plt.ylabel('y (m)')
    plt.axis('equal')
    plt.legend()
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    plot_orbit_with_planet(theta_planet=np.pi*3/4) # plantet at 45 degrees

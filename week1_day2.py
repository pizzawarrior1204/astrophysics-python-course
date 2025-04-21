def gravitational_force(m1, m2, r, G=6.67430e-11):
    """
    Calculate the gravitational force between two masses.
    m1, m2: masses in kg
    r: distance between the masses in meters
    G: gravitational constant
    Returns force in Newtons
    """
    return G * (m1 * m2) / r**2

if __name__ == "__main__":
    # Example usage
    mass1 = 5.972e24  # Mass of Earth in kg
    mass2 = 7.348e22  # Mass of Moon in kg
    distance = 3.844e8  # Distance between Earth and Moon in meters

    force = gravitational_force(mass1, mass2, distance)
    print(f"The gravitational force between the Earth and the Moon is {force:.2e} N")


import numpy as np
import matplotlib.pyplot as plt

r = np.linspace(1e6, 1e9, 100)
force = gravitational_force(5.972e24, 7.347e22, r)
plt.plot(r, force)
plt.title('Gravitational Force vs. Distance')
plt.xlabel('Distance (m)')
plt.ylabel('Force (N)')
plt.grid(True)
plt.show()
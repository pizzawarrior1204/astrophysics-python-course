import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Specify FFmpeg path (adjust based on `which ffmpeg`)
import matplotlib
# matplotlib.rcParams['animation.ffmpeg_path'] = '/opt/anaconda3/bin/ffmpeg'  # Or /usr/local/bin/ffmpeg
matplotlib.rcParams['animation.ffmpeg_path'] = '/opt/homebrew/bin/ffmpeg'  # Or /usr/local/bin/ffmpeg

def animate_orbit(a=1.496e9, e=0.0167, n=100, frames=100):
    """
    Animate planet on elliptical orbit and save as MP4.
    a: semi-major axis (meters)
    e: eccentricity
    n: points for orbit
    frames: animation frames
    """
    fig, ax = plt.subplots()
    theta = np.linspace(0, 2 * np.pi, n)
    r = a * (1 - e**2) / (1 + e * np.cos(theta))
    x = r * np.cos(theta)
    y = r * np.sin(theta)
    ax.plot(x, y, 'b-', label='Orbit')
    ax.scatter([0], [0], color='red', label='Star')
    planet, = ax.plot([], [], 'go', label='Planet')
    ax.set_title('Animated Elliptical Orbit')
    ax.set_xlabel('x (m)')
    ax.set_ylabel('y (m)')
    ax.axis('equal')
    ax.legend()
    ax.grid(True)

    def update(frame):
        theta_p = 2 * np.pi * frame / frames
        r_p = a * (1 - e**2) / (1 + e * np.cos(theta_p))
        x_p = r_p * np.cos(theta_p)
        y_p = r_p * np.sin(theta_p)
        planet.set_data([x_p], [y_p])
        return planet,

    ani = FuncAnimation(fig, update, frames=frames, interval=50, blit=False)
    ani.save('orbit.mp4', writer='ffmpeg', fps=20, dpi=100)
    # plt.show()  # Uncomment to display

if __name__ == "__main__":
    animate_orbit(e=0.5)
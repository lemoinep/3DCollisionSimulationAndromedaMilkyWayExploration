import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

G = 1.0  # Gravitational constant in arbitrary units
# Note: Units are normalized for simplicity; real astrophysical units require proper scaling.

class Galaxy:
    def __init__(self, center_mass, center_pos, center_vel, n_stars, radius):
        """
        Initialize a galaxy model consisting of a massive central core and orbiting stars.

        Parameters:
        - center_mass: Mass of the galactic core (dominant gravitational source).
        - center_pos: Initial 3D position vector of the galactic center.
        - center_vel: Initial 3D velocity vector of the galactic center.
        - n_stars: Number of stars orbiting the galactic center.
        - radius: Approximate radius of the stellar disk (defines spatial extent).

        The galaxy is simplified as a point mass core plus stars distributed in a disk.
        """
        self.center_mass = center_mass
        self.center_pos = np.array(center_pos, dtype=float)
        self.center_vel = np.array(center_vel, dtype=float)
        self.n_stars = n_stars
        self.radius = radius
        
        # Assign equal mass to each star for simplicity.
        # This assumes the total stellar mass is roughly equal to the central mass, 
        # which is a simplification.
        self.star_mass = center_mass / n_stars
        
        self.stars_pos = []
        self.stars_vel = []
        self.init_stars()

    def init_stars(self):
        """
        Initialize star positions and velocities distributed in a thin disk.
        
        Positions:
        - Stars are placed in a disk with radius 'self.radius'.
        - Radial distribution uses sqrt(random) to ensure uniform surface density.
        - Vertical distribution is Gaussian with small scale height (~5% of radius) 
          to simulate the disk thickness.
        
        Velocities:
        - Stars are given circular orbital velocities assuming a Keplerian potential 
          dominated by the central mass.
        - Velocity magnitude v = sqrt(G * M / r), where M is central mass and r is radius.
        - Velocity direction is tangential to the radius vector in the disk plane (xy-plane).
        - Stars inherit the velocity of the galactic center to maintain relative motion.
        """
        for _ in range(self.n_stars):
            # Radial distance with sqrt for uniform surface density in disk
            r = self.radius * np.sqrt(np.random.rand())
            theta = 2 * np.pi * np.random.rand()
            # Small vertical displacement to simulate disk thickness
            z = np.random.normal(0, self.radius * 0.05)

            # Convert polar to Cartesian coordinates in disk plane
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            pos = self.center_pos + np.array([x, y, z])

            # Calculate circular orbital velocity magnitude (Keplerian approximation)
            if r > 1e-5:
                v_mag = np.sqrt(G * self.center_mass / r)
            else:
                v_mag = 0.0

            # Velocity direction is perpendicular to radius vector in xy-plane
            vel_dir = np.array([-y, x, 0])
            norm = np.linalg.norm(vel_dir)
            if norm > 0:
                vel_dir /= norm
            else:
                vel_dir = np.zeros(3)

            # Total velocity is center velocity plus orbital velocity
            vel = self.center_vel + v_mag * vel_dir

            self.stars_pos.append(pos)
            self.stars_vel.append(vel)

        self.stars_pos = np.array(self.stars_pos)
        self.stars_vel = np.array(self.stars_vel)

def compute_gravity(pos1, mass1, pos2, mass2):
    """
    Compute the gravitational force exerted on mass1 at pos1 by mass2 at pos2.

    Uses Newton's law of universal gravitation:
    F = G * m1 * m2 / r^2 * r_hat

    Parameters:
    - pos1, pos2: Position vectors of the two masses.
    - mass1, mass2: Scalar masses.

    Returns:
    - Force vector acting on mass1 due to mass2.

    Notes:
    - A small distance threshold avoids singularities when r ~ 0.
    """
    r_vec = pos2 - pos1
    r = np.linalg.norm(r_vec)
    if r < 1e-5:
        return np.zeros(3)  # Avoid infinite force at zero distance
    force_mag = G * mass1 * mass2 / r**2
    force_dir = r_vec / r
    return force_mag * force_dir

def update_simulation(gal1, gal2, dt):
    """
    Update positions and velocities of the galaxies and their stars over a timestep dt.

    Procedure:
    1. Compute gravitational forces between the two galactic centers.
    2. Update the velocities and positions of the galactic centers (Euler integration).
    3. For each star in each galaxy:
       - Compute gravitational forces from both galactic centers.
       - Update star velocities and positions accordingly.

    Assumptions:
    - Stars do not gravitationally interact with each other (collisionless approximation).
    - Dominant gravitational influence on stars comes from the two galactic centers.
    - This reduces computational complexity while capturing main dynamics.
    """
    # Compute forces between galactic centers (point masses)
    force_on_1 = compute_gravity(gal1.center_pos, gal1.center_mass, gal2.center_pos, gal2.center_mass)
    force_on_2 = -force_on_1  # Newton's third law

    # Update galactic center velocities and positions (simple Euler method)
    gal1.center_vel += force_on_1 / gal1.center_mass * dt
    gal2.center_vel += force_on_2 / gal2.center_mass * dt
    gal1.center_pos += gal1.center_vel * dt
    gal2.center_pos += gal2.center_vel * dt

    # Update stars in galaxy 1
    for i in range(gal1.n_stars):
        f1 = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal1.center_pos, gal1.center_mass)
        f2 = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal2.center_pos, gal2.center_mass)
        total_force = f1 + f2
        a = total_force / gal1.star_mass  # acceleration = force / mass
        gal1.stars_vel[i] += a * dt
        gal1.stars_pos[i] += gal1.stars_vel[i] * dt

    # Update stars in galaxy 2
    for i in range(gal2.n_stars):
        f1 = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal1.center_pos, gal1.center_mass)
        f2 = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal2.center_pos, gal2.center_mass)
        total_force = f1 + f2
        a = total_force / gal2.star_mass
        gal2.stars_vel[i] += a * dt
        gal2.stars_pos[i] += gal2.stars_vel[i] * dt

# Initialize galaxies with parameters reflecting Andromeda larger than Milky Way

# Approximate radii (arbitrary units):
# Andromeda disk radius ~44 units (~220,000 ly scaled)
# Milky Way disk radius ~20 units (~100,000 ly scaled)

# Masses in arbitrary units (adjusted for demonstration):
# Andromeda mass ~4.5e5
# Milky Way mass ~1e6
# Note: Masses can be tuned; here Milky Way is more massive but smaller radius to illustrate dynamics.

galaxy1 = Galaxy(center_mass=4.5e5, center_pos=[-150, 0, 0], center_vel=[0, 0.4, 0], n_stars=2000, radius=44)  # Andromeda
galaxy2 = Galaxy(center_mass=1e6, center_pos=[150, 0, 0], center_vel=[0, -0.4, 0], n_stars=2000, radius=20)   # Milky Way

dt = 0.01  # Time step for numerical integration
steps = 200  # Number of simulation steps

# Setup 3D plot for visualization
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-200, 200)
ax.set_ylim(-200, 200)
ax.set_zlim(-60, 60)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_title("3D Collision Simulation: Andromeda (larger) - Milky Way")

# Initial scatter plots for stars and galactic centers
stars1_scatter = ax.scatter(galaxy1.stars_pos[:, 0], galaxy1.stars_pos[:, 1], galaxy1.stars_pos[:, 2], s=1, color='blue', label='Andromeda')
stars2_scatter = ax.scatter(galaxy2.stars_pos[:, 0], galaxy2.stars_pos[:, 1], galaxy2.stars_pos[:, 2], s=1, color='red', label='Milky Way')
center1_scatter = ax.scatter(galaxy1.center_pos[0], galaxy1.center_pos[1], galaxy1.center_pos[2], s=50, color='cyan', label='Center Andromeda')
center2_scatter = ax.scatter(galaxy2.center_pos[0], galaxy2.center_pos[1], galaxy2.center_pos[2], s=50, color='magenta', label='Center Milky Way')
ax.legend()

def animate(frame):
    """
    Animation update function called at each frame.
    Advances the simulation by one timestep and updates scatter plot data.
    """
    update_simulation(galaxy1, galaxy2, dt)
    stars1_scatter._offsets3d = (galaxy1.stars_pos[:, 0], galaxy1.stars_pos[:, 1], galaxy1.stars_pos[:, 2])
    stars2_scatter._offsets3d = (galaxy2.stars_pos[:, 0], galaxy2.stars_pos[:, 1], galaxy2.stars_pos[:, 2])
    center1_scatter._offsets3d = (np.array([galaxy1.center_pos[0]]), np.array([galaxy1.center_pos[1]]), np.array([galaxy1.center_pos[2]]))
    center2_scatter._offsets3d = (np.array([galaxy2.center_pos[0]]), np.array([galaxy2.center_pos[1]]), np.array([galaxy2.center_pos[2]]))
    return stars1_scatter, stars2_scatter, center1_scatter, center2_scatter

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=30, blit=False)

plt.show()

# Author(s): Dr. Patrick Lemoine

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# --- Physical constants and unit conversions ---
G_phys = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
Msun = 1.989e30       # Solar mass in kg
kpc = 3.086e19        # Kiloparsec in meters
Myr = 3.154e13        # Megayear in seconds

# Gravitational constant in kpc^3 / (Msun * Myr^2)
G = G_phys * (Msun**-1) * (kpc**3) * (Myr**-2)

def rotate_3d(vecs, axis, angle_deg):
    """
    Rotate a set of 3D vectors around an arbitrary axis by a given angle (degrees).
    Used to set the inclination of galactic disks.
    """
    angle = np.deg2rad(angle_deg)
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)
    ux, uy, uz = axis
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    R = np.array([
        [cos_a + ux**2*(1 - cos_a), ux*uy*(1-cos_a) - uz*sin_a, ux*uz*(1-cos_a) + uy*sin_a],
        [uy*ux*(1-cos_a) + uz*sin_a, cos_a + uy**2*(1-cos_a), uy*uz*(1-cos_a) - ux*sin_a],
        [uz*ux*(1-cos_a) - uy*sin_a, uz*uy*(1-cos_a) + ux*sin_a, cos_a + uz**2*(1-cos_a)]
    ])
    return vecs @ R.T

class Galaxy:
    """
    Represents a galaxy with a stellar disk and dark matter halo.
    The disk is initialized with stars in a quasi-equilibrium configuration.
    """
    def __init__(self, center_mass, dark_mass, center_pos, center_vel, 
                 n_stars, radius, inclination_deg=0, inclination_axis=[1,0,0]):
        self.center_mass = center_mass      # Visible (stellar) mass in solar masses
        self.dark_mass = dark_mass          # Dark matter halo mass in solar masses
        self.total_mass = center_mass + dark_mass
        self.center_pos = np.array(center_pos, dtype=float)    # Galactic center position [kpc]
        self.center_vel = np.array(center_vel, dtype=float)    # Galactic center velocity [kpc/Myr]
        self.n_stars = n_stars              # Number of stars in the disk (simulation particles)
        self.radius = radius                # Disk radius [kpc]
        self.star_mass = center_mass / n_stars
        self.inclination_deg = inclination_deg
        self.inclination_axis = inclination_axis
        self.stars_pos = []
        self.stars_vel = []
        self.init_stars()

    def init_stars(self):
        """
        Initialize the positions and velocities of stars in the disk.
        The distribution is exponential in radius and Gaussian in height.
        Velocities are set for approximate circular orbits.
        """
        pos_list = []
        vel_list = []
        for _ in range(self.n_stars):
            # Exponential disk in radius, Gaussian in z
            r = self.radius * (-0.5 * np.log(1 - np.random.rand()))**0.5
            theta = 2 * np.pi * np.random.rand()
            z = np.random.normal(0, self.radius * 0.05)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            pos = np.array([x, y, z])
            # Circular velocity (approximate, includes dark matter halo)
            v_mag = np.sqrt(G * self.total_mass * r**2 / (r**2 + (self.radius*0.2)**2)**1.5)
            vel_dir = np.array([-y, x, 0])  # Tangential direction in the disk plane
            norm = np.linalg.norm(vel_dir)
            vel_dir = vel_dir / norm if norm > 0 else np.zeros(3)
            vel = v_mag * vel_dir
            pos_list.append(pos)
            vel_list.append(vel)
        pos_arr = np.array(pos_list)
        vel_arr = np.array(vel_list)
        # Apply disk inclination
        if self.inclination_deg != 0:
            pos_arr = rotate_3d(pos_arr, self.inclination_axis, self.inclination_deg)
            vel_arr = rotate_3d(vel_arr, self.inclination_axis, self.inclination_deg)
        self.stars_pos = self.center_pos + pos_arr
        self.stars_vel = self.center_vel + vel_arr

def compute_gravity(pos1, mass1, pos2, mass2, softening=0.5):
    """
    Compute the gravitational force vector exerted by mass2 at pos2 on mass1 at pos1.
    Softening prevents singularities at small distances (Plummer softening).
    """
    r_vec = pos2 - pos1
    r = np.linalg.norm(r_vec)
    r_soft = np.sqrt(r**2 + softening**2)
    if r_soft < 1e-5:
        return np.zeros(3)
    force_mag = G * mass1 * mass2 / r_soft**2
    force_dir = r_vec / (r + 1e-10)
    return force_mag * force_dir

def update_all_galaxies(galaxies, dt):
    """
    Update the positions and velocities of all galaxies and their stars.
    Only the centers interact with each other (N-body), stars feel the centers.
    """
    # Update galaxy centers (mutual gravity)
    for i, gal1 in enumerate(galaxies):
        total_force = np.zeros(3)
        for j, gal2 in enumerate(galaxies):
            if i != j:
                total_force += compute_gravity(gal1.center_pos, gal1.total_mass, gal2.center_pos, gal2.total_mass)
        gal1.center_vel += total_force / gal1.total_mass * dt
        gal1.center_pos += gal1.center_vel * dt
    # Update stars (feel all galaxy centers)
    for gal in galaxies:
        for i in range(gal.n_stars):
            total_force = compute_gravity(gal.stars_pos[i], gal.star_mass, gal.center_pos, gal.total_mass)
            for other_gal in galaxies:
                if other_gal is not gal:
                    total_force += compute_gravity(gal.stars_pos[i], gal.star_mass, other_gal.center_pos, other_gal.total_mass)
            a = total_force / gal.star_mass
            gal.stars_vel[i] += a * dt
            gal.stars_pos[i] += gal.stars_vel[i] * dt

# --- Astrophysical parameters for major galaxies and satellites ---
m31_visible = 1.5e12     # Andromeda visible mass [Msun]
m31_dark = 1.5e13        # Andromeda dark matter halo [Msun]
mw_visible = 1.2e12      # Milky Way visible mass [Msun]
mw_dark = 1.2e13         # Milky Way dark matter halo [Msun]
sep = 800                # Initial separation between Milky Way and Andromeda [kpc]

# Andromeda (M31)
galaxy1 = Galaxy(
    center_mass=m31_visible,
    dark_mass=m31_dark,
    center_pos=[-sep/2, 50, -20],
    center_vel=[0.1226, 0.05, 0],
    n_stars=1000,
    radius=30,
    inclination_deg=77,
    inclination_axis=[1,0,0]
)

# Milky Way
galaxy2 = Galaxy(
    center_mass=mw_visible,
    dark_mass=mw_dark,
    center_pos=[sep/2, -30, 10],
    center_vel=[-0.118, -0.03, 0],
    n_stars=1000,
    radius=15,
    inclination_deg=45,
    inclination_axis=[0,1,0]
)

# Large Magellanic Cloud (LMC)
lmc = Galaxy(
    center_mass=1.4e11,
    dark_mass=1.2e12,
    center_pos=[-50, -40, 0],
    center_vel=[0, 0.3, 0],
    n_stars=300,
    radius=5,
    inclination_deg=35,
    inclination_axis=[0,1,0]
)

# Small Magellanic Cloud (SMC)
smc = Galaxy(
    center_mass=6.5e9,
    dark_mass=6e10,
    center_pos=[-63, -55, 0],
    center_vel=[0, 0.45, 0],
    n_stars=150,
    radius=3,
    inclination_deg=40,
    inclination_axis=[0,1,0]
)

# Triangulum Galaxy (M33)
m33 = Galaxy(
    center_mass=5e10,
    dark_mass=4.5e11,
    center_pos=[-850, 200, 50],
    center_vel=[0.05, 0.02, 0],
    n_stars=300,
    radius=7,
    inclination_deg=54,
    inclination_axis=[1,0,0]
)

# --- Major Milky Way satellites ---
# Each tuple: (name, visible mass, dark mass, pos_x, pos_y, pos_z, vel_x, vel_y, vel_z, n_stars, radius, color, inclination_deg, inclination_axis)
satellite_params = [
    ("SagDEG",   2.1e8,   2e9,     16, -2, -6,    0, 0.13, 0,     80,     2.5,   'brown',  20, [0,1,0]),
    ("CanMaj",   1e8,     1e9,    -13, -8, -2,    0, 0.12, 0,     50,     2.0,   'gold',   10, [1,0,0]),
    ("UrsaMin",  2e7,     2e8,    -50, 55, 50,    0, 0.10, 0,     40,     1.2,   'gray',   15, [0,1,0]),
    ("Sculptor", 3.6e7,   3e8,    -87, -13, -40,  0, 0.11, 0,     40,     1.5,   'olive',  20, [1,0,0]),
    ("Sextans",  5e6,     5e7,     95,  35,  50,  0, 0.09, 0,     30,     1.0,   'cyan',   10, [0,1,0]),
    ("Fornax",   7e7,     7e8,    -140, -90, 50,  0, 0.08, 0,     50,     2.0,   'pink',   25, [1,0,0]),
    ("LeoI",     5e7,     5e8,     250, 120, 70,  0, 0.07, 0,     40,     1.2,   'navy',   12, [0,1,0]),
    ("LeoII",    1.5e7,   1.5e8,   210, 140, 100, 0, 0.06, 0,     30,     1.0,   'purple', 10, [1,0,0]),
    ("Carina",   4e6,     4e7,    -110, 20, -20,  0, 0.08, 0,     30,     0.8,   'lime',   10, [0,1,0]),
    ("Draco",    2.9e7,   3e8,     80,  65,  50,  0, 0.09, 0,     30,     1.0,   'magenta',10, [1,0,0])
]

# Create Galaxy objects for each satellite
satellites = []
for sat in satellite_params:
    satellites.append(
        Galaxy(
            center_mass=sat[1],
            dark_mass=sat[2],
            center_pos=[sat[3], sat[4], sat[5]],
            center_vel=[sat[6], sat[7], sat[8]],
            n_stars=sat[9],
            radius=sat[10],
            inclination_deg=sat[12],
            inclination_axis=sat[13]
        )
    )

# Combine all galaxies into a single list for the simulation
galaxies = [galaxy1, galaxy2, lmc, smc, m33] + satellites

# --- Simulation parameters ---
dt = 0.01       # Time step in Myr
steps = 10000   # Number of simulation steps

# --- Visualization setup ---
fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1200, 1200)
ax.set_ylim(-1200, 1200)
ax.set_zlim(-1200, 1200)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.set_title("Simulation: Local Group and Milky Way Satellites")

# Prepare scatter plots for each galaxy
main_colors = ['blue', 'red', 'orange', 'green', 'purple']
sat_colors = [sat[11] for sat in satellite_params]
scatters = [
    ax.scatter([], [], [], s=1, color=col, alpha=0.5, label=lab)
    for col, lab in zip(main_colors, ['Andromeda', 'Milky Way', 'LMC', 'SMC', 'M33'])
]
scatters += [
    ax.scatter([], [], [], s=4, color=col, alpha=0.8, label=satellite_params[i][0])
    for i, col in enumerate(sat_colors)
]
ax.legend(loc='upper left', fontsize='small', markerscale=3)

def animate(frame):
    """
    Animation step: update all galaxies and plot their stars.
    """
    update_all_galaxies(galaxies, dt)
    for scatter, gal in zip(scatters, galaxies):
        scatter._offsets3d = (gal.stars_pos[:,0], gal.stars_pos[:,1], gal.stars_pos[:,2])
    return scatters

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=20, blit=False)
plt.show()

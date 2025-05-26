import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# --- Constantes physiques et conversions d'unités ---
G_phys = 6.67430e-11  # m^3 kg^-1 s^-2
Msun = 1.989e30       # kg
kpc = 3.086e19        # m
Myr = 3.154e13        # s

# Constante gravitationnelle en kpc³/(Msun·Myr²)
G = G_phys * (Msun**-1) * (kpc**3) * (Myr**-2)

def rotate_3d(vecs, axis, angle_deg):
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
    def __init__(self, center_mass, dark_mass, center_pos, center_vel, 
                 n_stars, radius, inclination_deg=0, inclination_axis=[1,0,0]):
        self.center_mass = center_mass
        self.dark_mass = dark_mass
        self.total_mass = center_mass + dark_mass
        self.center_pos = np.array(center_pos, dtype=float)
        self.center_vel = np.array(center_vel, dtype=float)
        self.n_stars = n_stars
        self.radius = radius
        self.star_mass = center_mass / n_stars
        self.inclination_deg = inclination_deg
        self.inclination_axis = inclination_axis
        self.stars_pos = []
        self.stars_vel = []
        self.init_stars()

    def init_stars(self):
        pos_list = []
        vel_list = []
        for _ in range(self.n_stars):
            r = self.radius * (-0.5 * np.log(1 - np.random.rand()))**0.5
            theta = 2 * np.pi * np.random.rand()
            z = np.random.normal(0, self.radius * 0.05)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            pos = np.array([x, y, z])
            v_mag = np.sqrt(G * self.total_mass * r**2 / (r**2 + (self.radius*0.2)**2)**1.5)
            vel_dir = np.array([-y, x, 0])
            norm = np.linalg.norm(vel_dir)
            vel_dir = vel_dir / norm if norm > 0 else np.zeros(3)
            vel = v_mag * vel_dir
            pos_list.append(pos)
            vel_list.append(vel)
        pos_arr = np.array(pos_list)
        vel_arr = np.array(vel_list)
        if self.inclination_deg != 0:
            pos_arr = rotate_3d(pos_arr, self.inclination_axis, self.inclination_deg)
            vel_arr = rotate_3d(vel_arr, self.inclination_axis, self.inclination_deg)
        self.stars_pos = self.center_pos + pos_arr
        self.stars_vel = self.center_vel + vel_arr

def compute_gravity(pos1, mass1, pos2, mass2, softening=0.5):
    r_vec = pos2 - pos1
    r = np.linalg.norm(r_vec)
    r_soft = np.sqrt(r**2 + softening**2)
    if r_soft < 1e-5:
        return np.zeros(3)
    force_mag = G * mass1 * mass2 / r_soft**2
    force_dir = r_vec / (r + 1e-10)
    return force_mag * force_dir

def update_all_galaxies(galaxies, dt):
    # Mise à jour des centres galactiques (forces principales entre galaxies)
    for i, gal1 in enumerate(galaxies):
        total_force = np.zeros(3)
        for j, gal2 in enumerate(galaxies):
            if i != j:
                total_force += compute_gravity(gal1.center_pos, gal1.total_mass, gal2.center_pos, gal2.total_mass)
        gal1.center_vel += total_force / gal1.total_mass * dt
        gal1.center_pos += gal1.center_vel * dt
    # Mise à jour des étoiles (influence des centres uniquement)
    for gal in galaxies:
        for i in range(gal.n_stars):
            total_force = compute_gravity(gal.stars_pos[i], gal.star_mass, gal.center_pos, gal.total_mass)
            for other_gal in galaxies:
                if other_gal is not gal:
                    total_force += compute_gravity(gal.stars_pos[i], gal.star_mass, other_gal.center_pos, other_gal.total_mass)
            a = total_force / gal.star_mass
            gal.stars_vel[i] += a * dt
            gal.stars_pos[i] += gal.stars_vel[i] * dt

# --- Paramètres astrophysiques réalistes ---
m31_visible = 1.5e12
m31_dark = 1.5e13
mw_visible = 1.2e12
mw_dark = 1.2e13

sep = 800  # kpc

# Andromède (M31)
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

# Voie Lactée
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

# Grand Nuage de Magellan (LMC)
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

# Petit Nuage de Magellan (SMC)
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

# Galaxie du Triangle (M33)
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

galaxies = [galaxy1, galaxy2, lmc, smc, m33]

dt = 0.01
steps = 10000

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1200, 1200)
ax.set_ylim(-1200, 1200)
ax.set_zlim(-1200, 1200)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.set_title("Simulation : Groupe Local avec LMC, SMC et M33")

# Prépare les scatters pour chaque galaxie
scatters = [
    ax.scatter([], [], [], s=1, color='blue', alpha=0.5, label='Andromède'),
    ax.scatter([], [], [], s=1, color='red', alpha=0.5, label='Voie Lactée'),
    ax.scatter([], [], [], s=2, color='orange', alpha=0.8, label='LMC'),
    ax.scatter([], [], [], s=2, color='green', alpha=0.8, label='SMC'),
    ax.scatter([], [], [], s=2, color='purple', alpha=0.8, label='M33')
]
ax.legend()

def animate(frame):
    update_all_galaxies(galaxies, dt)
    for scatter, gal in zip(scatters, galaxies):
        scatter._offsets3d = (gal.stars_pos[:,0], gal.stars_pos[:,1], gal.stars_pos[:,2])
    return scatters

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=20, blit=False)
plt.show()

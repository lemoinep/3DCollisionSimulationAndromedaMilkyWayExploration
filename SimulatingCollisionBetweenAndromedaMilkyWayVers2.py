# Author(s): Dr. Patrick Lemoine
# Correction of the inclination of the Andromede and Milky Way galaxy.

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

# --- Physical constants and unit conversions ---
G_phys = 6.67430e-11  # m^3 kg^-1 s^-2
Msun = 1.989e30       # kg
kpc = 3.086e19        # m
Myr = 3.154e13        # s

# Gravitational constant in kpc^3 / (Msun * Myr^2)
G = G_phys * (Msun**-1) * (kpc**3) * (Myr**-2)

def rotate_3d(vecs, axis, angle_deg):
    angle = np.deg2rad(angle_deg)
    axis = np.asarray(axis)
    axis = axis / np.linalg.norm(axis)
    ux, uy, uz = axis
    cos_a = np.cos(angle)
    sin_a = np.sin(angle)
    R = np.array([
        [cos_a + ux**2 * (1 - cos_a),      ux*uy*(1-cos_a) - uz*sin_a, ux*uz*(1-cos_a) + uy*sin_a],
        [uy*ux*(1-cos_a) + uz*sin_a, cos_a + uy**2*(1-cos_a),      uy*uz*(1-cos_a) - ux*sin_a],
        [uz*ux*(1-cos_a) - uy*sin_a, uz*uy*(1-cos_a) + ux*sin_a, cos_a + uz**2*(1-cos_a)]
    ])
    return vecs @ R.T

class Galaxy:
    def __init__(self, center_mass, center_pos, center_vel, n_stars, radius, inclination_deg=0, inclination_axis=[1,0,0]):
        self.center_mass = center_mass
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
            r = self.radius * np.sqrt(np.random.rand())
            theta = 2 * np.pi * np.random.rand()
            z = np.random.normal(0, self.radius * 0.05)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            pos = np.array([x, y, z])
            v_mag = np.sqrt(G * self.center_mass / (r + 0.1)) if r > 1e-5 else 0.0  # add 0.1 kpc softening
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

def compute_gravity(pos1, mass1, pos2, mass2, softening=1.0):
    r_vec = pos2 - pos1
    r = np.linalg.norm(r_vec)
    r_soft = np.sqrt(r**2 + softening**2)
    if r_soft < 1e-5:
        return np.zeros(3)
    force_mag = G * mass1 * mass2 / r_soft**2
    force_dir = r_vec / (r + 1e-10)
    return force_mag * force_dir

def update_simulation(gal1, gal2, dt):
    # Use a softening parameter to avoid infinite forces at close range
    softening = 2.0  # kpc
    force_on_1 = compute_gravity(gal1.center_pos, gal1.center_mass, gal2.center_pos, gal2.center_mass, softening)
    force_on_2 = -force_on_1
    gal1.center_vel += force_on_1 / gal1.center_mass * dt
    gal2.center_vel += force_on_2 / gal2.center_mass * dt
    gal1.center_pos += gal1.center_vel * dt
    gal2.center_pos += gal2.center_vel * dt
    for i in range(gal1.n_stars):
        f1 = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal1.center_pos, gal1.center_mass, softening)
        f2 = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal2.center_pos, gal2.center_mass, softening)
        a = (f1 + f2) / gal1.star_mass
        gal1.stars_vel[i] += a * dt
        gal1.stars_pos[i] += gal1.stars_vel[i] * dt
    for i in range(gal2.n_stars):
        f1 = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal1.center_pos, gal1.center_mass, softening)
        f2 = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal2.center_pos, gal2.center_mass, softening)
        a = (f1 + f2) / gal2.star_mass
        gal2.stars_vel[i] += a * dt
        gal2.stars_pos[i] += gal2.stars_vel[i] * dt

# --- Galaxy parameters in physical units ---
m31_mass = 1.5e12      # Msun
m31_radius = 30        # kpc
mw_mass = 1.5e12       # Msun
mw_radius = 15         # kpc

sep = 800  # kpc

galaxy1 = Galaxy(
    center_mass=m31_mass,
    center_pos=[-sep/2, 0, 0],
    center_vel=[0, 0, 0],
    n_stars=1000,  # fewer stars for speed and clarity
    radius=m31_radius,
    inclination_deg=77,
    inclination_axis=[1,0,0]
)
galaxy2 = Galaxy(
    center_mass=mw_mass,
    center_pos=[sep/2, 0, 0],
    center_vel=[0, 0, 0],
    n_stars=1000,
    radius=mw_radius,
    inclination_deg=0,
    inclination_axis=[1,0,0]
)

dt = 0.05   # Myr per step (smaller for stability)
steps = 600  # Number of simulation steps

fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-500, 500)
ax.set_ylim(-500, 500)
ax.set_zlim(-500, 500)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.set_title("3D Galaxy Collision: Andromeda (inclined) vs Milky Way [physical units]")

stars1_scatter = ax.scatter(galaxy1.stars_pos[:, 0], galaxy1.stars_pos[:, 1], galaxy1.stars_pos[:, 2], s=1, color='blue', label='Andromeda')
stars2_scatter = ax.scatter(galaxy2.stars_pos[:, 0], galaxy2.stars_pos[:, 1], galaxy2.stars_pos[:, 2], s=1, color='red', label='Milky Way')
center1_scatter = ax.scatter(galaxy1.center_pos[0], galaxy1.center_pos[1], galaxy1.center_pos[2], s=50, color='cyan', label='Center Andromeda')
center2_scatter = ax.scatter(galaxy2.center_pos[0], galaxy2.center_pos[1], galaxy2.center_pos[2], s=50, color='magenta', label='Center Milky Way')
ax.legend()

def animate(frame):
    update_simulation(galaxy1, galaxy2, dt)
    stars1_scatter._offsets3d = (galaxy1.stars_pos[:, 0], galaxy1.stars_pos[:, 1], galaxy1.stars_pos[:, 2])
    stars2_scatter._offsets3d = (galaxy2.stars_pos[:, 0], galaxy2.stars_pos[:, 1], galaxy2.stars_pos[:, 2])
    center1_scatter._offsets3d = (np.array([galaxy1.center_pos[0]]), np.array([galaxy1.center_pos[1]]), np.array([galaxy1.center_pos[2]]))
    center2_scatter._offsets3d = (np.array([galaxy2.center_pos[0]]), np.array([galaxy2.center_pos[1]]), np.array([galaxy2.center_pos[2]]))
    return stars1_scatter, stars2_scatter, center1_scatter, center2_scatter

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=30, blit=False)
plt.show()

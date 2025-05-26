import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

G_phys = 6.67430e-11
Msun = 1.989e30
kpc = 3.086e19
Myr = 3.154e13
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

def v_circ(M_enc, r):
    return np.sqrt(G * M_enc / (r + 1e-5))

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

    def mass_enclosed(self, r):
        Rd = self.radius / 3
        M_disk = self.center_mass
        M_halo = self.dark_mass
        M_enc_disk = M_disk * (1 - np.exp(-r/Rd)*(1 + r/Rd))
        Rh = self.radius
        M_enc_halo = M_halo * (r/(r+Rh))**3
        return M_enc_disk + M_enc_halo

    def init_stars(self):
        pos_list = []
        vel_list = []
        for _ in range(self.n_stars):
            r = -self.radius * np.log(1 - np.random.rand())
            theta = 2 * np.pi * np.random.rand()
            z = np.random.normal(0, self.radius * 0.05)
            x = r * np.cos(theta)
            y = r * np.sin(theta)
            pos = np.array([x, y, z])
            M_enc = self.mass_enclosed(np.abs(r))
            v_c = v_circ(M_enc, np.abs(r))
            # Ajout d'une dispersion de vitesse (stabilisation)
            sigma = 0.05 * v_c
            v_tan = v_c + np.random.normal(0, sigma)
            vel_dir = np.array([-y, x, 0])
            norm = np.linalg.norm(vel_dir)
            vel_dir = vel_dir / norm if norm > 0 else np.zeros(3)
            vel = v_tan * vel_dir
            # Dispersion radiale et verticale
            vel += np.random.normal(0, sigma, 3)
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

def update_simulation(gal1, gal2, dt):
    softening = 0.5
    friction_factor = 0.005
    force_on_1 = compute_gravity(gal1.center_pos, gal1.total_mass, gal2.center_pos, gal2.total_mass, softening)
    force_on_2 = -force_on_1
    gal1.center_vel += force_on_1 / gal1.total_mass * dt
    gal2.center_vel += force_on_2 / gal2.total_mass * dt
    gal1.center_vel *= (1 - friction_factor * dt)
    gal2.center_vel *= (1 - friction_factor * dt)
    gal1.center_pos += gal1.center_vel * dt
    gal2.center_pos += gal2.center_vel * dt
    for i in range(gal1.n_stars):
        f_self = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal1.center_pos, gal1.total_mass, softening)
        f_other = compute_gravity(gal1.stars_pos[i], gal1.star_mass, gal2.center_pos, gal2.total_mass, softening)
        a = (f_self + f_other) / gal1.star_mass
        gal1.stars_vel[i] += a * dt
        gal1.stars_pos[i] += gal1.stars_vel[i] * dt
    for i in range(gal2.n_stars):
        f_self = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal2.center_pos, gal2.total_mass, softening)
        f_other = compute_gravity(gal2.stars_pos[i], gal2.star_mass, gal1.center_pos, gal1.total_mass, softening)
        a = (f_self + f_other) / gal2.star_mass
        gal2.stars_vel[i] += a * dt
        gal2.stars_pos[i] += gal2.stars_vel[i] * dt

# Paramètres astrophysiques réalistes
m31_visible = 1.5e12
m31_dark = 1.5e13
mw_visible = 1.2e12
mw_dark = 1.2e13
sep = 800

galaxy1 = Galaxy(
    center_mass=m31_visible,
    dark_mass=m31_dark,
    center_pos=[-sep/2, 50, -20],
    #center_vel=[0.1226, 0.05, 0],
    center_vel=[0.25, 0.05, 0],
    n_stars=1000,
    radius=30,
    inclination_deg=77,
    inclination_axis=[1,0,0]
)
galaxy2 = Galaxy(
    center_mass=mw_visible,
    dark_mass=mw_dark,
    center_pos=[sep/2, -30, 10],
    #center_vel=[-0.118, -0.03, 0],
    center_vel=[-0.25, -0.03, 0],
    n_stars=1000,
    radius=15,
    inclination_deg=45,
    inclination_axis=[0,1,0]
)

dt = 0.001  # PAS DE TEMPS RÉDUIT
steps = 100

fig = plt.figure(figsize=(12, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim(-1200, 1200)
ax.set_ylim(-1200, 1200)
ax.set_zlim(-1200, 1200)
ax.set_xlabel('X [kpc]')
ax.set_ylabel('Y [kpc]')
ax.set_zlabel('Z [kpc]')
ax.set_title("Simulation réaliste : Collision Andromède-Voie Lactée")

stars1_scatter = ax.scatter([], [], [], s=1, color='blue', alpha=0.5, label='Andromède')
stars2_scatter = ax.scatter([], [], [], s=1, color='red', alpha=0.5, label='Voie Lactée')
center1_line, = ax.plot([], [], [], color='cyan', lw=2, label='Trajectoire Andromède')
center2_line, = ax.plot([], [], [], color='magenta', lw=2, label='Trajectoire Voie Lactée')
ax.legend()

center1_traj = [galaxy1.center_pos.copy()]
center2_traj = [galaxy2.center_pos.copy()]

def animate(frame):
    update_simulation(galaxy1, galaxy2, dt)
    stars1_scatter._offsets3d = (galaxy1.stars_pos[:,0], galaxy1.stars_pos[:,1], galaxy1.stars_pos[:,2])
    stars2_scatter._offsets3d = (galaxy2.stars_pos[:,0], galaxy2.stars_pos[:,1], galaxy2.stars_pos[:,2])
    center1_traj.append(galaxy1.center_pos.copy())
    center2_traj.append(galaxy2.center_pos.copy())
    x1, y1, z1 = np.array(center1_traj).T
    x2, y2, z2 = np.array(center2_traj).T
    center1_line.set_data_3d(x1, y1, z1)
    center2_line.set_data_3d(x2, y2, z2)
    return stars1_scatter, stars2_scatter, center1_line, center2_line

ani = animation.FuncAnimation(fig, animate, frames=steps, interval=20, blit=False)
plt.show()

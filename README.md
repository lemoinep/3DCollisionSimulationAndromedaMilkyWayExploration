# Simulation of the Andromeda-Milky Way Collision

## Objective

My objective is to explore various concepts in physics by formulating problems mathematically and conducting simulations. This approach will yield many interesting and valuable insights.

## Introduction

The Andromeda–Milky Way collision is a monumental cosmic event predicted to occur in about 4 to 4.5 billion years, involving the two largest galaxies in our Local Group: the Milky Way, which hosts our Solar System, and the neighboring Andromeda Galaxy. Driven by their mutual gravitational attraction, Andromeda is currently approaching the Milky Way at roughly 110 kilometers per second. Recent precise measurements from the Hubble Space Telescope have confirmed that their relative sideways motion is small enough to make a future collision inevitable.

Despite the immense scale of this encounter, individual stars within the galaxies are so widely spaced that direct stellar collisions are extremely unlikely. Instead, the gravitational forces will distort the structures of both galaxies over billions of years, triggering bursts of star formation and ultimately leading to their merger into a single, larger elliptical galaxy often nicknamed "Milkdromeda."

Simulations of this collision—like the one presented here—model the complex gravitational interplay between the galaxies’ massive cores and their vast populations of stars. By tracking the positions and velocities of stars orbiting each galactic center, these simulations reveal how the galaxies’ shapes warp, how stars are flung into new orbits, and how the two supermassive black holes at their centers may eventually merge. Such models help astronomers understand not only the fate of our cosmic neighborhood but also the general processes of galaxy formation and evolution throughout the universe.

## Summary of Scientific and Physical Concepts of the Simulation

- **Galaxy Modeling**: Simplified as a massive central core plus many stars orbiting in a disk.
- **Star Distribution**: Uniform surface density in the disk plane with small vertical thickness to mimic real spiral galaxies.
- **Orbital Velocities**: Calculated assuming stars orbit a point mass (central core) with circular velocity.

## Physical and Mechanical Equations of the Andromeda-Milky Way Collision

### Gravitational Force

The gravitational force between two point masses, such as the central cores of the galaxies, can be described by Newton's law of universal gravitation:

$$
F = G \\frac{m_1 m_2}{r^2}
$$

where:
- \( F \) is the gravitational force between the masses,
- \( G \) is the gravitational constant (\(6.67430 \\times 10^{-11} \\, \\text{m}^3 \\text{kg}^{-1} \\text{s}^{-2}\)),
- \( m_1 \) and \( m_2 \) are the masses of the two objects,
- \( r \) is the distance between the centers of the two masses.

### Orbital Velocity

The orbital velocity of stars within a galaxy can be approximated using the circular velocity formula derived from the balance between gravitational force and centripetal acceleration:

$$
v = \\sqrt{\\frac{GM}{r}}
$$

where:
- \( v \) is the orbital velocity of a star,
- \( G \) is the gravitational constant,
- \( M \) is the mass of the galactic core,
- \( r \) is the distance from the star to the galactic center.

### Potential Energy

The gravitational potential energy between two galaxies can be expressed as:

$$
U = -G \\frac{m_1 m_2}{r}
$$

This equation gives the potential energy of the system, which is crucial for understanding the energy dynamics during the collision.

### Dynamics of the Collision

The dynamics of the collision can be modeled using the principles of conservation of momentum and energy. For a simplified two-body problem:

- **Conservation of Momentum:**

$$
m_1 \\vec{v}_1 + m_2 \\vec{v}_2 = \\text{constant}
$$

- **Conservation of Energy:**

$$
\\frac{1}{2} m_1 v_1^2 + \\frac{1}{2} m_2 v_2^2 - G \\frac{m_1 m_2}{r} = \\text{constant}
$$

### Tidal Forces

Tidal forces play a significant role in distorting the shapes of the galaxies. The tidal force on a star in one galaxy due to the gravitational pull of the other galaxy can be approximated by:

$$
F_{\\text{tidal}} \\approx \\frac{2GMd}{r^3}
$$

where:
- \( d \) is the diameter of the star's orbit around its galaxy,
- \( r \) is the distance between the two galaxies.

### Simulation and Numerical Methods

In simulations, the Euler method is often used for numerical integration to update the positions and velocities of stars over time. The basic Euler integration steps are:

- **Position Update:**

$$
\\vec{r}(t + \\Delta t) = \\vec{r}(t) + \\vec{v}(t) \\Delta t
$$

- **Velocity Update:**

$$
\\vec{v}(t + \\Delta t) = \\vec{v}(t) + \\vec{a}(t) \\Delta t
$$

where:
- \(\vec{r}\) is the position vector,
- \(\vec{v}\) is the velocity vector,
- \(\vec{a}\) is the acceleration vector,
- \(\Delta t\) is the time step.


## Computational Aspects of the Simulation

### Initial Conditions

To accurately simulate the collision between the Andromeda Galaxy and the Milky Way, it is essential to set appropriate initial conditions that reflect current astronomical observations:

- **Mass Distribution**: Both galaxies are modeled with a central bulge and a disk. The mass of Andromeda is estimated to be slightly higher than that of the Milky Way.
- **Initial Velocities**: The relative velocity of Andromeda with respect to the Milky Way is set based on observational data, approximately 110 km/s.
- **Positions**: The initial separation between the two galaxies is set according to the current distance, approximately 2.5 million light-years.

### Simulation Parameters

- **Time Step**: The simulation uses a small time step to ensure stability and accuracy, typically in the range of thousands to millions of years per step.
- **Total Simulation Time**: The simulation runs for a total time sufficient to observe the complete merger, usually several billion years.
- **Resolution**: High resolution is crucial for capturing detailed interactions, especially in regions of high density such as galactic centers.

### Expected Outcomes

The simulation aims to provide insights into several key aspects of the collision:

- **Galactic Distortion**: As the galaxies approach, tidal forces will distort their spiral structures, creating long tidal tails and bridges of stars.
- **Star Formation**: The collision is expected to trigger bursts of star formation due to the compression of interstellar gas.
- **Merging of Supermassive Black Holes**: The central supermassive black holes of both galaxies are expected to eventually merge, emitting gravitational waves.

### Visualization

Visualization tools are employed to render the simulation data, providing intuitive insights into the complex dynamics:

- **Density Plots**: Visual representations of star density to show the evolution of galactic structures.
- **Velocity Fields**: Vector fields illustrating the velocity and direction of stars, highlighting regions of high dynamical activity.
- **3D Rendering**: Three-dimensional models to visualize the collision from multiple perspectives, enhancing the understanding of spatial relationships and movements.

## Implications of the Collision

### Astrophysical Insights

The simulation of the Andromeda-Milky Way collision offers profound insights into galaxy formation and evolution:

- **Galaxy Merger Dynamics**: Understanding the processes governing galaxy mergers, which are fundamental to the hierarchical model of galaxy formation.
- **Star and Gas Dynamics**: Insights into how stars and gas are redistributed during such events, influencing future star formation and galactic morphology.
- **Dark Matter Interaction**: Although not directly observable, the role of dark matter in mediating the collision dynamics can be inferred from the simulation.

### Future of the Solar System

While the collision will dramatically alter the structure of both galaxies, the fate of the Solar System is of particular interest:

- **Low Probability of Direct Impact**: Due to the vast distances between stars, direct collisions are unlikely, but the Solar System may be relocated within the new galaxy.
- **Changes in Orbit**: The gravitational perturbations may alter the orbit of the Solar System, potentially placing it in a different region of the merged galaxy.

### Conclusion

The simulation of the Andromeda-Milky Way collision is not only a fascinating exploration of our cosmic future but also a valuable tool for astrophysical research. By modeling this inevitable event, scientists can refine their understanding of galactic interactions, star formation, and the long-term evolution of the universe.

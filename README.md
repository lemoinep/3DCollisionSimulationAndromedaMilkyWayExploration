My objective is to explore various concepts in physics by formulating problems mathematically and conducting simulations. 
This approach will yield many interesting and valuable insights.


The Andromeda–Milky Way collision is a monumental cosmic event predicted to occur in about 4 to 4.5 billion years, involving the two largest galaxies in our Local Group: the Milky Way, which hosts our Solar System, and the neighboring Andromeda Galaxy. Driven by their mutual gravitational attraction, Andromeda is currently approaching the Milky Way at roughly 110 kilometers per second. Recent precise measurements from the Hubble Space Telescope have confirmed that their relative sideways motion is small enough to make a future collision inevitable.

Despite the immense scale of this encounter, individual stars within the galaxies are so widely spaced that direct stellar collisions are extremely unlikely. Instead, the gravitational forces will distort the structures of both galaxies over billions of years, triggering bursts of star formation and ultimately leading to their merger into a single, larger elliptical galaxy often nicknamed "Milkdromeda."

Simulations of this collision—like the one presented here—model the complex gravitational interplay between the galaxies’ massive cores and their vast populations of stars. By tracking the positions and velocities of stars orbiting each galactic center, these simulations reveal how the galaxies’ shapes warp, how stars are flung into new orbits, and how the two supermassive black holes at their centers may eventually merge. Such models help astronomers understand not only the fate of our cosmic neighborhood but also the general processes of galaxy formation and evolution throughout the universe.


### Summary of scientific and physical concepts of the simulation:

- **Galaxy modeling**: Simplified as a massive central core plus many stars orbiting in a disk.
- **Star distribution**: Uniform surface density in disk plane with small vertical thickness to mimic real spiral galaxies.
- **Orbital velocities**: Calculated assuming stars orbit a point mass (central core) with circular velocity $$ v = \sqrt{GM/r} $$.
- **Gravitational interactions**: Stars feel gravity from both galactic centers; stars do not interact mutually (collisionless approximation).
- **Numerical integration**: Euler method updates positions and velocities each timestep; simple but may require small dt for stability.
- **Units**: Arbitrary units used for demonstration; real simulations require astrophysical units and scaling.
- **Galaxy sizes and masses**: Andromeda is modeled with a larger disk radius than the Milky Way to reflect astrophysical observations.


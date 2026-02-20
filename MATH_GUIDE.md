# Math Guide: What the Simulation Is Computing

This document explains the math and physical meaning behind the MATLAB simulation files in this project. It is written for a college sophomore with calculus and basic physics.

## 1) Big Picture: What Are We Simulating?

We want the **intensity pattern** in a focal region after a beam passes through a circular aperture (or an annular mask) and a lens. The pattern depends on:

- The **beam** (wavelength, waist, Gaussian profile)
- The **aperture** (size, masked inner/outer radius)
- The **lens** (focal length, refractive index)
- **Aberrations** (spherical aberration modeled by a phase term)
- The **observation region** (a 2D grid in axial distance $z$ and radial distance $r$)

Mathematically, each point in the observation plane is computed by a **diffraction integral** that sums contributions from every ring on the aperture.

## 2) Coordinates and Normalization

The code uses normalized coordinates for stability and easier computation.

### Physical coordinates

- $z$: axial distance along the optical axis (meters)
- $r$: radial distance from the axis (meters)
- $\rho$: normalized pupil radius, $\rho \in [0,1]$

### Scaling factors

Two scaling constants convert between physical and normalized coordinates:

- Axial scale:
  $$z_{\text{scale}} = \frac{f^2}{k a^2}$$
- Radial scale:
  $$r_{\text{scale}} = \frac{f}{k a}$$

where:

- $f$ is the focal length
- $a$ is the aperture radius
- $k = \frac{2\pi}{\lambda}$ is the wave number
- $\lambda$ is the wavelength

The normalized coordinates are:

- $u = \frac{z}{z_{\text{scale}}}$
- $v = \frac{r}{r_{\text{scale}}}$

In code, these normalized grids are created by:

```matlab
u_values = linspace(z_min/z_scale, z_max/z_scale, grid_points_z);
v_values = linspace(r_min/r_scale, r_max/r_scale, grid_points_r);
```

## 3) The Core Diffraction Integral

The field at each point $(u, v)$ is computed using a radial diffraction integral:

$$
E(u, v) = \int_0^1 A(\rho)\,\exp\left(-i\left[k\,\phi(\rho) - \frac{u\rho^2}{2}\right]\right)\,J_0(v\rho)\,\rho\,d\rho
$$

The **intensity** is then:

$$
I(u, v) = |E(u, v)|^2
$$

In full form (Shvedov's diffraction integral), the intensity is:

$$
I(u, v) = \frac{S\pi a^4 P}{\lambda^2 f^2 w_0^2}\left|\int_0^1 e^{\frac{\rho^2}{w_0^2}} e^{-j(k\rho - \frac{u\rho^2}{2})} J_0(v\rho)\rho d\rho\right|^2
$$

where:
- $S$ is a source intensity factor
- $P$ is the incident power
- The exponential $e^{\frac{\rho^2}{w_0^2}}$ is the Gaussian apodization of the beam intensity
- The phase term includes both the spherical aberration $\phi(\rho)$ and the defocus term $\frac{u\rho^2}{2}$
- The Bessel function $J_0(v\rho)$ arises from the radial symmetry of the system

Each term has a physical meaning:

- $A(\rho)$: amplitude across the pupil (mask + Gaussian beam)
- $\phi(\rho)$: spherical aberration phase
- $\exp\left(-i\frac{u\rho^2}{2}\right)$: defocus phase term
- $J_0(v\rho)$: Bessel function from radial symmetry
- $\rho$: extra factor from polar coordinates

This integral is what every simulation function is evaluating, using different performance strategies.

## 4) What Each Variable Represents

All these values come from the parameter struct created in `create_simulation_params`.

### Physical parameters

- `lambda` ($\lambda$): wavelength of the light (meters)
- `a` ($a$): aperture radius (meters)
- `w0` ($w_0$): Gaussian beam waist (meters)
- `f` ($f$): lens focal length (meters)
- `n` ($n$): refractive index of the lens material
- `p` ($p$): position factor of the lens (dimensionless, describes lens position relative to object/image planes)
- `q` ($q$): shape factor of the lens (dimensionless, describes lens curvature)
- `k` ($k$): wave number $2\pi/\lambda$

### Grid parameters

- `z_min`, `z_max`: axial range to simulate (meters)
- `r_min`, `r_max`: radial range to simulate (meters)
- `grid_points_z`, `grid_points_r`: grid resolution in each direction

### Integration parameter

- `Nrho`: number of sample points for the integral over $\rho$

### Mask parameters

- `mask_radius_min`: inner radius of the mask (normalized)
- `mask_radius_max`: outer radius of the mask (normalized)

A full circular aperture uses `mask_radius_min = 0` and `mask_radius_max = 1`.

## 5) Amplitude Across the Aperture

The amplitude in the aperture plane is a product of three factors:

1. **Mask (annulus):** blocks inner/outer regions
2. **Gaussian beam profile:** reduces edges if the beam is smaller than the aperture
3. **Jacobian factor $\rho$:** comes from polar coordinates

In code:

```matlab
mask_vals = double((rho >= mask_radius_min) & (rho <= mask_radius_max));
gauss_vals = exp(-rho.^2 / (w0/a)^2);
base_vals = mask_vals .* gauss_vals .* rho;
```

Physical meaning:

- `mask_vals`: 1 where light passes, 0 where blocked
- `gauss_vals`: beam power falls off radially like $\exp(-(\rho a)^2 / w_0^2)$
- `rho`: gives more weight to outer rings because in polar coordinates, each ring has length proportional to $\rho$

## 6) Spherical Aberration: Shape and Position Factors

### The Shape Factor (q)

The shape factor describes the curvature of the lens surfaces:

$$q = \frac{R_2 + R_1}{R_2 - R_1}$$

where $R_1$ and $R_2$ are the radii of curvature of the two lens surfaces.

- For a plano-convex (PCX) lens with a flat side towards the source and collimated incident rays, this yields $q = 1$
- **Important:** Flipping the lens to have $q = -1$ produces much better light-dark contrast in the focal region
- Different values of $q$ produce different aberration signatures

### The Position Factor (p)

The position factor describes where the lens is positioned relative to the object and image planes:

$$p = 1 - \frac{2f}{S_2}$$

where $f$ is the focal length and $S_2$ is the distance from the lens to the image.

- For collimated input rays (object at infinity, image at $f$), this simplifies to $p = s = -1$
- The position factor affects the type and magnitude of spherical aberration introduced

## 7) Spherical Aberration Phase

The code uses a 4th-order phase term (Seidel aberration model that depends on both $p$ and $q$):

$$
\phi(\rho) = -\frac{(a\rho)^4}{32 n f^3 (n-1)}\left[\frac{n+2}{n-1}q^2 + 4\frac{n+1}{1}pq + (3n+2)(n-1)p^2 + \frac{n^3}{n-1}\right]
$$

In code, this is computed as:

```matlab
phi_vals = -((a*rho).^4 / (32*n*f^3*(n-1))) * ...
           (((n+2)/(n-1))*q^2 + (4*(n+1)*p*q) + ...
            ((3*n+2)*(n-1)*p^2) + ((n^3)/(n-1)));
```

Physical meaning:

- A perfect lens would have $\phi(\rho) = 0$ (flat phase)
- This term adds a phase error that grows as $\rho^4$ (worse at the edges)
- The parameters `p` and `q` control how strong and what type of aberration is modeled
- The coefficient structure shows how shape ($q$) and position ($p$) combine to produce the total spherical aberration

## 8) The Defocus Term

The term $\exp\left(-i\frac{u\rho^2}{2}\right)$ represents **defocus** along the axial direction. Larger $u$ means the observation point is farther from focus, causing more phase variation across the pupil and changing the intensity pattern.

In code:

```matlab
phase_u = exp(-1i * (k * phi_vals - u_val * rho.^2 / 2));
```

This combines both the aberration phase and the defocus phase.

## 9) The Bessel Function Term

The Bessel function $J_0(v\rho)$ appears because the system is radially symmetric and we integrate in polar coordinates. It describes how ring-shaped contributions interfere at radius $v$.

In code:

```matlab
J0_vals = besselj(0, v_val * rho);
```

This is the key term that shapes the radial structure of the focal spot.

## 10) Numerical Integration Over $\rho$

The integral is evaluated numerically using the trapezoidal rule:

```matlab
field = trapz(rho, common .* J0_vals);
```

- `rho` is a uniform grid between 0 and 1
- `common` combines the amplitude and phase terms
- `trapz` approximates the integral by summing trapezoids

The result is a complex field, and the intensity is its magnitude squared:

```matlab
col_intensity(iv) = abs(field)^2;
```

## 11) Normalization

Most CPU methods normalize the intensity to a maximum of 1:

```matlab
intensity = intensity / max(intensity(:));
```

This makes plots comparable across runs. The GPU ultra method intentionally leaves raw values (normalization is commented out) so you can decide how to scale afterward.

## 12) How Each File Uses the Math

### run_simulation.m (reference, slow)

- Loops over every $u,v$ point
- Evaluates the integral using MATLAB `integral` for each point
- Best for correctness checks, not speed

### run_simulation_fast.m (vectorized + parallel)

- Precomputes $\rho$ terms once
- Uses `parfor` over $u$ to compute columns in parallel
- Uses `trapz` for fast numerical integration

### run_simulation_gpu_ultra.m (full GPU)

- Builds a full 3D array with dimensions $(r, z, \rho)$ on the GPU
- Computes all Bessel and phase terms in parallel
- Integrates over $\rho$ using a vectorized trapezoid sum
- Extremely fast but memory heavy

## 13) Why the Parameters Matter

- Increasing `grid_points_z` or `grid_points_r` gives smoother plots but increases runtime
- Increasing `Nrho` improves integral accuracy but increases runtime
- Larger `a` or smaller `w0` increases edge weighting and can strengthen rings
- Changing `lambda` changes the scale of diffraction and the size of the focal spot
- Spherical aberration parameters (`p`, `q`) can blur or distort the focus
  - When using a PCX lens with collimated input, try `q = -1` instead of `q = 1` for better light-dark contrast in the focal plane

## 14) Typical Workflow (Math Perspective)

1. **Define physical parameters** (wavelength, aperture, lens)
2. **Define observation region** (ranges in $z$ and $r$)
3. **Compute normalized grids** ($u, v$)
4. **For each $u, v$** evaluate the diffraction integral
5. **Square magnitude** to get intensity
6. **Normalize and plot**

## 15) Quick Glossary

- **Pupil function:** The amplitude and phase pattern across the aperture
- **Diffraction integral:** The sum of contributions from all aperture points
- **Bessel function $J_0$:** Arises from circular symmetry
- **Spherical aberration:** Phase error that grows with radius, controlled by shape factor $q$ and position factor $p$
- **Defocus:** How far you are from the focal plane
- **Shape factor ($q$):** Describes the lens geometry; $q = -1$ provides better contrast than $q = 1$ for PCX lenses
- **Position factor ($p$):** Describes the lens position in the optical system

## 16) If You Want to Go Deeper

Good keywords to search in textbooks or online:

- “Debye integral”
- “Richards-Wolf diffraction”
- “Hankel transform in optics”
- “Seidel aberrations”

These are the theoretical foundations of the integral used here.

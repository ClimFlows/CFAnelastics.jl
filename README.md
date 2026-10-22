# CFAnelastics

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://dubosipsl.github.io/CFAnelastics.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://dubosipsl.github.io/CFAnelastics.jl/dev/)
[![Build Status](https://github.com/dubosipsl/CFAnelastics.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/dubosipsl/CFAnelastics.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/dubosipsl/CFAnelastics.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/dubosipsl/CFAnelastics.jl)

Anelastic binary fluid (a fluid with a conservative temperature-like variable `consvar` and a composition `q`) in a
2D box (x, z), discretized with energy-consistent finite volumes on a C-grid. Part of the
[ClimFlows](https://github.com/ClimFlows) ecosystem; sibling of
[CFCompressible](https://github.com/ClimFlows/CFCompressible.jl), whose box model (FC2D) it mirrors.

## Model

`AN2D` advances the tracers `m = ρ_R (consvar, q)` at cell centers and the velocity `(u, w)` on the faces, with
`ρ_R(z)` a hydrostatic reference density, under the constraint `∇·(ρ_R u) = 0`. The momentum equation has
- a vortex force (`VortexEnergyCons`, `VortexEnstrophyCons` or `VortexEnergyEnstrophyCons`, from CFDiffusionSchemes),
- a thermodynamic or buoyancy forcing, selected by a `BuoyancyScheme`: `ThermodynamicForcing` (default, energetically
  consistent, has a horizontal component), `ConsistentBuoyancy` and `LegacyBuoyancy`,
- viscosity (constant dynamic viscosity, Stokes hypothesis),
- the anelastic pressure, solved at every tendency evaluation by a fast Poisson solver (`PoissonSolver`: Fourier or
  cosine transforms horizontally, tridiagonal solves vertically).

Tracers diffuse by the flux laws of CFDiffusionSchemes, including barodiffusion. Walls (free- or no-slip, with
prescribed fluxes or values of the tracers) and periodic directions are set with `BoxBoundaryConditions`.

## Usage

```julia
using CFAnelastics, CFBoxes, CFDiffusionSchemes, CFTimeSchemes, ClimFluids
using CFPlanets: Tank2D
using LoopManagers: MultiThread, VectorizedCPU

space = Tank2D(Lx = 1.0, Lz = 0.5, g = 10.0)
fluid = NonlinearBinaryFluid(:potential_temperature, (; p0 = 1e5, T0 = 280.0, q0 = 0.035, Cp0 = 4000.0, v0 = 1e-3,
    α_T = 2e-4, α_q = 0.8, α_p = 1e-6, M_s = 3.14038218e-2, R = 8.31446261815324, α_TT = 1e-5, γ = 1e-7))
domain = Box2D(; Mx = 32, Mz = 16, boundary = (Periodic, Bounded))
wall = Wall(tangential = FreeSlip(), tracers = (Flux(0.0), Flux(0.0)))     # closed, insulating walls
model = AN2D((; mgr = MultiThread(VectorizedCPU()), domain, space, fluid,
    spatial_scheme = Centred2(), vortex_scheme = VortexEnergyCons(), buoyancy_scheme = ThermodynamicForcing(),
    thermoflux_scheme = ThermoFluxConsistent(k_T = 1e-4, k_q = 1e-5), viscosity_scheme = ViscosityScheme(dyn_visc = 1e-3),
    anelastic_reference = AnelasticReference(domain, space, 1000.0; pref = 1e5),
    boundary_conditions = BoxBoundaryConditions(domain; bottom = wall, top = wall), gauge_offset = 0.0))

state = CFAnelastics.initialize(model; consvar = (x, z) -> 280 + z / 0.5 + 0.3 * sin(2π * x), q = (x, z) -> 0.035)
state = CFAnelastics.loop(state, model, CFTimeSchemes.RungeKutta4(model),
                          (; Nslice = 2, slice_size = 1.0, dt_max = 0.1, cfl = 0.8))
```

`CFTimeSchemes.tendencies!(dstate, scratch, model, state, t)` returns `(dstate, scratch)` and, once `dstate` and
`scratch` are allocated (pass `void` the first time), does not allocate.

## Tests

```julia
using Pkg; Pkg.test("CFAnelastics")
```

The tests check the Poisson solver, then the tendencies: no dependence on uninitialised memory or halo values,
impermeable walls, conservation of mass and tracers, semi-discrete energy conservation to round-off (dE/dt equals the
wall input with forcing walls), absence of allocations, and the anelastic constraint.

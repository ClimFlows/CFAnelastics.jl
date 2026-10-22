```@meta
CurrentModule = CFAnelastics
```

# CFAnelastics

Anelastic binary fluid in a 2D box, discretized with energy-consistent finite volumes on a C-grid, in the
[ClimFlows](https://github.com/ClimFlows) ecosystem. The source is at
[CFAnelastics.jl](https://github.com/dubosipsl/CFAnelastics.jl).

## The model

[`AN2D`](@ref) advances the tracers ``m = \rho_R (\mathrm{consvar}, q)`` at cell centers and the velocity
``(u, w)`` on the faces, where ``\rho_R(z)`` is a hydrostatic reference density
([`AnelasticReference`](@ref)), under the constraint ``\nabla\cdot(\rho_R \mathbf{u}) = 0``. The pressure
``\varphi = (p - p_R)/\rho_R`` that enforces it is found at each evaluation of the tendencies by the
[`PoissonSolver`](@ref), with impermeable walls and periodic or bounded horizontal directions.

The forcing of the momentum is chosen by a [`BuoyancyScheme`](@ref):

| scheme                  | momentum forcing                                                                   |
|:------------------------|:-----------------------------------------------------------------------------------|
| `ThermodynamicForcing`  | ``-\mathrm{consvar}\,\nabla\pi - q\,\nabla\mu - \nabla B - g\hat z``, with a horizontal part |
| `ConsistentBuoyancy`    | along ``z`` only, ``b = g(\rho_R v - 1)``                                          |
| `LegacyBuoyancy`        | along ``z`` only, ``b = g(1 - 1/(\rho_R v))``                                      |

## Usage

```julia
using CFAnelastics, CFBoxes, CFDiffusionSchemes, CFTimeSchemes, ClimFluids
using CFPlanets: Tank2D
using LoopManagers: MultiThread, VectorizedCPU

space = Tank2D(Lx = 1.0, Lz = 0.5, g = 10.0)
fluid = NonlinearBinaryFluid(:potential_temperature, (; p0 = 1e5, T0 = 280.0, q0 = 0.035, Cp0 = 4000.0, v0 = 1e-3,
    α_T = 2e-4, α_q = 0.8, α_p = 1e-6, M_s = 3.14038218e-2, R = 8.31446261815324, α_TT = 1e-5, γ = 1e-7))
domain = Box2D(; Mx = 32, Mz = 16, boundary = (Periodic, Bounded))
wall = Wall(tangential = FreeSlip(), tracers = (Flux(0.0), Flux(0.0)))
model = AN2D((; mgr = MultiThread(VectorizedCPU()), domain, space, fluid,
    spatial_scheme = Centred2(), vortex_scheme = VortexEnergyCons(), buoyancy_scheme = ThermodynamicForcing(),
    thermoflux_scheme = ThermoFluxConsistent(k_T = 1e-4, k_q = 1e-5), viscosity_scheme = ViscosityScheme(dyn_visc = 1e-3),
    anelastic_reference = AnelasticReference(domain, space, 1000.0; pref = 1e5),
    boundary_conditions = BoxBoundaryConditions(domain; bottom = wall, top = wall), gauge_offset = 0.0))

state = CFAnelastics.initialize(model; consvar = (x, z) -> 280 + z / 0.5 + 0.3 * sin(2π * x), q = (x, z) -> 0.035)
state = CFAnelastics.loop(state, model, CFTimeSchemes.RungeKutta4(model),
                          (; Nslice = 2, slice_size = 1.0, dt_max = 0.1, cfl = 0.8))
```

The tendencies are `CFTimeSchemes.tendencies!(dstate, scratch, model, state, t)`, which returns
`(dstate, scratch)`; both are allocated when `void`, and the call does not allocate once they are.

## API

```@index
```

```@autodocs
Modules = [CFAnelastics]
```

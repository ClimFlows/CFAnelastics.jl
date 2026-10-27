using CFAnelastics
using CFTimeSchemes: tendencies!

# AN2D configuration and states

# reference density increasing with depth Lz - z
an_density(z) = ρ0 * (1 + 0.1 * (Lz - z) / Lz)

#=
32 × 16 box with every physical term on: buoyancy, vortex force, thermodynamic diffusion and viscosity, using the
energy-conserving schemes
=#
function an_model(; mgr = mgr, gauge_offset = 0.0, boundary = (Periodic, Bounded), walls = (; bottom = closed_wall(), top = closed_wall()))
    domain = Box2D(; Mx = 32, Mz = 16, boundary)
    anelastic_reference = AnelasticReference(domain, space, an_density; pref = 1e5)
    boundary_conditions = BoxBoundaryConditions(domain; walls...)
    return AN2D((; mgr, domain, space, fluid, spatial_scheme = Centred2(), vortex_scheme = VortexEnergyCons(),
        buoyancy_scheme = ThermodynamicForcing(), thermoflux_scheme = thermoflux, viscosity_scheme = ViscosityScheme(dyn_visc = 1e-3),
        anelastic_reference, boundary_conditions, gauge_offset))
end

# the model's tendencies!, in the argument order of the generic checks (helpers.jl)
an_tend!(dstate, scratch, model, state) = tendencies!(dstate, scratch, model, state, 0.0)

#=
active state (m = ρ_R(θ, q), u, w), with no symmetry, so that every term of the tendencies is non-zero; velocities
from a mass streamfunction ψ vanishing on the walls, so that ∇·(ρ_R u) = 0 discretely
=#
function an_state(model)
    # with walls in x, ψ also vanishes on x = 0, Lx
    wall_x(X) = boundary_topology(model.domain)[1] isa Bounded ? sin(π * X / Lx) : 1.0
    ψ(X, Z) = ρ0 * U * Lz * wall_x(X) * (sin(π * Z / Lz) * (1 + 0.5 * cos(2π * X / Lx + 0.3)) +
                                        0.3 * sin(2π * Z / Lz) * sin(4π * X / Lx))
    return CFAnelastics.initialize(model; consvar = θ_active, q = q_active,
                                   velocity = CFAnelastics.divergence_free_velocity(model, ψ))
end

# per-cell energy ΔxΔz [ρ_R K + ρ_R (h(p_R, s, q) + Φ) - p_R],  K = (ℑ(u², (1,), (Cell(),)) + ℑ(w², (2,), (Cell(),))) / 2
function an_cell_energies(model, (; m, u, w))
    (; domain, dx, dz, fluid, reference) = model
    (; ρ_R, p_R) = reference
    ic, jc = interior_xz(domain)
    return [begin
        K = (u[i, j]^2 + u[east(domain, i), j]^2 + w[i, j]^2 + w[i, j+1]^2) / 4
        h = fluid(:p, :consvar, :q).specific_enthalpy(p_R[j], m[i, j, 1] / ρ_R[j], m[i, j, 2] / ρ_R[j])
        Φ = g * (cellpoint(domain, dz, 2, j) - Lz)
        (ρ_R[j] * (K + h + Φ) - p_R[j]) * dx * dz
    end for i in ic, j in jc]
end

# walls exercising every kind of condition (consvar = θ for this fluid)
const forced_box = (;
    left   = Wall(tangential = TangentialVelocity(z -> 0.02 * z), tracers = (Value(T0 + 0.5), Value(q0))),
    right  = Wall(tangential = NoSlip(), tracers = (Value(z -> T0 - 0.5 + 0.1 * z), Flux(1e-6))),
    bottom = Wall(tangential = TangentialVelocity(0.05), tracers = (HeatFlux(x -> 40.0 + 10.0 * cos(2π * x / Lx)), Flux(0.0))),
    top    = Wall(tangential = Stress(x -> 2e-4 * sin(2π * x / Lx)), tracers = (Flux(1e-3), Flux(-1e-6))))

# discrete anelastic divergence ∇·(ρ_R u) at interior cells, and the size of its terms
function an_divergence(model, u, w)
    (; domain, dx, dz) = model
    (; ρ_R) = model.reference
    ic, jc = interior_xz(domain)
    div = [ρ_R[j] * (u[east(domain, i), j] - u[i, j]) / dx + (ℑ(ρ_R, (2,), (Face(),), j+1) * w[i, j+1] - ℑ(ρ_R, (2,), (Face(),), j) * w[i, j]) / dz for i in ic, j in jc]
    scale = max(maximum(abs, ρ_R[j] * u[i, j] / dx for i in ic, j in jc), maximum(abs, ℑ(ρ_R, (2,), (Face(),), j) * w[i, j] / dz for i in ic, j in jc), floatmin())
    return div, scale
end

# Tests

@testset "AN2D" begin
    model = an_model(mgr = VectorizedCPU())         # closed walls, periodic in x; serial, for the allocation test
    active = an_state(model)
    box = an_model(boundary = (Bounded, Bounded), walls = forced_box)    # threaded: memory safety also catches races

    # energy-conserving configuration: every flux through the walls vanishes, so mass, tracers and energy are conserved
    @testset "closed walls, periodic in x" begin
        div, scale = an_divergence(model, active.u, active.w)
        @test maximum(abs, div) / scale < 1e-13                              # the test state is divergence-free
        test_memory_safety(an_tend!, model, active)
        test_walls(an_tend!, model, active)
        test_conservation(an_tend!, model, active, (; mass = nothing, consvar = 1, q = 2))   # m = (ρ_R θ, ρ_R q)
        r = test_energy(an_tend!, an_cell_energies, model, active)
        @info "AN2D relative dE/dt residual" r
        ds, _ = evaluate(an_tend!, model, active)
        div, scale = an_divergence(model, ds.u, ds.w)
        @test maximum(abs, div) / scale < 1e-11                              # anelastic constraint ∇·(ρ_R ∂ₜu) = 0
        test_no_allocations(an_tend!, model, active)
    end

    # a short run of the time loop, 3 steps per slice
    @testset "loop" begin
        test_loop(CFAnelastics.loop, CFTimeSchemes.RungeKutta4(model), model, active, (; mass = nothing, q = 2),
                  (; Nslice = 2, slice_size = 0.2, dt_max = 1.0, cfl = 0.9))
    end

    # every kind of wall condition: the wall fluxes and the work of the wall stresses enter the budgets exactly
    @testset "forced walls, closed box, threaded" begin
        s = an_state(box)
        test_memory_safety(an_tend!, box, s)
        test_walls(an_tend!, box, s)
        test_conservation(an_tend!, box, s, (; mass = nothing, consvar = 1, q = 2))
        r = test_energy_budget(an_tend!, an_cell_energies, box, s)           # dE/dt = wall input
        @test abs(r.input) > 1e-8                                            # the forcing is felt
        ds, scratch = evaluate(an_tend!, box, s)
        @test wall_input(box, scratch, 2, 1).heat ≈ 40.0 * Lx rtol = 1e-12     # bottom HeatFlux delivers exactly ∫F dx
        div, scale = an_divergence(box, ds.u, ds.w)
        @test maximum(abs, div) / scale < 1e-11                              # anelastic constraint
    end
end

# the reference pressure is hydrostatic and takes the value `pref` at the height z0 where the geopotential vanishes, inside
# or outside the box; the gauge of the pressure changes the pressure and no tendency
@testset "AN2D reference pressure and gauge" begin
    domain = Box2D(; Mx = 8, Mz = 16, boundary = (Periodic, Bounded))
    dz, Hz, Mz = Lz / 16, halo_size(domain)[2], 16
    for z0 in (0.0, 0.2, Lz, -0.3, Lz + 0.4)
        tank = Tank2D(Lx = Lx, Lz = Lz, g = g, z0 = z0)
        (; ρ_R, p_R) = AnelasticReference(domain, tank, an_density; pref = 1e5)
        j = clamp(floor(Int, (z0 - cellpoint(domain, dz, 2, Hz + 1)) / dz) + Hz + 1, Hz + 1, Hz + Mz - 1)
        @test p_R[j] + (z0 - cellpoint(domain, dz, 2, j)) / dz * (p_R[j+1] - p_R[j]) ≈ 1e5 rtol = 1e-14
        @test all(k -> p_R[k] - p_R[k+1] ≈ g * dz * (ρ_R[k] + ρ_R[k+1]) / 2, 1:Mz+1)       # ∂p_R/∂z = -g ℑz(ρ_R)
    end
    c = 37.5
    for walls in ((; bottom = closed_wall(), top = closed_wall()), forced_box)
        boundary = walls === forced_box ? (Bounded, Bounded) : (Periodic, Bounded)
        model, shifted = an_model(; boundary, walls), an_model(; boundary, walls, gauge_offset = c)
        state = an_state(model)
        ds, scratch = evaluate(an_tend!, model, state)
        ds_c, scratch_c = evaluate(an_tend!, shifted, state)
        for name in (:u, :w, :m)
            @test maximum(abs, getproperty(ds_c, name) .- getproperty(ds, name)) < 1e-12 * maximum(abs, getproperty(ds, name))
        end
        (; domain) = model
        @test interior(domain, scratch_c.φ) .- interior(domain, scratch.φ) ≈ fill(c, size(interior(domain, scratch.φ))) rtol = 1e-12
    end
end

# the model in Float32 (every complex level of the Poisson solver is not 16-byte aligned there): it runs, keeps its float
# type, satisfies the anelastic constraint and conserves energy to Float32 round-off, and allocates nothing
@testset "AN2D in Float32" begin
    domain, space32 = Box2D(; Mx = 32, Mz = 16, boundary = (Periodic, Bounded)), Tank2D(Lx = Float32(Lx), Lz = Float32(Lz), g = Float32(g))
    wall = Wall(; tangential = FreeSlip(), tracers = (Flux(0f0), Flux(0f0)))
    model = AN2D((; mgr = VectorizedCPU(), domain, space = space32, fluid = NonlinearBinaryFluid(:potential_temperature, map(Float32, fluid_params)),
        spatial_scheme = Centred2(), vortex_scheme = VortexEnergyCons(), buoyancy_scheme = ThermodynamicForcing(),
        thermoflux_scheme = ThermoFluxConsistent(k_T = 1f-4, k_q = 1f-5), viscosity_scheme = ViscosityScheme(dyn_visc = 1f-3),
        anelastic_reference = AnelasticReference(domain, space32, z -> Float32(an_density(z)); pref = 1f5),
        boundary_conditions = BoxBoundaryConditions(domain; bottom = wall, top = wall), gauge_offset = 0f0))
    tend32!(dstate, scratch, model, state) = tendencies!(dstate, scratch, model, state, 0f0)
    state = an_state(model)
    ds, scratch = evaluate(tend32!, model, state)
    @test all(a -> eltype(a) == Float32, (values(state)..., values(ds)..., scratch.φ)) && all_finite(ds)
    div, scale = an_divergence(model, ds.u, ds.w)
    @test maximum(abs, div) / scale < 1e-4
    test_energy(tend32!, an_cell_energies, model, state; tol = 1e-5)
    test_no_allocations(tend32!, model, state)
end

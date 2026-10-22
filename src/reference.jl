# Anelastic reference state

"""
    reference = AnelasticReference(domain, space, ρ; pref)

Reference state of the anelastic model on the grid of `domain` spanning `space`, along the full z axis (halo
included): the density `ρ_R` and the hydrostatic pressure `p_R` at cell centers, with ∂p_R/∂z = -g ℑz(ρ_R) discretely
(the density of the w-face between two cells is their average, as in the model). The pressure is anchored at the height
`space.z0` where the geopotential vanishes: `pref` [Pa] is the value at z0 of the piecewise-linear interpolant of the
cell-center pressures (extrapolated linearly from the two nearest interior cells if z0 is outside the box). The density
`ρ` [kg m⁻³] is a number, a function ρ(z) of the height z above the bottom, or a vector of its values in the Mz cells,
from bottom to top (extended by its end values into the halo). Face values and the specific volume are computed where
needed: ρ_R on w-faces is the average of the two cells.
"""
struct AnelasticReference{F, Vec}
    "pressure at the height z0 of the geopotential origin [Pa]"
    pref::F
    "reference density at the cell centers of the full z axis, halo included [kg m⁻³]"
    ρ_R::Vec
    "hydrostatic reference pressure at the cell centers of the full z axis, halo included [Pa]"
    p_R::Vec
end
function AnelasticReference(domain, space, ρ; pref)
    (; Lz, g, z0) = space
    Mz, Hz = dims(domain)[end], halo_size(domain)[end]
    dz = Lz / Mz
    check_length(ρ, Mz)
    ρ_R = [cell_value(ρ, j - Hz, cellpoint(domain, dz, 2, j)) for j in 1:Mz+2Hz]
    p_R = anchored_pressure(hydrostatic_pressure(ρ_R, dz, g), domain, dz, z0, pref)
    return AnelasticReference(pref, ρ_R, p_R)
end

"""
    value = cell_value(f, k, z)

Value of a profile in cell k (k = 1, …, Mz in the interior) at height z: a number, a function of z, or a vector of
the Mz cell values, extended by its end values.
"""
cell_value(f::Number, k, z) = f
cell_value(f, k, z) = f(z)
cell_value(f::AbstractVector, k, z) = f[clamp(k, 1, length(f))]

"""
    check_length(ρ, Mz)

Throw an error if the density profile `ρ` is a vector without one value per interior cell.
"""
check_length(ρ, Mz) = nothing
check_length(ρ::AbstractVector, Mz) = length(ρ) == Mz || error("AnelasticReference: $(length(ρ)) density values for $Mz cells.")

"""
    p = hydrostatic_pressure(ρ, dz, g)

Hydrostatic pressure at the cell centers up to an additive constant (zero in the first cell): p_j - p_{j-1} =
-g Δz (ρ_{j-1} + ρ_j)/2, i.e. ∂p/∂z = -g ℑz(ρ) between centers.
"""
function hydrostatic_pressure(ρ, dz, g)
    p = zeros(eltype(ρ), length(ρ))
    for j in 2:length(ρ)
        p[j] = p[j-1] - g * dz * (ρ[j-1] + ρ[j]) / 2
    end
    return p
end

"""
    p = anchored_pressure(p, domain, dz, z0, pref)

`p` plus the constant such that its piecewise-linear interpolant between cell centers has the value `pref` at the height
z0. The segment is that between two interior cells containing z0, or the nearest one (linear extrapolation) outside the
box.
"""
function anchored_pressure(p, domain, dz, z0, pref)
    Mz, Hz = dims(domain)[end], halo_size(domain)[end]
    z_first = cellpoint(domain, dz, 2, Hz + 1)
    j = clamp(floor(Int, (z0 - z_first) / dz) + Hz + 1, Hz + 1, Hz + Mz - 1)     # lower cell of the segment
    p_z0 = p[j] + (z0 - cellpoint(domain, dz, 2, j)) / dz * (p[j+1] - p[j])
    return p .+ (pref - p_z0)
end

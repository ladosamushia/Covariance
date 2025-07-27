include("distance.jl")
include("bin_index.jl")

"""
    neighbor_counts!(
        count_histogram, X, Y, Z, L, rmax, Nr
    ) -> Nothing

Computes pairwise neighbor counts between particles in a 3D periodic box.
For each pair `(i, j)` with separation less than `rmax`, the function contributes to:
- A **count histogram** indexed by the radial bin `r_bin` for galaxy `i`
- A **Legendre-weighted histogram** (P₂) indexed by `Nr + r_bin` for galaxy `j`

This version is **serial** (not multithreaded).

# Arguments
- `count_histogram::Matrix{Float64}`:
    Preallocated array of size `(Ngal, 2 * Nr)`, where:
    - Columns `1:Nr` hold the pair counts.
    - Columns `Nr+1:2Nr` hold the Legendre-weighted (P₂) sums.
- `X, Y, Z::Vector{Float64}`:
    Coordinates of the particles; all vectors must be of length `Ngal`.
- `L::Float64`: Size of the cubic periodic box (applies in all three directions).
- `rmax::Float64`: Maximum radial distance for a pair to be counted.
- `Nr::Int`: Number of radial bins. Bins are uniform in linear distance.

# Binning Scheme
- Radial bin width is `Δr = rmax / Nr`
- Radial bin edges are `[(i * Δr)^2 for i in 0:Nr]` in squared distance.
- A pair is assigned to bin `r_bin` satisfying:
  r_bin = floor(Int, sqrt(r2) / Δr) + 1
"""
function neighbor_counts!(
    count_histogram::AbstractMatrix{Float64},
    X::AbstractVector{Float64},
    Y::AbstractVector{Float64},
    Z::AbstractVector{Float64},
    L::Float64,
    rmax::Float64,
    Nr::Int
)
    Ngal = length(X)
    Δr = rmax / Nr
    r2max = rmax * rmax

    @inbounds for i in 1:Ngal
        xi, yi, zi = X[i], Y[i], Z[i]
        for j in i+1:Ngal
            r2, μ2 = distance_squared(xi, yi, zi, X[j], Y[j], Z[j], L)

            if r2 < r2max
                r_bin = bin_index_squared(r2, Δr)
                count_histogram[i, r_bin] += 1
                count_histogram[j, Nr + r_bin] += 5.0 * (3.0 * μ2 - 1) / 2.0
            end
        end
    end

    return nothing
end
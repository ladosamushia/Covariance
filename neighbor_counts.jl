include("distance.jl")
include("bin_index.jl")

"""
    neighbor_counts!(count_histogram, X, Y, Z, L, dmax, Δr, Δμ)

Computes the pairwise neighbor counts between particles in 3D periodic space.

Each unique pair `(i, j)` is examined, and if the distance between particles `i` and `j`
is less than `dmax`, it contributes to a 2D histogram indexed by radial and angular bins.

The histogram counts are stored symmetrically:
- Both `i` and `j` increment the same bin.

# Arguments
- `count_histogram::Matrix{Int}`:
    A 2D array of size `(Ngal, Nr × Nμ)` where each row corresponds to a galaxy and each column to a (r, μ) bin.
- `X, Y, Z::Vector{Float64}`:
    Coordinates of the galaxies (same length `Ngal`).
- `L::Float64`: Box size (assumed cubic periodic).
- `dmax::Float64`: Maximum distance threshold (radial cutoff).
- `Δr::Float64`: Width of radial bins.
- `Δμ::Float64`: Width of μ² (cosine squared angle) bins.

# Notes
- Assumes bin edges for r and μ² are multiples of `Δr` and `Δμ`, respectively.
- Binning uses `bin_index_squared(d2, Δ)` assuming uniform bin spacing in squared units.
"""
function neighbor_counts!(
    count_histogram::AbstractMatrix{Int64},
    X::AbstractVector{Float64},
    Y::AbstractVector{Float64},
    Z::AbstractVector{Float64},
    L::Float64,
    dmax::Float64,
    Δr::Float64,
    Δμ::Float64
)
    Ngal = length(X)
    Nμ = Int(floor(1.0 / Δμ))  # μ² ∈ [0, 1)

    d2max = dmax^2

    @inbounds for i in 1:Ngal
        xi, yi, zi = X[i], Y[i], Z[i]

        for j in i+1:Ngal
            d2, μ2 = distance_squared(xi, yi, zi, X[j], Y[j], Z[j], L)

            if d2 < d2max
                r_bin = bin_index_squared(d2, Δr)  # no offset
                μ_bin = bin_index_squared(μ2, Δμ)

                rμ_bin = (r_bin - 1) * Nμ + μ_bin
                count_histogram[i, rμ_bin] += 1
                count_histogram[j, rμ_bin] += 1
            end
        end
    end

    return nothing
end

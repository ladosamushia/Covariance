using Base.Threads: @threads, threadid, nthreads

include("distance.jl")
include("bin_index.jl")

"""
    neighbor_counts!(
        count_histogram, X, Y, Z, L, dmax, Δr, Δμ
    ) -> Nothing

Computes the pairwise neighbor counts between particles in 3D periodic space
and fills a histogram indexed by radial and angular bins.

Each unique pair `(i, j)` is examined. If the distance between particles `i` and `j`
is less than `dmax`, the pair contributes to a 2D histogram indexed by radial bin `r_bin`
and cosine-squared angle bin `μ_bin`. The histogram bin index is flattened as:
`(r_bin - 1) * Nμ + μ_bin`.

The histogram counts are stored symmetrically:
- Both `i` and `j` increment the same bin.

This version is multithreaded using `Threads.@threads`, with thread-local histograms
to ensure race-free parallel updates. Results are summed into the final histogram at the end.

# Arguments
- `count_histogram::Matrix{Int}`:
    A 2D array of size `(Ngal, Nr × Nμ)` where each row corresponds to a galaxy
    and each column to a `(r, μ)` bin. The output is written in-place.
- `X, Y, Z::Vector{Float64}`:
    Coordinates of the galaxies (same length `Ngal`).
- `L::Float64`: Box size (assumed cubic and periodic).
- `dmax::Float64`: Maximum distance threshold (radial cutoff).
- `Δr::Float64`: Width of radial bins (in linear `r`).
- `Δμ::Float64`: Width of μ² (cosine squared angle) bins.

# Notes
- Assumes bin edges for `r²` and `μ²` are `[(iΔ)^2 for i=0:Nr]` and `[iΔμ for i=0:Nμ]`, respectively.
- Binning uses `bin_index_squared(value, Δ)` which computes:
  `i = floor(sqrt(value)/Δ) + 1`
- μ² is used instead of μ to avoid computing square roots.
- The function assumes `count_histogram` is preallocated and correctly sized.
"""
function neighbor_counts!(
    count_histogram::AbstractMatrix{Int},
    X::AbstractVector{Float64},
    Y::AbstractVector{Float64},
    Z::AbstractVector{Float64},
    L::Float64,
    dmax::Float64,
    Δr::Float64,
    Δμ::Float64
)
    Ngal = length(X)
    @assert size(count_histogram, 1) == Ngal

    Nμ = Int(floor(1.0 / Δμ))                # Number of μ² bins
    Nbins = size(count_histogram, 2)         # Total number of (r, μ) bins
    d2max = dmax * dmax                      # Precompute squared distance cutoff

    nthread = nthreads()

    # Thread-local histograms: one per thread to avoid race conditions
    local_histograms = [zeros(Int, Ngal, Nbins) for _ in 1:nthread]

    @threads for i in 1:Ngal
        tid = threadid()
        xi, yi, zi = X[i], Y[i], Z[i]
        H = local_histograms[tid]

        @inbounds for j in i+1:Ngal
            d2, μ2 = distance_squared(xi, yi, zi, X[j], Y[j], Z[j], L)

            if d2 < d2max
                r_bin = bin_index_squared(d2, Δr)
                μ_bin = bin_index_squared(μ2, Δμ)
                rμ_bin = (r_bin - 1) * Nμ + μ_bin

                H[i, rμ_bin] += 1
                H[j, rμ_bin] += 1
            end
        end
    end

    # Merge all local histograms into the final one
    @inbounds for tid in 1:nthread
        count_histogram .+= local_histograms[tid]
    end

    return nothing
end

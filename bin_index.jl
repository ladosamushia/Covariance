"""
    bin_index_squared(value, a, Δ) -> Int

Returns the bin index `i` such that the value falls into the bin
`[(a + (i - 1)Δ)^2, (a + iΔ)^2)`. Assumes value is in range.

# Arguments
- `value` : The value to bin (assumed ≥ 0 and in range)
- `a`     : Starting point of the linear grid before squaring (e.g. 0.0)
- `Δ`     : Step size of the linear grid (before squaring)

# Returns
- `i` : index such that `bin_edges[i] ≤ value < bin_edges[i+1]`
"""
@inline function bin_index_squared(value::Float64, a::Float64, Δ::Float64)
    s = sqrt(value)
    return Int(floor((s - a)/Δ)) + 1
end

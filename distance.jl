"""
    distance_squared(x1, y1, z1, x2, y2, z2, L) -> Float64

Computes the squared Euclidean distance between two points in 3D space
under periodic boundary conditions for a box of size L.

# Arguments
- `x1`, `y1`, `z1` : Coordinates of the first point (Float64)
- `x2`, `y2`, `z2` : Coordinates of the second point (Float64)
- `L`              : Box size (Float64), assumed to be the same in all dimensions

# Returns
- The **squared wrapped distance** between the two points as a `Float64`

# Notes
- This function avoids computing `sqrt`, making it suitable for performance-critical
  code where you only need relative distances or comparisons.
- Assumes a cubic periodic box with side length `L`.

# Example
```julia
julia> distance_squared(1.0, 1.0, 1.0, 9.5, 9.5, 9.5, 10.0)
6.75
"""
@inline function distance_squared(x1, y1, z1, x2, y2, z2, L)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    # Apply minimum image convention
    dx = dx >  L/2 ? dx - L : dx < -L/2 ? dx + L : dx
    dy = dy >  L/2 ? dy - L : dy < -L/2 ? dy + L : dy
    dz = dz >  L/2 ? dz - L : dz < -L/2 ? dz + L : dz

    return dx*dx + dy*dy + dz*dz
end
"""
    distance_squared(x1, y1, z1, x2, y2, z2, L) -> (Float64, Float64)

Computes the **squared Euclidean distance** between two 3D points under periodic
boundary conditions, along with the **cosine squared of the angle** with respect
to the z-axis.

# Arguments
- `x1`, `y1`, `z1` : Coordinates of the first point (`Float64`)
- `x2`, `y2`, `z2` : Coordinates of the second point (`Float64`)
- `L`              : Side length of the periodic cubic box (`Float64`)

# Returns
- `d²` : The **squared wrapped distance** between the two points
- `cos²θ` : The squared cosine of the angle between the vector and the z-axis (`dz² / d²`)

# Notes
- The minimum image convention is applied to compute the shortest displacement under periodic boundaries.
- This function avoids computing `sqrt`, making it suitable for performance-critical code.
- If `distance_squared == 0.0`, the return value of `cos²θ` will be `NaN`.

# Example
```julia
julia> d2, cos2θ = distance_squared(1.0, 1.0, 1.0, 9.5, 9.5, 9.5, 10.0)
(6.75, 0.4444)
"""
@inline function distance_squared(x1, y1, z1, x2, y2, z2, L)
    dx = x1 - x2
    dy = y1 - y2
    dz = z1 - z2

    # Apply minimum image convention
    dx = dx > L/2 ? dx - L : dx < -L/2 ? dx + L : dx
    dy = dy > L/2 ? dy - L : dy < -L/2 ? dy + L : dy
    dz = dz > L/2 ? dz - L : dz < -L/2 ? dz + L : dz

    @fastmath distance_squared = dx*dx + dy*dy + dz*dz
    angle_cosine_squared = dz*dz/distance_squared

    return distance_squared, angle_cosine_squared

end
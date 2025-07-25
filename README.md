# Covariance

# Data Structures

X, Y, Z: coordinates of galaxies in a periodic box. Vectors of length Ngal.

count_histogram: number of neighbors each galaxy has in r-μ space. 2D vector of size (Ngal, NrxNμ).

# Functions

distance_squared(x1, y1, z1, x2, y2, z2, L): compute distance assuming wrapping with a box size of L

bin_index_squared(value, Δ): bin edges are equally spaced from zero in steps of Δ. Find to which bin the squared value belongs.

neighbor_counts!(count_histogram, X, Y, Z, L, dmax, Δr, Δμ): Compute number of neighbours in r-μ space that every galaxy has individually.
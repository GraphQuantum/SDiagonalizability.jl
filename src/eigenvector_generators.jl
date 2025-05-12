# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    _potential_kernel_eigvecs(n)

Return a lazy generator of potential kernel `{-1,0,1}`-eigenvectors of a Laplacian.

Each vector is normalized so that its first nonzero entry is `1`, enforcing pairwise linear
independence between all generated vectors.

# Arguments
- `n::Int`: the order of the Laplacian matrix of some undirected graph for which to find
    potential kernel eigenvectors.

# Returns
- `eigvec_generator::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over
    all `{-1,0,1}`-vectors in `ℝⁿ`, unique up to span.

# Examples
Generate all potential kernel eigenvectors for an order `2` Laplacian matrix:
```jldoctest
julia> collect(Vector{Int}, SDiagonalizability._potential_kernel_eigvecs(2))
4-element Vector{Vector{Int64}}:
 [1, -1]
 [1, 0]
 [1, 1]
 [0, 1]
```
"""
function _potential_kernel_eigvecs(n::Int)
    # Cache to avoid redundant recomputations of the `leading` vector
    leading_cache = Dict{Int,Vector{Int}}()

    return (
        vcat(
            get!(leading_cache, k) do
                leading = Vector{Int}(undef, k)
                leading[1:(k - 1)] .= 0
                leading[k] = 1 # Normalize the leading entry to 1
                leading
            end,
            body,
        ) # Append the permuted entries to the leading [0, ..., 0, 1] vector
        for k in 1:n # Iterate over possible indices of first nonzero entry
        # Iterate over permutations taken from the `r`th Cartesian power of {-1,0,1}
        for body in multiset_permutations(
            let r = n - k
                entries = Vector{Int}(undef, 3r)
                entries[1:r] .= -1
                entries[(r + 1):2r] .= 0
                entries[(2r + 1):3r] .= 1
                entries
            end,
            n - k,
        )
    )
end

"""
    _potential_nonkernel_eigvecs(n)

Return a lazy generator of potential non-kernel `{-1,0,1}`-eigenvectors of a Laplacian.

Each vector is normalized so that its first nonzero entry is `1`, enforcing pairwise linear
independence between all generated vectors. Since all Laplacian matrices have pairwise
orthogonal eigenspaces and the all-ones vector is always in the kernel, every non-kernel
`{-1,0,1}`-eigenvector must also have an equal number of `-1`'s and `1`'s.

# Arguments
- `n::Int`: the order of the Laplacian matrix of some undirected graph for which to find
    potential non-kernel eigenvectors.

# Returns
- `eigvec_generator::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over
    all `{-1,0,1}`-vectors in `ℝⁿ` orthogonal to the all-ones kernel vector, unique up to
    span.

# Examples
Generate all potential non-kernel eigenvectors of an order `4` Laplacian matrix:
```jldoctest
julia> collect(Vector{Int}, SDiagonalizability._potential_nonkernel_eigvecs(4))
9-element Vector{Vector{Int64}}:
 [1, -1, 0, 0]
 [1, 0, -1, 0]
 [1, 0, 0, -1]
 [1, -1, -1, 1]
 [1, -1, 1, -1]
 [1, 1, -1, -1]
 [0, 1, -1, 0]
 [0, 1, 0, -1]
 [0, 0, 1, -1]
```
"""
function _potential_nonkernel_eigvecs(n::Int)
    # Caches to avoid redundant recomputations of the `leading` and `entries` vectors
    leading_cache = Dict{Int,Vector{Int}}()
    entries_cache = Dict{Int,Vector{Int}}()
    entries_key = 0

    return (
        vcat(
            get!(leading_cache, k) do
                leading = Vector{Int}(undef, k)
                leading[1:(k - 1)] .= 0
                leading[k] = 1 # Normalize the leading entry to 1
                leading
            end,
            body,
        ) # Append the permuted entries to the leading [0, ..., 0, 1] vector
        for k in 1:n # Iterate over indices of the first nonzero entry
        # Iterate over combinations of the remaining `r` entries, varying the number of -1's
        for j in 1:ceil(Int, (r = n - k) / 2)
        # Iterate over permutations of the `entries` combination with `j` -1's
        for body in multiset_permutations(
            get!(entries_cache, entries_key += 1) do
                entries = Vector{Int}(undef, r)
                entries[1:j] .= -1
                entries[(j + 1):(r - j + 1)] .= 0
                entries[(r - j + 2):r] .= 1
                entries
            end,
            r,
        )
    )
end

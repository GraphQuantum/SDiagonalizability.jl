# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    pot_kernel_s_eigvecs(n, S) -> Iterators.Flatten{<:Base.Generator}

[TODO: Write here]

# Arguments
[TODO: Write here]

# Returns
[TODO: Write here]

# Throws
- `DomainError`: if `n` is negative.

# Examples
[TODO: Write here]
"""
function pot_kernel_s_eigvecs(n::Integer, S::Tuple{Vararg{Integer}})
    if n < 0
        throw(DomainError(n, "Laplacian order must be non-negative"))
    end

    S = _sort_tuple(S)

    if S == (-1, 0, 1)
        gen = _pot_kernel_01neg_eigvecs(n)
    elseif S == (-1, 1)
        gen = _pot_kernel_1neg_eigvecs(n)
    else
        throw(ArgumentError("Unsupported entry set S: $S"))
    end

    return gen
end

"""
    pot_nonkernel_s_eigvecs(n, S) -> Iterators.Flatten{<:Base.Generator}

[TODO: Write here]

# Arguments
[TODO: Write here]

# Returns
[TODO: Write here]

# Throws
- `DomainError`: if `n` is negative.

# Examples
[TODO: Write here]
"""
function pot_nonkernel_s_eigvecs(n::Integer, S::Tuple{Vararg{Integer}})
    if n < 0
        throw(DomainError(n, "Laplacian order must be non-negative"))
    end

    S = _sort_tuple(S)

    if S == (-1, 0, 1)
        gen = _pot_nonkernel_01neg_eigvecs(n)
    elseif S == (-1, 1)
        gen = _pot_nonkernel_1neg_eigvecs(n)
    else
        throw(ArgumentError("Unsupported entry set S: $S"))
    end

    return gen
end

"""
    _pot_kernel_01neg_eigvecs(n) -> Iterators.Flatten{<:Base.Generator}

Lazily compute all potential kernel ``\\{-1, 0 ,1\\}``-eigenvectors of an ``n×n`` Laplacian.

Each vector is normalized so that its first nonzero entry is ``1``, enforcing pairwise
linear independence between all generated vectors.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential kernel eigenvectors.

# Returns
- `eigvec_generator::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over
    all ``\\{-1, 0, 1\\}``-vectors in ``ℝⁿ``, unique up to span.

# Examples
Generate all potential kernel eigenvectors for an order ``3`` Laplacian matrix:
```jldoctest
julia> hcat(SDiagonalizability._pot_kernel_01neg_eigvecs(3)...)
3×13 Matrix{Int64}:
  1   1   1   1  1  1   1  1  1   0  0  0  0
 -1  -1  -1   0  0  0   1  1  1   1  1  1  0
 -1   0   1  -1  0  1  -1  0  1  -1  0  1  1
```

# Notes
The number of potential kernel eigenvectors (unique up to span) for an order ``n`` Laplacian
matrix is given by ``(3ⁿ - 1) / 2``. See also the relevant OEIS sequence [Slo25].

Regrettably, the implementation here is rather clunky and unidiomatic, but it is worth
noting that eigenvector generation is one of two major bottlenecks in the overall
*S*-bandwidth minimization algorithm. Given how much potential there is for optimization in
this piece of code, we thus prioritize performance over readability in this particular case,
making every effort to include inline comments wherever clarification may be needed.

# References

- [Slo25](@cite): N. J. Sloane, *a(n) = (3^n - 1)/2*. Entry A003462 (2025). Accessed:
    2025-05-22. https://oeis.org/A003462.
"""
function _pot_kernel_01neg_eigvecs(n::Integer)
    # Cache to avoid redundant recomputations of the `leading` vector
    leading_cache = Dict{Int,Vector{Int}}()

    return (
        vcat(
            get!(leading_cache, k) do
                leading = Vector{Int}(undef, k)
                leading[1:(k - 1)] .= 0
                leading[k] = 1 # Normalize the leading entry to 1
                return leading
            end,
            body,
        ) # Append the permuted entries to the leading [0, ..., 0, 1] vector
        for k in 1:n # Iterate over possible indices of first nonzero entry
        # Iterate over permutations taken from the `r`-th Cartesian power of {-1, 0, 1}
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

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_kernel_01neg_eigvecs(0))) = Vector{Int}

"""
    _pot_kernel_1neg_eigvecs(n) -> Iterators.Flatten{<:Base.Generator}

[TODO: Write here. Also, comment inline]
"""
function _pot_kernel_1neg_eigvecs(n::Integer)
    if n == 0
        gen = (Int[] for _ in 1:0)
    else
        entries = Vector{Int}(undef, 2(n - 1))
        entries[1:(n - 1)] .= -1
        entries[n:(2(n - 1))] .= 1
        gen = (vcat(1, body) for body in multiset_permutations(entries, n - 1))
    end

    return gen
end

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_kernel_1neg_eigvecs(0))) = Vector{Int}

"""
    _pot_nonkernel_01neg_eigvecs(n) -> Iterators.Flatten{<:Base.Generator}

Lazily compute all potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors of an ``n×n`` Laplacian.

Each vector is normalized so that its first nonzero entry is ``1``, enforcing pairwise
linear independence between all generated vectors. Since all Laplacian matrices have
pairwise orthogonal eigenspaces and the all-ones vector is always in the kernel, every
non-kernel ``\\{-1, 0, 1\\}``-eigenvector must also have an equal number of ``-1``'s and
``1``'s.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential non-kernel eigenvectors.

# Returns
- `eigvec_generator::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over
    all ``\\{-1, 0, 1\\}``-vectors in ``ℝⁿ`` orthogonal to the all-ones kernel vector,
    unique up to span.

# Examples
Generate all potential non-kernel eigenvectors of an order ``4`` Laplacian matrix:
```jldoctest
julia> hcat(SDiagonalizability._pot_nonkernel_01neg_eigvecs(4)...)
4×9 Matrix{Int64}:
  1   1   1   1   1   1   0   0   0
 -1   0   0  -1  -1   1   1   1   0
  0  -1   0  -1   1  -1  -1   0   1
  0   0  -1   1  -1  -1   0  -1  -1
```

# Notes
The number of potential non-kernel eigenvectors (unique up to span) for an order ``n``
Laplacian matrix is, by non-trivial combinatorial arguments, equal to the number of humps in
all Motzkin paths of length ``n``. See also the relevant OEIS sequence [Deu25].

Regrettably, the implementation here is rather clunky and unidiomatic, but it is worth
noting that eigenvector generation is one of two major bottlenecks in the overall
*S*-bandwidth minimization algorithm. Given how much potential there is for optimization in
this piece of code, we thus prioritize performance over readability in this particular case,
making every effort to include inline comments wherever clarification may be needed.

# References

- [Deu25](@cite): E. Deutsch. *Number of humps in all Motzkin paths of length n*. Entry
    A097861 (2025). Accessed: 2025-05-22. https://oeis.org/A097861.
"""
function _pot_nonkernel_01neg_eigvecs(n::Integer)
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
                return leading
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
                return entries
            end,
            r,
        )
    )
end

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_nonkernel_01neg_eigvecs(0))) = Vector{Int}

"""
    _pot_nonkernel_1neg_eigvecs(n)

[TODO: Write here. Also, comment inline]
"""
function _pot_nonkernel_1neg_eigvecs(n::Integer)
    if n == 0 || isodd(n)
        gen = (Int[] for _ in 1:0)
    else
        entries = Vector{Int}(undef, n - 1)
        mid = div(n, 2)
        entries[1:mid] .= -1
        entries[(mid + 1):(n - 1)] .= 1
        gen = (vcat(1, body) for body in multiset_permutations(entries, n - 1))
    end

    return gen
end

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_nonkernel_1neg_eigvecs(0))) = Vector{Int}

# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    pot_kernel_s_eigvecs(n, S) -> Union{Base.Generator,Base.Iterators.Flatten{<:Base.Generator}}

Lazily compute all potential kernel `S`-eigenvectors of an `n×n` Laplacian.

The only accepted values of `S` are `(-1, 0, 1)`, `(-1, 1)`, and permutations thereof (e.g.,
`(1, -1, 0)`, which will be sorted internally to `(-1, 0, 1)` anyway), in line with the sets
studied by [JP25].

Each vector is normalized so that its first nonzero entry is ``1``, enforcing pairwise
linear independence between all generated vectors.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential kernel `S`-eigenvectors.
- `S::Tuple{Vararg{Integer}}`: the set of possible entries of the potential eigenvectors,
    provided in the form of a tuple of unique `Integer`'s. Must be either `(-1, 0, 1)`,
    `(-1, 1)`, or a permutation thereof.

# Returns
- `gen::Union{Base.Generator,Base.Iterators.Flatten{<:Base.Generator}}`: a lazily evaluated
    iterator over all `S`-vectors in ``ℝⁿ``, unique up to span. Eltype is `Vector{Int}`.

# Throws
- `DomainError`: if `n` is negative.
- `ArgumentError`: if `S` is not `(-1, 0, 1)`, `(-1, 1)`, or a permutation thereof.

# Examples
Generate all potential kernel ``\\{-1, 0, 1\\}``-eigenvectors of a ``4×4`` Laplacian matrix:
```jldoctest
julia> stack(SDiagonalizability.pot_kernel_s_eigvecs(4, (-1, 0, 1)))
4×40 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1  …   0   0  0  0   0  0  0   0  0  0  0
 -1  -1  -1  -1  -1  -1  -1  -1  -1      1   1  1  1   1  1  1   0  0  0  0
 -1  -1  -1   0   0   0   1   1   1     -1   0  0  0   1  1  1   1  1  1  0
 -1   0   1  -1   0   1  -1   0   1      1  -1  0  1  -1  0  1  -1  0  1  1
```

Generate all potential kernel ``\\{-1, 1\\}``-eigenvectors of a ``5×5`` Laplacian matrix,
with the `S` set provided in a different order:
```jldoctest
julia> stack(SDiagonalizability.pot_kernel_s_eigvecs(5, (1, -1)))
5×16 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1  1
 -1  -1  -1  -1  -1  -1  -1  -1   1   1   1   1   1   1   1  1
 -1  -1  -1  -1   1   1   1   1  -1  -1  -1  -1   1   1   1  1
 -1  -1   1   1  -1  -1   1   1  -1  -1   1   1  -1  -1   1  1
 -1   1  -1   1  -1   1  -1   1  -1   1  -1   1  -1   1  -1  1
```

Only the first two entry sets are supported, so the following command throws an error:
```jldoctest
julia> stack(SDiagonalizability.pot_kernel_s_eigvecs(6, (1, 2)))
ERROR: ArgumentError: Unsupported entry set S: (1, 2)
[...]
```

# References

- [JP25](@cite): N. Johnston and S. Plosker. *Laplacian {−1,0,1}- and {−1,1}-diagonalizable
    graphs*. Linear Algebra and its Applications **704**, 309–39 (2025).
    https://doi.org/10.1016/j.laa.2024.10.016.
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
    pot_nonkernel_s_eigvecs(n, S) -> Union{Base.Generator,Base.Iterators.Flatten{<:Base.Generator}}

Lazily compute all potential non-kernel `S`-eigenvectors of an `n×n` Laplacian.

The only accepted values of `S` are `(-1, 0, 1)`, `(-1, 1)`, and permutations thereof (e.g.,
`(1, -1, 0)`, which will be sorted internally to `(-1, 0, 1)` anyway), in line with the sets
studied by [JP25].

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential non-kernel `S`-eigenvectors.
- `S::Tuple{Vararg{Integer}}`: the set of possible entries of the potential eigenvectors,
    provided in the form of a tuple of unique `Integer`'s. Must be either `(-1, 0, 1)`,
    `(-1, 1)`, or a permutation thereof.

# Returns
- `gen::Union{Base.Generator,Base.Iterators.Flatten{<:Base.Generator}}`: a lazily evaluated
    iterator over all `S`-vectors in ``ℝⁿ`` orthogonal to the all-ones kernel vector, unique
    up to span. Eltype is `Vector{Int}`.

# Throws
- `DomainError`: if `n` is negative.
- `ArgumentError`: if `S` is not `(-1, 0, 1)`, `(-1, 1)`, or a permutation thereof.

# Examples
Generate all potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors of a ``5×5`` Laplacian
matrix, with the `S` set provided in a different order:
```jldoctest
julia> stack(SDiagonalizability.pot_nonkernel_s_eigvecs(5, (1, -1, 0)))
5×25 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1  …   0   0   0   0   0   0   0   0   0
 -1   0   0   0  -1  -1  -1  -1  -1      1   1   1   1   1   1   0   0   0
  0  -1   0   0  -1  -1   0   0   1     -1   0   0  -1  -1   1   1   1   0
  0   0  -1   0   0   1  -1   1  -1      0  -1   0  -1   1  -1  -1   0   1
  0   0   0  -1   1   0   1  -1   0      0   0  -1   1  -1  -1   0  -1  -1
```

Generate all potential non-kernel ``\\{-1, 1\\}``-eigenvectors of a ``6×6`` Laplacian
matrix:
```jldoctest
julia> stack(SDiagonalizability.pot_nonkernel_s_eigvecs(6, (-1, 1)))
6×10 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1   1
 -1  -1  -1  -1  -1  -1   1   1   1   1
 -1  -1  -1   1   1   1  -1  -1  -1   1
 -1   1   1  -1  -1   1  -1  -1   1  -1
  1  -1   1  -1   1  -1  -1   1  -1  -1
  1   1  -1   1  -1  -1   1  -1  -1  -1
```

Only the first two entry sets are supported, so the following command throws an error:
```jldoctest
julia> stack(SDiagonalizability.pot_nonkernel_s_eigvecs(7, (0, 3)))
ERROR: ArgumentError: Unsupported entry set S: (0, 3)
[...]
```

# References

- [JP25](@cite): N. Johnston and S. Plosker. *Laplacian {−1,0,1}- and {−1,1}-diagonalizable
    graphs*. Linear Algebra and its Applications **704**, 309–39 (2025).
    https://doi.org/10.1016/j.laa.2024.10.016.
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

Lazily compute all potential kernel ``\\{-1, 0 ,1\\}``-eigenvectors of an `n×n` Laplacian.

Each vector is normalized so that its first nonzero entry is ``1``, enforcing pairwise
linear independence between all generated vectors.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential kernel ``\\{-1, 0, 1\\}``-eigenvectors.

# Returns
- `::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over all
    ``\\{-1, 0, 1\\}``-vectors in ``ℝⁿ``, unique up to span. Eltype is `Vector{Int}`.

# Examples
Generate all potential kernel ``\\{-1, 0, 1\\}``-eigenvectors of a ``3×3`` Laplacian matrix:
```jldoctest
julia> stack(SDiagonalizability._pot_kernel_01neg_eigvecs(3))
3×13 Matrix{Int64}:
  1   1   1   1  1  1   1  1  1   0  0  0  0
 -1  -1  -1   0  0  0   1  1  1   1  1  1  0
 -1   0   1  -1  0  1  -1  0  1  -1  0  1  1
```

# Notes
The number of potential kernel ``\\{-1, 0, 1\\}``-eigenvectors (unique up to span) of an
``n×n`` Laplacian matrix is equal to ``(3ⁿ - 1) / 2``. See also the relevant OEIS sequence
[Slo10].

Regrettably, the implementation here is rather clunky and unidiomatic, but it is worth
noting that eigenvector generation is one of two major bottlenecks in the overall
*S*-bandwidth minimization algorithm. Given how much potential there is for optimization in
this piece of code, we thus prioritize performance over readability in this particular case,
making every effort to include inline comments wherever clarification may be needed.

# References

- [Slo10](@cite): N. J. Sloane, *a(n) = (3^n - 1)/2*. Entry A003462 (2010). Accessed:
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
    _pot_kernel_1neg_eigvecs(n) -> Base.Generator

Lazily compute all potential kernel ``\\{-1, 1\\}``-eigenvectors of an `n×n` Laplacian.

Each vector is normalized so that its first entry is ``1``, enforcing pairwise linear
independence between all generated vectors.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential kernel ``\\{-1, 1\\}``-eigenvectors.

# Returns
- `::Base.Generator`: a lazily evaluated iterator over all ``\\{-1, 1\\}``-vectors in
    ``ℝⁿ``, unique up to span. Eltype is `Vector{Int}`.

# Examples
Generate all potential kernel ``\\{-1, 1\\}``-eigenvectors of a ``5×5`` Laplacian matrix:
```jldoctest
julia> stack(SDiagonalizability._pot_kernel_1neg_eigvecs(5))
5×16 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1   1   1   1   1   1   1  1
 -1  -1  -1  -1  -1  -1  -1  -1   1   1   1   1   1   1   1  1
 -1  -1  -1  -1   1   1   1   1  -1  -1  -1  -1   1   1   1  1
 -1  -1   1   1  -1  -1   1   1  -1  -1   1   1  -1  -1   1  1
 -1   1  -1   1  -1   1  -1   1  -1   1  -1   1  -1   1  -1  1
```

# Notes
The number of potential kernel ``\\{-1, 1\\}``-eigenvectors (unique up to span) of an
``n×n`` Laplacian matrix is equal to ``0`` for ``n = 0`` and ``2ⁿ`` for ``n > 0``. See also
the relevant OEIS sequence [Slo14] for the ``n > 0`` case.

# References

- [Slo14](@cite): N. J. Sloane, *Powers of 2: a(n) = 2^n*. Entry A000079 (2014). Accessed:
    2025-08-13. https://oeis.org/A000079.
"""
function _pot_kernel_1neg_eigvecs(n::Integer)
    if n == 0
        entries = Int[]
    else
        entries = Vector{Int}(undef, 2(n - 1))
        entries[1:(n - 1)] .= -1
        entries[n:(2(n - 1))] .= 1
    end

    return (vcat(1, body) for body in multiset_permutations(entries, n - 1))
end

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_kernel_1neg_eigvecs(0))) = Vector{Int}

"""
    _pot_nonkernel_01neg_eigvecs(n) -> Iterators.Flatten{<:Base.Generator}

Lazily compute all potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors of an `n×n` Laplacian.

Each vector is normalized so that its first nonzero entry is ``1``, enforcing pairwise
linear independence between all generated vectors. Since all Laplacian matrices have
pairwise orthogonal eigenspaces and the all-ones vector is always in the kernel, every
non-kernel ``\\{-1, 0, 1\\}``-eigenvector must also have an equal number of ``-1``'s and
``1``'s.

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors.

# Returns
- `::Iterators.Flatten{<:Base.Generator}`: a lazily evaluated iterator over all
    ``\\{-1, 0, 1\\}``-vectors in ``ℝⁿ`` orthogonal to the all-ones kernel vector, unique up
    to span. Eltype is `Vector{Int}`.

# Examples
Generate all potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors of a ``4×4`` Laplacian
matrix:
```jldoctest
julia> stack(SDiagonalizability._pot_nonkernel_01neg_eigvecs(4))
4×9 Matrix{Int64}:
  1   1   1   1   1   1   0   0   0
 -1   0   0  -1  -1   1   1   1   0
  0  -1   0  -1   1  -1  -1   0   1
  0   0  -1   1  -1  -1   0  -1  -1
```

# Notes
The number of potential non-kernel ``\\{-1, 0, 1\\}``-eigenvectors (unique up to span) of an
``n×n`` Laplacian matrix is, by non-trivial combinatorial arguments, equal to the number of
humps in all Motzkin paths of length ``n``. See also the relevant OEIS sequence [Deu21].

Regrettably, the implementation here is rather clunky and unidiomatic, but it is worth
noting that eigenvector generation is one of two major bottlenecks in the overall
*S*-bandwidth minimization algorithm. Given how much potential there is for optimization in
this piece of code, we thus prioritize performance over readability in this particular case,
making every effort to include inline comments wherever clarification may be needed.

# References

- [Deu21](@cite): E. Deutsch. *Number of humps in all Motzkin paths of length n*. Entry
    A097861 (2021). Accessed: 2025-05-22. https://oeis.org/A097861.
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
    _pot_nonkernel_1neg_eigvecs(n) -> Base.Generator

Lazily compute all potential non-kernel ``\\{-1, 1\\}``-eigenvectors of an `n×n` Laplacian.

Each vector is normalized so that its first entry is ``1``, enforcing pairwise linear
independence between all generated vectors. Since all Laplacian matrices have pairwise
orthogonal eigenspaces and the all-ones vector is always in the kernel, every non-kernel
``\\{-1, 1\\}``-eigenvector must half exactly `n / 2` ``-1``'s and `n / 2` ``1``'s. (As a
direct corollary, if `n` is odd, an empty iterator is returned.)

# Arguments
- `n::Integer`: the order of the Laplacian matrix of some undirected graph for which to find
    potential non-kernel ``\\{-1, 1\\}``-eigenvectors.

# Returns
- `::Base.Generator``: a lazily evaluated iterator over all ``\\{-1, 1\\}``-vectors in
    ``ℝⁿ`` orthogonal to the all-ones kernel vector, unique up to span. Eltype is
    `Vector{Int}`.

# Examples
Generate all potential non-kernel ``\\{-1, 1\\}``-eigenvectors of a ``6×6`` Laplacian
matrix:
```jldoctest
julia> stack(SDiagonalizability._pot_nonkernel_1neg_eigvecs(6))
6×10 Matrix{Int64}:
  1   1   1   1   1   1   1   1   1   1
 -1  -1  -1  -1  -1  -1   1   1   1   1
 -1  -1  -1   1   1   1  -1  -1  -1   1
 -1   1   1  -1  -1   1  -1  -1   1  -1
  1  -1   1  -1   1  -1  -1   1  -1  -1
  1   1  -1   1  -1  -1   1  -1  -1  -1
```

# Notes

The number of potential non-kernel ``\\{-1, 1\\}``-eigenvectors (unique up to span) of an
``n×n`` Laplacian matrix is equal to ``n`` for ``n ≤ 1`` and, by non-trivial combinatorial
arguments, the number of Motzkin ``(n + 1)``-paths with exactly one flat step for ``n > 1``.
See also the relevant OEIS sequence [Sut14] for the ``n > 1`` case.

# References

- [Sut14](@cite): A. V. Sutherland. *The number of Motzkin n-paths with exactly one flat
    step*. Entry A138364 (2014). Accessed: 2025-08-13.
"""
function _pot_nonkernel_1neg_eigvecs(n::Integer)
    if n == 0 || isodd(n)
        entries = Int[]
    else
        entries = Vector{Int}(undef, n - 1)
        mid = div(n, 2)
        entries[1:mid] .= -1
        entries[(mid + 1):(n - 1)] .= 1
    end

    return (vcat(1, body) for body in multiset_permutations(entries, n - 1))
end

# Specify the return type of the generator for type inference and stability
Base.eltype(::typeof(_pot_nonkernel_1neg_eigvecs(0))) = Vector{Int}

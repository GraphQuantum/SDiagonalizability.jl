# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    KOrthogonality

An abstract type representing the property of `k`-orthogonality of a collection of vectors.

Recall that an (ordered) collection of vectors ``v₁, v₂, ..., vₙ`` is said to be
*``k``-orthogonal* if we have the inner product ``⟨`vᵢ, vⱼ⟩ = 0`` whenever ``|i - j| ≥ k``
(i.e., if every pair of vectors at least ``k`` indices apart is orthogonal). This is
equivalent to the vectors' Gram matrix having bandwidth at most `k`, where we define the
bandwidth of a matrix ``A`` to be the minimum integer ``k ∈ \\{1, 2, …, n\\}`` such that
``Aᵢⱼ = 0`` whenever ``|i - j| ≥ k`` [JP25; p. 313](@cite). (Note that many texts instead
define matrix bandwidth using zero-based indexing—that is, with the condition
``|i - j| > k`` [Maf14; p. 186](@cite).)

This type is used as a template for concretely defined properties corresponding to specific
values of ``k``. In the context of the overarching *S*-bandwidth algorithm, we perform a
different depth-first search for each family of values of ``k`` on our "tree" of
*S*-eigenvectors to determine whether there exists a ``k``-orthogonal collection of them.

# Interface
Concrete subtypes of `KOrthogonality` **must** implement the following fields:
- `k::Int`: the ``k``-orthogonality parameter. Must be a positive integer.
"""
abstract type KOrthogonality end

"""
    Orthogonality <: KOrthogonality

The property of pairwise orthogonality for a collection of vectors.

Recall that a collection of vectors ``v₁, v₂, ..., vₙ`` is said to be pairwise *orthogonal*
if we have the inner product ``⟨vᵢ, vⱼ⟩ = 0`` whenever ``|i - j| ≥ 1`` (i.e., if every pair
of vectors is orthogonal). This is equivalent to the vectors' Gram matrix being diagonal.

# Fields
- `k::Int`: the ``k``-orthogonality parameter; always necessarily ``1``.

# Supertype Hierarchy
`Orthogonality` <: [`KOrthogonality`](@ref)

# Constructors
- `Orthogonality()`: constructs a new `Orthogonality` object with `k = 1`.
"""
struct Orthogonality <: KOrthogonality
    k::Int

    Orthogonality() = new(1)
end

"""
    QuasiOrthogonality <: KOrthogonality

The property of quasi-orthogonality for a collection of vectors.

Recall that an (ordered) collection of vectors ``v₁, v₂, ..., vₙ`` is said to be
*quasi-orthogonal* if we have the inner product ``⟨vᵢ, vⱼ⟩ = 0`` whenever ``|i - j| ≥ 2``
(i.e., if every pair of vectors at least ``2`` indices apart is orthogonal). This is
equivalent to the vectors' Gram matrix being tridiagonal [JP25; p. 313](@cite).

# Fields
- `k::Int`: the ``k``-orthogonality parameter; always necessarily ``2``.

# Supertype Hierarchy
`QuasiOrthogonality` <: [`KOrthogonality`](@ref)

# Constructors
- `QuasiOrthogonality()`: constructs a new `QuasiOrthogonality` object with `k = 2`.
"""
struct QuasiOrthogonality <: KOrthogonality
    k::Int

    QuasiOrthogonality() = new(2)
end

"""
    WeakOrthogonality <: KOrthogonality

The property of "weak orthogonality" for a collection of vectors.

In particular, "weak orthogonality" is an *ad hoc* term used to refer to the property of
``k``-orthogonality for ``k > 2``. Recall that an (ordered) collection of vectors
``v₁, v₂, ..., vₙ`` is said to be *`k`-orthogonal* if we have the inner product
``⟨vᵢ, vⱼ⟩ = 0`` whenever ``|i - j| ≥ k`` (i.e., if every pair of vectors at least ``k``
indices apart is orthogonal). This is equivalent to the vectors' Gram matrix having
bandwidth at most ``k``.

# Fields
- `k::Int`: the ``k``-orthogonality parameter. Must be greater than ``2``.

# Supertype Hierarchy
`WeakOrthogonality` <: [`KOrthogonality`](@ref)

# Notes
The term "weak orthogonality" is not standard terminology in the literature, but it is used
here to emphasize the weaker nature of this property compared to orthogonality and
quasi-orthogonality. It is an *ad hoc* term coined for this module and is not intended to
be formally introduced in the broader literature.
"""
struct WeakOrthogonality <: KOrthogonality
    k::Int

    function WeakOrthogonality(k::Integer)
        if k <= 2
            throw(DomainError(k, "k-orthogonality parameter must be greater than 2"))
        end

        return new(Int(k))
    end
end

"""
    classify_k_orthogonality(k)

Classifies the `k`-orthogonality property based on the given `k` parameter.

When searching for a ``k``-orthogonal *S*-basis of a given Laplacian eigenspace, the family
of values to which our ``k`` parameter belongs informs our choice of algorithm.

# Arguments
- `k::Integer`: the `k`-orthogonality parameter to classify. Must be a positive integer.

# Returns
- `::KOrthogonality`: An instance of a concrete [`KOrthogonality`](@ref) subtype
    corresponding to `k`.

# Throws
- `DomainError`: if ``k`` is not a positive integer.

# Examples
```jldoctest
julia> SDiagonalizability.classify_k_orthogonality(1)
SDiagonalizability.Orthogonality(1)

julia> SDiagonalizability.classify_k_orthogonality(2)
SDiagonalizability.QuasiOrthogonality(2)

julia> SDiagonalizability.classify_k_orthogonality(3)
SDiagonalizability.WeakOrthogonality(3)

julia> SDiagonalizability.classify_k_orthogonality(13)
SDiagonalizability.WeakOrthogonality(13)
```
"""
function classify_k_orthogonality(k::Integer)
    if k < 1
        throw(DomainError(k, "k-orthogonality parameter must be a positive integer"))
    end

    if k == 1
        prop = Orthogonality()
    elseif k == 2
        prop = QuasiOrthogonality()
    else
        prop = WeakOrthogonality(k)
    end

    return prop
end

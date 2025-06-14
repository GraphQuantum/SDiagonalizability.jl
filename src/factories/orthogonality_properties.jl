# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    abstract type _KOrthogonality

An abstract type representing the property of `k`-orthogonality of a collection of vectors.

Recall that an (ordered) collection of vectors `v₁, v₂, ..., vₙ` is said to be
*`k`-orthogonal* if we have the inner product `⟨vᵢ, vⱼ⟩ = 0` whenever `|i - j| ≥ k` (i.e.,
if every pair of vectors at least `k` indices apart is orthogonal). This is equivalent to
the vectors' Gram matrix having bandwidth at most `k`, where we define the bandwidth of a
matrix `A` to be the minimum integer `k ≥ 1` such that `aᵢⱼ = 0` whenever `|i - j| ≥ k`
[JP25; p. 313](@cite). (Note that many texts instead define matrix bandwidth using
zero-based indexing—that is, with the condition `|i - j| > k` [Maf14; p. 186](@cite).)

This type is used as a template for concretely defined properties corresponding to specific
values of `k`. For each family of values, we apply a separate algorithm to determine whether
there exists a permutation of a given collection of vectors that induces `k`-orthogonality.

# Interface
Concrete subtypes of `_KOrthogonality` **must** implement the following properties:
- `k::Int`: the `k`-orthogonality parameter.

`k` need not necessarily be a field of the subtype (especially in the case of single-value
types like [`_Orthogonality`](@ref) and [`_QuasiOrthogonality`](@ref), in which case making
`k` a field may erroneously imply that it is variable), but it must be accessible via
`Base.getproperty`.

# Subtypes
- [`_Orthogonality`](@ref): `k`-orthogonality for `k = 1`.
- [`_QuasiOrthogonality`](@ref): `k`-orthogonality for `k = 2`.
- [`_WeakOrthogonality`](@ref): `k`-orthogonality for `k > 2`.

# Notes
Given the domain restriction of `k` to the positive integers, combined with the relatively
small nature (`≤ 25`) of bandwidths examined by this package's principal *S*-bandwidth
minimization algorithm, it may seem prudent to ask why we do not type the `k` field as
`UInt8` instead of `Int`. An equally compelling argument can be made to type it as `Integer`
for genericity as well. However, many methods (both within this module and externally) often
expect `Int`-typed parameters, so we settle on `Int` as our field type. However, in any
constructors/factories of `_KOrthogonality` subtypes, we allow the `k`-orthogonality input
parameter to be any generic `Integer`, simply casting it to an `Int` later on.

On that note, the stipulation in the interface contract that `k` be an accessible property
even for semantically unique types like [`_Orthogonality`](@ref) and
[`_QuasiOrthogonality`](@ref) may seem somewhat odd, as it will never accessed in practice
by [`_find_k_orthogonal_basis`](@ref). (Rather, the type itself—which already implicitly
encodes the value of `k`—directly informs the choice of algorithm in said function.)
However, we require it nonetheless for the sake of consistency future extensibility, just in
case the API is later modified to require more generic handling of `_KOrthogonality`
instances whose specific subtypes are not known based on their `k`-orthogonality parameter.

Another design choice worth scrutinizing is the "`_WeakOrthogonality`" name. "Orthogonality"
is a term as old as time, and "quasi-orthogonality" too has risen to prominence in recent
years (e.g., [JP25; p. 313](@cite)), but "weak orthogonality," strictly speaking, is not
standard terminology. We have taken the liberty of coining it here for `k > 2` to emphasize
its weaker nature compared to orthogonality and quasi-orthogonality. That said, we make no
attempts to formally introduce as a new definition in the broader literature beyond its *ad
hoc* use (and but internally, at that) in this module.

Finally, we note that it may be slightly misleading to state that this type is used to
inform our choice of algorithm in determining whether a collection of vectors is
*`k`-orthogonalizable* (i.e., whether there exists a `k`-orthogonal permutation of said
collection). Rather, we use it when determining whether a large, linearly dependent set of
vectors contains an `k`-orthogonal independent spanning subset. This is still related to the
long-standing matrix bandwidth minimization problem, but it differs in that we can take
advantage of dynamic programming techniques when recursively searching for such a subset.
(See the [`_find_k_orthogonal_basis`](@ref) documentation for more details.)
"""
abstract type _KOrthogonality end

"""
    struct _Orthogonality

The property of pairwise orthogonality for a collection of vectors.

Recall that a collection of vectors `v₁, v₂, ..., vₙ` is said to be pairwise *orthogonal* if
we have the inner product `⟨vᵢ, vⱼ⟩ = 0` whenever `|i - j| ≥ 1` (i.e., if every pair of
vectors is orthogonal). This is equivalent to the vectors' Gram matrix being diagonal.

# Properties
- `k::Int`: the `k`-orthogonality parameter. (Necessarily `1`.)

# Supertype Hierarchy
`_Orthogonality` <: [`_KOrthogonality`](@ref) <: Any

# Notes
Since this type is semantically unique inasmuch as it simply encodes the information
`k = 1`, we define the singleton instance [`_ORTHOGONALITY`](@ref) to always be used when
needed in favor of instantiating an identical object every time.
"""
struct _Orthogonality <: _KOrthogonality end

function Base.getproperty(::_Orthogonality, name::Symbol)
    if name == :k
        return 1
    else
        # This should throw an error, unless we later add fields to this type
        return getfield(_Orthogonality, name)
    end
end

"""
    _ORTHOGONALITY::_Orthogonality

A singleton instance of the semantically unique [`_Orthogonality`](@ref) type.

This constant is used to represent the property of pairwise orthogonality for a collection
of vectors and should always be used in favor of instantiating a new object with
`_Orthogonality()` every time.

# Notes
Simply by definition, we are certain that the `k` property of any [`_Orthogonality`](@ref)
instance will always be `1`, so given that the struct is immutable and we do not intend to
add any other distinguishing fields/properties, it is idiomatic to define this sort of
singleton instance for universal use. This design choice also avoids unnecessary allocations
(although any overhead will be negligible in practice).
"""
const _ORTHOGONALITY = _Orthogonality()

"""
    struct _QuasiOrthogonality

The property of quasi-orthogonality for a collection of vectors.

Recall that an (ordered) collection of vectors `v₁, v₂, ..., vₙ` is said to be
*quasi-orthogonal* if we have the inner product `⟨vᵢ, vⱼ⟩ = 0` whenever `|i - j| ≥ 2` (i.e.,
if every pair of vectors at least `2` indices apart is orthogonal). This is equivalent to
the vectors' Gram matrix being tridiagonal [JP25; p. 313](@cite).

# Properties
- `k::Int`: the `k`-orthogonality parameter. (Necessarily `2`.)

# Supertype Hierarchy
`_QuasiOrthogonality` <: [`_KOrthogonality`](@ref) <: Any

# Notes
Since this type is semantically unique inasmuch as it simply encodes the information
`k = 2`, we define the singleton instance [`_QUASI_ORTHOGONALITY`](@ref) to always be used
when needed in favor of instantiating an identical object every time.
"""
struct _QuasiOrthogonality <: _KOrthogonality end

function Base.getproperty(::_QuasiOrthogonality, name::Symbol)
    if name == :k
        return 2
    else
        # This should throw an error, unless we later add fields to this type
        return getfield(_QuasiOrthogonality, name)
    end
end

"""
    _QUASI_ORTHOGONALITY::_QuasiOrthogonality

A singleton instance of the semantically unique [`_QuasiOrthogonality`](@ref) type.

This constant is used to represent the property of quasi-orthogonality for a collection of
vectors and should always be used in favor of instantiating a new object with
`_QuasiOrthogonality()` every time.

# Notes
Simply by definition, we are certain that the `k` property of any
[`_QuasiOrthogonality`](@ref) instance will always be `2`, so given that the struct is
immutable and we do not intend to add any other distinguishing fields/properties, it is
idiomatic to define this sort of singleton instance for universal use. This design choice
also avoids unnecessary allocations (although any overhead will be negligible in practice).
"""
const _QUASI_ORTHOGONALITY = _QuasiOrthogonality()

"""
    struct _WeakOrthogonality

The property of "weak orthogonality" for a collection of vectors.

In particular, "weak orthogonality" is an *ad hoc* term used to refer to the property of
`k`-orthogonality for `k > 2`. Recall that an (ordered) collection of vectors
`v₁, v₂, ..., vₙ` is said to be *`k`-orthogonal* if we have the inner product `⟨vᵢ, vⱼ⟩ = 0`
whenever `|i - j| ≥ k` (i.e., if every pair of vectors at least `k` indices apart is
orthogonal). This is equivalent to the vectors' Gram matrix having bandwidth at most `k`.

# Fields
- `k::Int`: the `k`-orthogonality parameter. (Necessarily a positive integer greater than
    `2`.)

# Supertype Hierarchy
`_WeakOrthogonality` <: [`_KOrthogonality`](@ref) <: Any

# Constructors
- `_WeakOrthogonality(::Integer)`: constructs a `_WeakOrthogonality` object with the given
    `k`-orthogonality parameter, casting it to an `Int` for type stability. Throws a
    `DomainError` if the argument is less than or equal to `2`.

# Notes
The term "weak orthogonality" is not standard terminology in the literature, but it is used
here to emphasize the weaker nature of this property compared to orthogonality and
quasi-orthogonality. It is an *ad hoc* term coined for this module and is not intended to
be formally introduced in the broader literature. See also the documentation for parent type
[`_KOrthogonality`](@ref) for further discussion of this topic, as well as justification for
enforcing the `k` field to be an `Int` rather than the more generic `Integer`.
"""
struct _WeakOrthogonality <: _KOrthogonality
    k::Int

    function _WeakOrthogonality(k::Integer)
        if k <= 2
            throw(
                DomainError(
                    k, "k-orthogonality parameter must be a positive integer greater than 2"
                ),
            )
        end

        return new(Int(k))
    end
end

"""
    _classify_orthogonality_property(k)

Classifies the `k`-orthogonality property based on the given `k` parameter.

When searching for a `k`-orthogonal *S*-basis of a given Laplacian eigenspace, the family of
values to which our `k` parameter belongs informs our choice of algorithm.

# Arguments
- `k::Integer`: the `k`-orthogonality parameter to classify. Must be a positive integer.

# Returns
- `::_KOrthogonality`: An instance of a concrete [`_KOrthogonality`](@ref) subtype
    corresponding to `k`.

# Throws
- `DomainError`: if `k` is not a positive integer.

# Examples
```jldoctest; setup = :(using SDiagonalizability)
julia> SDiagonalizability._classify_orthogonality_property(1)
SDiagonalizability._Orthogonality()

julia> SDiagonalizability._classify_orthogonality_property(2)
SDiagonalizability._QuasiOrthogonality()

julia> SDiagonalizability._classify_orthogonality_property(3)
SDiagonalizability._WeakOrthogonality(3)

julia> SDiagonalizability._classify_orthogonality_property(13)
SDiagonalizability._WeakOrthogonality(13)
```

# Notes
See the [`_find_k_orthogonal_basis`](@ref) documentation for more details on how this
function is used in the context of finding a `k`-orthogonal basis.
"""
function _classify_orthogonality_property(k::Integer)
    if k < 1
        throw(DomainError(k, "k-orthogonality parameter must be a positive integer"))
    end

    if k == 1
        prop = _ORTHOGONALITY
    elseif k == 2
        prop = _QUASI_ORTHOGONALITY
    else
        prop = _WeakOrthogonality(k)
    end

    return prop
end

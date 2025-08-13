# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractSDiagResult

Abstract base type for all ``S``-diagonalizability and ``S``-bandwidth problem results.

# Interface
Concrete subtypes of `AbstractSDiagResult` *must* implement parametric types
- `A<:Union{AbstractGraph,AbstractMatrix{<:Integer}}`;
- `B<:Tuple{Vararg{Integer}}`; and
- `C<:Union{Nothing,Eigen}`,

alongside the following fields:
- `network::A`: the network whose `S`-bandwidth is investigated;
- `S::B`: the set from whose entries we are allowed to construct eigenvectors;
- `diagonalization::C`: an `S`-diagonalization of the matrix representation of the network,
    if it satisfies the specified `S`-bandwidth constraints; otherwise, `nothing`.
"""
abstract type AbstractSDiagResult end

"""
    SBandMinimizationResult{A,B,C,D} <: AbstractSDiagResult

[TODO: Write here]

# Supertype Hierarchy
`SBandMinimizationResult` <: [`AbstractSDiagResult`](@ref)
"""
struct SBandMinimizationResult{
    A<:Union{AbstractGraph,AbstractMatrix{<:Integer}},
    B<:Tuple{Vararg{Integer}},
    C<:Union{Nothing,Eigen},
    D<:Union{Int,Float64},
} <: AbstractSDiagResult
    network::A
    S::B
    s_diagonalization::C
    s_bandwidth::D
end

#= The `Base.show` override here takes heavy inspiration from the `MatrixBandwidth.jl`
package (written by myself), which in turn took inspiration from `Optim.jl`. =#
function Base.show(io::IO, res::SBandMinimizationResult{A}) where {A}
    if A<:AbstractGraph
        n = nv(res.network)
    else
        n = size(res.network, 1)
    end

    println(io, "Results of S-Bandwidth Minimization")
    println(io, " * S: $(Int.(res.S))")
    println(io, " * S-Bandwidth: $(res.s_bandwidth)")
    print(io, " * Graph Order: $n")

    return nothing
end

"""
    SBandRecognitionResult{A,B,C} <: AbstractSDiagResult

[TODO: Write here]

# Supertype Hierarchy
`SBandRecognitionResult` <: [`AbstractSDiagResult`](@ref)
"""
struct SBandRecognitionResult{
    A<:Union{AbstractGraph,AbstractMatrix{<:Integer}},
    B<:Tuple{Vararg{Integer}},
    C<:Union{Nothing,Eigen},
} <: AbstractSDiagResult
    network::A
    S::B
    s_diagonalization::C
    k::Integer
    s_band_at_most_k::Bool
end

#= The `Base.show` override here takes heavy inspiration from the `MatrixBandwidth.jl`
package (written by myself), which in turn took inspiration from `Optim.jl`. =#
function Base.show(io::IO, res::SBandRecognitionResult{A}) where {A}
    if A<:AbstractGraph
        n = nv(res.network)
    else
        n = size(res.network, 1)
    end

    println(io, "Results of S-Bandwidth Recognition")
    println(io, " * S: $(Int.(res.S))")
    println(io, " * S-Bandwidth Threshold k: $(res.k)")
    println(io, " * Has S-Bandwidth ≤ k: $(res.s_band_at_most_k)")
    print(io, " * Graph Order: $n")

    return nothing
end

"""
    SDiagonalizabilityResult{A,B,C} <: AbstractSDiagResult

[TODO: Write here]

# Supertype Hierarchy
`SDiagonalizabilityResult` <: [`AbstractSDiagResult`](@ref)
"""
struct SDiagonalizabilityResult{
    A<:Union{AbstractGraph,AbstractMatrix{<:Integer}},
    B<:Tuple{Vararg{Integer}},
    C<:Union{Nothing,Eigen},
} <: AbstractSDiagResult
    network::A
    S::B
    s_diagonalization::C
    has_s_diagonalization::Bool
end

#= The `Base.show` override here takes heavy inspiration from the `MatrixBandwidth.jl`
package (written by myself), which in turn took inspiration from `Optim.jl`. =#
function Base.show(io::IO, res::SDiagonalizabilityResult{A}) where {A}
    if A<:AbstractGraph
        n = nv(res.network)
    else
        n = size(res.network, 1)
    end

    println(io, "Results of S-Diagonalizability Check")
    println(io, " * S: $(Int.(res.S))")
    println(io, " * S-Diagonalizable: $(res.has_s_diagonalization)")
    print(io, " * Graph Order: $n")

    return nothing
end

"""
    SSpectra{A,B,C,D,E}

[TODO: Write here]
"""
struct SSpectra{
    A<:AbstractMatrix{<:Integer},
    B<:Tuple{Vararg{Integer}},
    C<:Union{Nothing,OrderedDict{Int,Int}},
    D<:Union{Nothing,OrderedDict{Int,<:AbstractMatrix{<:Integer}}},
    E<:Union{Nothing,OrderedDict{Int,<:Union{Nothing,AbstractMatrix{<:Integer}}}},
}
    matrix::A
    S::B
    multiplicities::C
    s_eigenspaces::D
    s_eigenbases::E
    s_diagonalizable::Bool
end

"""
    struct SpectrumIntegralResult

Data on whether a matrix is spectrum integral (i.e., whether its eigenvalues are integers).

This struct also contains a map from each eigenvalue to its multiplicity, provided that the
eigenvalues are indeed all integers. (Otherwise, the associated field is simply `nothing`.)

# Fields
- `matrix::Matrix{Int}`: the matrix whose eigenvalues and their integrality are of interest.
- `spectrum_integral::Bool`: whether the eigenvalues of `matrix` are all integers.
- `multiplicities::Union{Nothing,OrderedDict{Int,Int}}`: a map from each eigenvalue to its
    multiplicity, sorted first by ascending multiplicity then by ascending eigenvalue. (This
    field is `nothing` if and only if `spectrum_integral` is false, since we cannot map
    non-integer eigenvalues to data.)

# Notes
If an undirected graph with integer edge weights is ``\\{-1, 0, 1\\}``-diagonalizable (or,
more restrictively, ``\\{-1, 1\\}``-diagonalizable), then its Laplacian matrix has integer
eigenvalues [JP25, p. 312]. Hence, validating Laplacian integrality serves as a useful
screening step in this package's principal ``S``-bandwidth minimization algorithm.

# References
- [JP25](@cite): N. Johnston and S. Plosker. *Laplacian {−1,0,1}- and {−1,1}-diagonalizable
    graphs*. Linear Algebra and its Applications **704**, 309–39 (2025).
    https://doi.org/10.1016/j.laa.2024.10.016.
"""
struct SpectrumIntegralResult{T<:Union{Nothing,OrderedDict{Int,Int}}}
    matrix::AbstractMatrix{<:Integer}
    spectrum_integral::Bool
    multiplicities::T
end

"""
    NotImplementedError{Nothing}(f, subtype, abstracttype)
    NotImplementedError{Symbol}(f, arg, subtype, abstracttype)

An exception indicating that a function lacks dispatch to handle a specific argument type.

Semantically, this differs from `MethodError` in that it connotes a developer-side failure
to implement a method rather than erroneous user input. Throughout this package, it is often
used to warn when an existing function with multiple dispatch on some abstract type is
called on a newly created subtype for which no method has been defined.

# Fields
- `f::Function`: the function called.
- `arg::Symbol`: the name of the argument with the unsupported type, if the function has
    multiple arguments. If the function has only one argument, this field should be set to
    `nothing`.
- `subtype::Type`: the type of the argument. May be the actual concrete type or some
    intermediate supertype. (For instance, if the relevant input has concrete type `A` with
    hierarchy `A <: B <: C` and the `abstracttype` field is `C`, then both `A` and `B` are
    perfectly valid choices for `subtype`.)
- `abstracttype::Type`: the abstract type under which the argument is meant to fall.

# Supertype Hierarchy
`NotImplementedError` <: `Exception`

# Constructors
- `NotImplementedError(::Function, ::Type, ::Type)`: constructs a new `NotImplementedError`
    instance for a single-argument function. Throws an error if the second type is not
    abstract or the first type is not a subtype of the second.
- `NotImplementedError(::Function, ::Symbol, ::Type, ::Type)`: constructs a new
    `NotImplementedError` instance for a multi-argument function. Throws an error if the
    second type is not abstract or the first type is not a subtype of the second.
"""
struct NotImplementedError{T<:Union{Nothing,Symbol}} <: Exception
    f::Function
    arg::T
    subtype::Type
    abstracttype::Type

    function NotImplementedError(f::Function, subtype::Type, abstracttype::Type)
        return NotImplementedError(f, nothing, subtype, abstracttype)
    end

    function NotImplementedError(
        f::Function, arg::T, subtype::Type, abstracttype::Type
    ) where {T<:Union{Nothing,Symbol}}
        if !isabstracttype(abstracttype)
            throw(ArgumentError("Expected an abstract type, got $abstracttype"))
        end

        if !(subtype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $subtype"))
        end

        return new{T}(f, arg, subtype, abstracttype)
    end
end

function Base.showerror(io::IO, e::NotImplementedError{Nothing})
    print(
        io,
        """NotImplementedError with $(e.subtype):
        $(e.f) is not yet implemented for this subtype of $(e.abstracttype).
        Try defining method dispatch manually if this is a newly created subtype.""",
    )
    return nothing
end

function Base.showerror(io::IO, e::NotImplementedError{Symbol})
    print(
        io,
        """NotImplementedError with argument $(e.arg)::$(e.subtype):
        $(e.f) is not yet implemented for this subtype of $(e.abstracttype).
        Try defining method dispatch manually if this is a newly created subtype.""",
    )
    return nothing
end

"""
    EfficiencyWarning(msg)

[TODO: Write here]
"""
struct EfficiencyWarning
    msg::String
end

function Base.show(io::IO, w::EfficiencyWarning)
    print(io, "EfficiencyWarning: $(w.msg)")
    return nothing
end

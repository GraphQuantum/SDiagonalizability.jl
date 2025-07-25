# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    AbstractSBandResult

Abstract base type for all *S*-bandwidth problem results.

# Interface
Concrete subtypes of `AbstractSBandResult` *must* implement parametric types
- `A<:Union{AbstractGraph,AbstractMatrix{<:Integer}}`;
- `B<:Tuple`; and
- `C<:Union{Nothing,Eigen}`,

alongside the following fields:
- `network::A`: the network whose *S*-bandwidth is investigated;
- `S::B`: the set *S* from whose entries we are allowed to construct eigenvectors;
- `diagonalization::C`: an *S*-diagonalization of the matrix representation of the network,
    if it satisfies the specified *S*-bandwidth constraints; otherwise, `nothing`.
"""
abstract type AbstractSBandResult end

"""
    SBandMinimizationResult{A,B,C,D} <: AbstractSBandResult

[TODO: Write here]
"""
struct SBandMinimizationResult{
    A<:Union{AbstractGraph,AbstractMatrix{<:Integer}},
    B<:Tuple,
    C<:Union{Nothing,Eigen},
    D<:Union{Int,Float64},
} <: AbstractSBandResult
    network::A
    S::B
    diagonalization::C
    band::D
end

"""
    SBandRecognitionResult{A,B,C} <: AbstractSBandResult

[TODO: Write here]
"""
struct SBandRecognitionResult{
    A<:Union{AbstractGraph,AbstractMatrix{<:Integer}},B<:Tuple,C<:Union{Nothing,Eigen}
} <: AbstractSBandResult
    network::A
    S::B
    diagonalization::C
    k::Int
    has_band_k_diag::Bool
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
If an undirected graph with integer edge weights is `{-1,0,1}`-diagonalizable (or, more
restrictively, `{-1,1}`-diagonalizable), then its Laplacian matrix has integer eigenvalues
[JP25; p. 312](@cite). Hence, validating Laplacian integrality serves as a useful screening
step in this package's principal *S*-bandwidth minimization algorithm.
"""
struct SpectrumIntegralResult{T<:Union{Nothing,OrderedDict{Int,Int}}}
    matrix::AbstractMatrix{<:Integer}
    spectrum_integral::Bool
    multiplicities::T
end

"""
    _SSpectra{A,B,C,D}

[TODO: Write here]
"""
struct _SSpectra{A<:Tuple,B<:Union{Nothing,OrderedDict{Int}},C<:B,D<:B}
    matrix::AbstractMatrix{<:Integer}
    S::A
    multiplicities::B
    s_eigenspaces::C
    s_eigenbases::D
    s_diagonalizable::Bool
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

# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

# TODO: Add comments and docstrings to this file

const _LAPLACIAN_ERR_MSGS = Dict(
    "asymmetric" => "Matrix is not symmetric; cannot be an (undirected) Laplacian",
    "row_sums" => "Matrix has nonzero row sums; cannot be a Laplacian",
)

# TODO: Extend functionality to multi-argument functions
struct NotImplementedError <: Exception
    f::Function
    concretetype::Type
    abstracttype::Type

    function NotImplementedError(f::Function, concretetype::Type, abstracttype::Type)
        if !isconcretetype(concretetype)
            throw(ArgumentError("Expected a concrete type, got $concretetype"))
        end

        if !isabstracttype(abstracttype)
            throw(ArgumentError("Expected an abstract type, got $abstracttype"))
        end

        if !(concretetype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $concretetype"))
        end

        return new(f, concretetype, abstracttype)
    end
end

function Base.showerror(io::IO, e::NotImplementedError)
    print(
        io,
        """NotImplementedError with `$(e.concretetype)`:
        The function `$(e.f)` is not yet implemented for this subtype of `$(e.abstracttype)`.""",
    )
end

function assert_matrix_is_undirected_laplacian(L::AbstractMatrix{<:Integer})
    issymmetric(L) || throw(DomainError(L, _LAPLACIAN_ERR_MSGS["asymmetric"]))
    # Symmetricity implies the column sums (faster to compute) are equal to the row sums
    iszero(sum(L; dims=1)) || throw(DomainError(L, _LAPLACIAN_ERR_MSGS["row_sums"]))
    return nothing
end

# TODO: Again, add more docstrings

#= This tolerance computation is adapted from NumPy's `numpy.linalg.matrix_rank`, but with
an additional square root operation to provide even more robustness. (See
`https://numpy.org/doc/2.1/reference/generated/numpy.linalg.matrix_rank.html`.) =#
function rank_rtol(A::AbstractMatrix{<:Real})
    return sqrt(maximum(size(A)) * eps())
end

function rank_rtol(A::AbstractMatrix{T}) where {T<:AbstractFloat}
    return sqrt(maximum(size(A)) * Float64(eps(T)))
end

# TODO: Almost certainly add more utilities as we develop the rest of the library

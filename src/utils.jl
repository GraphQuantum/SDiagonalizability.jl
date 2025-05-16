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
            throw(ArgumentError("$concretetype must be a concrete type"))
        end

        if !isabstracttype(abstracttype)
            throw(ArgumentError("$abstracttype must be an abstract type"))
        end

        if !(concretetype <: abstracttype)
            throw(ArgumentError("Expected a subtype of $abstracttype, got $concretetype"))
        end

        return new(f, concretetype, abstracttype)
    end
end

function Base.show(io::IO, e::NotImplementedError)
    # TODO: Implement
    println(
        "$(e.f) is not yet implemented for subtype $(e.concretetype) <: $(e.abstracttype)"
    )
end

function assert_matrix_is_undirected_laplacian(L::AbstractMatrix{<:Integer})
    issymmetric(L) || throw(DomainError(L, _LAPLACIAN_ERR_MSGS["asymmetric"]))
    # Symmetricity implies the column sums (faster to compute) are equal to the row sums
    iszero(sum(L; dims=1)) || throw(DomainError(L, _LAPLACIAN_ERR_MSGS["row_sums"]))
    return nothing
end

# TODO: Almost certainly add more utilities as we develop the rest of the library

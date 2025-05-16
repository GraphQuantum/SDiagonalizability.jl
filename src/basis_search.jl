# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

# TODO: Extract this typing logic out into a separate file

abstract type _KOrthogonality end

struct _Orthogonality <: _KOrthogonality end
struct _QuasiOrthogonality <: _KOrthogonality end
struct _WeakOrthogonality <: _KOrthogonality
    k::Int
end

function Base.getproperty(obj::_Orthogonality, name::Symbol)
    name == :k && return 1
    error("type $(typeof(obj)) has no field $name")
end

function Base.getproperty(obj::_QuasiOrthogonality, name::Symbol)
    name == :k && return 2
    error("type $(typeof(obj)) has no field $name")
end

function _classify_orthogonality_property(k::Int)
    if k < 1
        throw(DomainError(k, "k-orthogonality parameter must be a positive integer"))
    end

    if k == 1
        prop = _Orthogonality()
    elseif k == 2
        prop = _QuasiOrthogonality()
    else
        prop = _WeakOrthogonality(k)
    end

    return prop
end

# TODO: Some impl's?

function _has_k_orthogonal_basis(column_space::AbstractMatrix{Int}, k::Int, rank::Int)
    prop = _classify_orthogonality_property(k)
    return _has_basis_with_property(column_space, prop, rank)
end

function _has_basis_with_property(_::AbstractMatrix{Int}, prop::_KOrthogonality, rank::Int)
    throw(NotImplementedError(_has_basis_with_property, typeof(prop), _KOrthogonality))
end

# TODO: Make concrete implementations

function _has_basis_with_property(
    column_space::AbstractMatrix{Int}, prop::_WeakOrthogonality, rank::Int
)
    k = prop.k
    # TODO: Implement
end

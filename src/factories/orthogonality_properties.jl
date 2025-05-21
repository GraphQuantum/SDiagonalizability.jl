# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

# TODO: Add comments and docstrings to this file

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
    k < 1 && throw(DomainError(k, "k-orthogonality parameter must be a positive integer"))

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

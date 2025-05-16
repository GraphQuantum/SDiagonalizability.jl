# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

# TODO: Add docstrings to this file, and some comments to the structs and associated methods

struct SpectrumIntegralResult
    matrix::Matrix{Int}
    spectrum_integral::Bool
    multiplicities::Union{Nothing,OrderedDict{Int,Int}}
end

struct _Eigenspace01Neg
    dimension::Int
    eigvecs_01neg::AbstractMatrix{Int}
    indices_1neg::BitVector
    basis_01neg::Union{Nothing,Matrix{Int}}
    basis_1neg::Union{Nothing,Matrix{Int}}
end

struct _LaplacianSpectrum01Neg
    laplacian_matrix::Matrix{Int}
    diagonalizable_01neg::Bool
    diagonalizable_1neg::Bool
    eigspaces_01neg::Union{Nothing,OrderedDict{Int,_Eigenspace01Neg}}
end

function Base.getproperty(obj::_LaplacianSpectrum01Neg, name::Symbol)
    if name == :dimension || !(name in fieldnames(_Eigenspace01Neg))
        error("type $(typeof(obj)) has no field $name")
    end

    if name != :eigspaces_01neg && name in fieldnames(_LaplacianSpectrum01Neg)
        value = getfield(obj, name)
    elseif obj.diagonalizable_01neg
        eigspaces_01neg = getfield(obj, :eigspaces_01neg)
        name == :multiplicities && (name = :dimension)
        value = OrderedDict(
            eigval => getfield(eigspace, name) for (eigval, eigspace) in eigspaces_01neg
        )
    else
        value = nothing
    end

    return value
end

function check_spectrum_integrality(X::AbstractMatrix{<:Integer})
    X_copy = Matrix{Int}(X) # Avoid shared mutability and cast to `Matrix{Int}`
    eigvals_float = eigvals(X_copy)
    eigvals_int = round.(Int, real.(eigvals_float))
    spectrum_integral = isapprox(eigvals_float, eigvals_int)

    # Sort first by ascending multiplicity then by ascending eigenvalue
    if spectrum_integral
        multiplicities = OrderedDict(sort!(collect(counter(eigvals_int)); by=reverse))
    else
        multiplicities = OrderedDict{Int,Int}()
    end

    SpectrumIntegralResult(X_copy, spectrum_integral, multiplicities)
end

function _extract_independent_cols(E::AbstractMatrix{<:Integer})
    # Factorize `E = QR` for some orthogonal matrix `Q` and upper triangular matrix `R`
    F = qr(E, ColumnNorm())
    pivots = F.p # The first `rank(E)` pivots correspond to independent columns of `E`
    ortho_coefs = diag(F.R) # Scaling factors of the columns of the `Q` matrix

    #= The tolerance is taken from NumPy's `numpy.linalg.matrix_rank` function, which uses
    singular values. It is easily adapted to scaling coefficients from a QR factorization.
    (See `https://numpy.org/doc/2.1/reference/generated/numpy.linalg.matrix_rank.html`.) =#
    tol = maximum(abs, ortho_coefs) * maximum(size(E)) * eps()
    r = count(x -> abs(x) > tol, ortho_coefs) # The rank of `E` within floating-point error
    return E[:, pivots[1:r]] # A largest independent subset of the columns of `E`
end

function _find_indices_1neg(eigvecs_01neg::AbstractDict{Int,<:AbstractMatrix{Int}})
    #= The order of the Laplacian matrix. Accessing the zero key is safe, as every
    undirected graph has a rank-deficient Laplacian. =#
    n = size(eigvecs_01neg[0], 1)

    # Only even-ordered Laplacians have non-kernel {-1,1}-eigenvectors
    if n % 2 == 0 # Check all eigenspaces for {-1,1}-eigenvectors
        indices_1neg = Dict(
            eigval => findall(v -> all(v .!= 0), eachcol(eigvecs)) for
            (eigval, eigvecs) in eigvecs_01neg
        )
    else # Only check the kernel for {-1,1}-eigenvectors
        indices_1neg = Dict(eigval => Int[] for eigval in keys(eigvecs_01neg))
        indices_1neg[0] = findall(v -> all(v .!= 0), eachcol(eigvecs_01neg[0]))
    end

    return indices_1neg # Indices of the {-1,0,1}-eigenvectors without 0's
end

function _laplacian_spectra_01neg(L::AbstractMatrix{<:Integer})
    TL = _cast_to_typed_laplacian(L)
    return _typed_laplacian_spectra_01neg(TL)
end

function _typed_laplacian_spectra_01neg(TL::_TypedLaplacian)
    throw(NotImplementedError(_typed_laplacian_spectra_01neg, typeof(TL), _TypedLaplacian))
end

function _typed_laplacian_spectra_01neg(TL::_NullGraphLaplacian)
    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()
    return _LaplacianSpectrum01Neg(TL.matrix, true, true, eigspaces_01neg)
end

function _typed_laplacian_spectra_01neg(TL::_EmptyGraphLaplacian)
    L = TL.matrix
    n = size(L, 1)

    kernel_eigvecs_01neg = hcat(_pot_kernel_eigvecs_01neg(n)...)
    eigvecs_01neg = Dict(0 => kernel_eigvecs_01neg)
    indices_1neg = _find_indices_1neg(eigvecs_01neg)
    pot_kernel_basis_1neg = _extract_independent_cols(kernel_eigvecs_01neg)

    if size(pot_kernel_basis_1neg, 2) == n # The potential {-1,1}-basis spans the kernel
        kernel_basis_01neg = pot_kernel_basis_1neg
        kernel_basis_1neg = copy(pot_kernel_basis_1neg) # Avoid shared mutability
        diagonalizable_1neg = true
    else # No {-1,1}-basis exists, so compute a {-1,0,1}-eigenbasis (guaranteed to exist)
        kernel_basis_01neg = _extract_independent_cols(kernel_eigvecs_01neg)
        kernel_basis_1neg = nothing
        diagonalizable_1neg = false
    end

    eigspaces_01neg = OrderedDict(
        0 => _Eigenspace01Neg(
            n,
            kernel_eigvecs_01neg,
            indices_1neg[0],
            kernel_basis_01neg,
            kernel_basis_1neg,
        ),
    )
    return _LaplacianSpectrum01Neg(L, true, diagonalizable_1neg, eigspaces_01neg)
end

function _typed_laplacian_spectra_01neg(TL::_CompleteGraphLaplacian)
    L = TL.matrix
    n = size(L, 1)

    #= For every `n ≥ 2`, `Kₙ` (assuming a uniform edge weight of `w`) has eigenvalues 0
    with multiplicity 1 and `n * w` with multiplicity `n - 1`. (`n` is guaranteed to be at
    least 2 in this context due to the Laplacian type casting logic.) =#
    eigval_nonzero = n * TL.weight

    kernel_eigvecs_01neg = ones(Int, n, 1) # The all-ones vector spans the kernel
    # The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel
    nonkernel_eigvecs_01neg = hcat(_pot_nonkernel_eigvecs_01neg(n)...)

    # All the {-1,0,1}-eigenvectors corresponding to each eigenvalue, as columns of a matrix
    eigvecs_01neg = Dict(
        0 => kernel_eigvecs_01neg, eigval_nonzero => nonkernel_eigvecs_01neg
    )
    indices_1neg = _find_indices_1neg(eigvecs_01neg)

    # Again, the set containing the all-ones vector forms a basis for the kernel
    kernel_basis_01neg = copy(kernel_eigvecs_01neg) # Avoid shared mutability
    kernel_basis_1neg = copy(kernel_eigvecs_01neg) # Avoid shared mutability
    pot_nonkernel_basis_1neg = _extract_independent_cols(
        nonkernel_eigvecs_01neg[:, indices_1neg[eigval_nonzero]]
    )

    #= Find a {-1,0,1}- and, if one exists, a {-1,1}-basis for the remaining eigenspace.
    {-1,1}-bases are preferable, so they take priority whenever they exist. =#
    if size(pot_nonkernel_basis_1neg, 2) == n - 1 # The potential {-1,1}-basis spans the space
        nonkernel_basis_01neg = pot_nonkernel_basis_1neg
        nonkernel_basis_1neg = copy(pot_nonkernel_basis_1neg) # Avoid shared mutability
        diagonalizable_1neg = true
    else # No {-1,1}-basis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist). =#
        nonkernel_basis_01neg = _extract_independent_cols(nonkernel_eigvecs_01neg)
        nonkernel_basis_1neg = nothing
        diagonalizable_1neg = false
    end

    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()
    kernel_01neg = _Eigenspace01Neg(
        n, kernel_eigvecs_01neg, indices_1neg[0], kernel_basis_01neg, kernel_basis_1neg
    )
    nonkernel_01neg = _Eigenspace01Neg(
        n - 1,
        nonkernel_eigvecs_01neg,
        indices_1neg[eigval_nonzero],
        nonkernel_basis_01neg,
        nonkernel_basis_1neg,
    )

    # Sort first by ascending multiplicity then by ascending eigenvalue
    if n == 2 && eigval_nonzero < 0
        eigspaces_01neg[eigval_nonzero] = nonkernel_01neg
        eigspaces_01neg[0] = kernel_01neg
    else
        eigspaces_01neg[0] = kernel_01neg
        eigspaces_01neg[eigval_nonzero] = nonkernel_01neg
    end

    return _LaplacianSpectrum01Neg(L, true, diagonalizable_1neg, eigspaces_01neg)
end

function _typed_laplacian_spectra_01neg(TL::_ArbitraryGraphLaplacian)
    L = TL.matrix
    res = check_spectrum_integrality(L)

    if !res.spectrum_integral
        # For integer edge weights, {-1,0,1}-diagonalizability implies Laplacian integrality
        return _LaplacianSpectrum01Neg(L, false, false, nothing)
    end

    n = size(L, 1)
    multiplicities = res.multiplicities
    eigvals_nonzero = filter(!iszero, keys(multiplicities))
    eigvecs_01neg = Dict{Int,AbstractMatrix{Int}}(
        # Initialize resizeable matrices to store an unknown number of eigenvectors
        eigval => ElasticMatrix{Int}(undef, n, 0) for eigval in eigvals_nonzero
    )

    # The kernel is simpler and handled differently from the other eigenspaces
    if multiplicities[0] == 1 # The all-ones vector spans the kernel
        eigvecs_01neg[0] = ones(Int, n, 1)
    else # Check all {-1,0,1}-vectors in ℝⁿ, unique up to span 
        eigvecs_01neg[0] = hcat(
            Iterators.filter(v -> iszero(L * v), _pot_kernel_eigvecs_01neg(n))...
        )
    end

    # TODO: Multithread this; lazy evaluation means speed will become an issue before memory
    # Now fill up the remaining (non-kernel) eigenspaces
    for v in _pot_kernel_eigvecs_01neg(n)
        for eigval in eigvals_nonzero
            if L * v == eigval * v
                append!(eigvecs_01neg[eigval], v)
                break
            end
        end
    end

    indices_1neg = _find_indices_1neg(eigvecs_01neg)

    bases_01neg = Dict{Int,Union{Nothing,Matrix{Int}}}(
        eigval => _extract_independent_cols(vecs) for (eigval, vecs) in eigvecs_01neg
    )
    bases_1neg = Dict(
        eigval => _extract_independent_cols(vecs[:, indices_1neg[eigval]]) for
        (eigval, vecs) in eigvecs_01neg
    )
    diagonalizable_01neg = true
    diagonalizable_1neg = true

    #= If they exist, find a {-1,0,1}- and {-1,1}-basis for each non-kernel eigenspace.
    {-1,1}-bases are preferable, so they take priority whenever they exist. =#
    for eigval in eigvals_nonzero
        multiplicity = multiplicities[eigval]
        pot_basis_01neg = bases_01neg[eigval]
        pot_basis_1neg = bases_1neg[eigval]

        if size(pot_basis_01neg, 2) < multiplicity # The eigenspace has no {-1,0,1}-basis
            bases_01neg[eigval] = bases_1neg[eigval] = nothing
            diagonalizable_01neg = diagonalizable_1neg = false
            # There exists a {-1,0,1}- but not a {-1,1}-eigenbasis
        elseif size(pot_basis_1neg, 2) < multiplicity
            bases_1neg[eigval] = nothing
            diagonalizable_1neg = false
        else # There exists a {-1,1}-eigenbasis, so use it in both dictionaries
            bases_01neg[eigval] = copy(pot_basis_1neg) # Avoid shared mutability
        end
    end

    # Create and populate our usual wrapper for the eigenspace data
    eigspaces_01neg = OrderedDict{Int,_Eigenspace01Neg}()

    for (eigval, multiplicity) in multiplicities
        eigspaces_01neg[eigval] = _Eigenspace01Neg(
            multiplicity,
            eigvecs_01neg[eigval],
            indices_1neg[eigval],
            bases_01neg[eigval],
            bases_1neg[eigval],
        )
    end

    return _LaplacianSpectrum01Neg(
        L, res.spectrum_integral, diagonalizable_1neg, eigspaces_01neg
    )
end

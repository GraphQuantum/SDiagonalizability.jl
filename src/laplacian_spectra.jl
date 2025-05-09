struct LaplacianIntegralResult
    laplacian_matrix::Symmetric{Int,Matrix{Int}}
    laplacian_integral::Bool
    eigvals_sorted::Vector{Int}
    eigval_counts::OrderedDict{Int,Int}
end

struct LaplacianSpectra
    laplacian_matrix::Symmetric{Int,Matrix{Int}}
    eigvals_sorted::Vector{Int}
    eigval_counts::Dict{Int,Int}
    eigspaces_01neg::Dict{Int,Matrix{Int}}
    idxs_1neg::Dict{Int,Vector{Int}}
    bases_01neg::Dict{Int,Matrix{Int}}
    bases_1neg::Dict{Int,Matrix{Int}}

    function LaplacianSpectra(
        res::LaplacianIntegralResult,
        eigspaces_01neg::Dict{Int,Matrix{Int}},
        idxs_1neg::Dict{Int,Vector{Int}},
        bases_01neg::Dict{Int,Matrix{Int}},
        bases_1neg::Dict{Int,Matrix{Int}},
    )
        return new(
            res.laplacian_matrix,
            res.eigvals_sorted,
            Dict(res.eigval_counts),
            eigspaces_01neg,
            idxs_1neg,
            bases_01neg,
            bases_1neg,
        )
    end
end


function laplacian_integrality(L::Symmetric{Int,Matrix{Int}})
    eigvals_float = eigvals(L)
    eigvals_int = round.(Int, eigvals_float) # Necessarily real by symmetricity of L
    laplacian_integral = isapprox(eigvals_float, eigvals_int)

    if !laplacian_integral
        eigvals_sorted = Int[]
        eigval_counts = OrderedDict{Int,Int}()
    else
        # Sort unique eigenvalues by ascending multiplicity then b y ascending value
        eigval_counts = OrderedDict(sort!(collect(counter(eigvals_int)); by=reverse))
        # All eigenvalues (including repetitions) in the same order as `eigval_counts`
        eigvals_sorted = collect(Iterators.flatmap(pair -> fill(pair...), eigval_counts))
    end

    return LaplacianIntegralResult(L, laplacian_integral, eigvals_sorted, eigval_counts)
end

function laplacian_spectra(L::Symmetric{Int,Matrix{Int}})
    n = size(L, 1)

    if n == 0
        res = LaplacianIntegralResult(L, true, Int[], OrderedDict{Int,Int}())
        return _null_graph_laplacian_spectra(res)
    end

    if iszero(L)
        res = LaplacianIntegralResult(L, true, zeros(Int, n), OrderedDict(0 => n))
        return _empty_graph_laplacian_spectra(res)
    end

    res = laplacian_integrality(L)

    if !res.spectrum_integral
        return LaplacianSpectra(
            res,
            Dict{Int,Matrix{Int}}(),
            Dict{Int,Vector{Int}}(),
            Dict{Int,Matrix{Int}}(),
            Dict{Int,Matrix{Int}}(),
        )
    end

    if _is_complete_graph_res(res)
        return _complete_graph_laplacian_spectra(res)
    end

    return _arbitrary_graph_laplacian_spectra(res)
end


@inline function _extract_independent_cols(E::Matrix{Int})
    F = qr(E, ColumnNorm()) # Factorize E = QR for some orthogonal Q and upper triangular R
    pivots = F.p # The first `rank(E)` pivots correspond to independent columns of E
    ortho_coefs = diag(F.R) # Scaling factors of the columns of the Q matrix

    # The tolerance is taken from NumPy's `numpy.linalg.matrix_rank` function, which uses
    # singular values. It is easily adapted to scaling coefficients from a QR factorization.
    # (See `https://numpy.org/doc/2.1/reference/generated/numpy.linalg.matrix_rank.html`.)
    tol = maximum(abs, ortho_coefs) * maximum(size(E)) * eps()
    r = count(x -> abs(x) > tol, ortho_coefs) # The rank of E within floating-point error
    return E[:, pivots[1:r]] # A largest independent subset of the columns of E
end

@inline function _find_idxs_1neg(eigspaces_01neg::Dict{Int,Matrix{Int}}, n::Int)
    # Only even-ordered Laplacians have non-kernel {-1,1}-eigenvectors
    if n % 2 == 0 # Check all eigenspaces for {-1,1}-eigenvectors
        idxs_1neg = Dict(
            eigval => findall(v -> all(v .!= 0), eachcol(eigvecs)) for
            (eigval, eigvecs) in eigspaces_01neg
        )
    else # Only check the kernel for {-1,1}-eigenvectors
        idxs_1neg = Dict(eigval => Int[] for eigval in keys(eigspaces_01neg))
        idxs_1neg[0] = findall(v -> all(v .!= 0), eachcol(eigspaces_01neg[0]))
    end

    return idxs_1neg # Indices of the {-1,0,1}-eigenvectors without 0's
end


@inline function _null_graph_laplacian_spectra(res::LaplacianIntegralResult)
    eigspaces_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    idxs_1neg = Dict(0 => Int[])
    bases_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    bases_1neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

@inline function _empty_graph_laplacian_spectra(res::LaplacianIntegralResult)
    n = size(res.laplacian_matrix, 1)
    eigspaces_01neg = Dict(0 => Matrix{Int}(undef, n, 0))
    idxs_1neg = Dict(0 => Int[])
    eigspaces_01neg = Dict(0 => hcat(_potential_kernel_eigvecs(n)...))
    idxs_1neg = _find_idxs_1neg(eigspaces_01neg, n)
    pot_basis_1neg = _extract_independent_cols(eigspaces_01neg[0])

    if size(pot_basis_1neg, 2) == n # The potential {-1,1}-basis spans the kernel
        bases_1neg = Dict(0 => pot_basis_1neg)
        bases_01neg = Dict(0 => copy(pot_basis_1neg)) # Avoid shared mutability
    else # No {-1,1}-basis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        bases_1neg = Dict(0 => Matrix{Int}(undef, n, 0)) # Use an empty matrix as a sentinel
        # Compute a basis with 0's
        bases_01neg = Dict(0 => _extract_independent_cols(eigspaces_01neg[0]))
    end

    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _complete_graph_laplacian_spectra(res::LaplacianIntegralResult)
    #= For every `n ≥ 2`, `Kₙ` (assuming uniform edge weights) has one zero eigenvalue with
    multiplicity 1 and one nonzero eigenvalue with multiplicity `n - 1`. (`n` is always
    guaranteed to be at least 2 in whatever context this helper function is called.) =#
    n = size(res.laplacian_matrix, 1)
    val1, val2 = keys(res.eigval_counts) # There are precisely two eigenvalues, so this is safe
    eigval_nonzero = val1 == 0 ? val2 : val1

    kernel_01neg = ones(Int, n, 1) # The all-ones vector spans the (one-dimensional) kernel
    # The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel
    nonkernel_01neg = hcat(_potential_nonkernel_eigvecs(n)...)

    # All the {-1,0,1}-eigenvectors corresponding to each eigenvalue, as columns of a matrix
    eigspaces_01neg = Dict(0 => kernel_01neg, eigval_nonzero => nonkernel_01neg)
    idxs_1neg = _find_idxs_1neg(eigspaces_01neg, n) # The indices of the {-1,1}-columns

    #= Find {-1,0,1}- and, if they exist, {-1,1}-eigenbases. {-1,1}-bases are preferable, so
    use them wherever applicable. =#
    null_basis_1neg = copy(kernel_01neg) # Avoid shared mutability
    null_basis_01neg = copy(kernel_01neg) # Avoid shared mutability
    nonnull_idxs_1neg = idxs_1neg[eigval_nonzero]
    nonnull_pot_basis_1neg = _extract_independent_cols(nonkernel_01neg[:, nonnull_idxs_1neg])

    if size(nonnull_pot_basis_1neg, 2) == n - 1 # The potential {-1,1}-basis spans the space
        nonnull_basis_1neg = nonnull_pot_basis_1neg
        nonnull_basis_01neg = copy(nonnull_pot_basis_1neg) # Avoid shared mutability
    else # No {-1,1}-eigenbasis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        nonnull_basis_1neg = Matrix{Int}(undef, n, 0) # Use an empty matrix as a sentinel
        nonnull_basis_01neg = _extract_independent_cols(nonkernel_01neg) # Compute a basis with 0's
    end

    #= A {-1,0,1}- and, if it exists, a {-1,1}-basis for each eigenspace. If no {-1,1}-basis
    exists for the non-kernel eigenspace, an empty matrix is used as a placeholder. (The
    kernel always has one, as it is spanned by the all-ones vector.) =#
    bases_01neg = Dict(0 => null_basis_01neg, eigval_nonzero => nonnull_basis_01neg)
    bases_1neg = Dict(0 => null_basis_1neg, eigval_nonzero => nonnull_basis_1neg)

    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _arbitrary_graph_laplacian_spectra(res::LaplacianIntegralResult)
    # TODO: Implement
    return LaplacianSpectra(res, Dict(), Dict(), Dict(), Dict())
end

function _is_complete_graph_res(res::LaplacianIntegralResult)
    n = size(res.laplacian_matrix, 1)
    return length(res.eigval_counts) == 2 && all(values(res.eigval_counts) .== [1, n - 1])
end

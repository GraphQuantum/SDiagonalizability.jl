# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

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
    eigvals_int = round.(Int, eigvals_float) # Necessarily real by symmetricity of `L`
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

    if n == 0 # `L` represents the null graph on zero vertices
        res = LaplacianIntegralResult(L, true, Int[], OrderedDict{Int,Int}())
        return _null_graph_laplacian_spectra(res)
    end

    if iszero(L) # `L` represents an empty graph (a graph with no edges)
        res = LaplacianIntegralResult(L, true, zeros(Int, n), OrderedDict(0 => n))
        return _empty_graph_laplacian_spectra(res)
    end

    res = laplacian_integrality(L)

    #= We are only interested in Laplacian integral graphs, a precondition for
    S-diagonalizability in the case of integer edge weights. =#
    if !res.laplacian_integral
        return LaplacianSpectra(
            res,
            Dict{Int,Matrix{Int}}(),
            Dict{Int,Vector{Int}}(),
            Dict{Int,Matrix{Int}}(),
            Dict{Int,Matrix{Int}}(),
        )
    end

    if _is_complete_graph_res(res) # `L` represents a complete graph on `n ≥ 2` vertices
        return _complete_graph_laplacian_spectra(res)
    end

    # `L` represents some Laplacian integral graph for which no special properties are known
    return _arbitrary_graph_laplacian_spectra(res)
end

function _extract_independent_cols(E::Matrix{Int})
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

function _find_idxs_1neg(eigspaces_01neg::Dict{Int,Matrix{Int}})
    #= The order of the Laplacian matrix. Accessing the zero key is safe, as every
    undirected graph has a rank-deficient Laplacian. =#
    n = size(eigspaces_01neg[0], 1)

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

function _null_graph_laplacian_spectra(res::LaplacianIntegralResult)
    eigspaces_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    idxs_1neg = Dict(0 => Int[])
    bases_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    bases_1neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _empty_graph_laplacian_spectra(res::LaplacianIntegralResult)
    n = size(res.laplacian_matrix, 1)
    eigspaces_01neg = Dict(0 => Matrix{Int}(undef, n, 0))
    idxs_1neg = Dict(0 => Int[])
    eigspaces_01neg = Dict(0 => hcat(_potential_kernel_eigvecs(n)...))
    idxs_1neg = _find_idxs_1neg(eigspaces_01neg)
    pot_basis_1neg = _extract_independent_cols(eigspaces_01neg[0])

    if size(pot_basis_1neg, 2) == n # The potential {-1,1}-basis spans the kernel
        bases_1neg = Dict(0 => pot_basis_1neg)
        bases_01neg = Dict(0 => copy(pot_basis_1neg)) # Avoid shared mutability
    else # No {-1,1}-basis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        #= Use an empty matrix as a sentinel for the {-1,1}-basis, and compute a new basis
        allowing for 0's for the {-1,0,1}-basis (guaranteed to exist). =#
        bases_1neg = Dict(0 => Matrix{Int}(undef, n, 0))
        bases_01neg = Dict(0 => _extract_independent_cols(eigspaces_01neg[0]))
    end

    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _complete_graph_laplacian_spectra(res::LaplacianIntegralResult)
    #= For every `n ≥ 2`, `Kₙ` (assuming uniform edge weights) has one zero eigenvalue with
    multiplicity 1 and one nonzero eigenvalue with multiplicity `n - 1`. (`n` is always
    guaranteed to be at least 2 in whatever context this helper function is called.) =#
    n = size(res.laplacian_matrix, 1)
    val1, val2 = keys(res.eigval_counts) # This is safe--there are precisely two eigenvalues
    eigval_nonzero = val1 == 0 ? val2 : val1

    kernel_01neg = ones(Int, n, 1) # The all-ones vector spans the (one-dimensional) kernel
    # The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel
    nonkernel_01neg = hcat(_potential_nonkernel_eigvecs(n)...)

    # All the {-1,0,1}-eigenvectors corresponding to each eigenvalue, as columns of a matrix
    eigspaces_01neg = Dict(0 => kernel_01neg, eigval_nonzero => nonkernel_01neg)
    idxs_1neg = _find_idxs_1neg(eigspaces_01neg) # The indices of the {-1,1}-columns

    #= Find {-1,0,1}- and, if they exist, {-1,1}-eigenbases. {-1,1}-bases are preferable, so
    use them wherever applicable. =#
    null_basis_1neg = copy(kernel_01neg) # Avoid shared mutability
    null_basis_01neg = copy(kernel_01neg) # Avoid shared mutability
    nonnull_idxs_1neg = idxs_1neg[eigval_nonzero]
    nonnull_pot_basis_1neg = _extract_independent_cols(
        nonkernel_01neg[:, nonnull_idxs_1neg]
    )

    if size(nonnull_pot_basis_1neg, 2) == n - 1 # The potential {-1,1}-basis spans the space
        nonnull_basis_1neg = nonnull_pot_basis_1neg
        nonnull_basis_01neg = copy(nonnull_pot_basis_1neg) # Avoid shared mutability
    else # No {-1,1}-eigenbasis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        #= Use an empty matrix as a sentinel for the {-1,1}-basis, and compute a new basis
        allowing for 0's for the {-1,0,1}-basis (guaranteed to exist). =#
        nonnull_basis_1neg = Matrix{Int}(undef, n, 0)
        nonnull_basis_01neg = _extract_independent_cols(nonkernel_01neg)
    end

    #= A {-1,0,1}- and, if it exists, a {-1,1}-basis for each eigenspace. If no {-1,1}-basis
    exists for the non-kernel eigenspace, an empty matrix is used as a placeholder. (The
    kernel always has one, as it is spanned by the all-ones vector.) =#
    bases_01neg = Dict(0 => null_basis_01neg, eigval_nonzero => nonnull_basis_01neg)
    bases_1neg = Dict(0 => null_basis_1neg, eigval_nonzero => nonnull_basis_1neg)

    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _arbitrary_graph_laplacian_spectra(res::LaplacianIntegralResult)
    L = res.laplacian_matrix
    n = size(L, 1)
    eigvals_nonzero = filter(!iszero, keys(res.eigval_counts))

    eigspaces_01neg_elastic = Dict(
        # Initialize resizeable matrices to store an unknown number of eigenvectors
        eigval => ElasticMatrix{Int}(undef, n, 0) for eigval in eigvals_nonzero
    )

    # The kernel is simpler and handled differently from the other eigenspaces
    if res.eigval_counts[0] == 1 # The all-ones vector spans the kernel
        eigspaces_01neg_elastic[0] = ElasticMatrix(ones(Int, n, 1))
    else # Check all {-1,0,1}-vectors in ℝⁿ, unique up to span
        eigspaces_01_neg_elastic[0] = ElasticMatrix(
            hcat(Iterators.filter(v -> iszero(L * v), _potential_kernel_eigvecs(n))...)
        )
    end

    # TODO: Multithread this; lazy evaluation means speed will become an issue before memory
    # Now fill up the remaining (non-kernel) eigenspaces
    for v in _potential_kernel_eigvecs(n)
        for eigval in eigvals_nonzero
            if L * v == eigval * v
                append!(eigspaces_01neg_elastic[eigval], v)
                break
            end
        end
    end

    # Convert the eigenspace data back to statically sized matrices for API compatibility
    eigspaces_01neg = Dict(
        eigval => Matrix(eigspace_elastic) for
        (eigval, eigspace_elastic) in eigspaces_01neg_elastic
    )
    idxs_1neg = _find_idxs_1neg(eigspaces_01neg) # The indices of the {-1,1}-columns

    # Find largest independent subsets of eigenvectors in each eigenspace
    bases_01neg = Dict(
        eigval => _extract_independent_cols(eigspace) for
        (eigval, eigspace) in eigspaces_01neg
    )
    bases_1neg = Dict(
        eigval => _extract_independent_cols(eigspace[:, idxs_1neg[eigval]]) for
        (eigval, eigspace) in eigspaces_01neg
    )

    # {-1,1}-bases are preferable, so use them in lieu of {-1,0,1}-bases wherever applicable
    for (eigval, basis_01neg) in bases_01neg
        basis_1neg = bases_1neg[eigval]

        if size(basis_01neg) == size(basis_1neg)
            bases_1neg[eigval] = copy(basis_01neg) # Avoid shared mutability
        end
    end

    return LaplacianSpectra(res, eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg)
end

function _is_complete_graph_res(res::LaplacianIntegralResult)
    n = size(res.laplacian_matrix, 1)
    return length(res.eigval_counts) == 2 && all(values(res.eigval_counts) .== [1, n - 1])
end

struct LaplacianIntegralResult
    spectrum_integral::Bool
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
end

function is_laplacian_integral(L::Symmetric{Int,Matrix{Int}})
    eigvals_float = eigvals(L)
    eigvals_int = round.(Int, eigvals_float) # Necessarily real by symmetricity of L
    laplacian_integral = isapprox(eigvals_float, eigvals_int)

    if !laplacian_integral
        eigvals_sorted = Int[]
        eigval_counts = OrderedDict{Int,Int}()
    else
        eigval_counts = OrderedDict(sort!(collect(counter(eigvals_int)); by=reverse))
        eigvals_sorted = collect(Iterators.flatmap(pair -> fill(pair...), eigval_counts))
    end

    return LaplacianIntegralResult(laplacian_integral, eigvals_sorted, eigval_counts)
end

@inline function _column_span_basis(E::Matrix{Int})
    rref, pivots = rref_with_pivots(E)
    mat_rank = rank(rref, 1e-5)

    # Iterate over potential bases until a basis is found (it is guaranteed that one exists)
    for comb in combinations(pivots, mat_rank)
        basis = E[:, comb]
        (rank(basis, 1e-5) == mat_rank) && return basis
    end
end

@inline function _get_idxs_1neg(eigspaces_01neg::Dict{Int,Matrix{Int}}, n::Int)
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

    return idxs_1neg
end

@inline function _null_graph_laplacian_spectra()
    eigspaces_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    idxs_1neg = Dict(0 => Int[])
    bases_01neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    bases_1neg = Dict(0 => Matrix{Int}(undef, 0, 0))
    return eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg
end

@inline function _empty_graph_laplacian_spectra(n::Int)
    eigspaces_01neg = Dict(0 => hcat(_potential_kernel_eigvecs(n)...))
    idxs_1neg = _get_idxs_1neg(eigspaces_01neg, n)
    pot_basis_1neg = _column_span_basis(eigspaces_01neg[0])

    if size(pot_basis_1neg, 2) == n # The potential {-1,1}-basis spans the kernel
        bases_1neg = Dict(0 => pot_basis_1neg)
        bases_01neg = Dict(0 => copy(pot_basis_1neg)) # Avoid shared mutability
    else # No {-1,1}-basis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        bases_1neg = Dict(0 => Matrix{Int}(undef, n, 0)) # Use an empty matrix as a sentinel
        # Compute a basis with 0's
        bases_01neg = Dict(0 => _column_span_basis(eigspaces_01neg[0]))
    end

    return eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg
end

function _complete_graph_laplacian_spectra(eigvals::Tuple{Int,Int}, n::Int)
    #= For every `n ≥ 2`, `Kₙ` (assuming uniform edge weights) has one zero eigenvalue with
    multiplicity 1 and one nonzero eigenvalue with multiplicity `n - 1`. (`n` is always
    guaranteed to be at least 2 in whatever context this helper function is called.) =#
    val1, val2 = eigvals
    eigval_nonzero = val1 == 0 ? val2 : val1

    kernel_01neg = ones(Int, n, 1) # The all-ones vector spans the (one-dimensional) kernel
    # The non-kernel eigenspace necessarily contains all vectors orthogonal to the kernel
    nonkernel_01neg = hcat(_potential_nonkernel_eigvecs(n)...)

    # All the {-1,0,1}-eigenvectors corresponding to each eigenvalue, as columns of a matrix
    eigspaces_01neg = Dict(0 => kernel_01neg, eigval_nonzero => nonkernel_01neg)
    idxs_1neg = _get_idxs_1neg(eigspaces_01neg, n) # The indices of the {-1,1}-columns

    #= Find {-1,0,1}- and, if they exist, {-1,1}-eigenbases. {-1,1}-bases are preferable, so
    use them wherever applicable. =#
    null_basis_1neg = copy(kernel_01neg) # Avoid shared mutability
    null_basis_01neg = copy(kernel_01neg) # Avoid shared mutability
    nonnull_idxs_1neg = idxs_1neg[eigval_nonzero]
    nonnull_pot_basis_1neg = _column_span_basis(nonkernel_01neg[:, nonnull_idxs_1neg])

    if size(nonnull_pot_basis_1neg, 2) == n - 1 # The potential {-1,1}-basis spans the space
        nonnull_basis_1neg = nonnull_pot_basis_1neg
        nonnull_basis_01neg = copy(nonnull_pot_basis_1neg) # Avoid shared mutability
    else # No {-1,1}-eigenbasis exists, so use a {-1,0,1}-eigenbasis (guaranteed to exist)
        nonnull_basis_1neg = Matrix{Int}(undef, n, 0) # Use an empty matrix as a sentinel
        nonnull_basis_01neg = _column_span_basis(nonkernel_01neg) # Compute a basis with 0's
    end

    #= A {-1,0,1}- and, if it exists, a {-1,1}-basis for each eigenspace. If no {-1,1}-basis
    exists for the non-kernel eigenspace, an empty matrix is used as a placeholder. (The
    kernel always has one, as it is spanned by the all-ones vector.) =#
    bases_01neg = Dict(0 => null_basis_01neg, eigval_nonzero => nonnull_basis_01neg)
    bases_1neg = Dict(0 => null_basis_1neg, eigval_nonzero => nonnull_basis_1neg)

    return eigspaces_01neg, idxs_1neg, bases_01neg, bases_1neg
end

struct _SpectrumIntegralResult
    spectrum_integral::Bool
    eigvals_sorted::Vector{Int}
    eigval_counts::OrderedDict{Int,Int}
end

struct LaplacianSpectra
    laplacian_matrix::Symmetric{Int,Matrix{Int}}
    eigvals_sorted::Vector{Int}
    eigval_counts::OrderedDict{Int,Int}
    bases_zerooneneg::OrderedDict{Int,Matrix{Int}}
    bases_oneneg::OrderedDict{Int,Matrix{Int}}
end

function _is_spectrum_integral(L::Symmetric{Int,Matrix{Int}})
    eigvals_float = eigvals(L)
    eigvals_int = round.(Int, eigvals_float) # Necessarily real by symmetricity of L
    spectrum_integral = isapprox(eigvals_float, eigvals_int)

    if !spectrum_integral
        eigvals_sorted = Int[]
        eigval_counts = OrderedDict{Int,Int}()
    else
        eigval_counts = OrderedDict(sort!(collect(counter(eigvals_int)); by=reverse))
        eigvals_sorted = collect(Iterators.flatmap(pair -> fill(pair...), eigval_counts))
    end

    return _SpectrumIntegralResult(spectrum_integral, eigvals_sorted, eigval_counts)
end

@inline function _column_span_basis(E::Matrix{Int})
    reduced_form = rref_with_pivots(E)
    reduced_rank = rank(reduced_form[1], 1e-5)

    # Iterate over potential bases until a basis is found (it is guaranteed that one exists)
    for comb in combinations(reduced_form[2], reduced_rank)
        basis = E[:, comb]
        (rank(basis, 1e-5) == reduced_rank) && return basis
    end
end

function _get_idxs_oneneg(eigspaces::OrderedDict{Int,Matrix{Int}}, idx_kernel::Int, n::Int)
    # Only even-ordered Laplacians have non-kernel {-1,1}-eigenvectors
    if n % 2 == 0 # Check all eigenspaces for {-1,1}-eigenvectors
        idxs_oneneg = OrderedDict(
            eigval => findall(v -> all(v .!= 0), eachcol(eigvecs)) for
            (eigval, eigvecs) in eigspaces
        )
    else # Only check the kernel for {-1,1}-eigenvectors
        idxs_oneneg = OrderedDict(eigval => Int[] for eigval in keys(eigspaces))
        idxs_oneneg[idx_kernel] = findall(v -> all(v .!= 0), eachcol(eigspaces[idx_kernel]))
    end

    return idxs_oneneg
end

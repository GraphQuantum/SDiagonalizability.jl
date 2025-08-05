# Copyright 2025 Luis M. B. Varona, Nathaniel Johnston, and Sarah Plosker
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    find_k_orthogonal_basis(
        col_space, col_rank, k, pre_basis=nothing
    ) -> Union{Nothing,AbstractMatrix{<:Integer}}

Find a `k`-orthogonal spanning subset of the column set of a matrix.

This function does not allow for any `k`-orthogonal basis of `col_space` (in which case it
would be trivial to construct a *pairwise* orthogonal basis via a QR decomposition).
Instead, it restricts its search to spanning subsets comprised exclusively of columns
already in `col_space`, returning `nothing` if no such basis exists.

# Arguments
- `col_space::T<:AbstractMatrix{<:Integer}`: the matrix whose column space is searched for a
    `k`-orthogonal basis.
- `col_rank::Integer`: the rank of the column space of `col_space`, pre-computed for
    efficiency.
- `k::Integer`: the minimum ``k``-orthogonality parameter of the desired basis. Must be a
    positive integer.
- `pre_basis::Union{Nothing,AbstractMatrix{<:Integer}}=nothing`: an optional precomputed
    submatrix of `col_space` whose columns form a spanning subset of the column space of
    `col_space`. In case `col_space` has at most `k` columns, `existing_basis` is guaranteed
    to already be a `k`-orthogonal basis and thus is used to skip unnecessary computations.

# Returns
- `::Union{Nothing,AbstractMatrix{<:Integer}}`: a `k`-orthogonal spanning subset of the
    columns of `col_space`, if one exists; otherwise, `nothing`.
"""
function find_k_orthogonal_basis(
    col_space::AbstractMatrix{<:Integer},
    col_rank::Integer,
    k::Integer,
    pre_basis::Union{Nothing,AbstractMatrix{<:Integer}}=nothing,
)
    if isnothing(pre_basis) && size(col_space, 2) <= k
        # Every collection of `k` or less vectors is `k`-orthogonal
        basis = _extract_independent_cols(col_space)
    #! format: off
    #= `MatrixBandwidth.jl` uses zero-based indexing for bandwidth, not one-based, so we use
    `k - 1` instead of `k`. =#
    #! format: on
    elseif !isnothing(pre_basis) && bandwidth(pre_basis'pre_basis) <= k - 1
        basis = pre_basis
    else
        prop = classify_k_orthogonality(k)
        basis = _find_basis_with_property(col_space, col_rank, prop)
    end

    return basis
end

"""
    _find_basis_with_property(col_space::T, col_rank, prop) -> Union{Nothing,T}

[TODO: Write here]
"""
function _find_basis_with_property(
    ::AbstractMatrix{<:Integer}, ::Integer, ::K
) where {K<:KOrthogonality}
    throw(NotImplementedError(_find_basis_with_property, K, KOrthogonality))
end

function _find_basis_with_property(
    col_space::AbstractMatrix{<:Integer}, col_rank::Integer, prop::Orthogonality
)
    for root in 1:(size(col_space, 2) - col_rank + 1)
        basis_idxs = _find_basis_idxs_with_prop([root], prop, col_space, col_rank, 1)

        if !(isnothing(basis_idxs))
            return col_space[:, basis_idxs]
        end
    end

    return nothing
end

function _find_basis_with_property(
    col_space::AbstractMatrix{<:Integer}, col_rank::Integer, prop::QuasiOrthogonality
)
    #= We can restrict each `root` to end no later than `size(col_space, 2) - col_rank + 2`,
    since we must fill `col_rank - 2` positions afterwards. =#
    for root in combinations(1:(size(col_space, 2) - col_rank + 2), 2)
        #= The current vertex degree and connected component of `root[1]` will always be 0
        and 1, respectively, before the DFS begins. =#
        nodes = [_QOBasisSearchNode(1, 0)]
        union_find = Dict{UInt16,Vector{UInt16}}(1 => [1])

        basis_idxs = _find_basis_idxs_with_prop(
            root, prop, col_space, col_rank, 2, nodes, union_find
        )

        if !(isnothing(basis_idxs))
            return col_space[:, basis_idxs]
        end
    end

    return nothing
end
#lalalala
function _find_basis_with_property(
    col_space::AbstractMatrix{<:Integer}, col_rank::Integer, prop::WeakOrthogonality
)
    k = prop.k
    num_columns = size(col_space, 2)

    gram_matrix = falses(col_rank, col_rank)

    #= Every collection of `k` or less vectors is `k`-orthogonal, so we could start
    searching from `k`-combinations of the columns. However, we start from 2-combinations to
    allow for more aggressive pruning of linearly dependent subsets. We can restrict each
    `root` to end no later than `num_columns - col_rank + 2`, since we must fill
    `col_rank - 2` positions afterwards. =#
    for root in combinations(1:(num_columns - col_rank + 2), 2)
        partial_basis = view(col_space, :, root)
        gram_submatrix = (!iszero).(partial_basis'partial_basis)
        gram_matrix[1:2, 1:2] .= gram_submatrix

        basis_idxs = _find_basis_idxs_with_prop(
            root, prop, col_space, col_rank, 2, gram_matrix
        )

        if !(isnothing(basis_idxs))
            return col_space[:, basis_idxs]
        end
    end

    return nothing
end

"""
    _find_basis_idxs_with_prop(
        curr_idxs, prop::Orthogonality, col_space, col_rank, depth
    ) -> Union{Nothing,AbstractVector{Int}}
    _find_basis_idxs_with_prop(
        curr_idxs, prop::QuasiOrthogonality, col_space, col_rank, depth, nodes, union_find
    ) -> Union{Nothing,AbstractVector{Int}}
    _find_basis_idxs_with_prop(
        curr_idxs, prop::WeakOrthogonality, col_space, col_rank, depth, gram_matrix
    ) -> Union{Nothing,AbstractVector{Int}}

[TODO: Write here]
"""
function _find_basis_idxs_with_prop(
    ::AbstractVector{Int}, ::K, ::AbstractMatrix{Int}, ::Integer, ::Int, args...
) where {K<:KOrthogonality}
    throw(NotImplementedError(_find_basis_idxs_with_prop, :prop, K, KOrthogonality))
end

function _find_basis_idxs_with_prop(
    curr_idxs::AbstractVector{Int},
    prop::Orthogonality,
    col_space::AbstractMatrix{Int},
    col_rank::Integer,
    depth::Int,
)
    num_columns = size(col_space, 2)
    parent_partial_basis = view(col_space, :, 1:(depth - 1))
    curr_column = view(col_space, :, depth)

    if !(iszero(parent_partial_basis'curr_column))
        return nothing
    end

    remaining_levels = col_rank - depth

    if remaining_levels == 0
        return curr_idxs
    end

    for i in (curr_idxs[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_idxs, i)
        basis_idxs = _find_basis_idxs_with_prop(
            curr_idxs, prop, col_space, col_rank, depth + 1
        )
        pop!(curr_idxs)

        if !isnothing(basis_idxs)
            return basis_idxs
        end
    end

    return nothing
end

"""
    _QOBasisSearchNode

[TODO: Write here]
"""
mutable struct _QOBasisSearchNode
    component::UInt16
    degree::UInt16
end

function _find_basis_idxs_with_prop(
    curr_idxs::AbstractVector{Int},
    prop::QuasiOrthogonality,
    col_space::AbstractMatrix{Int},
    col_rank::Integer,
    depth::Int,
    nodes::Vector{_QOBasisSearchNode},
    union_find::Dict{UInt16,Vector{UInt16}},
)
    num_columns = size(col_space, 2)
    partial_basis = view(col_space, :, curr_idxs)
    # Use a more robust tolerance (NumPy's and MATLAB's default) than Julia's default
    rtol = _rank_rtol(partial_basis)

    if rank(partial_basis; rtol=rtol) < depth
        return nothing
    end

    parent_partial_basis = view(partial_basis, :, 1:(depth - 1))
    curr_column = view(partial_basis, :, depth)
    neighbors = findall(!iszero, parent_partial_basis'curr_column)
    degree = length(neighbors)

    if degree == 0
        component = maximum(keys(union_find)) + 1
        union_find[component] = UInt16[]
    elseif degree == 1
        neighbor = nodes[neighbors[1]]

        if neighbor.degree == 2
            return nothing
        end

        neighbor.degree += 1
        component = neighbor.component
    elseif degree == 2
        neighbor1, neighbor2 = nodes[neighbors]

        if neighbor1.component == neighbor2.component
            return nothing
        end

        if neighbor1.degree == 2 || neighbor2.degree == 2
            return nothing
        end

        neighbor1.degree += 1
        neighbor2.degree += 1
        component = min(neighbor1.component, neighbor2.component)
        deleted_component = max(neighbor1.component, neighbor2.component)

        for i in union_find[deleted_component]
            nodes[i].component = component
            push!(union_find[component], i)
        end

        delete!(union_find, deleted_component)
    else
        return nothing
    end

    node = _QOBasisSearchNode(component, degree)
    push!(nodes, node)
    push!(union_find[component], depth)

    remaining_levels = col_rank - depth

    if remaining_levels == 0
        return curr_idxs[vcat(values(union_find)...)]
    end

    for i in (curr_idxs[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_idxs, i)
        nodes_new = deepcopy(nodes)
        union_find_new = deepcopy(union_find)

        basis_idxs = _find_basis_idxs_with_prop(
            curr_idxs, prop, col_space, col_rank, depth + 1, nodes_new, union_find_new
        )

        pop!(curr_idxs)

        if !(isnothing(basis_idxs))
            return basis_idxs
        end
    end

    return nothing
end

function _find_basis_idxs_with_prop(
    curr_idxs::AbstractVector{Int},
    prop::WeakOrthogonality,
    col_space::AbstractMatrix{Int},
    col_rank::Integer,
    depth::Int,
    gram_matrix::BitMatrix,
)
    k = prop.k
    num_columns = size(col_space, 2)
    partial_basis = view(col_space, :, curr_idxs)
    # Use a more robust tolerance (NumPy's and MATLAB's default) than Julia's default
    rtol = _rank_rtol(partial_basis)

    if rank(partial_basis; rtol=rtol) < depth
        return nothing
    end

    parent_partial_basis = view(partial_basis, :, 1:(depth - 1))
    curr_column = view(partial_basis, :, depth)
    neighbors = findall(!iszero, parent_partial_basis'curr_column)

    if length(neighbors) > 2k - 2
        return nothing
    end

    gram_matrix[:, depth] .= gram_matrix[depth, :] .= false # Reset due to shared mutability
    gram_matrix[depth, neighbors] .= gram_matrix[neighbors, depth] .= true

    gram_submatrix = view(gram_matrix, 1:depth, 1:depth)
    #= `MatrixBandwidth.jl` uses zero-based indexing for bandwidth, not one-based, so we use
    `k - 1` instead of `k`. =#
    res = has_bandwidth_k_ordering(gram_submatrix, k - 1, Recognition.DelCorsoManzini())

    if !res.has_ordering
        return nothing
    end

    remaining_levels = col_rank - depth

    if remaining_levels == 0
        return curr_idxs[res.ordering]
    end

    for i in (curr_idxs[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_idxs, i)
        basis_idxs = _find_basis_idxs_with_prop(
            curr_idxs, prop, col_space, col_rank, depth + 1, gram_matrix
        )
        pop!(curr_idxs)

        if !(isnothing(basis_idxs))
            return basis_idxs
        end
    end

    return nothing
end

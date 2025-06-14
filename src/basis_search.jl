# Copyright 2025 Luis M. B. Varona and Nathaniel Johnston
#
# Licensed under the MIT license <LICENSE or
# http://opensource.org/licenses/MIT>. This file may not be copied, modified, or
# distributed except according to those terms.

"""
    mutable struct _QOBasisSearchNodeData

TODO: Write here
"""
mutable struct _QOBasisSearchNodeData
    component::UInt8
    vertex_degree::UInt8
end

"""
    _find_k_orthogonal_basis(column_space, column_rank, k)

Find a `k`-orthogonal spanning subset of the column set of a matrix.

This function does not allow for any `k`-orthogonal basis of `column_space` (in which case
it would be trivial to construct a *pairwise* orthogonal basis via a QR decomposition).
Instead, it restricts its search to spanning subsets comprised exclusively of columns
already in `column_space`, returning `nothing` if no such basis exists.

# Arguments
- `column_space::T<:AbstractMatrix{Int}`: the matrix whose column space is search for a
    `k`-orthogonal basis.
- `column_rank::Int`: the rank of the column space of `column_space`, pre-computed for
    efficiency.
- `k::Integer`: the minimum `k`-orthogonality parameter of the desired basis. (Must be a
    positive integer.)

# Returns
- `::Union{Nothing,T}`: a `k`-orthogonal spanning subset of the columns of `column_space`,
    with the same type as `column_space`. If no such basis exists, this is simply `nothing`.

# Notes
TODO: Write here. In particular, comment on the "basis" naming and perhaps `k`. Also,
comment on typing decisions. (Should `k` just be `Int`?)
"""
function _find_k_orthogonal_basis(
    column_space::AbstractMatrix{Int}, column_rank::Int, k::Integer
)
    prop = _classify_orthogonality_property(k)
    return _find_basis_with_property(column_space, column_rank, prop)
end

"""
    _find_basis_with_property(column_space, column_rank, prop::_KOrthogonality)

TODO: Write here
"""
function _find_basis_with_property(::AbstractMatrix{Int}, ::Int, prop::_KOrthogonality)
    throw(NotImplementedError(_find_basis_with_property, typeof(prop), _KOrthogonality))
end

"""
    _find_basis_with_property(column_space, column_rank, prop:::_Orthogonality)

TODO: Write here, and add comments
"""
function _find_basis_with_property(
    column_space::AbstractMatrix{Int}, column_rank::Int, prop::_Orthogonality
)
    num_columns = size(column_space, 2)

    if num_columns == 1
        return column_space[:, 1]
    end

    for root_index in 1:(num_columns - column_rank + 1)
        basis_indices = _find_basis_indices([root_index], column_space, column_rank, prop)

        if !(isnothing(basis_indices))
            return column_space[:, basis_indices]
        end
    end

    return nothing
end

"""
    _find_basis_with_property(column_space, column_rank, prop::_QuasiOrthogonality)

TODO: Write here, and add comments
"""
function _find_basis_with_property(
    column_space::AbstractMatrix{Int}, column_rank::Int, prop::_QuasiOrthogonality
)
    num_columns = size(column_space, 2)

    if num_columns <= 2
        return column_space[:, 1:column_rank]
    end

    #= The current vertex degree and connected component of any root node (each
    `root_indices[1]` in the loop below) will always be 0 and 1 before the DFS begins. =#
    nodes = [_QOBasisSearchNodeData(1, 0)]
    union_find = Dict{UInt8,Vector{UInt8}}(1 => [1])

    for root_indices in combinations(1:num_columns, 2)
        basis_indices = _find_basis_indices(
            root_indices, column_space, column_rank, nodes, union_find, prop
        )

        if !(isnothing(basis_indices))
            return column_space[:, basis_indices]
        end
    end

    return nothing
end

# TODO: Note that we assume `k < column_rank` (else, just use QR results). Add assert or no?
"""
    _find_basis_with_property(column_space, column_rank, prop::_WeakOrthogonality)

TODO: Write here
"""
function _find_basis_with_property(
    column_space::AbstractMatrix{Int}, column_rank::Int, prop::_WeakOrthogonality
)
    k = prop.k
    num_columns = size(column_space, 2)

    if num_columns <= k
        return column_space[:, 1:column_rank]
    end

    #= Construct the orthogonality graph complement of the "least `k`-orthogonal" matrix in
    `Mₘ,ₙ` to apply a subgraph monomorphism search. (We later plan to replace this with a
    more efficient matrix bandwidth minimization algorithm.) =#
    G = Graph(column_rank) # First, initialize the empty graph on `column_rank` vertices

    #= Iterate over all possible distances between non-orthogonal columns in the matrix
    whose orthogonality graph complement `G` represents. =#
    for d in 1:(k - 1), i in 1:(column_rank - d)
        add_edge!(G, i, i + d)
    end

    #= Every collection of `k` or less vectors is `k`-orthogonal, so we can start searching
    from `k`-combinations of our columns. =#
    for root_indices in combinations(1:num_columns, k)
        # Avoid reallocations by taking a view
        partial_basis = view(column_space, :, root_indices)

        #= Bandwidth minimization of the Gram matrix of our candidate basis vectors is
        equivalent to testing the resulting graph for subgraph monomorphism to `G`. (We do,
        however, plan to replace this with a more specialized algorithm in the future.) =#
        ortho_graph_compl_adjacency = .!iszero.(partial_basis' * partial_basis)
        ortho_graph_compl_adjacency[diagind(ortho_graph_compl_adjacency)] .= false

        #= Before running a subgraph monomorphism search, we can filter out (partial) bases
        that lack a sufficient number of orthogonality relations to be `k`-orthogonal,
        regardless of the column ordering. =#
        if all(sum(ortho_graph_compl_adjacency; dims=1) .<= 2k - 2)
            H = Graph(column_rank) # Initialize the empty graph on `column_rank` vertices

            #= The aforementioned "resulting graph" from the Gram matrix is constructed by
            [TODO: Write here] =#
            for i in 1:(k - 1), j in (i + 1):k
                if ortho_graph_compl_adjacency[i, j]
                    add_edge!(H, i, j)
                end
            end

            basis_indices = _find_basis_indices(
                root_indices, column_space, column_rank, G, H, prop
            )

            # If there is no monomorphism, terminate searching this branch of the tree
            if !(isnothing(basis_indices))
                return column_space[:, basis_indices]
            end
        end
    end

    #= If all children are searched/pruned and basis is found, return `nothing` to signal
    that this branch should be pruned. =#
    return nothing
end

"""
    _find_basis_indices(curr_indices, column_space, column_rank, prop::_Orthogonality)

TODO: Write here
"""
function _find_basis_indices(
    curr_indices::AbstractVector{Int},
    column_space::AbstractMatrix{Int},
    column_rank::Int,
    prop::_Orthogonality,
)
    depth = length(curr_indices)
    num_columns = size(column_space, 2)
    parent_partial_basis = view(column_space, :, 1:(depth - 1))
    curr_column = view(column_space, :, depth)

    if !(iszero(parent_partial_basis' * curr_column))
        return nothing
    end

    if depth == column_rank
        return copy(curr_indices)
    end

    remaining_levels = column_rank - depth

    for i in (curr_indices[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_indices, i)
        basis_indices = _find_basis_indices(curr_indices, column_space, column_rank, prop)
        pop!(curr_indices)

        if !(isnothing(basis_indices))
            return basis_indices
        end
    end

    return nothing
end

"""
    _find_basis_indices(
        curr_indices,
        column_space,
        column_rank,
        nodes,
        union_find,
        prop::_QuasiOrthogonality,
    )

TODO: Write here
"""
function _find_basis_indices(
    curr_indices::AbstractVector{Int},
    column_space::AbstractMatrix{Int},
    column_rank::Int,
    nodes::Vector{_QOBasisSearchNodeData},
    union_find::Dict{UInt8,Vector{UInt8}},
    prop::_QuasiOrthogonality,
)
    depth = length(curr_indices)
    num_columns = size(column_space, 2)
    partial_basis = view(column_space, :, curr_indices)
    # Use a more robust tolerance (NumPy's and MATLAB's default) than Julia's default
    rtol = _rank_rtol(partial_basis)

    if rank(partial_basis; rtol=rtol) < depth
        return nothing
    end

    parent_partial_basis = view(partial_basis, :, 1:(depth - 1))
    curr_column = view(partial_basis, :, depth)
    neighbor_indices = findall(!iszero, parent_partial_basis' * curr_column)
    degree = length(neighbor_indices)

    if degree == 0
        component = maximum(keys(union_find)) + 1
        union_find[component] = UInt8[]
    elseif degree == 1
        neighbor = nodes[neighbor_indices[1]]

        if neighbor.degree == 2
            return nothing
        end

        neighbor.degree += 1
        component = neighbor.component
    elseif degree == 2
        neighbor1, neighbor2 = nodes[neighbor_indices]

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

    node = _QOBasisSearchNodeData(component, degree)
    push!(nodes, node)
    push!(union_find[component], depth)

    if depth == column_rank
        basis_indices = copy(curr_indices)
        #= Sorting within each connected component of the union-find structure
        yields a "nicer," more natural ordering. =#
        order = collect(Iterators.flatmap(sort, values(union_find)))
        return basis_indices[order]
    end

    remaining_levels = column_rank - depth

    for i in (curr_indices[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_indices, i)
        nodes = deepcopy(nodes)
        union_find = deepcopy(union_find)

        basis_indices = _find_basis_indices(
            curr_indices, column_space, column_rank, nodes, union_find, prop
        )

        pop!(curr_indices)

        if !(isnothing(basis_indices))
            return basis_indices
        end
    end

    return nothing
end

"""
    _find_basis_indices(
        curr_indices,
        column_space,
        column_rank,
        G::Graph,
        H::Graph,
        prop::_WeakOrthogonality,
    )

TODO: Write here
"""
function _find_basis_indices(
    curr_indices::AbstractVector{Int},
    column_space::AbstractMatrix{Int},
    column_rank::Int,
    G::Graph,
    H::Graph,
    prop::_WeakOrthogonality,
)
    k = prop.k
    depth = length(curr_indices)
    num_columns = size(column_space, 2)
    partial_basis = view(column_space, :, curr_indices)
    # Use a more robust tolerance (NumPy's and MATLAB's default) than Julia's default
    rtol = _rank_rtol(partial_basis)

    if rank(partial_basis; rtol=rtol) < depth
        return nothing
    end

    parent_partial_basis = view(partial_basis, :, 1:(depth - 1))
    curr_column = view(partial_basis, :, depth)
    neighbor_indices = findall(!iszero, parent_partial_basis' * curr_column)

    if length(neighbor_indices) > 2k - 2
        return nothing
    end

    for u in neighbor_indices
        add_edge!(H, depth, u)
    end

    monomorphism = iterate(subgraph_monomorphisms(G, H))

    if isnothing(monomorphism)
        return nothing
    end

    if depth == column_rank
        basis_indices = copy(curr_indices)
        # The monomorphism mapping gives a `k`-orthogonal column ordering
        order = last.(sort(collect(monomorphism[1])))
        return basis_indices[order]
    end

    remaining_levels = column_rank - depth

    for i in (curr_indices[depth] + 1):(num_columns - remaining_levels + 1)
        push!(curr_indices, i)
        basis_indices = _find_basis_indices(
            curr_indices, column_space, column_rank, G, H, prop
        )
        pop!(curr_indices)

        if !(isnothing(basis_indices))
            return basis_indices
        end
    end

    return nothing
end

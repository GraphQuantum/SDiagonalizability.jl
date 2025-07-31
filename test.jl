using SDiagonalizability, MatrixBandwidth, GraphIO.Graph6, Graphs

g6 = "J}ox~~}~~~_"
L = laplacian_matrix(Graph6._g6StringToGraph(g6))
res = SDiagonalizability.laplacian_s_spectra(L, (-1, 0, 1))

(n, k) = (11, 2)
out = SDiagonalizability.find_k_orthogonal_basis(hcat(values(res.s_eigenspaces)...), n, k) # `nothing`

C = hcat(values(res.s_eigenspaces)...)
E = [
    1 1 0 0 1 1 1 1 1 1 0;
    1 1 0 0 -1 -1 -1 1 1 1 0;
    1 0 1 1 0 0 0 0 -1 -1 0;
    1 0 -1 0 0 0 0 0 -1 -1 0;
    1 0 0 -1 0 0 0 0 -1 -1 0;
    1 -1 0 0 -1 0 1 1 1 1 0;
    1 -1 0 0 1 0 -1 1 1 1 0;
    1 0 0 0 0 0 -1 -1 0 0 1;
    1 0 0 0 0 0 1 -1 0 0 1;
    1 0 0 0 0 0 0 -1 0 -1 -1;
    1 0 0 0 0 0 0 -1 -1 0 -1
]
for col in eachcol(E)
    @assert col in eachcol(C)
end
bandwidth(E'E) <= k - 1

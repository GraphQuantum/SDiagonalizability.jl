# using .SDiagonalizability: _find_k_orthogonal_basis

# column_space = [1 0 0 0 3; 0 1 0 0 4; 0 0 1 0 15; 0 0 0 1 0; 0 0 0 0 0]

# _find_k_orthogonal_basis(column_space, 4, 1)
# _find_k_orthogonal_basis(column_space, 4, 2)
# _find_k_orthogonal_basis(column_space, 4, 3)
# _find_k_orthogonal_basis(column_space, 4, 4)
# _find_k_orthogonal_basis(column_space, 4, 5)

# using LinearAlgebra, Random

# function slow_decay_matrix(m, n, k; σ_max=100.0, σ_min=1e-15,
#                            noise=1e-12, decay=:geometric, tail_scale=1e-3,
#                            rng=MersenneTwister(87))
#     U, _ = qr!(randn(rng, m, m))
#     V, _ = qr!(randn(rng, n, n))

#     σ = Vector{Float64}(undef, m)
#     σ[1:k] .= 1 .+ (σ_max - 1) .* rand(rng, k)

#     tail_len = m - k
#     if tail_len > 0
#         tail_σ_max = minimum(σ[1:k]) * tail_scale
#         tail_σ_min = σ_min

#         if decay == :geometric
#             σ[k+1:end] = exp.(range(log(tail_σ_max), log(tail_σ_min), length=tail_len+1))[2:end]
#         elseif decay == :linear
#             σ[k+1:end] = range(tail_σ_max, tail_σ_min, length=tail_len+1)[2:end]
#         else
#             error("Unsupported decay type: $decay")
#         end
#     end

#     Σ = Diagonal(σ)
#     A = U * Σ * V[:,1:m]' + noise * randn(rng, m, n)
#     A = A[:, shuffle(rng, 1:n)]
#     return A
# end

# # Construct test matrix
# m, n, k = 12, 1000, 6
# A = slow_decay_matrix(m, n, k)

# # Compute relative tolerance based on largest singular value
# svd_vals = svdvals(A)
# σ₁ = maximum(svd_vals)
# rtol = maximum(size(A)) * eps(real(eltype(A))) * σ₁

# # SVD-based numerical rank
# svd_rank = count(>(rtol), svd_vals)

# # QRCP-based numerical rank
# F = qr(A, ColumnNorm())
# qr_rank = count(x -> abs(x) > rtol, diag(F.R))

# println("Consistent rtol: ", rtol)
# println("SVD-based rank: ", svd_rank)
# println("QRCP-based rank: ", qr_rank)

using LinearAlgebra, Random

# Generate k sorted points roughly uniform on log scale with jitter and sorting
function sorted_log_spaced(rng, k, log_lo, log_hi)
    if k == 1
        return [10 ^ (rand(rng) * (log_hi - log_lo) + log_lo)]
    else
        base = range(0, 1; length=k)
        jitter = 0.2 * randn(rng, k)
        pts = clamp.(base .+ jitter, 0, 1)
        sort!(pts)
        return 10 .^ (pts .* (log_hi - log_lo) .+ log_lo)
    end
end

function natural_decay_matrix(
    m, n, k; σ_max=100.0, σ_min=1e-15, noise=1e-12, tail_scale=1e-3, rng=MersenneTwister(87)
)
    U, _ = qr!(randn(rng, m, m))
    V, _ = qr!(randn(rng, n, n))

    σ = Vector{Float64}(undef, m)

    k_big = fld(k, 2)
    k_border = k - k_big
    tail_len = m - k

    # 1) Largest singular value close but not equal to σ_max:
    σ[1] = rand(rng)*(σ_max - 0.8*σ_max) + 0.8*σ_max  # Uniform in [0.9*σ_max, σ_max]

    # 2) Remaining big block (2..k_big): spread roughly 1–2 orders of magnitude below σ_max,
    #    with a higher lower bound so values don't get too small.
    if k_big > 1
        log_hi1 = log10(σ_max / 5)
        log_lo1 = log10(σ_max / 20)  # Increased lower bound here (was /500 before)
        σ[2:k_big] .= reverse(sorted_log_spaced(rng, k_big - 1, log_lo1, log_hi1))
    end

    # 3) Borderline block: 2–3 orders of magnitude, fixed range
    if k_border > 0
        border_hi = 1e-3
        border_lo = 1e-5
        log_lo2 = log10(border_lo)
        log_hi2 = log10(border_hi)
        σ[(k_big + 1):k] .= sorted_log_spaced(rng, k_border, log_lo2, log_hi2)
    end

    # 4) Tail block: very small singular values spanning several orders of magnitude
    if tail_len > 0
        tail_hi = 1e-6 * tail_scale
        tail_lo = σ_min
        log_lo3 = log10(tail_lo)
        log_hi3 = log10(tail_hi)
        σ[(k + 1):end] .= sorted_log_spaced(rng, tail_len, log_lo3, log_hi3)
    end

    # Add jitter ±50% to decorrelate singular values a bit
    σ .= σ .* (1 .+ 0.50 .* randn(rng, m))

    Σ = Diagonal(σ)
    A = U * Σ * V[:, 1:m]' + noise * randn(rng, m, n)
    A = A[:, shuffle(rng, 1:n)]
    return A, σ
end

# === Example usage ===

m, n, k = 12, 1000, 6
A, σ_true = natural_decay_matrix(m, n, k)

svd_vals = svdvals(A)
σ₁ = maximum(svd_vals)
rtol = maximum(size(A)) * eps(real(eltype(A))) * σ₁

svd_rank = count(>(rtol), svd_vals)
F = qr(A, ColumnNorm())
qr_rank = count(x -> abs(x) > rtol, diag(F.R))

println("True σ (before noise):")
println(round.(σ_true; sigdigits=5))

println("\nSVDvals of A:")
println(round.(svd_vals; sigdigits=5))

println("\nConsistent rtol: ", rtol)
println("SVD-based rank: ", svd_rank)
println("QRCP-based rank: ", qr_rank)

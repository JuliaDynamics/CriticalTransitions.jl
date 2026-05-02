# =====================================================================
# TPT-specific numerical helpers. Kept separate from the general
# diffusion-operator utilities since these expressions are specific to
# the reactive-trajectory decomposition.
# =====================================================================

"""
$(TYPEDSIGNATURES)

A → B reactive transition rate from a CTMC generator `G` (positive
off-diagonals = transition rates):

    k_{AB} = ∑_{i ∈ A, j ∉ A} ρ[i] · v[i] · q⁻[i] · G[i, j] · q⁺[j] .

Computed in `O(nnz(G))` by walking sparse columns once.
"""
function _reactive_rate(
        G::SparseMatrixCSC{Float64, Int},
        ρ::Vector{Float64},
        qplus::Vector{Float64},
        qminus::Vector{Float64},
        weights::Vector{Float64},
        A_mask::BitVector,
    )::Float64
    rv = rowvals(G)
    nz = nonzeros(G)
    rate = 0.0
    @inbounds for col in 1:size(G, 2)
        A_mask[col] && continue
        qpc = qplus[col]
        for p in nzrange(G, col)
            row = rv[p]
            row != col || continue
            A_mask[row] || continue
            rate += ρ[row] * weights[row] * qminus[row] * nz[p] * qpc
        end
    end
    return rate
end

# using LinearAlgebra

export vec_from_density_op
export thermal_state

function thermal_state(H, T)
    # kb = 1.380649 * 10^(-23)
    B = 1.0 / T
    out = exp(-B * H)
    return out / tr(out)
end

function vec_from_density_op(rho::Matrix{<:Number})
    u, s, v = svd(rho)
    result = findfirst(cumsum(s) .> rand())
    return u[:, result]
end
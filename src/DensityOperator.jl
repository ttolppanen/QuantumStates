# using LinearAlgebra

export vec_from_density_op
export thermal_state

# thermal states won't work properly if the temperature is too high
# when comparing the local system dimenstion d, e.g. if you have non-zero
# probability to occupy states with more than one excitation, then
# the qubit state won't be accurately sampled.

function thermal_state(H::Matrix{<:Number}, T)
    # kb = 1.380649 * 10^(-23)
    B = 1.0 / T
    out = exp(-B * H)
    return out / tr(out)
end
function thermal_state(H::AbstractVector, T)
    out = exp.(-H / T)
    return out / sum(out)
end

function vec_from_density_op(rho::Matrix{<:Number})
    u, s, v = svd(rho)
    result = findfirst(cumsum(s) .> rand())
    return u[:, result]
end
# for diagonal density operator
function vec_from_density_op(rho::AbstractVector)
    result = findfirst(cumsum(real.(rho)) .> rand())
    out = complex(spzeros(size(rho)))
    out[result] = 1.0
    return out
end
# using QuantumOperators
# using LinearAlgebra
# using Statistics

⊗(a, b) = kron(a, b)

function nop(d::Integer)
    out = spzeros(d, d)
    for i in 1:d
        out[i, i] = i - 1
    end
    return complex(out)
end

function singlesite(op::AbstractMatrix{<:Number}, L::Integer, target::Integer)
    d = size(op)[1]
    return sparse(I, d^(target-1), d^(target-1)) ⊗ op ⊗ sparse(I, d^(L-target), d^(L-target))
end

function nall(d::Integer, L::Integer)
    n = nop(d)
    out = singlesite(n, L, 1)
    for i in 2:L
        out += singlesite(n, L, i)
    end
    return out
end

function bosehubbard(d::Integer, L::Integer; w=1, U=1, J=1) # Order might be reversed here? The last site could be the first?
    if isa(w, Number) w = [w for _ in 1:L] end
    if isa(U, Number) U = [U for _ in 1:L] end
    if isa(J, Number) J = [J for _ in 1:L] end
    w = reverse(w)
    U = reverse(U)
    J = reverse(J)

    Is::Vector{Int64} = []
    Js::Vector{Int64} = []
    Vs::Vector{ComplexF64} = []
    n_s = zeros(L)
    for i in 1:d^L
        for site in 1:L
            n_s[site] = floor((i-1) / d^(site-1)) % d
        end
        for site in (L-1):-1:1
            n_next = n_s[site + 1]
            n = n_s[site]
            if n + 1 < d && n_next > 0 
                i_j = i + d^(site - 1) - d^site
                push!(Is, i)
                push!(Js, i_j)
                push!(Vs, J[site] * sqrt((n + 1) * n_next))
            end
        end
        push!(Is, i)
        push!(Js, i)
        H_i_i = 0
        for site in 1:L
            n = n_s[site]
            H_i_i += w[site] * n - U[site] / 2 * n * (n - 1)
        end
        push!(Vs, H_i_i)
        for site in 1:(L-1)
            n = n_s[site]
            n_next = n_s[site + 1]
            if n > 0 && n_next + 1 < d
                i_j = i - d^(site - 1) + d^site
                push!(Is, i)
                push!(Js, i_j)
                push!(Vs, J[site] * sqrt(n * (n_next + 1)))
            end
        end
    end
    return sparse(Is, Js, Vs)
end

function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number})
    out = dot(state, op, state)
    return real(out)
end
function expval(state::AbstractVector{<:Number}, op::AbstractMatrix{<:Number}, subspace_indeces::AbstractVector{<:Integer})
    @views out = dot(state[subspace_indeces], op[subspace_indeces, subspace_indeces], state[subspace_indeces])
    return real(out)
end
function expval(state::MPS, op::Union{Matrix{<:Number}, String}; kwargs...) # keyword arguments for ITensors.expect
    expect(state, op; kwargs...) # kwargs can be {sites} 
end

@testset "Thermal State" begin

@testset "Bose-Einstein Distribution" begin
    d = 3; L = 2
    w = 7.5 * 10^9 * 2 * pi
    J = 5 * 10^6 * 2 * pi
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34
    T = 0.1
    T_scale = kb / (hbar * J)

    H = bosehubbard(d, L; w = w / J, U = 0, J = 0)
    rho_T = thermal_state(diag(H), T * T_scale)

    N = 1000000
    n_op = nall(d, L)
    out = []
    for _ in 1:N
        v = vec_from_density_op(rho_T)
        push!(out, expval(v, n_op))
    end
    n_exp = mean(out)

    # error vs. bose-einstein distribution
    be_dist(e, T) = 1 / (exp(e / (kb * T)) - 1)
    dist_error = abs(n_exp - 2 * be_dist(w * hbar, T))
    @test dist_error / n_exp < 0.01

    # error vs. analytical solution
    e_term(n) = exp(-(n * hbar * w) / (kb * T))
    e_norm = sum([e_term(n) for n in 0:(d-1)])
    analytical_solution = (1 * e_term(1) + 2 * e_term(2)) / e_norm
    anal_error = abs(n_exp - 2 * analytical_solution)
    @test anal_error / n_exp < 0.01
end

@testset "Sample Thermal State Vector" begin
    d = 3; L = 2
    w = 7.5 * 10^9 * 2 * pi
    U = 0
    T = 0.1
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34

    N = 1000000
    n_op = nall(d, L)
    out = []
    for _ in 1:N
        v = sample_transmon_thermal_state(d, L, T, w, U)
        push!(out, expval(v, n_op))
    end
    n_exp = mean(out)

    # error vs. bose-einstein distribution
    be_dist(e, T) = 1 / (exp(e / (kb * T)) - 1)
    dist_error = abs(n_exp - 2 * be_dist(w * hbar, T))
    @test dist_error / n_exp < 0.01
end

end # testset
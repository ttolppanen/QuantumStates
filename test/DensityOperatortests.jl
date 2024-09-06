# using QuantumOperators
# using LinearAlgebra
# using Statistics

@testset "Density Operator" begin

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

    be_dist(e, T) = 1 / (exp(e / (kb * T)) - 1)
    @test abs(mean(out) - 2 * be_dist(w * hbar, T)) < 0.001
end

end # testset
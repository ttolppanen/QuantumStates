# using QuantumOperators
# using LinearAlgebra
# using Statistics

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
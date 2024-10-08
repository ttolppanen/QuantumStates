# using LinearAlgebra
# using SparseArrays

@testset "Complete Space" begin

function thermal_state(d, L)
    w = 7.5 * 10^9 * 2 * pi
    U = 50 / 1000 * w
    T = 0.1
    sample_transmon_thermal_state(d, L, T, w, U)
end

@testset "Orthogonality" begin
    @test QuantumStates.zerostate(3)' * QuantumStates.onestate(3) == 0.0
    @test zeroone(2, 2)' * onezero(2, 2) == 0.0
    @test onezero(3, 4)' * allone(3, 4) == 0.0
    @test bosonstack(4, 4, 2)' * bosonstack(4, 4, 1) == 0.0
    @test productstate(4, [0, 1, 3])' * productstate(4, [2, 1, 3]) == 0.0
end

@testset "Normalization" begin
    d = 5; L = 5
    @test norm(zeroone(d, L)) == 1.0
    @test norm(onezero(d, L)) == 1.0
    @test norm(allone(d, L)) == 1.0
    @test norm(allzero(d, L)) == 1.0
    @test norm(singleone(d, L, 2)) == 1.0
    @test norm(bosonstack(3, L, 1)) == 1.0
    @test norm(productstate(3, [1, 0, 2, 1, 0])) == 1.0
    @test norm(thermal_state(d, L)) == 1.0
end

@testset "Dimension" begin
    d = 3; L = 4
    @test length(zeroone(d, L)) == d^L
    @test length(onezero(d, L)) == d^L
    @test length(allone(d, L)) == d^L
    @test length(allzero(d, L)) == d^L
    @test length(singleone(d, L, 1)) == d^L
    @test length(bosonstack(d - 1, L, 1)) == d^L
    @test length(productstate(d, [1 for _ in 1:L])) == d^L
    @test length(thermal_state(d, L)) == d^L
end

function testtypes(state)
    @test isa(state, AbstractSparseArray{<:Number})
    @test typeof(state) == typeof(complex(state))
end
@testset "Types" begin
    d = 4; L = 4
    @testset "zeroone" begin testtypes(zeroone(d, L)) end
    @testset "onezero" begin testtypes(onezero(d, L)) end
    @testset "allone" begin testtypes(allone(d, L)) end
    @testset "allzero" begin testtypes(allzero(d, L)) end
    @testset "singleone" begin testtypes(singleone(d, L, 1)) end
    @testset "bosonstack" begin testtypes(bosonstack(d - 1, L, 1)) end
    @testset "productstate" begin testtypes(productstate(d, [0, 1, 0, 1])) end
    @testset "thermal state" begin testtypes(thermal_state(d, L)) end
end

@testset "Errors" begin
    @test_throws ArgumentError bosonstack(3, 3, 1; d=2)
end

end # testset
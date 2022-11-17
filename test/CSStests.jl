using LinearAlgebra
using SparseArrays

@testset "Complete Space" begin
    
@testset "Orthogonality" begin
    @test QuantumStates.zerostate(3)' * QuantumStates.onestate(3) == 0.0
    @test zeroone(2, 2)' * onezero(2, 2) == 0.0
    @test onezero(3, 4)' * allone(3, 4) == 0.0
    @test bosonstack(4, 4, 2)' * bosonstack(4, 4, 1) == 0.0
end

@testset "Normalization" begin
    d = 5; L = 5
    @test norm(zeroone(d, L)) == 1.0
    @test norm(onezero(d, L)) == 1.0
    @test norm(allone(d, L)) == 1.0
    @test norm(singleone(d, L, 2)) == 1.0
    @test norm(bosonstack(3, L, 1)) == 1.0
end

@testset "Dimension" begin
    d = 3; L = 4
    @test length(zeroone(d, L)) == d^L
    @test length(onezero(d, L)) == d^L
    @test length(allone(d, L)) == d^L
    @test length(singleone(d, L, 1)) == d^L
    @test length(bosonstack(d - 1, L, 1)) == d^L
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
    @testset "singleone" begin testtypes(singleone(d, L, 1)) end
    @testset "bosonstack" begin testtypes(bosonstack(d - 1, L, 1)) end
end

@testset "Errors" begin
    @test_throws ArgumentError bosonstack(3, 3, 1; d=2)
end

end #testset
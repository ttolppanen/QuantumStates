using ITensors

@testset "Matrix Product States" begin
    
@testset "Orthogonality" begin
    d = 4; L = 4
    inds = siteinds("Boson", L; dim = d)
    @test inner(zeroonemps(inds), zeroonemps(inds)) == 1.0
    @test inner(zeroonemps(inds), onezeromps(inds)) == 0.0
    @test inner(singleonemps(inds, 1), allonemps(inds)) == 0.0
end

@testset "Normalization" begin
    d = 4; L = 5
    @test norm(zeroonemps(d, L)) == 1.0
    @test norm(onezeromps(d, L)) == 1.0
    @test norm(allonemps(d, L)) == 1.0
    @test norm(singleonemps(d, L, 2)) == 1.0
end

@testset "Types" begin
    d = 3; L = 2
    @test isa(zeroonemps(d, L), MPS)
    @test isa(onezeromps(d, L), MPS)
    @test isa(allonemps(d, L), MPS)
    @test isa(singleonemps(d, L, 2), MPS)
end

end #testset

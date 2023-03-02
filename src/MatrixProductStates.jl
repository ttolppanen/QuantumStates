# using ITensors

# include("CompleteSpaceStates.jl")

export onezeromps
export zeroonemps
export allonemps
export allzeromps
export singleonemps
export bosonstackmps
export mps_to_array

# NOTE: MPS do not want to handle sparse vectors

# d : dimension; e.g. with qubits d = 2
# L : number of systems;
# indices : physical indices; refering to the physical indices of the MPS

function mps_to_array(state::MPS)
    out = prod(state) * combiner(siteinds(state))
    return Array(out, inds(out)...)
end

function onezeromps(indices::Vector{Index{Int64}})
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(onezero(d, L))
    return MPS(stateVector, indices)
end
function onezeromps(d::Integer, L::Integer)
     indices = siteinds("Boson", L; dim = d)
     return onezeromps(indices)
end

function zeroonemps(indices::Vector{Index{Int64}})
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(zeroone(d, L))
    return MPS(stateVector, indices)
end
function zeroonemps(d::Integer, L::Integer)
     indices = siteinds("Boson", L; dim = d)
     return zeroonemps(indices)
end

function allonemps(indices::Vector{Index{Int64}})
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(allone(d, L))
    return MPS(stateVector, indices)
end
function allonemps(d::Integer, L::Integer)
     indices = siteinds("Boson", L; dim = d)
     return allonemps(indices)
end

function allzeromps(indices::Vector{Index{Int64}})
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(allzero(d, L))
    return MPS(stateVector, indices)
end
function allzeromps(d::Integer, L::Integer)
     indices = siteinds("Boson", L; dim = d)
     return allzeromps(indices)
end

function singleonemps(indices::Vector{Index{Int64}}, i::Integer)
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(singleone(d, L, i))
    return MPS(stateVector, indices)
end
function singleonemps(d::Integer, L::Integer, i::Integer)
     indices = siteinds("Boson", L; dim = d)
     return singleonemps(indices, i)
end

function bosonstackmps(indices::Vector{Index{Int64}}, N::Integer, i::Integer)
    d = dim(indices[1])
    L = length(indices)
    stateVector = Vector(bosonstack(N, L, i; d))
    return MPS(stateVector, indices)
end
function bosonstackmps(N::Integer, L::Integer, i::Integer; d::Integer = N + 1)
     indices = siteinds("Boson", L; dim = d)
     return bosonstackmps(indices, N, i)
end
# using SparseArrays

export onezero
export zeroone
export allone
export allzero
export singleone
export bosonstack
export productstate

# d : dimension; e.g. with qubits d = 2
# L : number of systems;

function zerostate(d::Integer)
    out = spzeros(d)
    out[1] = 1.0
    return complex(out)
end

function onestate(d::Integer)
    out = spzeros(d)
    out[2] = 1.0
    return complex(out)
end

function onezero(d::Integer, L::Integer)
    zero = zerostate(d)
    one = onestate(d)
    sites = [i % 2 == 1 ? one : zero for i in 1:L] # Makes an array something like this [[1, 0], [0, 1], [1, 0]...]
    return kron(sites...)
end

function zeroone(d::Integer, L::Integer)
    zero = zerostate(d)
    one = onestate(d)
    sites = [i % 2 == 1 ? zero : one for i in 1:L] # Makes an array something like this [[1, 0], [0, 1], [1, 0]...]
    return kron(sites...)
end

function allone(d::Integer, L::Integer)
    one = onestate(d)
    sites = [one for _ in 1:L]
    return kron(sites...)
end

function allzero(d::Integer, L::Integer)
    zero = zerostate(d)
    sites = [zero for _ in 1:L]
    return kron(sites...)
end

function singleone(d::Integer, L::Integer, i::Integer)
    zero = zerostate(d)
    one = onestate(d)
    sites = [zero for _ in 1:L]
    sites[i] = one
    return kron(sites...)
end

function bosonstack(N::Integer, L::Integer, i::Integer; d::Integer = N + 1)
    d <= N ? throw(ArgumentError("d <= N")) : nothing
    zero = zerostate(d)
    stack = complex(spzeros(d))
    stack[N + 1] = 1.0
    sites = [zero for _ in 1:L]
    sites[i] = stack
    return kron(sites...)
end

function productstate(d::Integer, state::AbstractVector{<:Integer})
    sites = [complex(spzeros(d)) for _ in state]
    for (i, val) in enumerate(state)
        sites[i][val + 1] = 1.0
    end
    return kron(sites...)
end
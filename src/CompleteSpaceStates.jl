using SparseArrays

export onezero
export zeroone
export allone
export singleone
export bosonstack

#d : dimension; e.g. with qubits d = 2
#L : number of systems;

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
    sites = [i % 2 == 1 ? one : zero for i in 1:L] #Makes an array something like this [[1, 0], [0, 1], [1, 0]...]
    return kron(sites...)
end

function zeroone(d::Integer, L::Integer)
    zero = zerostate(d)
    one = onestate(d)
    sites = [i % 2 == 1 ? zero : one for i in 1:L] #Makes an array something like this [[1, 0], [0, 1], [1, 0]...]
    return kron(sites...)
end

function allone(d::Integer, L::Integer)
    one = onestate(d)
    sites = [one for _ in 1:L]
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
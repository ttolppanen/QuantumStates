# using SparseArrays

export onezero
export zeroone
export allone
export allzero
export singleone
export bosonstack
export productstate
export sample_transmon_thermal_state

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

function sample_transmon_thermal_state(d, T, w, U; err = 10^-10)
    if T == 0.0
        state = complex(spzeros(d))
        state[1] = 1.0
        return state
    end
    kb = 1.3806 * 10^-23
    hbar = 1.0546 * 10^-34
    exp_val(n) = exp(-hbar * (w * n - U / 2 * n * (n - 1)) / (kb * T))
    e_val = 1
    e_sum = 1
    n = 0
    while e_val / e_sum > err
        n += 1
        e_val = exp_val(n)
        e_sum += e_val
    end
    r = rand()
    p = 0.0
    for i in 1:(d-1)
        p += exp_val(i - 1) / e_sum
        if r < p
            state = complex(spzeros(d))
            state[i] = 1.0
            return state
        end
    end
    state = complex(spzeros(d))
    state[d] = 1.0
    return state
end

function sample_transmon_thermal_state(d, L, T, w, U; kwargs...)
    # sample a single site state
    s_s() = sample_transmon_thermal_state(d, T, w, U; kwargs...)
    if L == 1
        return s_s()
    end
    return kron([s_s() for _ in 1:L]...)
end
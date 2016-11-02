using Cliffords, BenchmarkTools

rand_phases = rand([1,-1,im,-im], 10_000)
function benchpauliphase(ns)
    r = Id
    XX = kron(X,X)
    for n in ns
        r = n ∘ XX
    end
    return r
end

rand_paulis = rand(allpaulis(2), 10_000)
benchpaulis(ps) = reduce(*, ps)

rand_c1seq = rand(1:24, 1_000)
benchc1(seq) = reduce(*, localclifford(n) for n in seq)

# create some of C2
C2 = Clifford[]
for _ in 1:100
    c = kron(localclifford(rand(1:24)), localclifford(rand(1:24)))
    if rand() > 0.5
        push!(C2, c)
    else
        push!(C2, c * CNOT * c)
    end
end

rand_c2seq = rand(1:100, 1_000)
benchc2(seq) = reduce(*, C2[n] for n in seq)

function benchpaulicliff()
    c = rand(C2)
    r = paulieye(2)
    IX = kron(Id, X)
    XI = kron(X, Id)
    for _ in 1:10_000
        r = c * IX
        r = c * XI
    end
    r
end

println("Pauli phase multiplication benchmarking:\n")
r0 = @benchmark benchpauliphase(rand_phases)
println(r0)

println("Pauli benchmarking:\n")
r1 = @benchmark benchpaulis(rand_paulis)
println(r1)

println("C1 benchmarking:\n")
r2 = @benchmark benchc1(rand_c1seq)
println(r2)

println("C2 benchmarking:\n")
r3 = @benchmark benchc2(rand_c2seq)
println(r3)

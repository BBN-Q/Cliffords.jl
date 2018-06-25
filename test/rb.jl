using Cliffords, Test
# randomized benchmarking

C2, C2dict = Cliffords.fullC2()

function rb_seq(n)
	seqn = rand(1:11520, N)
	seq = mapreduce(x -> C2[x], *, reverse(seqn))
	push!(seqn, C2dict[inv(seq)])
end

N = 20
for _ = 1:100
	@test mapreduce(x -> C2[x], *, reverse(rb_seq(N))) == cliffordeye(2)
end

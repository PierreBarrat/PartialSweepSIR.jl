################## Parameters ##################

@with_kw mutable struct Parameters
	# SIR
	N :: Int = 2 # Number of viruses
	α :: Float64 = 10
	γ :: Float64 = .01
	δ :: Float64 = 1

	# Geographic structure
	M :: Int = 1 # Number of geographic regions
	c :: Float64 = .1 # Strength of connection between region, in `[0,1]`
	C :: Matrix{Float64} = c * ones(M, M) .+ diagm(ones(M)*(1 - M*c))
end

################## Virus ##################

@with_kw mutable struct Virus
	S :: Float64
	I :: Float64
	C :: Float64
	R :: Float64 = 1 - S - I - C
	function Virus(S, I, C, R)
		@assert S <= 1 "S must be lower than 1; got $S"
		@assert I <= 1 "I must be lower than 1; got $I"
		return new(S, I, C, R)
	end
end

"""
	population(v::Virus)

Total host population when summing all compartments relative to `v`.
"""
population(v::Virus) = _population(v.S, v.I, v.C, v.R)
_population(S, I, C, R) = S + I + C + R

Base.vec(v::Virus) = [v.S, v.I, v.C, v.R]

################## Region ##################

Base.@kwdef mutable struct Region
	viruses :: Vector{Virus}
	K :: Matrix{Float64}

	function Region(viruses, K)
		@assert length(viruses) == size(K,1) == size(K,2) "Inconsistency between number\
			of viruses $(length(viruses)) and cross-immunity matrix $(size(K))"
		# @assert mapreduce(
		# 	v -> isapprox(population(v), 1, rtol = 1e-14),
		# 	*,
		# 	viruses,
		# 	init=true
		# ) "Host pop. is not one for some viruses."
		return new(viruses, K)
	end
end


"""
	Region(
		N::Int;
		I::AbstractVector = zeros(N),
		C::AbstractVector = zeros(N),
		R::AbstractVector = zeros(N),
		b=0.,
		f=0.,
		K = cross_immunity(N, b, f),
	)

	Region(N::Int; I = 0., C = 0., R = 0., b = 0., f = 0., K = cross_immunity(N, b, f))

`N`: number of viruses.
"""
function Region(
	N::Int;
	I::AbstractVector = zeros(N),
	C::AbstractVector = zeros(N),
	R::AbstractVector = zeros(N),
	b=0.,
	f=0.,
	K = cross_immunity(N, b, f),
)
	S = ones(N) - I - C - R
	viruses = [Virus(s, i, c, r) for (s, i, c, r) in zip(S, I, C, R)]
	return Region(; viruses, K)
end

viruses(r::Region) = r.viruses
Base.size(region::Region) = length(viruses(region))
Base.vec(region::Region) = @chain region PSS.viruses map(vec, _) vcat(_...)

function check_population_size(r::Region)
	mapreduce(*, viruses(region); init=true) do v
		@chain v population isapprox(1, rtol=1e-14)
	end
end

function Base.push!(
	region::Region, v::Virus;
	b::AbstractVector{Float64} = zeros(size(region)),
	f::AbstractVector{Float64} = zeros(size(region)),
	K = grow_cross_immunity(region.K, b, f),
)
	@assert length(b) == length(f) == size(K,1)-1 == size(K,1)-1 == size(region) "Size mismatch"
	region.K = K
	push!(region.viruses, virus)
end


################## SIRState ##################

"""
	SIRState(; regions, parameters)
"""
Base.@kwdef mutable struct SIRState
	parameters :: Parameters
	regions :: Vector{Region}

	function SIRState(parameters, regions)
		@assert !isempty(regions) "Must have at least one region"
		@assert length(regions) == parameters.M
		for (i,r) in enumerate(regions)
			@assert size(r) == parameters.N "Region $i has size $(size(r)) \
			instead of N=$(parameters.N)"
		end

		return new(parameters, regions)
	end
end

regions(X::SIRState) = X.regions
parameters(X::SIRState) = X.parameters

"""
	SIRState(parameters::Parameters; I0=0., R0=0., b=0., f=0.)
"""
function SIRState(parameters::Parameters; I0=0., C0=0., R0=0., b=0., f=0.)
	N = parameters.N
	regions = [
		Region(
			N;
			I = I0*ones(N), C = C0*ones(N), R = R0*ones(N), b, f,
		) for m in 1:parameters.M
	]
	return SIRState(parameters, regions)
end
"""
	SIRState(dat::Vector{Float64}, Ks, parameters::SIRParameters)

Create `SIRState` from a vector and an array of cross-immunity matrices `Ks`.
Such that `vec(out::SIRState) == dat`.
"""
function SIRState(dat::Vector{Float64}, Ks, parameters::Parameters)
	@assert length(dat) == parameters.M * parameters.N * 4
	regions = Region[]
	for i in 1:parameters.M
		viruses = Vector{Virus}(undef, parameters.N)
		for a in 1:parameters.N
			viruses[a] = Virus(;
				S = dat[sir_index(i, a, :S, parameters)],
				I = dat[sir_index(i, a, :I, parameters)],
				C = dat[sir_index(i, a, :C, parameters)],
				R = dat[sir_index(i, a, :R, parameters)],
			)
		end
		r = Region(; viruses, K = Ks[i])
		push!(regions, r)
	end
	return SIRState(; parameters, regions)
end

function Base.getindex(X::SIRState, i::Int, a::Int, g)
	g = compartment_to_symb(g)
	return getfield(X.regions[i].viruses[a], g)
end
Base.getindex(X::SIRState, i::Int, A::AbstractRange, g) = [X[i, a, g] for a in A]
Base.getindex(X::SIRState, i::Int, ::Colon, g) = X[i, 1:(X.parameters.N), g]
Base.getindex(X::SIRState, I::AbstractRange, a::Int, g) = [X[i, a, g] for i in I]
Base.getindex(X::SIRState, Colon, a::Int, g) = X[1:(X.parameters.M), a, g]
function Base.getindex(X::SIRState, ::Colon, ::Colon, g)
	return vcat([X[i, :, g] for i in 1:X.parameters.M]...)
end

Base.vec(X::SIRState) = @chain X regions map(vec, _) vcat(_...)

"""
"""
function Base.push!(
	X::SIRState,
	v::Virus,
	b::AbstractVector{<:AbstractVector},
	f::AbstractVector{<:AbstractVector}
)
end

################## Utils ##################

function symb_to_virus(a::Symbol)
	return if a == :wt
		1
	elseif a == :m
		2
	else
		error("Unknown symbol $a")
	end
end
symb_to_virus(a::Number) = a

function symb_to_compartment(g::Symbol)
	if g == :S
		return 1
	elseif g == :I
		return 2
	elseif g == :C
		return 3
	elseif g == :R
		return 4
	else
		error()
	end
end
symb_to_compartment(g::Number) = g
compartment_to_symb(g::Symbol) = g
function compartment_to_symb(g::Number)
	return if g == 1
		:S
	elseif g == 2
		:I
	elseif g == 3
		:C
	elseif g == 4
		:R
	else
		error()
	end
end


"""
	sir_index(i, a, g, N::Int)

Index of strain `a` in SIR class `g` in region `i`.

- `a` is a number in `[1, N]` (or `:wt`/`:m` if `N==2`)
- `g` can be any of `(:S, :I, :R)`
- `i` is a number in `[1, M]`: it's the region

## Drawing
```
[
	## a --> Virus
	[S, I, C, R]
	## i --> Region
	[V1, V2, ..., VN] # where Va ~ [Sa, Ia, Ca, Ra]
	## The whole state
	[R1, R2, ..., RM] # where Ri ~ [V1, V2, ..., VN]
]
```

- Linear size of each region: `4*N`
- Total size of array `4*N*M` if `M` is the number of regions
"""
function sir_index(i, a, g, N::Int)
	# @assert 1 <= i <= M "Error accessing region $i in size $M SIR model"
	# @assert 1 <= a <= N "Error accessing region $i in size $M SIR model"

	a = symb_to_virus.(a)
	g = symb_to_compartment.(g)
	vsize = 4 # number of compartments
	rsize = vsize*N

	region_offset = (i .- 1) * rsize
	virus_offset = (a .- 1) * vsize

	# @info "region_offset: $(region_offset), compartment_offset = $(compartment_offset)"

	return region_offset .+ virus_offset .+ g
end
sir_index(i, ::Colon, g, N::Int) = sir_index(i, 1:N, g, N)
sir_index(i, a, g, p::Parameters) = sir_index(i, a, g, p.N)
sir_index(i, a, g, s::SIRState) = sir_index(i, a, g, s.parameters)




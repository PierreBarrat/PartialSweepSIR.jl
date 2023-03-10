Base.@kwdef struct Solution
	sol
	tspan :: Tuple{Float64, Float64}
	parameters :: Parameters
	ϕ :: Vector{Vector{Float64}} # ϕ[j][a]
	K :: Vector{Matrix{Float64}}
end

function (sol::Solution)(t::Real)
	@assert t <= sol.tspan[2] "Time span of solution $(sol.tspan) - got $t"
	return SIRState(sol.sol(t), sol.K, sol.parameters, sol.ϕ)
end

function Base.getindex(sol::Solution, tvals::AbstractVector, i, a, g)
	idx = sir_index(i, a, g, sol.parameters)
	return map(t -> sol.sol(t)[idx], tvals)
end
function Base.getindex(sol::Solution, t::Number, i, g, a)
	idx = sir_index(i, g, a, sol.parameters)
	sol.sol(t)[idx]
end

"""
	simulate(X::SIRState, tspan)

Return a `Solution` object.
"""
function simulate(X::SIRState, tspan)
	u0 = vec(X)
	ϕ = map(regions(X)) do r
		map(v -> v.ϕ, r.viruses)
	end
	p = (
		params = X.parameters,
		K = [r.K for r in regions(X)],
		ϕ = ϕ,
	)
	return Solution(;
		sol = solve(ODEProblem(SIR!, u0, tspan, p), Tsit5()),
		tspan = tspan,
		parameters = parameters(X),
		ϕ = ϕ,
		K = p.K,
	)
end

function gradient(X::SIRState)
	u = vec(X)
	du = similar(u)
	p = (
		params = X.parameters,
		K = [r.K for r in regions(X)],
		ϕ = map(regions(X)) do r
			map(v -> v.ϕ, r.viruses)
		end
	)
	SIR!(du, u, p, 0)

	dX = deepcopy(X)
	(N, M) = (X.parameters.M, X.parameters.N)
	for i in 1:M, a in 1:N, g in (:I, :S, :C, :R)
		idx = sir_index(i, a, g, N)
		setfield!(dX.regions[i].viruses[a], g, du[idx])
	end
	return dX
end

function SIR!(du, u, p, t)
	#=
	Organization of u
	u ~ [R1, R2, ..., RM]
	Ri ~ [V1, V2, ..., VN]
	Va ~ [S, I, C, R]
	--> u is a vector of size 4*N*M

	p is a named tuple with fields
	- `params` which is the Params struct
	- `K` which is a vector of cross-immunities
	=#

	@unpack M, N, α, γ, δ, C = p.params # K is in p.K
	du .= 0.

	for i in 1:M, a in 1:N
		# S
		idx = sir_index(i, a, :S, N)
		for j in 1:M, b in 1:N
			# region j infecting region i
			# and virus b generating immunity to a
			du[idx] -=
				α * p.ϕ[j][b] * C[i,j] * u[sir_index(i,a,:S,N)] * p.K[i][a,b] * u[sir_index(j,b,:I,N)]
		end
		du[idx] += γ * u[sir_index(i, a, :R, N)]

		# I
		idx = sir_index(i, a, :I, N)
		for j in 1:M
			# region j infecting region i
			du[idx] += α * p.ϕ[j][a] * C[i,j] * u[sir_index(i, a, :S, N)] * u[sir_index(j, a, :I, N)]
		end
		du[idx] -= δ * u[idx]

		# C
		idx = sir_index(i, a, :C, N)
		for j in 1:M, b in Iterators.filter(!=(a), 1:N)
			du[idx] +=
				α * p.ϕ[j][b] * C[i,j] * u[sir_index(i,a,:S,N)] * p.K[i][a,b] * u[sir_index(j,b,:I,N)]
		end
		du[idx] -= δ * u[idx]

		# R
		idx = sir_index(i, a, :R, N)
		du[idx] += δ * (u[sir_index(i,a,:I,N)] + u[sir_index(i,a,:C,N)])
		du[idx] -= γ * u[idx]
	end

	@debug _conserved(du, p), u, p
	return nothing
end
function _conserved(du, parameters)
	return map(1:parameters.sir.M) do i
		map(1:parameters.sir.N) do a
			z = du[sir_index(i, a, :S, parameters.sir)]
			z += du[sir_index(i, a, :I, parameters.sir)]
			z += du[sir_index(i, a, :C, parameters.sir)]
			z += du[sir_index(i, a, :R, parameters.sir)]
		end
	end
end

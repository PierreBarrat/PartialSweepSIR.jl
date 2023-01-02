### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ a6e6086e-64cb-11ed-0ccb-d379d6ea917f
begin
	using Revise
	using Pkg; Pkg.activate()
	using Parameters
	using PartialSweepSIR
	using Plots
end

# ╔═╡ fdf3f050-56b0-445e-b657-b8e0ec3572d8
R1 = let
	S0 = .4
	I0 = 1e-6
	C0 = 0.
	R0 = 1 - S0 - I0 - C0
	v1 = PSS.Virus(S=S0, I=I0, C=C0, R=R0)
	R1 = PSS.Region(viruses = [v1], K = ones(1,1))
end;

# ╔═╡ a52a2bb7-67eb-4ec0-9025-3e4a7983f490
params = let
	N = 1
	M = 1
	PSS.Parameters(; N, M)
end

# ╔═╡ a114f1f2-4dbd-4d9b-a927-cfefa72bf55c
state = PSS.SIRState(; regions = [R1], parameters = params);

# ╔═╡ 655cff68-1d4a-40b6-9626-a92add581918
state.regions[1].K

# ╔═╡ 90e4f8bc-bfeb-480a-adb3-14aeb2a4d6e8
let
	u = vec(state)
	du = similar(u)
	p = (params = state.parameters, K = [r.K for r in PSS.regions(state)])
	PSS.SIR!(du, u, p, 1)
	du
end

# ╔═╡ 8d289bdb-4c16-4a2d-a193-ea3aa2bc8187
sol = PSS.simulate(state, (0, 100));

# ╔═╡ ac37cad1-8a24-4c7e-a7aa-b1e0fd5b9ca8
state_eq = PSS.equilibrium(state)

# ╔═╡ becfd96b-02a5-4923-8cf7-bad20d51e588
let
	tvals = range(sol.tspan..., length=100)
	@unpack α, γ, δ = params

	p = plot(ylim = (1e-9, 1), yscale=:log10, legend=:bottomright)
	for (i, g) in enumerate((:S, :I, :C, :R))
		X1 = sol[tvals, 1, 1, g]
		Xeq = state_eq[1, 1, g]

		plot!(tvals, X1, label="$g", color = i)
		hline!([Xeq], label="", color = i, linestyle=:dash)
	end
	p

end

# ╔═╡ Cell order:
# ╠═a6e6086e-64cb-11ed-0ccb-d379d6ea917f
# ╠═fdf3f050-56b0-445e-b657-b8e0ec3572d8
# ╠═a52a2bb7-67eb-4ec0-9025-3e4a7983f490
# ╠═a114f1f2-4dbd-4d9b-a927-cfefa72bf55c
# ╠═655cff68-1d4a-40b6-9626-a92add581918
# ╠═90e4f8bc-bfeb-480a-adb3-14aeb2a4d6e8
# ╠═8d289bdb-4c16-4a2d-a193-ea3aa2bc8187
# ╠═becfd96b-02a5-4923-8cf7-bad20d51e588
# ╠═ac37cad1-8a24-4c7e-a7aa-b1e0fd5b9ca8

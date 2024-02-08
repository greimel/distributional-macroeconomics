### A Pluto.jl notebook ###
# v0.19.38

using Markdown
using InteractiveUtils

# ╔═╡ 5d9a1abe-a51d-41c6-b86d-a1360140c05a
md"""
!!! danger "Under construction!"

	**This notebook is not ready for public consumption.** Use at your own risk.
"""

# ╔═╡ 9e4ec644-4566-4066-a958-ba7255fbe2f5
md"""
`continuous-time.jl` | **Version 0.1** | *last updated: Apr 26 2022*
"""

# ╔═╡ 36a22f3b-9f7d-4075-95d3-5787e40dbad7
md"""
# Heterogenous Agents in Continuous Time

*The model setup is adapted from [an example in `EconPDEs.jl`.](https://github.com/matthieugomez/EconPDEs.jl/blob/ecdefd0b1e52b2ec3bcebb9498d39057c26afd74/examples/ConsumptionProblem/AchdouHanLasryLionsMoll_TwoAssets.jl)*

"""

# ╔═╡ bdbf11b8-e6c1-4f9c-8ab4-cb714655fed0
function μa_from_va(va, (; y, a), (; r, γ))
 	c = va^(-1 / γ)
	μa = y + r * a - c
	
	(; c, μa)
end

# ╔═╡ 8a3e8d6a-23f4-41ba-93bb-235d10bdd534
begin
	Base.@kwdef struct AchdouHanLasryLionsMollModel
	    # income process parameters
	    κy::Float64 = 0.1
	    ybar::Float64 = 1.0
	    σy::Float64 = 0.07
	
	    r::Float64 = 0.03
	
	    # utility parameters
	    ρ::Float64 = 0.05
	    γ::Float64 = 2.0
	
	    amin::Float64 = 0.0
	    amax::Float64 = 500
	end
	
	function (m::AchdouHanLasryLionsMollModel)(state::NamedTuple, value::NamedTuple)
	    (; κy, σy, ybar, r, ρ, γ, amin, amax) = m    
	    (; y, a) = state
	    (; v, vy_up, vy_down, va_up, va_down, vyy, vya, vaa) = value
	    μy = κy * (ybar - y)
		
		vy = μy ≥ 0 ? vy_up : vy_down
	
	    va_up = max(va_up, eps())
	    va = va_up
	    (; c, μa) = μa_from_va(va, state, m)
		if μa ≤ 0
			va = va_down
			(; c, μa) = μa_from_va(va, state, m)
		end
		if (a ≈ amin) && (μa ≤ 0.0)
	        va = (y + r * amin)^(-γ)
	        c = y + r * amin
	        μa = 0.0
	    end
	    vt = - (c^(1 - γ) / (1 - γ) + μa * va + μy * vy + 0.5 * vyy * σy^2 - ρ * v)
	    return (; vt), (; μy, μa, c, va, vy, vyy)
	end
	
end

# ╔═╡ 23084da8-6d9d-40c4-a384-3d4066de62e5
@chain df_out begin
	#@groupby(:a)
	#@combine(:π = sum(:π))
	data(_) * mapping(:a, :π, color = :y => nonnumeric) * visual(Lines)
	AlgebraOfGraphics.draw
end

# ╔═╡ 212cca0b-314a-49df-9c20-057d70bfe277
m = AchdouHanLasryLionsMollModel(r=0.045, amax = 5)

# ╔═╡ f4da2066-e7fd-4641-96ae-22b2708e0399
distribution = Gamma(2 * m.κy * m.ybar / m.σy^2, m.σy^2 / (2 * m.κy))

# ╔═╡ 0708f272-da24-401e-b3c9-a212dd003a38
stategrid = OrderedDict(:y => range(quantile(distribution, 0.001), quantile(distribution, 0.9999), length = 10), 
                        :a =>  range(m.amin, m.amax, length = 200)
                        )

# ╔═╡ 9dd13311-3e5e-4f4b-998b-4b6f2bdd2bbb
states_df = let 
	ks = tuple(keys(stategrid)...)
	vs = values(stategrid)

	Iterators.product(vs...) .|> NamedTuple{ks} |> vec |> DataFrame
end

# ╔═╡ ae84468f-ee75-4763-812c-79106747908a
md"""
# Infinite horizon
"""

# ╔═╡ 920b0920-4b52-44e0-8779-ead83552d9cf
result_∞ = let
	yend = OrderedDict(:v => [log(y + max(a, 0.0)) for y in stategrid[:y], a in stategrid[:a]])
	result = pdesolve(m, stategrid, yend)
	@assert result.residual_norm <= 1e-5
	result
end

# ╔═╡ 0010335f-1657-4e2e-98f9-60d917f98f8c
df_out = let
	#A
	y_grid = stategrid[:y] |> collect
	a_grid = stategrid[:a] |> collect

	μa = result_∞.optional[:μa]
	va = result_∞.optional[:va]

	μy  = result_∞.optional[:μy]
	vy  = result_∞.optional[:vy]
	vyy = result_∞.optional[:vyy]
	(; σy) = m

	n_y = length(y_grid)
	n_a = length(a_grid)

	state_indices = Iterators.product(1:n_y, 1:n_a) |> collect
	lin_indices = LinearIndices(state_indices)
	@assert size(μa) == size(state_indices)

	A = zeros(n_y * n_a, n_y * n_a)
	for (i, (i_y, i_a)) ∈ enumerate(state_indices)
		i_to = nothing
		if μa[i] > 0 && i_a < n_a
			i_to = lin_indices[i_y, i_a + 1]
		elseif i_a > 1
			i_to = lin_indices[i_y, i_a - 1]
		end
		if !isnothing(i_to)
			A[i, i_to] = abs(μa[i]) * va[i]
		end
	end

	for (i, (i_y, i_a)) ∈ enumerate(state_indices)
		i_to = nothing
		if μy[i] > 0 && i_y < n_y
			i_to = lin_indices[i_y + 1, i_a]
		elseif i_y > 1
			i_to = lin_indices[i_y - 1, i_a]
		end
		if !isnothing(i_to)
			A[i, i_to] = abs(μy[i]) * vy[i] + 0.5 * vyy[i] * σy^2
		end
	end

	for i ∈ n_y * n_a
		A[i,i] = -sum(A[i,:])
	end

	df = copy(states_df)
	π = QuantEcon.gth_solve(A)
	df.π = vec(π)

	df
end

# ╔═╡ 2e3dac19-8b48-44f8-8c79-c8ef6d4fb54a
result_∞.optional[:va]

# ╔═╡ 7c8e2c71-47a5-49b6-811c-551a0f7ad756
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!.(Ref(ax), Ref(stategrid[:a]), eachrow(result_∞.optional[:μa])) 

	fig
end

# ╔═╡ 8cc86d11-e206-49ed-972f-583731d43bf1
let
	fig = Figure()
	ax = Axis(fig[1,1])

	lines!.(Ref(ax), eachrow(result_∞.optional[:c])) 

	fig
end

# ╔═╡ a719e2d7-1666-4d59-9c1e-43475e908e14
md"""
## Simulate
"""

# ╔═╡ eb9589e6-7310-416e-bd7f-d1a42880eed5
function simulate(m::AchdouHanLasryLionsMollModel, stategrid, result_ip, T)
   yT, aT, cT = zeros(T), zeros(T), zeros(T)
   y, a = 1.0, 1.0
   for t in 1:T
        y += m.κy * (m.ybar - y) + m.σy * rand(Normal())
        a += result_ip[:μa](y, a)
        yT[t] = y
        aT[t] = a
        cT[t] = result_ip[:c](y, a)
   end
    return (; t = 1:T, yT, aT, cT)
end

# ╔═╡ 63697439-bc45-4ef5-bdc7-9691b15f1332
sim_df = let
	yend = OrderedDict(:v => [log(y + max(a, 0.0)) for y in stategrid[:y], a in stategrid[:a]])
	(; zero, residual_norm, optional) = pdesolve(m, stategrid, yend)
	result = optional
	#(; y, result, distance)
	result_ip =  Dict(x => interpolate(tuple(values(stategrid)...), result[x], Gridded(Linear())) for x in keys(result))
	out = simulate(m, stategrid, result_ip, 90)
	
	DataFrame(out)
end

# ╔═╡ 3e4e0e34-8db8-4e0f-a593-5ac822bb6c36
let
	fig, ax, _ = lines(sim_df.yT, label = L"y_t")
	lines!(sim_df.aT, label = L"a_t")
	lines!(sim_df.cT, label = L"c_t")

	axislegend(ax)
	fig
end

# ╔═╡ 55ab0920-7ece-44e0-b937-7ee0e9d7cc3d
J = 20

# ╔═╡ 8ebbc2fd-6dc1-401e-a002-b2d7ad71dc31
md"""
## Finite horizon over $J years
"""

# ╔═╡ 8db3df4c-f2a2-40b7-82ee-91b9dc609e37
result_T = let
	yend = OrderedDict(:v => [max(a + y)^(1-m.γ)/(1-m.γ) for y in stategrid[:y], a in stategrid[:a]])

	τs = range(0, stop = J, step = 1)

	result  = pdesolve(m, stategrid, yend, τs)
	@assert maximum(result.residual_norm) <= 1e-5
	result
end

# ╔═╡ 22a8ac01-2af9-4b1f-a1c7-7a663c670828
md"""
# Neoclassical Growth Model and Real Business Cycle Model
"""

# ╔═╡ e2d7380e-19bf-445a-96ef-36613a42943e
md"""
## Neoclassical Growth Model

The HJB equation is
```math
\rho v(k) = \max_{c} u(c) + v'(k)(F(k) - c - \delta k)
```
The first order condition is ``u'(c) = v'(k) \iff c^* = (u')^{-1}(v'(k))``
or
```math
\rho v(k) = u(c^*) + v'(k)(F(k) - c^* - \delta k)
```

## Real Business Cycle Model

The HJB equation is
```math
\rho v(k, z) = \max_{c} u(c) + v_k(k, z)(F(k, z) - c - \delta k) + \underbrace{v_z(k, z) \mu(z) + \frac{1}{2} v_{zz}(k, z) \sigma^2(z)}_{\text{exo}}
```
The first order condition is ``u'(c) = v'(k) \iff c^* = (u')^{-1}(v'(k))``
or
```math
\rho v(k) = \underbrace{u(c^*) + v'(k)(F(k) - c^* - \delta k)}_{\text{endo}} + \underbrace{v_z(k, z) \mu(z) + \frac{1}{2} v_{zz}(k, z) \sigma^2(z)}_{\text{exo}}
```

"""

# ╔═╡ 5347b4ff-543a-42cf-8aea-a173e9b66f4d
u(c, (; γ)) = c > 0 ? γ == 1 ? log(c) : c^(1-γ)/(1-γ) : 10 * c - 100

# ╔═╡ 5c378bf5-f5d6-4e25-a1dc-3e5951fca01b
u_prime(c, (; γ)) = c^(-γ)

# ╔═╡ 2d563334-4032-45ce-b90b-65fe9b59efda
u_prime_inv(x, (; γ)) = x^(-1/γ)

# ╔═╡ 296a8245-a813-484c-86b8-506791d55ac7
map(c -> u_prime_inv(u_prime(c, (; γ=2.0)), (; γ=2.0)), 1:100)

# ╔═╡ fc8b67ef-faf4-4e8a-ad62-78c018d01747
NGM_steady_state((; α, ρ, δ)) = (α*1/(ρ+δ))^(1/(1-α))

# ╔═╡ 1f246a97-77b0-477f-b97d-7c0280f78310
(; out, k_grid) = let
	N = 300
	param = (α = 0.3, γ = 5.0, δ = 0.05, ρ = 0.05, κy = 0.1, ybar = 1.0, σy = 0.05) # √(eps()))
	kss = NGM_steady_state(param)

	k_min = 0.01 #0.0001 * kss
	k_max = 5.0 #1.5 * kss

	@info (; kss, k_min, k_max)
	k_grid = range(k_min, k_max, length = N)

	# exogenous state: productivity
	distribution = Gamma(2 * param.κy * param.ybar / param.σy^2, param.σy^2 / (2 * param.κy))
	y_grid = range(quantile(distribution, 0.01), quantile(distribution, 0.99), length = 10) 
	
	rbc = let
		stategrid = OrderedDict(:k => k_grid, :y => y_grid)
		solend = OrderedDict(:v => [u(y * k ^ param.α, param) / param.ρ for k ∈ k_grid, y ∈ y_grid])
		(; stategrid, solend)
	end

	ngm = let
		stategrid = OrderedDict(:k => k_grid)
		solend = OrderedDict(:v => [u(1 * k ^ param.α, param) / param.ρ for k ∈ k_grid])
		(; stategrid, solend)
	end

	function c_kdot_vk(vk, (; y, k), (; α, δ, γ))
		vk = max(vk, eps())
		c_star = u_prime_inv(vk, (; γ))
		k_dot = y * k^α - δ*k - c_star
		(; vk, c_star, k_dot)
	end

	function c₀_kdot_vk((; y, k), (; α, δ, γ))
		c_star = y * k^α - δ * k
        vk = u_prime(c_star, (; γ))
        k_dot = 0.0
		
		(; vk, c_star, k_dot)
	end
	
	function f(state::NamedTuple, sol::NamedTuple)
		(; α, γ, δ, ρ, κy, ybar, σy) = param
		
		# Only relevant for RBC
		if haskey(state, :y)
			(; y) = state
			(; vy_up, vy_down, vyy) = sol
			μy = κy * (ybar - y)
    		vy = (μy >= 0) ? vy_up : vy_down
			exo = μy * vy + 0.5 * vyy * σy^2
		else
			exo = 0
			state = (; y = 1, state...)
		end

		(; k) = state
        (; v, vk_up, vk_down) = sol

		(; vk, c_star, k_dot) = c_kdot_vk(vk_up, state, param)
		if k_dot ≤ 0.0
			(; vk, c_star, k_dot) = c_kdot_vk(vk_down, state, param)
		end
		if (k ≈ k_min) && (k_dot ≤ 0.0)
			(; vk, c_star, k_dot) = c₀_kdot_vk(state, param)
    	end
		endo = u(c_star, (; γ)) + vk * k_dot
        
		vt = ρ * v - (endo + exo)
		(vt = vt,), (; k_dot, c_star)
	end

	#ys, residual_norms = pdesolve(f, stategrid, solend, range(0, 1000, length = 100))

	out = pdesolve(f, rbc...)
	(; out, k_grid)
end

# ╔═╡ 73a5630f-e4ad-402c-85e4-409c1092a5e8
let
	(; zero, residual_norm, optional) = out
	fig = Figure()
	lines(fig[1,1], k_grid, zero[:v][:,1], axis=(; title = "value"))
	lines(fig[1,2], k_grid, optional[:c_star][:,1], axis=(; title = "consumption"))
	lines(fig[1,3], k_grid, optional[:k_dot][:,1], axis=(; title = "savings"))
	if length(size(zero[:v])) > 1
		surface(fig[2,1], zero[:v], axis = (type = Axis3, title = "value"))
		surface(fig[2,2], optional[:c_star], axis = (type = Axis3, title = "consumption"))
	end	
	fig
end

# ╔═╡ 103b4796-162b-4ac4-b34c-42fd584efb5c
md"""
# Consumption-saving
"""

# ╔═╡ 900cf795-953b-4437-970e-944b05c279d7
md"""
The HJB equation is
```math
\rho v(a, z) = \max_{c} u(c) + v'(a)(z w + r a - c) + \text{stuff related to incomes}
```
The first order condition is ``u'(c) = v_a(a,z) \iff c^* = (u')^{-1}(v_a(a,z))``
or
```math
\rho v(a, z) = u(c^*) + v_a(a,z)(z w + r a - c^*) + \text{stuff related to incomes}
```
"""

# ╔═╡ f4547de5-9254-4556-96ed-03e2221db3e8
md"""
## Useful functions
"""

# ╔═╡ ba08a29a-7480-4d2c-85ee-83ed67f3f5f5
function c_adot_va(va, (; z, a), (; r, w, γ))
	va = max(va, eps())
	cstar = u_prime_inv(va, (; γ))
	adot = z * w + a * r - cstar
	
	(; va, cstar, adot)
end

# ╔═╡ 40adb798-988e-4905-91c6-b0a14ceecacb
function c₀_adot_va((; z, a), (; r, w, γ))
	cstar = z * w + a * r
    va = u_prime(cstar, (; γ))
	
	(; va, cstar, adot=0.0)
end

# ╔═╡ 0bbf0fcf-9671-4535-8006-b345d3a43c69
function huggett_inner((; va_up, va_down), z, state, param)
	(; γ, amin) = param
	(; a) = state
	state = (; z, state...)
	
	(; va, cstar, adot) = c_adot_va(va_up, state, param)
	if adot ≤ 0.0
		(; va, cstar, adot) = c_adot_va(va_down, state, param)
	end
	if (a ≈ amin) && (adot ≤ 0.0)
		(; va, cstar, adot) = c₀_adot_va(state, param)
   	end
	endo = u(cstar, param) + va * adot

	(; endo, adot, va, cstar)
end

# ╔═╡ 381fdb71-38af-4570-b6e1-0bfe72c7042d
begin
	Base.@kwdef struct HuggettDiffusion
		γ = 2.0
		ρ = 0.05
		r = 0.04
		w = 1.0
		amin = 0.0
		amax = 15.0
		an = 100

		κz = 0.1
		zbar = 1.0
		σz = 0.1
	end

	function (m::HuggettDiffusion)(state::NamedTuple, sol::NamedTuple)
		(; ρ, κz, zbar, σz) = m

		begin
			(; z) = state
			(; vz_up, vz_down, vzz) = sol
			μz = κz * (zbar - z)
    		vz = (μz >= 0) ? vz_up : vz_down
			exo = μz * vz + 0.5 * vzz * σz^2
		end
		
		(; endo, adot, cstar, va) = huggett_inner(sol, state.z, state, m)

		(; v) = sol

		vt = ρ * v - (endo + exo)

		(; vt), (; adot, cstar, va, vz, vzz, μz)
	end	
end

# ╔═╡ 113206ad-361e-4214-9208-99fef39ea96a
out2, grids2, m2 = let
	m = HuggettDiffusion()
		
	agrid = range(m.amin, m.amax, m.an)
	zgrid = let
		distribution = Gamma(2 * m.κz * m.zbar / m.σz^2, m.σz^2 / (2 * m.κz))
		range(quantile(distribution, 0.01), quantile(distribution, 0.99), length = 10)
	end
	stategrid = OrderedDict(:a => agrid, :z => zgrid)
	
	solend = OrderedDict(
		:v => [log(a + z) for a ∈ agrid, z ∈ zgrid]
	)

	out = pdesolve(m, stategrid, solend)
	(; out, grids = (; agrid, zgrid), m)
end

# ╔═╡ 65ce5940-3f7c-4036-bf51-bc90e2f8b235
@chain begin
	stationary_distribution(m2, grids2, out2.optional)
	data(_) * mapping(:a, :z, :π) * visual(Surface)
	AlgebraOfGraphics.draw(axis = (type = Axis3, ))
end

# ╔═╡ fb9684f8-31e6-4da9-94fb-8cee9c9fe69e
function stationary_distribution(m, grids, optional)

	(; agrid, zgrid) = grids

	va = optional[:va]
	μa = optional[:adot]

	μz  = optional[:μz]
	vz  = optional[:vz]
	vzz = optional[:vzz]
	(; σz) = m

	nz = length(zgrid)
	na = length(agrid)

	state_indices = Iterators.product(1:na, 1:nz) |> collect
	lin_indices = LinearIndices(state_indices)
	
	@assert size(μa) == size(state_indices)

	A = zeros(nz * na, nz * na)
	for (i, (i_a, i_z)) ∈ enumerate(state_indices)
		i_to = nothing
		if μa[i] > 0 && i_a < na
			i_to = lin_indices[i_a + 1, i_z]
		elseif i_a > 1
			i_to = lin_indices[i_a - 1, i_z]
		end
		if !isnothing(i_to)
			A[i, i_to] = abs(μa[i]) * va[i]
		end
	end

	for (i, (i_a, i_z)) ∈ enumerate(state_indices)
		i_to = nothing
		if μz[i] > 0 && i_z < nz
			i_to = lin_indices[i_a, i_z + 1]
		elseif i_z > 1
			i_to = lin_indices[i_a, i_z - 1]
		end
		if !isnothing(i_to)
			A[i, i_to] = abs(μz[i]) * vz[i] + 0.5 * vzz[i] * σz^2
		end
	end

	for i ∈ nz * na
		A[i,i] = -sum(A[i,:])
	end

	df = @chain begin
		Iterators.product(agrid, zgrid)
		collect
		vec
		DataFrame
		rename([:a, :z])
	end
	# make sure its irreducible!!!
	π = QuantEcon.gth_solve(A)
	df.π = vec(π)

	df
end

# ╔═╡ b22052d8-e122-46d9-a00a-30b1d7d7555e
let
	out = out2
	(; agrid) = grids2
	(; zero, residual_norm, optional) = out
	fig = Figure()
	lines(fig[1,1], agrid, zero[:v][:,1], axis=(; title = "value"))
	lines(fig[1,2], agrid, optional[:cstar][:,1], axis=(; title = "consumption"))
	lines(fig[1,3], agrid, optional[:adot][:,1], axis=(; title = "savings"))
	surface(fig[2,1], zero[:v], axis = (type = Axis3, title = "value"))
	surface(fig[2,2], optional[:cstar], axis = (type = Axis3,))
	
	fig
end

# ╔═╡ 786c5de0-2639-47f8-8c12-a89813d57075
function clean_variables(nt, solname, statename, n)
	map(1:n) do i
		sol_key = Symbol(solname, i)
		up_key = Symbol(sol_key, statename, "_up") => Symbol(solname, statename, "_up")
		down_key = Symbol(sol_key, statename, "_down") => Symbol(solname, statename, "_down")
		
		(; solname => nt[sol_key], up_key[2] => nt[up_key[1]], down_key[2] => nt[down_key[1]])
	end |> DataFrame
end

# ╔═╡ 2f577cd6-2f4d-4998-99ae-dbf22d7a5a65
begin
	Base.@kwdef struct HuggettMC
		γ = 2.0
		ρ = 0.05
		r = 0.03
		w = 1.0
		amin = 0.0
		amax = 50.0
		an = 100
		zgrid = [.8, 2.0, 3.6]
		Λ = [-0.5 0.3 0.2;
			0.25 -0.5 0.25;
			0.2 0.3 -0.5]
	end
		
	function (m::HuggettMC)(state::NamedTuple, sol::NamedTuple)
		(; ρ, Λ, zgrid) = m
		nz = length(zgrid)
		
		sol_clean = clean_variables(sol, :v, :a, nz)
		
		nts = map(enumerate(eachrow(sol_clean))) do (i, row)
			(; i, huggett_inner(row, zgrid[i], state, m)...)
		end

		(; v) = sol_clean
		(; endo, adot, cstar) = DataFrame(nts)

		vt = ρ * v - (endo + Λ * v)

		(; (Symbol(:v, i, :t) => vt[i] for i ∈ 1:nz)...),
		(; 
			(Symbol(:c, i) => cstar[i] for i ∈ 1:nz)...,
			(Symbol(:s, i) => adot[i] for i ∈ 1:nz)...,
		)
	end	
end

# ╔═╡ 6ee07c2c-de88-4d38-9c4e-afb3e60b4e66
out3, agrid3 = let
	m = HuggettMC()
		
	agrid = range(m.amin, m.amax, m.an)
	stategrid = OrderedDict(:a => agrid)
	
	solend = OrderedDict(
		(Symbol(:v, i) => [log(a + m.zgrid[i]) for a ∈ agrid]) for i ∈ 1:length(m.zgrid)
	)

	out = pdesolve(m, stategrid, solend)
	(; out, agrid)
end

# ╔═╡ ae2a8f97-1d31-4b54-a0d3-31bfb87e4607
let
	using NamedTupleTools: namedtuple
	nt = namedtuple((:v1, :v1k_up, :v1k_down, :v2, :v2k_up, :v2k_down), rand(6))
	test_df = clean_variables(nt, :v, :k, 2)
	@test names(test_df) == string.([:v, :vk_up, :vk_down])
	@test size(test_df, 1) == 2
end

# ╔═╡ 002849f9-dba0-474b-b8bf-9a79627e56a8
let
	fig = Figure()
	ax = Axis(fig[1,1])
	
	map(collect(out3.zero)) do (k, v)
		lines!(ax, agrid3, v, label = string(k))
	end
	axislegend(ax)

	fig
end

# ╔═╡ 9fbffaa8-027d-4766-93df-8f2051f3e2cb
df3 = mapreduce((x,y) -> leftjoin(x, y, on = [:a, :i_z]), [:v => r"^v", :c => r"^c", :s => r"^s"]) do (v, rgx)
	@chain out3.optional begin
		DataFrame
		@transform(:a = @c agrid3)
		DataFrames.select(rgx, :a)
		stack(rgx, value_name = v)
		@transform(:i_z = parse(Int, :variable[end]))
		DataFrames.select(Not(:variable))
	end
end

# ╔═╡ c60a4fe6-e6a5-421f-a1f6-d05297f404b9
@chain df3 begin
	stack(Not([:a, :i_z]))
	data(_) * mapping(:a, :value, color = :i_z => nonnumeric, layout = :variable) * visual(Lines)
	AlgebraOfGraphics.draw(facet = (linkyaxes = false, ))
end

# ╔═╡ f2b06e57-7414-47e4-b0a3-ac94b740f850
md"""
# Appendix
"""

# ╔═╡ 095027d6-cf28-431b-84a9-447f750c1676
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
NamedTupleTools = "d9ec5142-1e00-5aa0-9d6a-321866360f50"

[compat]
NamedTupleTools = "~0.14.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "e69b688917c852f40a8f6e0753c1e4e03e4ba09d"

[[deps.NamedTupleTools]]
git-tree-sha1 = "90914795fc59df44120fe3fff6742bb0d7adb1d0"
uuid = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
version = "0.14.3"
"""

# ╔═╡ Cell order:
# ╟─5d9a1abe-a51d-41c6-b86d-a1360140c05a
# ╟─9e4ec644-4566-4066-a958-ba7255fbe2f5
# ╟─36a22f3b-9f7d-4075-95d3-5787e40dbad7
# ╠═8a3e8d6a-23f4-41ba-93bb-235d10bdd534
# ╠═bdbf11b8-e6c1-4f9c-8ab4-cb714655fed0
# ╠═0010335f-1657-4e2e-98f9-60d917f98f8c
# ╠═23084da8-6d9d-40c4-a384-3d4066de62e5
# ╠═9dd13311-3e5e-4f4b-998b-4b6f2bdd2bbb
# ╠═212cca0b-314a-49df-9c20-057d70bfe277
# ╠═f4da2066-e7fd-4641-96ae-22b2708e0399
# ╠═0708f272-da24-401e-b3c9-a212dd003a38
# ╟─ae84468f-ee75-4763-812c-79106747908a
# ╠═920b0920-4b52-44e0-8779-ead83552d9cf
# ╠═2e3dac19-8b48-44f8-8c79-c8ef6d4fb54a
# ╠═7c8e2c71-47a5-49b6-811c-551a0f7ad756
# ╠═8cc86d11-e206-49ed-972f-583731d43bf1
# ╟─a719e2d7-1666-4d59-9c1e-43475e908e14
# ╠═eb9589e6-7310-416e-bd7f-d1a42880eed5
# ╠═63697439-bc45-4ef5-bdc7-9691b15f1332
# ╠═3e4e0e34-8db8-4e0f-a593-5ac822bb6c36
# ╟─8ebbc2fd-6dc1-401e-a002-b2d7ad71dc31
# ╠═55ab0920-7ece-44e0-b937-7ee0e9d7cc3d
# ╠═8db3df4c-f2a2-40b7-82ee-91b9dc609e37
# ╟─22a8ac01-2af9-4b1f-a1c7-7a663c670828
# ╟─e2d7380e-19bf-445a-96ef-36613a42943e
# ╠═5347b4ff-543a-42cf-8aea-a173e9b66f4d
# ╠═5c378bf5-f5d6-4e25-a1dc-3e5951fca01b
# ╠═2d563334-4032-45ce-b90b-65fe9b59efda
# ╠═296a8245-a813-484c-86b8-506791d55ac7
# ╠═fc8b67ef-faf4-4e8a-ad62-78c018d01747
# ╠═73a5630f-e4ad-402c-85e4-409c1092a5e8
# ╠═1f246a97-77b0-477f-b97d-7c0280f78310
# ╟─103b4796-162b-4ac4-b34c-42fd584efb5c
# ╟─900cf795-953b-4437-970e-944b05c279d7
# ╠═381fdb71-38af-4570-b6e1-0bfe72c7042d
# ╟─f4547de5-9254-4556-96ed-03e2221db3e8
# ╠═ba08a29a-7480-4d2c-85ee-83ed67f3f5f5
# ╠═40adb798-988e-4905-91c6-b0a14ceecacb
# ╠═0bbf0fcf-9671-4535-8006-b345d3a43c69
# ╠═113206ad-361e-4214-9208-99fef39ea96a
# ╠═65ce5940-3f7c-4036-bf51-bc90e2f8b235
# ╠═fb9684f8-31e6-4da9-94fb-8cee9c9fe69e
# ╠═b22052d8-e122-46d9-a00a-30b1d7d7555e
# ╠═2f577cd6-2f4d-4998-99ae-dbf22d7a5a65
# ╠═6ee07c2c-de88-4d38-9c4e-afb3e60b4e66
# ╠═786c5de0-2639-47f8-8c12-a89813d57075
# ╠═ae2a8f97-1d31-4b54-a0d3-31bfb87e4607
# ╠═002849f9-dba0-474b-b8bf-9a79627e56a8
# ╠═9fbffaa8-027d-4766-93df-8f2051f3e2cb
# ╠═c60a4fe6-e6a5-421f-a1f6-d05297f404b9
# ╟─f2b06e57-7414-47e4-b0a3-ac94b740f850
# ╠═095027d6-cf28-431b-84a9-447f750c1676
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ f14f868b-d6ad-43b5-82f6-5eba76fc344b
using DINA

# ╔═╡ d8cfaacb-b3de-41de-a2d1-6d1eeabaa2c7
using StatsBase: weights

# ╔═╡ 90f6bf62-752d-4018-beec-fa2239fc2db1
using LinearAlgebra: I, dot

# ╔═╡ 16f6d020-561c-43ea-bd1f-606890b7e009
using Roots

# ╔═╡ fa55d156-c645-4b16-b253-05f7802cfc36
using Optimization

# ╔═╡ fda2e873-d974-42f3-97fe-94223ab0b077
using OptimizationMultistartOptimization

# ╔═╡ 4644fac4-a6ed-48f4-b822-40f1df4977f8
using OptimizationNLopt

# ╔═╡ 5be2de51-2391-457a-8b94-c07e364d3eef
using Sobol

# ╔═╡ 8ca25b51-89fa-47af-b2a0-e109e9b0f98a
using Statistics

# ╔═╡ 43174c3f-e313-46c0-ba79-33b5ace09c21
using DataFrames, DataFrameMacros, Chain

# ╔═╡ ce061fb9-1678-47c9-a13e-90cf90d0338c
using AlgebraOfGraphics, CairoMakie

# ╔═╡ 6819f7ad-e0b2-4ef4-be9e-665629bbc127
using HypertextLiteral

# ╔═╡ ab5b369e-5834-4154-ac75-fa10cc1d0d8b
using Symbolics

# ╔═╡ d2590b30-c26a-4502-a374-1e388894abc3
using CSV: CSV

# ╔═╡ 6a218989-d242-4318-9ceb-fceddbfeaabc
using PlutoUI: TableOfContents, Slider, Select

# ╔═╡ 00690911-2b64-4b26-aa58-23fc3e23bb48
using StatFiles

# ╔═╡ d58962fa-4073-467c-b3ec-418e4ecdce42
begin
	using StatsBase: StatsBase, AbstractWeights
	using CategoricalArrays
	formatter(from, to, i; kwargs...) = "$i"
	function CategoricalArrays.cut(x::AbstractArray, w::AbstractWeights,
		q; labels=formatter)
  		cut(x, StatsBase.quantile(x, w, q); extend=true, labels)
	end

	function CategoricalArrays.cut(x::AbstractArray, w::AbstractWeights,
		ngroups::Integer; labels=formatter)
  		cut(x, StatsBase.quantile(x, w, (1:ngroups-1)/ngroups); extend=true, labels)
	end
end

# ╔═╡ 9e825451-37d8-4812-86c4-609e9296a179
md"""
## Targets from DINA
"""

# ╔═╡ 15c84e2b-e828-4320-9140-09a4a6c88968
dina_80 = get_dina(1980) |> DataFrame

# ╔═╡ e5d0b7c5-9d96-4702-b7dd-70f9c9e69ef2
dina_80_groups = @chain dina_80 begin
	disallowmissing!
	@transform(:is_owner = :ownerhome > 0)
	@subset(:is_owner == true, :id > 0)
	@transform(:income_bin = @c cut(
		:peinc, weights(:dweght), [0.0, 0.5, 0.9, 1.0],
		labels = (x, y, i; kwargs...) -> ["bottom 50", "middle 40", "top 10"][i]
	))
	stack([:peinc, :ownerhome, :ownermort], [:income_bin, :dweght])
	@groupby(:variable, :income_bin)
	@combine(:value = mean(:value, weights(:dweght)), :weight = sum(:dweght))
	@groupby(:variable)
	@transform(:weight = @c :weight / sum(:weight))
	unstack(:variable, :value)
	@select(
		:income_bin,
		:weight,
		:income = @c(:peinc ./ sum(:peinc, weights(:weight))),
		:d2y = -:ownermort / :peinc, 
		:ph2y = :ownerhome / :peinc,
		:d2ph = -:ownermort / :ownerhome
	)
	@transform(:income_share = :income * :weight)
end

# ╔═╡ 673d0937-a670-45a0-a7a5-6b9d9d270610
dina_targets(; avg_sensitivity) = @chain dina_80 begin
	disallowmissing!
	@transform(:is_owner = :ownerhome > 0)
	@subset(:is_owner == true, :id > 0)
	stack([:peinc, :ownerhome, :ownermort], [:dweght])
	@groupby(:variable)
	@combine(:value = mean(:value, weights(:dweght)), :weight = sum(:dweght))
	@groupby(:variable)
	@transform(:weight = @c :weight / sum(:weight))
	unstack(:variable, :value)
	@select(
		:ph2y = :ownerhome / :peinc,
		:d2y = -:ownermort / :peinc, 
		:d2ph = -:ownermort / :ownerhome
	)
	only
	(; _..., avg_sensitivity, hx2y = 0.3, kind = :data)
end

# ╔═╡ fb00431c-4eee-4041-bba6-504cd00693d6
md"""
## Calibration targets
"""

# ╔═╡ 1d6866f4-af9c-444f-ba93-65b5900b9ee1
md"""
* ``\bar L``: The value of new land permits ``\bar L`` is set so that employment in the construction sector is 5% of total employment, consistent with Bureau of Labor Statistics data for 1998.
* ``\xi``
* ``\frac{1}{1-\varepsilon} \in [0.15, 1.25]``: micro vs macro estimates
* ``\frac{\alpha}{1-\alpha} = 1.5`` median value across MSAs in Saiz (2010)
* ``
* ``\gamma``: 1.5
* ``\delta``: 0.021
* ``\rho``
* ``\phi``
"""

# ╔═╡ 0b1a51ea-073f-11ed-1d50-b1e334a7bcb0
md"""
# Falling Behind - Tractable Model
"""

# ╔═╡ e8d0027f-d4e3-41eb-9181-821e3c586a22
md"""
## Solving the model
"""

# ╔═╡ cd60add2-cb80-431c-8119-def8a1b961ef
Base.@kwdef struct TractableModel
	ρ = 0.01 # discount rate
	r = ρ
	ξ = 0.3 # utility weight of consumption
	ela = 0.15
	ε = 1 - 1/ela
	δ = 0.022 #1
	ϕ = 0.7 # strength of the comparison motive
	G = [1//3  1//3  1//3;
		 0     1//2  1//2;
		 0        0  1//1;]
	p = 1.0
	HS_ela = 1.5
	α = HS_ela / (1 + HS_ela)
	L̄ = 0.5 # check this!
	groups = ["bottom 50", "middle 40", "top 10"]
	group_weights = weights([0.5, 0.4, 0.1])
end

# ╔═╡ 1555e6b0-2fbe-466c-884d-186f5b20e6e3
md"""
## Calibration
"""

# ╔═╡ dd95cda0-620e-4fee-8d35-e58ce9aebfd6
md"""
```math
\text{sensitivity} = - \frac{\text{ela}_{\tilde h}}{ela_h} = - \frac{\frac{\partial u}{\partial s}\frac{\partial s}{\partial \tilde h}\frac{\tilde h}{u}}{\frac{\partial u}{\partial s}\frac{\partial s}{\partial h}\frac{h}{u}} =
\frac{\phi}{1} \frac{\tilde h}{h} \stackrel{!}{=} 0.7
```
"""

# ╔═╡ a6f923d1-0292-4168-a40a-27929be0815c
bellet_sensitivity(h, h̃, (; ϕ)) = ϕ * h̃/h

# ╔═╡ 37618646-be2e-474d-972f-e2e84b914765
function choices(p, (; r, δ, ξ, ε, ela, ϕ, G, group_weights, groups, HS_ela, α, L̄), y)
	κ₀ = ((r + δ) * (1-ξ)/ξ * p)^(1/(1-ε))
	κ₁ = 1 / (p * (r + δ)/κ₀ + 1)
	κ₂ = κ₁/κ₀
	κ₃ = 1 / (1 + p*r/(δ * p + κ₀))

	if ε == 0.0
		@assert κ₀ ≈ p * (r + δ) * (1-ξ)/ξ
		@assert κ₁ ≈ 1 / (ξ/ (1-ξ) + 1)
		@assert κ₁ ≈ 1-ξ
		@assert κ₂ ≈ ξ / (p * (r + δ))
		@assert κ₃ ≈ 1 / (1 + p*r/(δ * p + κ₀))
	end
	
	a₀ = zeros(size(y))
	𝒴 = r .* a₀ + y
	h = (I - κ₁ * ϕ * G) \ (κ₂ * 𝒴)
	h̃ = G * h
	debt = y - κ₃ * 𝒴 + (1-κ₃) * ((I - κ₁ * ϕ * G) \ 𝒴 - 𝒴)

	sensitivities = bellet_sensitivity.(h, h̃, Ref((; ϕ)))
	avg_sensitivity = mean(sensitivities, group_weights)
	
	∑𝒴    = sum(𝒴, group_weights)
	∑h    = sum(h, group_weights)
	∑debt = sum(debt, group_weights)

	Iₕ_demand = ∑h * δ
	@assert HS_ela ≈ α/(1-α)
	Iₕ_supply = (α * p)^HS_ela * L̄
	ζₕ = Iₕ_demand - Iₕ_supply
	ζₕ_rel = ζₕ / maximum(abs, [Iₕ_demand, Iₕ_supply])

	group_tbl = (; h, h̃, debt, 𝒴, groups, group_weights, sensitivities)
	
	(; p, Iₕ_demand, Iₕ_supply, ζₕ, ζₕ_rel, ∑h, ∑debt,
		d2y = ∑debt / ∑𝒴, d2ph = ∑debt / (p * ∑h), ph2y = p * ∑h / ∑𝒴,
		hx2y = δ * p * ∑h + r * ∑debt,
		ϕ, ela, ε, group_tbl, avg_sensitivity
	)
end

# ╔═╡ e5abb6c6-7add-4f73-9ffd-26173c92b5a7
function market_price(par, 𝒴; price_bracket = (eps(), 100.0))
	p = find_zero(p -> choices(p, par, 𝒴).ζₕ_rel, price_bracket)
end

# ╔═╡ eea9988e-3fcf-48d8-93c3-67d94e1ab949
function general_equilibrium(par, 𝒴; kwargs...)
	p = market_price(par, 𝒴; kwargs...)
	choices(p, par, 𝒴)
end

# ╔═╡ 7b4f2179-15fe-4651-94d0-1765c6ed8b83
md"""
#### Specifying the loss function
"""

# ╔═╡ 6c8cf040-c64e-4a03-b42c-b1df688e3572
# using NamedTupleTools

# ╔═╡ 6ac092f8-8f03-4a1f-91f3-1177efd99745
model_statistics(out) = (;
	out.ph2y, out.d2y, out.d2ph, out.avg_sensitivity, out.hx2y,
	kind = :model
)

# ╔═╡ 279101cf-0d99-4446-b25a-b7e2dd7b75f4
md"""
### Trying a few Sobol numbers
"""

# ╔═╡ 2816d37a-7ea1-4c0a-9149-373c9a4361a6
md"""
### Trying `TikTak` with different local solvers
"""

# ╔═╡ 043624c1-0ee5-48db-80d2-82331ee5552e
elas = [0.5, 1.0, 1.25]

# ╔═╡ b7ff7f0f-a9e3-46c2-bf39-84c6ef67d272
𝒴 = dina_80_groups.income

# ╔═╡ 8174336f-bb98-4fb5-a275-53d16e4e1ab3
target_weights = (; 
	ph2y = 1//3,
	d2y = 1//3,
	d2ph = 0,
	avg_sensitivity = 1//3,
	hx2y = 0, #1//3,
	kind = :weights
)

# ╔═╡ 6e308459-f126-4742-8e6f-33a6420bda19
local_solver = NLopt.LN_BOBYQA()

# ╔═╡ e5d106d5-069b-464d-8555-8bd546a0a042
md"""
## Calibration targets
"""

# ╔═╡ a2d3f53b-2c15-4b17-83e9-25f6b975ee65
md"""
* ``ξ``: $(@bind _ξ_ Slider(0.0:0.01:1.0, default = 0.3, show_value = true))
* ``ϕ``: $(@bind _ϕ_ Slider(0.0:0.01:1.0, default = 0.2, show_value = true))
* ``ρ``: $(@bind _ρ_ Slider(0.0:0.001:0.2, default = 0.01, show_value = true))
"""

# ╔═╡ e0135a19-bf81-46c0-8b25-15e320a39e78
par = TractableModel(; ϕ = _ϕ_, ela = 1.0, ξ=_ξ_, ρ = _ρ_)

# ╔═╡ a6d14acc-8b2a-4736-a774-a03174f3ca42
out = general_equilibrium(par, 𝒴)

# ╔═╡ 5f8a5a6c-a628-46a6-a90c-8df775643e9c
bounds = [(eps(), 1-eps()), (0, 0.7), (eps(), 0.20)]

# ╔═╡ b4852dd9-97f5-402f-916e-e387e52c0694
moments(out, targets, target_weights) = let
	df = DataFrame([
		model_statistics(out), targets, target_weights
	])

	@chain df begin
		stack(Not(:kind), variable_name = :target)
		unstack(:kind, :value)
		disallowmissing!
		@transform(:abs_pc_dev = abs((:model - :data)/:data))
	end
	
end

# ╔═╡ 2eb5df98-f294-481d-b0b3-947f020b0d47
deviation(out, targets, target_weights) = @combine(moments(out, targets, target_weights), :loss = mean(:abs_pc_dev, weights(:weights))).loss |> only

# ╔═╡ 6d3a317e-8e14-457b-99a5-90baf85bfbaf
function loss(targets, 𝒴, other_params=(;))
	function (x, p; return_details=false, append = (;))
		named_x = NamedTuple{(:ξ, :ϕ, :ρ)}(x)
		#named_p = NamedTuple{(:ela,)}(p)
	
		model = TractableModel(; named_x..., #=named_p...,=# other_params...)
		out = general_equilibrium(model, 𝒴)
	
		#targets = dina_targets(; avg_sensitivity=0.7)
		loss = abs(deviation(out, targets, target_weights))
	
		if return_details
			
			return (; loss, named_x..., #=named_p..., =# other_params, append..., model_statistics(out)..., out.p)
		else
			return loss
		end
	end	
end

# ╔═╡ 1146e3f5-d701-42fa-929f-4699efb5a5f6
function calibrate(other_params, bounds, targets, target_weights, 𝒴; local_solver = NLopt.LN_NELDERMEAD())
	lb = first.(bounds)
	ub = last.(bounds)
	x0 = lb .+ ub ./ 2

	ℓ = loss(targets, 𝒴, other_params)
	f = OptimizationFunction(ℓ)
	prob = Optimization.OptimizationProblem(f, x0, []; lb, ub)

	sol = solve(prob, MultistartOptimization.TikTak(100), local_solver)
	ℓ(sol.u, prob.p, return_details=true, append = (; sol.retcode))
end

# ╔═╡ 7d6088f8-2ffe-4238-9cb9-c0b7e765cc4d


# ╔═╡ 082e201e-6732-432f-9323-162c84b68a25
md"""
## Analysis
"""

# ╔═╡ 62d8a4e2-1a24-48c8-b044-6c25147dcf51
df = let
	ϕ = 0.2
	pars = [
		"micro" => TractableModel(; ϕ = 0.0),
		"micro" => TractableModel(; ϕ),
		"CD" => TractableModel(; ϕ = 0.0, ela = 1.0),
		"CD" => TractableModel(; ϕ, ela = 1.0),
		"macro" => TractableModel(; ϕ = 0.0, ela = 1.25),
		"macro" => TractableModel(; ϕ, ela = 1.25)
	]


	𝒴₀ = [0.5, 1.0, 2.5]
	𝒴₁ = copy(𝒴₀)
	𝒴₁[3] *= 2.0

	𝒴s = ["1980" => 𝒴₀, "2007" => 𝒴₁]

	mapreduce(vcat, pars) do (label, par)
		# initial equilibrium
		p₈₀ = market_price(par, 𝒴₀)
		out_80 = merge(choices(p₈₀, par, 𝒴₀), (; label, time = "early", 𝒴_label="1980"))
		out_PE = merge(choices(p₈₀, par, 𝒴₁), (; label, time = "late", 𝒴_label="2007 PE"))

		p₀₇ = market_price(par, 𝒴₁)
		out_GE = merge(choices(p₀₇, par, 𝒴₁), (; label, time = "late", 𝒴_label="2007 GE"))

		[out_80, out_PE, out_GE]
	end |> DataFrame
end
	

# ╔═╡ d0925b77-4f4f-4a97-a6bf-982b318c1f2a
Makie.current_default_theme().axis

# ╔═╡ 616bc951-ad6c-4475-8394-3104a7d336b8
font_theme = Theme(
	font = "CMU",
	Axis = (titlefont = "CMU", ),
	Legend = (framevisible = false, position = :top, titleposition = :left, )
)

# ╔═╡ 2db97ad1-7662-48f6-b546-ff11b267013b
@chain df begin
	@subset(:ela ≈ 1.0)
	select(:label, :ϕ, :𝒴_label, :group_tbl => AsTable)
	flatten([:h, :debt, :𝒴, :groups, :group_weights])
	stack([:h, :debt, :𝒴])
	#data(_) * mapping(:𝒴_label, :value, row = :variable, col = :label, stack = :groups, color = :groups) * visual(BarPlot)
	#draw(facet = (linkyaxes = false, ))
	@groupby(:𝒴_label, :variable, :label, :ϕ)
	@transform(:total = @c mean(:value, weights(:group_weights)))
	@groupby(:label, :ϕ, :variable, :groups)
	@transform(
		:Δ = @c(:value .- first(:value)),
		:total_Δ = @c(:total .- first(:total)),
		:total_growth = @c(:total ./ first(:total) .- 1)
	)
	@transform(:wtd_Δ = :Δ * :group_weights)
	@subset(:𝒴_label != "1980")
	@transform(:xxx = :wtd_Δ / :total_Δ * :total_growth)
	data(_) * mapping(:𝒴_label, :xxx, row = :variable, col = :ϕ => nonnumeric, stack = :groups, color = :groups) * visual(BarPlot)
	draw(facet = (linkyaxes = true, ))
	#		:total = fill(total, size(:value))
	#		:share = :group_weights .* :value ./ total
	#	end
	#)
end

# ╔═╡ 0d6b97b7-a1ef-4cd5-9868-4d9ea74591e4
vars = [:p, :∑debt, :∑h, :d2y, :d2ph, :ph2y]

# ╔═╡ f8888d57-2e06-4dcc-b7b6-db830840d426
@chain df begin
	select(vars..., :label, :𝒴_label)
	stack(vars)
	@groupby(:label, :variable)
	@transform(:value = @c :value ./ first(:value) .- 1)
	@subset(:𝒴_label != "1980")
	data(_) * mapping(:label => "year", :value, layout = :variable, dodge = :𝒴_label, color = :𝒴_label) * visual(BarPlot)
	draw()#(facet = (linkyaxes = false, ))
end

# ╔═╡ b29b41b6-cc32-46a8-8c0c-dd6ccc990f71
md"""
# Tractable model -- Enhanced version
"""

# ╔═╡ b4c18055-3168-474d-8dfe-6ea5105e0e80
md"""
## Partial Equilibrium
"""

# ╔═╡ fcb8db73-9f24-4814-98c4-b08049d8029a
md"""
## General equilibrium effects (on the housing market)
"""

# ╔═╡ 04153526-eacc-4a6f-a175-d7f9ed44ce7e
md"""
### House prices: Closed forms with Cobb-Douglas
"""

# ╔═╡ b4aa9303-f440-4cea-82b4-943a3c14bfac
md"""
### House prices: Numerical analysis with general CES
"""

# ╔═╡ fcabed90-a37b-4f1d-aea9-4914a4f95a1c
vars2 = [:ϕ, :ε, :ela]

# ╔═╡ bd7a67c1-e465-41c4-8cb3-dfa58e79386b
vars3 = [:∑debt => "aggregate debt", :p => "house price", :d2y => "debt to income", :d2ph => "d2ph"]

# ╔═╡ bb63db54-1065-46e7-989c-50d080230938
@chain df begin
	select(vars3..., :𝒴_label, vars2..., :time)
	@transform(:𝒴_label = :𝒴_label[6:end])
	#@subset(:ela == 1.0)
	stack(last.(vars3))
	@aside early = @chain _ begin
		@subset(:time == "early")
		select(Not([:time, :𝒴_label]))
		rename(:value => :early)
	end
	@subset(:time == "late")
	rename(:value => :late)
	leftjoin!(early, on = [vars2; :variable])
	@transform(:change = :late - :early)
	@transform(:pc_change = :change / :early)
	@subset(:𝒴_label == "GE")
	data(_) * mapping(
		:ela => nonnumeric => L"elasticity of substitution $\frac{1}{1-\varepsilon}$",
		:pc_change => "percentage change",
		layout = :variable,
		color = :ϕ => nonnumeric => L"comparison motive $\phi$",
		dodge = :ϕ => nonnumeric
	) * visual(BarPlot)
	with_theme(font_theme) do 
		draw(_, 
			legend = (; position = :top),
		)
	end
end

# ╔═╡ 5fc360fe-62d3-446e-91d4-88bb5bdf3d9f
md"""
#### Analysis by income type
"""

# ╔═╡ a128a97c-0192-4cb7-8021-e5abb9f152da
md"""
### Mortgage debt: Closed form with Cobb-Douglas
"""

# ╔═╡ aec3912e-5847-4a46-a658-bce97258c54c
md"""
### Mortgage debt: Numerical analysis with general CES
"""

# ╔═╡ 41dd7623-d71e-46ca-8580-c44c2acca2cb
details_prices = md"""
## Housing demand and house prices

When we assume Cobb-Douglas flow utility function (``\varepsilon = 0``), we can solve for the equilibrium price in closed form.

The ``\kappa_i`` parameters simplify to
```math
	κ₁ = 1-ξ, \qquad κ₂ = \frac{ξ}{p \cdot (r + δ)}.
```
This allows us to write housing demand
```math
\begin{align}
h(p) &= (I - (1-ξ)ϕG)^{-1} \mathcal{Y} \frac{ξ}{p \cdot (r + δ)} \\
\implies H &= \sum_i ωᵢ hᵢ(p) = \underbrace{\omega^T (I - (1-ξ)ϕG)^{-1} \mathcal{Y}}_{=: θ \in ℝ} \frac{ξ}{p \cdot (r + δ)} \\
&=  \frac{θ ξ}{p (r + δ)}.
\end{align}
```
From the optimality condition of the construction sector we get
```math
I_h = (α p)^\frac{α}{1-α} \bar L.
```
The market clearing condition is ``I_h = δ H``. From there we can derive the equilibrium price
```math
\begin{align}
(α p)^\frac{α}{1-α} \bar L &= δ \frac{θ ξ}{p (r + δ)} \\
α^\frac{α}{1-α} p^\frac{1}{1-α} &= \frac{δ θ ξ}{\bar L (r + δ)}  \\
p &= α^{-α} \Bigl(\frac{δ θ ξ}{\bar L (r + δ)}\Big)^{1-\alpha}.
\end{align}
```

Let's look at the ``\theta`` term.
```math
\begin{align}
\theta &= \omega^T \Bigl(\sum_{i=0}^\infty ((1-ξ)ϕG)^i \Bigr) \mathcal{Y} \\
 &= \omega^T \Bigl(I + \sum_{i=1}^\infty ((1-ξ)ϕG)^i \Bigr) \mathcal{Y}
\end{align}
```
We can see that rising top incomes drive up house prices even in the absence of social comparison motive. This is because the price reacts to the aggregate housing demand, which is a function of total aggregate income.
""";

# ╔═╡ 3895355e-c93c-4225-89e3-3c907b6d11de
@htl("""
<details> <summary> Proof </summary>

$(details_prices)

</details>
""")

# ╔═╡ 73f900a1-c1da-4efa-847a-48c5e8cb0717
md"""
## Debt
"""

# ╔═╡ fe128ce0-87a3-4c13-b596-0e2c34650058
function price((; ξ, r, δ, L̄, α, ϕ, G, group_weights), y)
	a₀ = zeros(size(y))
	𝒴 = r .* a₀ + y

	#κ₀(p) = p * (r + δ) * (1-ξ)/ξ
	κ₁ = 1-ξ
	κ₂(p) = ξ / (p * (r + δ))
	
	θ = dot(group_weights, (I - ϕ * κ₁ * G) \ 𝒴)
	@assert θ ≈ dot(group_weights, inv(I - ϕ * (1-ξ) * G), 𝒴)
	
	p = (δ * θ * ξ / ((r+δ)*L̄))^(1-α) * α^(-α)
	


	∑h = θ * ξ / (r+δ) / p
	let
		h = (I - κ₁ * ϕ * G) \ (κ₂(p) * 𝒴)
		@assert ∑h ≈ sum(h, group_weights)
	end

	(; p, ∑h)
end

# ╔═╡ 468ca176-8079-423e-bad4-465d12d1b78a
let
	par = TractableModel(ϕ = 0.7, ela = 0.0)

	y = [0.5, 1.0, 2.5]
	CD = price(par, y)

	out = choices(CD.p, par, y)
	@info (; out.ζₕ_rel, out.ζₕ)
end

# ╔═╡ 9eb5034f-af9b-4c46-b480-8b03c977757c
md"""
## CRRA
"""

# ╔═╡ dad754a5-4b2f-4b3e-9910-7a7caea711a9
@variables δ ξ p r ε

# ╔═╡ e6ffda4f-a1a0-42ad-bb42-e755076b1fb9
κ₀ = ((r + δ) * (1-ξ)/ξ)^(1/(1-ε)) * p^(1/(1-ε)) # p^(1/(1-ε)) * (r + δ) * (1-ξ)/ξ

# ╔═╡ 4981e7b4-d43b-479a-bd35-752462cf84d2
κ₀_alt = ((r + δ) * (1-ξ)/ξ)^(1/(1-ε)) * p * p^(ε/(1-ε))

# ╔═╡ 2bda56a6-e49f-45c6-a814-016ce5a4cc4c
κ₀ / κ₀_alt |> simplify

# ╔═╡ dd774ed7-6928-45c9-bede-ecfeb3884335
κ₁ = 1 / ((r+δ)* p / κ₀ + 1) #|> simplify

# ╔═╡ c3d5e417-f8e0-40ff-8fcf-a8b4111f4ec3
κ₂ = κ₁ / κ₀ |> simplify

# ╔═╡ 054a70a2-5f63-4742-8782-102f229d47ff
κ₃ = (δ * p + κ₀) * κ₂ |> simplify

# ╔═╡ 17b1c24a-7337-4917-9a0f-5c2defa1bf96
κ₃_num = p^(1/(1-ε)) * ((r + δ) * (1-ξ)/ξ)^(1/(1-ε)) + p * δ

# ╔═╡ 0860e535-33f1-4e22-aa46-eba71dd95769
κ₃_den = p^(1/(1-ε)) * ((r + δ) * (1-ξ)/ξ)^(1/(1-ε)) + p * (r + δ)

# ╔═╡ b9fee53f-9e16-484b-9514-ba821ce8d85e
κ₃ / κ₃_num * κ₃_den |> simplify

# ╔═╡ 65e3be85-d508-4817-a632-4b32d75bffcf
@variables G[1:3, 1:3]

# ╔═╡ fb0e25b0-8066-4370-8f5b-527a35091017
κ₁ * G

# ╔═╡ 900491e1-1962-4543-ab5c-04ee8b95bc5c
h = (I - κ₁ * ϕ * G) \ (κ₂ * 𝒴)

# ╔═╡ bb9fa247-42c8-4ab4-ae36-2892fa5f1651
md"""
# Appendix
"""

# ╔═╡ 293fa922-ee60-4b7d-8085-8d41d5d2363a
function theorem_header(type, number, title)
	out = type
	if !isnothing(number) && length(number) > 0
		out = out * " " * string(number)
	end
	if !isnothing(title) && length(title) > 0
		out = out * ": " * string(title)
	end
	out
end

# ╔═╡ ec31d541-2087-4cd5-93c5-6a718d61ce4b
begin
	admonition(kind, title, text) = Markdown.MD(Markdown.Admonition(kind, title, [text]))
	proposition(text; number = nothing, title = nothing) = 
		admonition("correct", theorem_header("Proposition", number, title), text)
	corollary(text; number = nothing, title = nothing) =
		admonition("note", theorem_header("Corollary", number, title), text)
#	danger(text, title="Danger")   = admonition("danger",  title, text)
#	correct(text, title="Correct") = admonition("hint", title, text)

end

# ╔═╡ 89960fdf-82f3-49df-855e-35d21fb8b24d
proposition(md"Closed forms for ``h`` and ``-ra`` in partial equilibrium. (``\kappa_i`` are functions of ``p``.)", number=1, title="Closed forms")

# ╔═╡ b3ed4900-2279-49af-891d-0e52a29c4b08
proposition(md"``h_i`` is increasing and ``a_i`` is decreasing in ``\mathcal{Y}_j`` whenever ``i`` is linked (_cares about_) ``j``, directly or indirectly.", number=2, title="Housing, debt and others' incomes")

# ╔═╡ 730a6ef9-619e-4008-b03a-120cb5313850
proposition(md"Impact of change in ``j``'s income on total housing and debt is proportional to its popularity.", number=3, title="Popularity and total impact")

# ╔═╡ 9d96e880-41f7-41e4-a1c6-3e39f5cbc398
corollary(md"When houses are non-durable, debt does not depend on others' incomes.", number="3 → 2", title="Effect on debt depends on durability")

# ╔═╡ 8d76aad5-59b6-4ffe-83da-ad80bcb9a893
proposition(md"""
Under Cobb-Douglas aggregation, the equilibrium house price is
```math
p = α^{-α} \Bigl(\frac{δ θ ξ}{\bar L (r + δ)}\Big)^{1-\alpha}
```
where
```math
\theta = \omega^T \Bigl(I + \sum_{i=1}^\infty ((1-ξ)ϕG)^i \Bigr) \mathcal{Y}.
```
That is, house prices rise even if there are no social comparisons. But house prices rise _more strongly with social comparisons_.
""",
	number=4, title="Decomposing rising house prices (Cobb-Douglas)"
)

# ╔═╡ 787ca466-6d72-45f1-b6fc-878ddd2ab75b
corollary(md"Under Cobb-Douglas: Effects on debt independent on house prices", number="2 → 3", title="Cobb-Douglas, debt and GE")

# ╔═╡ 30ac4eda-2829-4879-854a-36a08acebb36
md"""
## Packages
"""

# ╔═╡ 89668f54-762b-4cf3-bf7c-44f43b796b29
TableOfContents()

# ╔═╡ a2c782a8-aa86-49f9-9b96-2ced327a876f
md"""
## DataDeps
"""

# ╔═╡ 4ae24902-ddb0-4abf-b2e5-ecbad83db5a0
md"""
## SCF
"""

# ╔═╡ e6ebdc70-fcb2-4d6e-874b-951671480e33
const SCF_YEARS = 1989:3:2019

# ╔═╡ 6821ead6-0cd5-4d92-8fab-0f3356de5637
scf_checksums = Dict(
	1989 => "3600d39fa908f2b6b32a63518fe38d2b7ca7d9c8047e947b054ea2e7d36f7b9e",
	1992 => "d86da7dc07819adb1f08683d7b04f782d2fb1c135ea94880ce4b6550db3c0ccd",
	1995 => "c28169ca73855a1b1f22a999655ea573fc5391b11884b9272c9cf2bf1ee5c442",
	1998 => "69fda43abc88df203f03b9a9e8cf5bf5bc0a626da4b99d44f2fd11dee2c0b11e",
	2001 => "28f5548f91d5f851ad643d9c172e00ecfeb1d6fe47126fb2496af5b980f75ffd",
	2004 => "bb08a6122a25348f6507fcd9377511820972108e8642d6ec1ab3ddf262c21071",
	2007 => "85324789b2ab6f5e5dfc05a8de294e4fd837f7e2174f41e5be5e146435a68aa9",
	2010 => "a85ee57748ec28b3366a4f0b9446ec8b0c34710f14097b6a6e03089c9ad8823a",
	2013 => "f13ed12756798c7e696dbf56ec26438cc2e0c46de4c3343afd3fd05ddfb9e6e8",
	2016 => "11e92c267f333fe10678c9cbb9752c10290085c35abd7c52f7f21c8df45dc468",
	2019 => "87766da9024f7b6742d277c955234cf6bee439248fbc5f689c81d91880fd1b05",
)

# ╔═╡ e303a91c-33a8-41d1-96b4-8a213dab5698
begin
	using DataDeps; ENV["DATADEPS_ALWAYS_ACCEPT"] = true
	
	register(DataDep(
   		"MacroHistory",
   		"",
   			["http://data.macrohistory.net/JST/JSTdatasetR5.dta"],				
 			["e4691b8ac90e8d6947c65ae197841862ded27a00c167a578c38eccf4f73043ce"]
	))

	for year in SCF_YEARS
		register(DataDep(
    		"SCF$(year)",
    		"",
    		"https://www.federalreserve.gov/econres/files/scfp$(year)excel.zip",
   			scf_checksums[year];# [checksum::Union{String,Vector{String}...},]; # Optional, if not provided will generate
    		post_fetch_method=unpack
		))
	end

	
end

# ╔═╡ 4fc16124-5f72-4c1f-ab18-a99f6bcb027e
function get_scf(year)
	@assert year ∈ SCF_YEARS

	str = "SCF$(year)"
	path = @datadep_str str

	CSV.File(joinpath(path, "SCFP$(year).csv")) |> DataFrame
end

# ╔═╡ d0da815f-65ea-40d7-ba77-d93afcb9a889
data_raw = @chain get_scf(1989) begin
	@subset(:HOUSES > 0)
	stack([:INCOME, :NETWORTH, :ASSET, :DEBT, :NH_MORT, :HOUSES], [:WGT])
	@groupby(:variable)
	@combine(:value = mean(:value, weights(:WGT)))
	NamedTuple{tuple(Symbol.(_.variable)...)}(tuple(_.value...))
end

# ╔═╡ 26cb90ad-e3ce-4a9f-9d9b-36e43a103057
scf_targets(; avg_sensitivity) = (; 
	ph2y = data_raw.HOUSES / data_raw.INCOME, 
	d2y = data_raw.NH_MORT / data_raw.INCOME,
	d2ph = data_raw.NH_MORT / data_raw.HOUSES,
	avg_sensitivity, hx2y = 0.3,
	kind = :data
)

# ╔═╡ a44b1c72-3bb5-4ff7-9378-275dffba2358
scf_targets(avg_sensitivity=0.7)

# ╔═╡ e518738f-6c0b-450d-a119-b17f1f71e5f0
t = scf_targets(avg_sensitivity=0.7)

# ╔═╡ 7fda0a57-1148-4f0d-b913-c148a7e6236f
targets = (; t.ph2y, t.d2y, t.d2ph, t.avg_sensitivity, t.hx2y, t.kind)

# ╔═╡ 6f3e7b25-95f2-4f77-b087-14c6543117c8
let
	n₀ = 1
	ela = 0.7
	
	s = SobolSeq(first.(bounds), last.(bounds))
	s = skip(s, n₀-1, exact = true)
	
	map(enumerate(first(s, 1_000))) do (i, x)
		ℓ = loss(targets, 𝒴, (; ela))
		(; i = i + n₀ - 1, ℓ(x, p, return_details=true)...)
	end |> DataFrame |> x -> sort(x, :loss)
end

# ╔═╡ b052b26f-9ac2-49b0-a713-3a81169e7dcb
df_out = map(elas) do ela
	calibrate((; ela, δ = 0.05), bounds, targets, target_weights, 𝒴; local_solver)
end |> DataFrame

# ╔═╡ 309a39ef-a3a8-45b0-8891-a8cd0aebfe5a
calibrate((; ela=0.5), bounds, targets, target_weights, 𝒴; local_solver)

# ╔═╡ eb45f538-6384-475e-9de5-6bb6152a3cc1
calibrate((; ela=1.0), bounds, targets, target_weights, 𝒴; local_solver)

# ╔═╡ 6f14d6c3-ba2c-41f1-acbf-8125fec69d7c
calibrate((; ela=1.25), bounds, targets, target_weights, 𝒴; local_solver)

# ╔═╡ 121b5939-216b-4864-9845-30ec53989f56
md"""
## Macro History
"""

# ╔═╡ 7bbe5a83-6d2d-4339-a3dd-3af947ee2172
function get_macro_history()
	DataFrame(StatFiles.load(joinpath(datadep"MacroHistory", "JSTdatasetR5.dta")))
end

# ╔═╡ fa04eb87-3159-4505-80d4-300f91962ab8
@chain begin
	get_macro_history()
	@transform(:rhpi = :hpnom / :cpi)
	@subset(:year > 1950)
	@subset(!ismissing(:rhpi))
	data(_) * mapping(:year, :rhpi, color = :country) * visual(Lines)
	draw
end

# ╔═╡ ce84af63-f8fc-47e6-8e0b-a3ac18283534


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DINA = "2e9d5fb1-712f-4dff-9379-c18f12b3746d"
DataDeps = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Optimization = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
OptimizationMultistartOptimization = "e4316d97-8bbb-4fd3-a7d8-3851d2a72823"
OptimizationNLopt = "4e6fcdb7-1186-4e1f-a706-475e75c168bb"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
Sobol = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
StatFiles = "1463e38c-9381-5320-bcd4-4134955f093a"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"

[compat]
AlgebraOfGraphics = "~0.6.9"
CSV = "~0.10.4"
CairoMakie = "~0.8.9"
CategoricalArrays = "~0.10.6"
Chain = "~0.4.10"
DINA = "~0.1.3"
DataDeps = "~0.7.9"
DataFrameMacros = "~0.2.1"
DataFrames = "~1.3.4"
HypertextLiteral = "~0.9.4"
Optimization = "~3.8.2"
OptimizationMultistartOptimization = "~0.1.0"
OptimizationNLopt = "~0.1.0"
PlutoUI = "~0.7.39"
Roots = "~2.0.1"
Sobol = "~1.5.0"
StatFiles = "~0.8.0"
StatsBase = "~0.33.19"
Symbolics = "~4.10.3"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

[[deps.AbstractAlgebra]]
deps = ["GroupsCore", "InteractiveUtils", "LinearAlgebra", "MacroTools", "Markdown", "Random", "RandomExtensions", "SparseArrays", "Test"]
git-tree-sha1 = "774a46f9b06d1667a295e58e6e3cf0838f5dbb4d"
uuid = "c3fe647b-3220-5bb0-a1ea-a7954cac585d"
version = "0.27.2"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "5c0b629df8a5566a06f5fef5100b53ea56e465a0"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.2"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "78fba56df46874fa7a4f4376b324207c0719079d"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.10"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "0582b5976fc76523f77056e888e454f0f7732596"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.22"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "40debc9f72d0511e12d817c7ca06a721b6423ba3"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.17"

[[deps.ArrayInterfaceStaticArrays]]
deps = ["Adapt", "ArrayInterface", "ArrayInterfaceStaticArraysCore", "LinearAlgebra", "Static", "StaticArrays"]
git-tree-sha1 = "efb000a9f643f018d5154e56814e338b5746c560"
uuid = "b0d46f97-bff5-4637-a19a-dd75974142cd"
version = "0.1.4"

[[deps.ArrayInterfaceStaticArraysCore]]
deps = ["Adapt", "ArrayInterfaceCore", "LinearAlgebra", "StaticArraysCore"]
git-tree-sha1 = "a1e2cf6ced6505cbad2490532388683f1e88c3ed"
uuid = "dd5226c6-a4d4-4bc7-8575-46859f9c95b9"
version = "0.1.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AutoHashEquals]]
git-tree-sha1 = "45bb6705d93be619b81451bb2006b7ee5d4e4453"
uuid = "15f4f7f2-30c1-5605-9d31-71845cf9641f"
version = "0.2.0"

[[deps.Automa]]
deps = ["Printf", "ScanByte", "TranscodingStreams"]
git-tree-sha1 = "d50976f217489ce799e366d9561d56a98a30d7fe"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "0.8.2"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "Future", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables", "ZygoteRules"]
git-tree-sha1 = "b15a6bc52594f5e4a3b825858d1089618871bf9d"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.36"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bijections]]
git-tree-sha1 = "fe4f8c5ee7f76f2198d5c2a06d3961c249cce7bd"
uuid = "e2ed5e7c-b2de-5872-ae92-c73ca462fb04"
version = "0.1.4"

[[deps.BinaryProvider]]
deps = ["Libdl", "Logging", "SHA"]
git-tree-sha1 = "ecdec412a9abc8db54c0efc5548c64dfce072058"
uuid = "b99e7846-7c00-51b0-8f62-c81ae34c0232"
version = "0.5.10"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "Tables", "Unicode", "WeakRefStrings"]
git-tree-sha1 = "873fb188a4b9d76549b81465b1f75c82aaf59238"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.4"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA"]
git-tree-sha1 = "0483b660cfc7bff57e986c6d5af5179bcccdd8dc"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.12"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5f5a975d996026a8dd877c35fe26a7b8179c02ba"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.6"

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "2e62a725210ce3c3c2e1a3080190e7ca491f18d7"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.7.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonSolve]]
git-tree-sha1 = "332a332c97c7071600984b3c31d9067e1a4e6e25"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.1"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "59d00b3139a9de4eb961057eabb65ac6522be954"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.0"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DINA]]
deps = ["CategoricalArrays", "Chain", "DataDeps", "DataFrames", "LinearAlgebra", "ReadableRegex", "StatFiles", "Statistics", "StatsBase", "TableOperations"]
git-tree-sha1 = "8e8af8b3028642ebf499bdc0f1ea6d8591543d5c"
uuid = "2e9d5fb1-712f-4dff-9379-c18f12b3746d"
version = "0.1.3"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataDeps]]
deps = ["BinaryProvider", "HTTP", "Libdl", "Reexport", "SHA", "p7zip_jll"]
git-tree-sha1 = "8e2713f4c3394fba835c0babd0bfd2b3cd007fd9"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.9"

[[deps.DataFrameMacros]]
deps = ["DataFrames"]
git-tree-sha1 = "cff70817ef73acb9882b6c9b163914e19fad84a9"
uuid = "75880514-38bc-4a95-a458-c2aea5a3a702"
version = "0.2.1"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Reexport", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "daa21eb85147f72e41f6352a57fccea377e310a9"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.3.4"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.DataValues]]
deps = ["DataValueInterfaces", "Dates"]
git-tree-sha1 = "d88a19299eba280a6d062e135a43f00323ae70bf"
uuid = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"
version = "0.4.13"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random"]
git-tree-sha1 = "36bc84c68847edd2a3f97f32839fa484d1e1bce7"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.22"

[[deps.DiffResults]]
deps = ["StaticArrays"]
git-tree-sha1 = "c18e98cba888c6c25d1c3b048e4b3380ca956805"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.0.3"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "28d605d9a0ac17118fe2c5e9ce0fbb76c3ceb120"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.11.0"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "3258d0659f812acde79e8a74b11f17ac06d0ca04"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.7"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "aafa0665e3db0d3d0890cdc8191ea03dc279b042"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.66"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "StaticArrays", "Statistics"]
git-tree-sha1 = "ac425eea956013b51e7891bef3c33684b7d37029"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.11"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.DynamicPolynomials]]
deps = ["DataStructures", "Future", "LinearAlgebra", "MultivariatePolynomials", "MutableArithmetics", "Pkg", "Reexport", "Test"]
git-tree-sha1 = "d0fa82f39c2a5cdb3ee385ad52bc05c42cb4b9f0"
uuid = "7c1d4256-1411-5781-91ec-d7bc3513ac07"
version = "0.4.5"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "90630efff0894f8142308e334473eba54c433549"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.5.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "129b104185df66e408edd6625d480b7f9e9823a0"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.18"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "246621d23d1f43e3b9c368bf3b72b2331a27c286"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.2"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.ForwardDiff]]
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions", "StaticArrays"]
git-tree-sha1 = "2f18915445b248731ec5db4e4a17e451020bf21e"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.30"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "cabd77ab6a6fdff49bfd24af2ebe76e6e018a2b4"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.0.0"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "b5c7fe9cea653443736d264b85466bad8c574f4a"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.9.9"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "039118892476c2bf045a43b88fcb75ed566000ff"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "d88b17a38322e153c519f5a9ed8d91e9baa03d8f"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.1"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "53c7e69a6ffeb26bd594f5a1421b889e7219eeaa"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.Groebner]]
deps = ["AbstractAlgebra", "Combinatorics", "Logging", "MultivariatePolynomials", "Primes", "Random"]
git-tree-sha1 = "f52b2f66320da340bba1f749d43261b151124600"
uuid = "0b43b601-686d-58a3-8a1c-6623616c7cd4"
version = "0.2.9"

[[deps.GroupsCore]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "9e1a5e9f3b81ad6a5c613d181664a0efc6fe6dd7"
uuid = "d5909c97-4eac-4ecc-a3dc-fdd0858a4120"
version = "0.4.0"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "ed47af35905b7cc8f1a522ca684b35a212269bd8"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.2.0"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "d19f9edd8c34760dca2de2b503f969d8700ed288"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.4"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "f366daebdfb079fd1fe4e3d560f99a0c892e15bc"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.0"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d979e54b71da82f3a65b62553da4fc3d18c9004c"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2018.0.3+2"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "23e651bbb8d00e9971015d0dd306b780edbdb6b9"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.3"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "57af5939800bce15980bddd2426912c4f83012d8"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.1"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IterableTables]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Requires", "TableTraits", "TableTraitsUtils"]
git-tree-sha1 = "70300b876b2cebde43ebc0df42bc8c94a144e1b4"
uuid = "1c8ee90f-4401-5389-894e-7a04a3dc0f4d"
version = "1.0.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "9816b296736292a80b9a3200eb7fbb57aaa3917a"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.LabelledArrays]]
deps = ["ArrayInterfaceCore", "ArrayInterfaceStaticArrays", "ChainRulesCore", "LinearAlgebra", "MacroTools", "PreallocationTools", "RecursiveArrayTools", "StaticArrays"]
git-tree-sha1 = "b86aeff13358dfef82efd9f66a9d44705c9a4746"
uuid = "2ee39098-c373-598a-b85f-a56591580800"
version = "1.11.1"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics"]
git-tree-sha1 = "46efcea75c890e5d820e670516dc156689851722"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.5.4"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "361c2b088575b07946508f135ac556751240091c"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.17"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "e595b205efd49508358f7dc670a940c790204629"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "c27ed640732b1e9bd7bb8f40d987873d8f5b4bca"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.12"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9bd42b962d5c6182fa0d74b1970edb075fe313e5"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.3.6"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "e652a21eb0b38849ad84843a50dcbab93313e537"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.6.1"

[[deps.MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "114ef48a73aea632b8aebcb84f796afcc510ac7c"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.3"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "14cb991ee7ccc6dabda93d310400575c3cae435b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.2"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Metatheory]]
deps = ["AutoHashEquals", "DataStructures", "Dates", "DocStringExtensions", "Parameters", "Reexport", "TermInterface", "ThreadsX", "TimerOutputs"]
git-tree-sha1 = "a160e323d3684889e6026914576f1f4288de131d"
uuid = "e9d8d322-4543-424a-9be4-0cc815abe26c"
version = "1.3.4"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "6bb7786e4f24d44b4e29df03c69add1b63d88f01"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"

[[deps.MultistartOptimization]]
deps = ["ArgCheck", "DocStringExtensions", "NLopt", "Parameters", "Requires", "Sobol"]
git-tree-sha1 = "54b819498e65b7024ce5216aa7a1c40f5a00e7a9"
uuid = "3933049c-43be-478e-a8bb-6e0f7fd53575"
version = "0.1.3"

[[deps.MultivariatePolynomials]]
deps = ["ChainRulesCore", "DataStructures", "LinearAlgebra", "MutableArithmetics"]
git-tree-sha1 = "393fc4d82a73c6fe0e2963dd7c882b09257be537"
uuid = "102ac46a-7ee4-5c85-9060-abc95bfdeaa3"
version = "0.4.6"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4e675d6e9ec02061800d6cfb695812becbd03cdf"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.4"

[[deps.NLopt]]
deps = ["MathOptInterface", "MathProgBase", "NLopt_jll"]
git-tree-sha1 = "5a7e32c569200a8a03c3d55d286254b0321cd262"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "0.6.5"

[[deps.NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9b1f15a08f9d00cdb2761dcfa6f453f5d0d6f973"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.7.1+0"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"

[[deps.Observables]]
git-tree-sha1 = "dfd8d34871bc3ad08cd16026c1828e271d554db9"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optimization]]
deps = ["ArrayInterfaceCore", "ConsoleProgressMonitor", "DiffResults", "DocStringExtensions", "Logging", "LoggingExtras", "Pkg", "Printf", "ProgressLogging", "Reexport", "Requires", "SciMLBase", "SparseArrays", "TerminalLoggers"]
git-tree-sha1 = "f2dbe632d3aad1fb1e5ee7dbbeb4896aabb39da1"
uuid = "7f7a1694-90dd-40f0-9382-eb1efda571ba"
version = "3.8.2"

[[deps.OptimizationMultistartOptimization]]
deps = ["MultistartOptimization", "Optimization", "Reexport"]
git-tree-sha1 = "cd4b39bf624e4bb5964685ff9dbdf05968d5806f"
uuid = "e4316d97-8bbb-4fd3-a7d8-3851d2a72823"
version = "0.1.0"

[[deps.OptimizationNLopt]]
deps = ["NLopt", "Optimization", "Reexport"]
git-tree-sha1 = "690b3b3a1a6e6a8b527e330e6d37d40ab745a6ed"
uuid = "4e6fcdb7-1186-4e1f-a706-475e75c168bb"
version = "0.1.0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "1155f6f937fa2b94104162f01fa400e192e4272f"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.4.2"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a121dfbba67c94a5bec9dde613c3d0cbcf3a12b"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.3+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "0044b23da09b5608b4ecacb4e5e6c6332f833a7e"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.3.2"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PreallocationTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ForwardDiff"]
git-tree-sha1 = "ba66bf03b84ca3bd0a26aa2bbe96cd9df2f4f9b9"
uuid = "d236fae5-4411-538c-8e31-a6e3d9e00b46"
version = "0.4.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.PrettyTables]]
deps = ["Crayons", "Formatting", "Markdown", "Reexport", "Tables"]
git-tree-sha1 = "dfb54c4e414caa595a1f2ed759b160f5a3ddcba5"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "1.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "311a2aa90a64076ea0fac2ad7492e914e6feeb81"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "78aadffb3efd2155af139781b8a8df1ef279ea39"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.4.2"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RandomExtensions]]
deps = ["Random", "SparseArrays"]
git-tree-sha1 = "062986376ce6d394b23d5d90f01d81426113a3c9"
uuid = "fb686558-2515-59ef-acaa-46db3789a887"
version = "0.4.3"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.ReadStat]]
deps = ["DataValues", "Dates", "ReadStat_jll"]
git-tree-sha1 = "f8652515b68572d3362ee38e32245249413fb2d7"
uuid = "d71aba96-b539-5138-91ee-935c3ee1374c"
version = "1.1.1"

[[deps.ReadStat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "afd287b1031406b3ec5d835a60b388ceb041bb63"
uuid = "a4dc8951-f1cc-5499-9034-9ec1c3e64557"
version = "1.1.5+0"

[[deps.ReadableRegex]]
git-tree-sha1 = "befcfa33f50688319571a770be4a55114b71d70a"
uuid = "cbbcb084-453d-4c4c-b292-e315607ba6a4"
version = "0.3.2"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecursiveArrayTools]]
deps = ["Adapt", "ArrayInterfaceCore", "ArrayInterfaceStaticArraysCore", "ChainRulesCore", "DocStringExtensions", "FillArrays", "GPUArraysCore", "LinearAlgebra", "RecipesBase", "StaticArraysCore", "Statistics", "ZygoteRules"]
git-tree-sha1 = "4ce7584604489e537b2ab84ed92b4107d03377f0"
uuid = "731186ca-8d62-57ce-b412-fbd966d074cd"
version = "2.31.2"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Referenceables]]
deps = ["Adapt"]
git-tree-sha1 = "e681d3bfa49cd46c3c161505caddf20f0e62aaa9"
uuid = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"
version = "0.1.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "bf3188feca147ce108c76ad82c2792c57abe7b1f"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.0"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "68db32dff12bb6127bac73c209881191bf0efbb7"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.3.0+0"

[[deps.Roots]]
deps = ["CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "50f945fb7d7fdece03bbc76ff1ab96170f64a892"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.2"

[[deps.RuntimeGeneratedFunctions]]
deps = ["ExprTools", "SHA", "Serialization"]
git-tree-sha1 = "cdc1e4278e91a6ad530770ebb327f9ed83cf10c4"
uuid = "7e49a35a-f44a-4d26-94aa-eba1b4ca6b47"
version = "0.5.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.SciMLBase]]
deps = ["ArrayInterfaceCore", "CommonSolve", "ConstructionBase", "Distributed", "DocStringExtensions", "IteratorInterfaceExtensions", "LinearAlgebra", "Logging", "Markdown", "RecipesBase", "RecursiveArrayTools", "StaticArraysCore", "Statistics", "Tables"]
git-tree-sha1 = "3c955ccc10d4ca8910fe6c39ffcf05f27412b975"
uuid = "0bca4576-84f4-4d90-8ffe-ffa030f20462"
version = "1.46.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "Requires"]
git-tree-sha1 = "38d88503f695eb0301479bc9b0d4320b378bafe5"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "0.8.2"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "22395afdcf37d6709a5a0766cc4a5ca52cb85ea0"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "1.0.0"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SignedDistanceFields]]
deps = ["Random", "Statistics", "Test"]
git-tree-sha1 = "d263a08ec505853a5ff1c1ebde2070419e3f28e9"
uuid = "73760f76-fbc4-59ce-8f25-708e95d2df96"
version = "0.4.0"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.Sobol]]
deps = ["DelimitedFiles", "Random"]
git-tree-sha1 = "5a74ac22a9daef23705f010f72c81d6925b19df8"
uuid = "ed01d8cd-4d21-5b2a-85b4-cc3bdc58bad4"
version = "1.5.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "39c9f91521de844bad65049efd4f9223e7ed43f9"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.14"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StatFiles]]
deps = ["DataValues", "FileIO", "IterableTables", "IteratorInterfaceExtensions", "ReadStat", "TableShowUtils", "TableTraits", "TableTraitsUtils", "Test"]
git-tree-sha1 = "28466ea10caec61c476a262172319d2edf248187"
uuid = "1463e38c-9381-5320-bcd4-4134955f093a"
version = "0.8.0"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f94f9d627ba3f91e41a815b9f9f977d729e2e06f"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.7.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "23368a3313d12a2326ad0035f0db0c0966f438ef"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.2"

[[deps.StaticArraysCore]]
git-tree-sha1 = "66fe9eb253f910fe8cf161953880cfdaef01cdf0"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.0.1"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "2c11d7290036fe7aac9038ff312d3b3a2a5bf89e"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.4.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "0005d75f43ff23688914536c5e9d5ac94f8077f7"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.20"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "4352d5badd1bc8bf0a8c825e886fa1eda4f0f967"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.30"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SymbolicUtils]]
deps = ["AbstractTrees", "Bijections", "ChainRulesCore", "Combinatorics", "ConstructionBase", "DataStructures", "DocStringExtensions", "DynamicPolynomials", "IfElse", "LabelledArrays", "LinearAlgebra", "Metatheory", "MultivariatePolynomials", "NaNMath", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "TermInterface", "TimerOutputs"]
git-tree-sha1 = "027b43d312f6d52187bb16c2d4f0588ddb8c4bb2"
uuid = "d1185830-fcd6-423d-90d6-eec64667417b"
version = "0.19.11"

[[deps.Symbolics]]
deps = ["ArrayInterfaceCore", "ConstructionBase", "DataStructures", "DiffRules", "Distributions", "DocStringExtensions", "DomainSets", "Groebner", "IfElse", "Latexify", "Libdl", "LinearAlgebra", "MacroTools", "Markdown", "Metatheory", "NaNMath", "RecipesBase", "Reexport", "Requires", "RuntimeGeneratedFunctions", "SciMLBase", "Setfield", "SparseArrays", "SpecialFunctions", "StaticArrays", "SymbolicUtils", "TermInterface", "TreeViews"]
git-tree-sha1 = "4072e46467cfcaca1f7fe2a14f9b060da5edf7d2"
uuid = "0c5d862f-8b57-4792-8d23-62f2024744c7"
version = "4.10.3"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableShowUtils]]
deps = ["DataValues", "Dates", "JSON", "Markdown", "Test"]
git-tree-sha1 = "14c54e1e96431fb87f0d2f5983f090f1b9d06457"
uuid = "5e66a065-1f0a-5976-b372-e0b8c017ca10"
version = "0.2.5"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.TableTraitsUtils]]
deps = ["DataValues", "IteratorInterfaceExtensions", "Missings", "TableTraits"]
git-tree-sha1 = "78fecfe140d7abb480b53a44f3f85b6aa373c293"
uuid = "382cd787-c1b6-5bf2-a167-d5b971a19bda"
version = "1.0.2"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TermInterface]]
git-tree-sha1 = "7aa601f12708243987b88d1b453541a75e3d8c7a"
uuid = "8ea1fca8-c5ef-4a55-8b96-4e9afe9c9a3c"
version = "0.2.3"

[[deps.TerminalLoggers]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "987a3ebb20307530775f4def7eb9109cfa881748"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.ThreadsX]]
deps = ["ArgCheck", "BangBang", "ConstructionBase", "InitialValues", "MicroCollections", "Referenceables", "Setfield", "SplittablesBase", "Transducers"]
git-tree-sha1 = "d223de97c948636a4f34d1f84d92fd7602dc555b"
uuid = "ac1d9e8a-700a-412c-b207-f0111f4b6c0d"
version = "0.1.10"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "fcf41697256f2b759de9380a7e8196d6516f0310"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.0"

[[deps.TimerOutputs]]
deps = ["ExprTools", "Printf"]
git-tree-sha1 = "464d64b2510a25e6efe410e7edab14fffdc333df"
uuid = "a759f4b9-e2f1-59dc-863e-4aeb61b1ea8f"
version = "0.5.20"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "c76399a3bbe6f5a88faa33c8f8a65aa631d95013"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.73"

[[deps.TreeViews]]
deps = ["Test"]
git-tree-sha1 = "8d0d7a3fe2f30d6a7f833a5f19f7c7a5b396eae6"
uuid = "a2a6695c-b41b-5b7d-aed9-dbfdeacea5d7"
version = "0.3.0"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.UnPack]]
git-tree-sha1 = "387c1f73762231e86e0c9c5443ce3b4a0a9a0c2b"
uuid = "3a884ed6-31ef-47d7-9d2a-63182c4928ed"
version = "1.0.2"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.WeakRefStrings]]
deps = ["DataAPI", "InlineStrings", "Parsers"]
git-tree-sha1 = "b1be2855ed9ed8eac54e5caff2afcdb442d52c23"
uuid = "ea10d353-3f73-51f8-a26c-33c1cb351aa5"
version = "1.4.2"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"

[[deps.ZygoteRules]]
deps = ["MacroTools"]
git-tree-sha1 = "8c1a8e4dfacb1fd631745552c8db35d0deb09ea0"
uuid = "700de1a5-db45-46bc-99cf-38207098b444"
version = "0.2.2"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"
"""

# ╔═╡ Cell order:
# ╟─9e825451-37d8-4812-86c4-609e9296a179
# ╠═f14f868b-d6ad-43b5-82f6-5eba76fc344b
# ╠═15c84e2b-e828-4320-9140-09a4a6c88968
# ╠═e5d0b7c5-9d96-4702-b7dd-70f9c9e69ef2
# ╠═673d0937-a670-45a0-a7a5-6b9d9d270610
# ╟─fb00431c-4eee-4041-bba6-504cd00693d6
# ╠═1d6866f4-af9c-444f-ba93-65b5900b9ee1
# ╟─0b1a51ea-073f-11ed-1d50-b1e334a7bcb0
# ╠═d8cfaacb-b3de-41de-a2d1-6d1eeabaa2c7
# ╠═90f6bf62-752d-4018-beec-fa2239fc2db1
# ╟─e8d0027f-d4e3-41eb-9181-821e3c586a22
# ╠═cd60add2-cb80-431c-8119-def8a1b961ef
# ╠═37618646-be2e-474d-972f-e2e84b914765
# ╠═e5abb6c6-7add-4f73-9ffd-26173c92b5a7
# ╠═eea9988e-3fcf-48d8-93c3-67d94e1ab949
# ╠═16f6d020-561c-43ea-bd1f-606890b7e009
# ╟─1555e6b0-2fbe-466c-884d-186f5b20e6e3
# ╟─dd95cda0-620e-4fee-8d35-e58ce9aebfd6
# ╠═a6f923d1-0292-4168-a40a-27929be0815c
# ╠═e0135a19-bf81-46c0-8b25-15e320a39e78
# ╠═a6d14acc-8b2a-4736-a774-a03174f3ca42
# ╟─7b4f2179-15fe-4651-94d0-1765c6ed8b83
# ╠═6c8cf040-c64e-4a03-b42c-b1df688e3572
# ╠═6d3a317e-8e14-457b-99a5-90baf85bfbaf
# ╠═6ac092f8-8f03-4a1f-91f3-1177efd99745
# ╟─279101cf-0d99-4446-b25a-b7e2dd7b75f4
# ╠═6f3e7b25-95f2-4f77-b087-14c6543117c8
# ╟─2816d37a-7ea1-4c0a-9149-373c9a4361a6
# ╠═fa55d156-c645-4b16-b253-05f7802cfc36
# ╠═fda2e873-d974-42f3-97fe-94223ab0b077
# ╠═4644fac4-a6ed-48f4-b822-40f1df4977f8
# ╠═043624c1-0ee5-48db-80d2-82331ee5552e
# ╠═b7ff7f0f-a9e3-46c2-bf39-84c6ef67d272
# ╠═a44b1c72-3bb5-4ff7-9378-275dffba2358
# ╠═e518738f-6c0b-450d-a119-b17f1f71e5f0
# ╠═7fda0a57-1148-4f0d-b913-c148a7e6236f
# ╠═8174336f-bb98-4fb5-a275-53d16e4e1ab3
# ╠═b052b26f-9ac2-49b0-a713-3a81169e7dcb
# ╠═1146e3f5-d701-42fa-929f-4699efb5a5f6
# ╠═6e308459-f126-4742-8e6f-33a6420bda19
# ╠═309a39ef-a3a8-45b0-8891-a8cd0aebfe5a
# ╠═eb45f538-6384-475e-9de5-6bb6152a3cc1
# ╠═6f14d6c3-ba2c-41f1-acbf-8125fec69d7c
# ╟─e5d106d5-069b-464d-8555-8bd546a0a042
# ╟─a2d3f53b-2c15-4b17-83e9-25f6b975ee65
# ╠═5be2de51-2391-457a-8b94-c07e364d3eef
# ╠═5f8a5a6c-a628-46a6-a90c-8df775643e9c
# ╠═2eb5df98-f294-481d-b0b3-947f020b0d47
# ╠═b4852dd9-97f5-402f-916e-e387e52c0694
# ╠═26cb90ad-e3ce-4a9f-9d9b-36e43a103057
# ╠═7d6088f8-2ffe-4238-9cb9-c0b7e765cc4d
# ╠═d0da815f-65ea-40d7-ba77-d93afcb9a889
# ╟─082e201e-6732-432f-9323-162c84b68a25
# ╠═62d8a4e2-1a24-48c8-b044-6c25147dcf51
# ╠═8ca25b51-89fa-47af-b2a0-e109e9b0f98a
# ╠═d0925b77-4f4f-4a97-a6bf-982b318c1f2a
# ╠═616bc951-ad6c-4475-8394-3104a7d336b8
# ╠═2db97ad1-7662-48f6-b546-ff11b267013b
# ╠═f8888d57-2e06-4dcc-b7b6-db830840d426
# ╠═0d6b97b7-a1ef-4cd5-9868-4d9ea74591e4
# ╠═43174c3f-e313-46c0-ba79-33b5ace09c21
# ╠═ce061fb9-1678-47c9-a13e-90cf90d0338c
# ╠═6819f7ad-e0b2-4ef4-be9e-665629bbc127
# ╟─b29b41b6-cc32-46a8-8c0c-dd6ccc990f71
# ╟─b4c18055-3168-474d-8dfe-6ea5105e0e80
# ╟─89960fdf-82f3-49df-855e-35d21fb8b24d
# ╟─b3ed4900-2279-49af-891d-0e52a29c4b08
# ╟─730a6ef9-619e-4008-b03a-120cb5313850
# ╟─9d96e880-41f7-41e4-a1c6-3e39f5cbc398
# ╟─fcb8db73-9f24-4814-98c4-b08049d8029a
# ╟─04153526-eacc-4a6f-a175-d7f9ed44ce7e
# ╟─8d76aad5-59b6-4ffe-83da-ad80bcb9a893
# ╟─3895355e-c93c-4225-89e3-3c907b6d11de
# ╟─b4aa9303-f440-4cea-82b4-943a3c14bfac
# ╠═fcabed90-a37b-4f1d-aea9-4914a4f95a1c
# ╠═bd7a67c1-e465-41c4-8cb3-dfa58e79386b
# ╟─bb63db54-1065-46e7-989c-50d080230938
# ╟─5fc360fe-62d3-446e-91d4-88bb5bdf3d9f
# ╟─a128a97c-0192-4cb7-8021-e5abb9f152da
# ╟─787ca466-6d72-45f1-b6fc-878ddd2ab75b
# ╟─aec3912e-5847-4a46-a658-bce97258c54c
# ╟─41dd7623-d71e-46ca-8580-c44c2acca2cb
# ╠═73f900a1-c1da-4efa-847a-48c5e8cb0717
# ╠═fe128ce0-87a3-4c13-b596-0e2c34650058
# ╠═468ca176-8079-423e-bad4-465d12d1b78a
# ╟─9eb5034f-af9b-4c46-b480-8b03c977757c
# ╠═ab5b369e-5834-4154-ac75-fa10cc1d0d8b
# ╠═dad754a5-4b2f-4b3e-9910-7a7caea711a9
# ╠═e6ffda4f-a1a0-42ad-bb42-e755076b1fb9
# ╠═4981e7b4-d43b-479a-bd35-752462cf84d2
# ╠═2bda56a6-e49f-45c6-a814-016ce5a4cc4c
# ╠═dd774ed7-6928-45c9-bede-ecfeb3884335
# ╠═c3d5e417-f8e0-40ff-8fcf-a8b4111f4ec3
# ╠═054a70a2-5f63-4742-8782-102f229d47ff
# ╠═17b1c24a-7337-4917-9a0f-5c2defa1bf96
# ╠═0860e535-33f1-4e22-aa46-eba71dd95769
# ╠═b9fee53f-9e16-484b-9514-ba821ce8d85e
# ╠═65e3be85-d508-4817-a632-4b32d75bffcf
# ╠═fb0e25b0-8066-4370-8f5b-527a35091017
# ╠═900491e1-1962-4543-ab5c-04ee8b95bc5c
# ╟─bb9fa247-42c8-4ab4-ae36-2892fa5f1651
# ╠═ec31d541-2087-4cd5-93c5-6a718d61ce4b
# ╠═293fa922-ee60-4b7d-8085-8d41d5d2363a
# ╟─30ac4eda-2829-4879-854a-36a08acebb36
# ╠═d2590b30-c26a-4502-a374-1e388894abc3
# ╠═6a218989-d242-4318-9ceb-fceddbfeaabc
# ╠═89668f54-762b-4cf3-bf7c-44f43b796b29
# ╠═00690911-2b64-4b26-aa58-23fc3e23bb48
# ╠═d58962fa-4073-467c-b3ec-418e4ecdce42
# ╟─a2c782a8-aa86-49f9-9b96-2ced327a876f
# ╠═e303a91c-33a8-41d1-96b4-8a213dab5698
# ╟─4ae24902-ddb0-4abf-b2e5-ecbad83db5a0
# ╠═e6ebdc70-fcb2-4d6e-874b-951671480e33
# ╠═4fc16124-5f72-4c1f-ab18-a99f6bcb027e
# ╠═6821ead6-0cd5-4d92-8fab-0f3356de5637
# ╟─121b5939-216b-4864-9845-30ec53989f56
# ╠═7bbe5a83-6d2d-4339-a3dd-3af947ee2172
# ╠═fa04eb87-3159-4505-80d4-300f91962ab8
# ╠═ce84af63-f8fc-47e6-8e0b-a3ac18283534
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

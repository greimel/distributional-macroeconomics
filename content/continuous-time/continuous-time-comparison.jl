### A Pluto.jl notebook ###
# v0.19.38

#> [frontmatter]
#> chapter = 6
#> section = 3
#> order = 3
#> title = "Comparing different solution methods"
#> layout = "layout.jlhtml"
#> tags = ["continuous-time"]
#> description = ""

using Markdown
using InteractiveUtils

# ╔═╡ 5cbd313e-b357-4b00-bb7d-4c87f8522cdd
using EconPDEs

# ╔═╡ 02c84fbb-f2cb-4051-8bda-62f5e2a4c6d1
using PlutoTest

# ╔═╡ d704e278-6ae8-4166-bade-599e71e54336
using DataFrameMacros

# ╔═╡ d3612723-d4eb-4887-b902-fa47e6a26226
using InfinitesimalGenerators

# ╔═╡ 5bc621c2-1d87-49a9-8dc6-d81d5f998640
using QuantEcon: gth_solve

# ╔═╡ abff0506-1792-420d-80d8-cb3645483443
using Arpack

# ╔═╡ 70930f84-fe69-4fa3-a37c-eec9a377526c
using StructArrays

# ╔═╡ 6f92c837-1760-423e-9777-9db9ad758475
using PlutoUI

# ╔═╡ 563c267e-0c6a-4b90-81a7-fc8ff4c73c75
using AlgebraOfGraphics, CairoMakie

# ╔═╡ 9fa33ac0-bb2f-4055-b124-8c7c19621226
using LinearAlgebra

# ╔═╡ 276e171b-271d-4610-b9bb-01b212193123
using SparseArrays

# ╔═╡ aecea3fe-f0ee-4953-857b-29d6a5640530
using DataFrames

# ╔═╡ 9353566b-a70e-4dbe-8f02-b16d2c0570d2
using Chain: @chain

# ╔═╡ 0b9c3978-0eb9-4f80-b51d-540d4ed88c3e
using Roots: find_zero, Brent

# ╔═╡ 639eeead-5c9e-445b-8e40-c811634c998a
using MarkdownLiteral: @markdown

# ╔═╡ 674afcf7-925e-4cc1-ac4e-9dba4cdc4693
using HypertextLiteral

# ╔═╡ cdd24518-1cd3-43af-9a88-f0a8cd3cd6ed
md"""
`continuous-time-comparison.jl` | **Version 1.0** | *last updated: May 31, 2023* | 
"""

# ╔═╡ 0edf9214-7213-4b1d-8fa6-54a365525d29
@markdown("""
# Huggett model in continuous time

In this notebook, we consider a Huggett economy in continuous time. Income is a $(@htl("<s>Poisson process</s>")) continuous time Markov Chain with two states. 

We compare the implementation of 

* Achdou et al (2021): (algorithm described in [their online appendix](https://benjaminmoll.com/wp-content/uploads/2020/02/HACT_Numerical_Appendix.pdf), code follows the [Matlab code snippets on Ben Moll's website](https://benjaminmoll.com/codes/)
* EconPDEs.jl by Mattieu Gomez: code adapted from package tests
""")

# ╔═╡ 1c23af35-64f5-4b96-9eab-98ee9b1ecb4c
md"""
## Model
"""

# ╔═╡ 64d44910-988c-4e24-baf6-6408da65bd21
@kwdef struct Moll
	
	σ::Float64 = 2.   # risk aversion coefficient u'(c) = c^(-σ)
	ρ::Float64 = 0.05 # rate of time preference

	z::Matrix{Float64} = [0.1 0.2]  # income state (row vector)
	λ::Matrix{Float64} = [0.02 0.03] # intensities (row vector)

	# asset grid parameters
	N_a::Int64 = 500
	aₘᵢₙ::Float64 = - 0.1
	aₘₐₓ::Float64 = 1.0

	Δa::Float64 = (aₘₐₓ - aₘᵢₙ)/(N_a - 1)
end

# ╔═╡ 95764a2a-0a92-447c-98a8-df22167abda6
function generator((; λ))
	λ₁₂, λ₂₁ = λ
	Λ = [-λ₁₂ λ₁₂;
	    λ₂₁ -λ₂₁]
end

# ╔═╡ a5233ee4-ba60-4428-8de5-4b012ed43408
function construct_A_switch((; λ, N_a))
	id = I(N_a)
	λ₁₂, λ₂₁ = λ
	A_switch_1 = hcat(-λ₁₂ * id,  λ₁₂ * id)
	A_switch_2 = hcat( λ₂₁ * id, -λ₂₁ * id)

	[-λ₁₂ * id   λ₁₂ * id;
	  λ₂₁ * id  -λ₂₁ * id]

	#A_switch   = vcat(A_switch_1, A_switch_2)
end

# ╔═╡ ac27dd69-1216-4a02-ba62-ebb098b33fa1
md"""
We work with an equi-spaced asset grid with $N_a$ grid points. The difference between two grid points is denoted by $\Delta a$. (Section 7 in the online appendix explains how to deal with non-uniform grids.)
"""

# ╔═╡ 3dd0c186-04f0-425c-93dc-825c7d4b606e
function construct_a(m)
	(; N_a, aₘᵢₙ, aₘₐₓ) = m
	range(aₘᵢₙ, aₘₐₓ, N_a) # asset grid (column vector)
end

# ╔═╡ 23a991be-7c8b-45e2-bd75-af5e146fc6b0
m = Moll();

# ╔═╡ 2289c7a7-3493-4bfb-8ffe-31b074d17b14
md"""
## HJB equation (implicit method)
"""

# ╔═╡ 0caee0bf-77c0-4542-a68e-9052d230ca74
md"""
$$\rho v_1(a) = \max_c u(c) + v_1'(a) (z_1 + ra - c) + \lambda_1(v_2(a) - v_1(a))$$
$$\rho v_2(a) = \max_c u(c) + v_2'(a) (z_2 + ra - c) + \lambda_2(v_1(a) - v_2(a))$$

```math
\rho \pmatrix{v_1(a) \\ v_2(a)} = \pmatrix{u(c^*_1(a)) \\ u(c^*_2(a))} + \pmatrix{v'_1(a)(z_1 + ra - c_1^*(a)) \\ v'_1(a)(z_1 + ra - c_2^*(a))} + \pmatrix{-\lambda_1 & \lambda_1 \\ \lambda_2 & -\lambda_2} \pmatrix{v_1(a) \\ v_2(a)} 
```
"""

# ╔═╡ 38555a5e-5d39-4c61-9381-94c8fb67a257
md"""
**Algorithm**

(Notation: $v_{i,j}$ is short-hand notation for $v_j(a_i)$.)

Start with an initial guess for the value function $v_{i,j}^0$. A natural choice is
$$v_{i,j}^0 = \frac{u(z_j + ra_i)}{\rho}$$.

For i = 1, ... maxit 

(1) Approximate $(v_{i,j}^n)'$, $j=1,2$ using a finite difference method:

(Notation: Superscript $n$ is omitted.) 

- forward difference $v_{i,j,F}' = \frac{v_{i+1,j} - v_{i,j}}{\Delta a}$
- backward diffrence $v_{i,j,B}' = \frac{v_{i,j} - v_{i-1,j}}{\Delta a}$

The state constraint $a \ge a_\text{min}$ needs to be enforced by setting $v'_{1,j,B} = u'(z_j + ra_\text{min})$.

- savings according to forward difference $s_{i,j,F} = z_j + ra_i - (u')^{-1}(v'_{i,j,F})$
- savings according to backward difference $s_{i,j,B} = z_j + ra_i - (u')^{-1}(v'_{i,j,B})$

The finite difference approximation of $(v_{i,j}^n)'$ is

$$v'_{i,j} = v'_{i,j,F} 1_{s_{i,j,F}>0} + v'_{i,j,B} 1_{s_{i,j,B}<0} + \bar{v}_{i,j} 1_{s_{i,j,F} \le 0 \le s_{i,j,B}}$$

where $\bar{v}_{i,j} = u'(s_j + r a_i)$.

(We assume concavity of the value function here so that the case $s_{i,j,F}>0$ and $s_{i,j,B}<0$ cannot occur.)

(2) Compute the consumption policy implied by the value function $c_{i,j}^n = (u')^{-1}[(v_{i,j}^n)']$

(3) Find updated value function $v^{n+1}$:

```math
\begin{align}
&\frac{v_{i,j}^{n+1} - v_{i,j}^{n}}{\Delta} + \rho v_{i,j}^{n+1} = \\
&u(c_{i,j}^n) + (v_{i,j,F}^{n+1})'[z_j + ra_i - c_{i,j,F}^n]^+ + (v_{i,j,B}^{n+1})'[z_j + ra_i - c_{i,j,B}^n]^- + \lambda_j (v_{i,-j}^{n+1} - v_{i,j}^{n+1})
\end{align}
```

where $\Delta$ is the step size.

This is a system of $2N_a$ linear equations. Since $v^{n+1}$ is implicitly defined by the equations above, this approach is referred to as the implicit method.

The system of equations can be written in matrix notation as 

$$\frac{1}{\Delta} (v^{n+1} - v^n) + \rho v^{n+1} = u^n + A^n v^{n+1}$$

The $2N_a \times 2N_a$ matrix $A^n$ can be written as a sum of two matrices $\bar{A}^n$ and $A_\text{switch}$:

$A^n = \bar{A}^n + A_\text{switch} = \begin{pmatrix} \bar{A}_{11}^n & 0 \\ 0 & \bar{A}_{22}^n \end{pmatrix} + \begin{pmatrix} -\lambda_1 I & \lambda_1 I \\ \lambda_2 I & -\lambda_2 I \end{pmatrix}$

where $I$ is a $N_a \times N_a$ identity matrix. Since $A_\text{switch}$ stays unchanged, it can be pre-computed outside the for-loop.

The $N_a \times N_a$ submatrices $\bar{A}_{11}^n$ and $\bar{A}_{22}^n$ are tri-diagonal:
- The -1 diagonal is filled with $x_{i,j} = - \frac{(s^n_{i,j,B})^-}{\Delta a}$, $i=2, \dots N_a$
- The main diagonal is filled with $y_{i,j} = - \frac{(s^n_{i,j,F})^+}{\Delta a} + \frac{(s^n_{i,j,B})^-}{\Delta a}$, $i=1, \dots N_a$
- The +1 diagonal is filled with $z_{i,j} = \frac{(s^n_{i,j,F})^+}{\Delta a}$, $i=1, \dots N_a - 1$

Since $A^n$ is a sparse matrix, computers can solve the system of linear equations quickly even for large $N_a$.

(4) Stop if $v^{n+1}$ is close enough to $v^n$. 

"""

# ╔═╡ fb88f286-0bce-4cb0-9d52-83b373f6fbcf
md"""
## Solving HJB with `EconPDEs.jl`
"""

# ╔═╡ e507dfb0-08a4-45a8-b342-09ca989e05f5
md"""
## Test from EconPDEs
"""

# ╔═╡ 69583d36-6190-4483-87fb-291d9cf82c77
states(a, z) = (; a, z)

# ╔═╡ db13a6b9-2de2-47ce-8af0-99d868183354
function dvs_v((; v1, v1a_up, v1a_down, v2, v2a_up, v2a_down))
	dvf = [v1a_up, v2a_up]
	dvb = [v1a_down, v2a_down]
	v = [v1, v2]

	(; dvf, dvb, v)
end

# ╔═╡ 56234834-6de9-4e44-95af-c5e3fb828a9d
u_prime_inv(x, (; σ)) = x^(-1/σ)

# ╔═╡ a10ef170-045c-439e-a23b-122e812853aa
u_prime(c, (; σ)) = c^(-σ)

# ╔═╡ 68aeb6cf-cff8-49e9-8c39-d94aadd5a444
u(c, (; σ)) = c > 0 ? σ == 1 ? log(c) : c^(1-σ)/(1-σ) : 10.0 * c - 100.0

# ╔═╡ e197b24a-7e51-45e4-8186-11f88bf48de6
# ╠═╡ disabled = true
#=╠═╡
function clean_variables(nt, solname, statename, n)
	map(1:n) do i
		sol_key = Symbol(solname, i)
		up_key = Symbol(sol_key, statename, "_up") => Symbol(solname, statename, "_up")
		down_key = Symbol(sol_key, statename, "_down") => Symbol(solname, statename, "_down")
		
		(; solname => nt[sol_key], up_key[2] => nt[up_key[1]], down_key[2] => nt[down_key[1]])
	end |> DataFrame
end
  ╠═╡ =#

# ╔═╡ e57a2dfa-0937-49cf-9161-fa1e39fb5e80
function consumption_and_drift((; a, z), dv, (; r, σ))
	#dv = max(dv, eps(0.0))
	c = u_prime_inv(dv, (; σ))
	ȧ = z + r * a - c

	(; c, ȧ, dv)
end

# ╔═╡ 91c8dce8-07b8-46d2-8c61-e6bccced64e4
function consumption_and_drift₀((; a, z), (; r, σ))
	ȧ = zero(a)
	c = z + r * a

	dv = u_prime(c, (; σ))

	(; c, ȧ, dv)
end

# ╔═╡ 7734c75e-5f2b-4807-8b9a-82e7e010edac
function consumption_and_drift_upwind(state, dvf, dvb, (; σ, r, aₘᵢₙ, aₘₐₓ))
	# consumption and savings with forward difference
	(; dv, ȧ, c) =  consumption_and_drift(state, dvf, (; σ, r))
	if ȧ > 0 && state.a < aₘₐₓ
		return (; dv, ȧf = ȧ, ȧb = 0.0, c, ȧ)
	end
	# consumption and savings with backward difference
	(; dv, ȧ, c) = consumption_and_drift(state, dvb, (; σ, r))
	if ȧ < 0 && state.a > aₘᵢₙ
		return (; dv, ȧf = 0.0, ȧb = ȧ, c, ȧ)
	end
	# consumption and derivate of value function at steady state
	(; dv, ȧ, c) = consumption_and_drift₀(state, (; σ, r))
	return (; dv, ȧf = 0.0, ȧb = 0.0, c, ȧ)
end

# ╔═╡ ad4b46b0-b6fb-4711-9501-0912eda3d6c3
function optional_to_df(optional, agrid, zgrid)
	@chain optional begin
		DataFrame
		@transform!(:a = @bycol collect(agrid))
		stack(Not(:a))
		@transform!(:i_z = parse(Int, :variable[end]))
		@transform!(:variable = :variable[1:end-1])
		unstack(:variable, :value)
		@transform!(:z = zgrid[:i_z])
		select(Not(:i_z))		
		rename!(:s => :ȧ)
	end
end

# ╔═╡ aeefc434-3375-400e-a360-4dc263b96252
# ╠═╡ disabled = true
#=╠═╡
test = let
	r = 0.03
 	M = EconPDEsFast(m, r)
	
	@elapsed (; residual_norm, optional, agrid, zgrid) = _solve_HJB_econpdes(M, r)
end
  ╠═╡ =#

# ╔═╡ 16b48db6-ef8d-4d7b-9e9f-b281143d0a90
md"""
### Check results
"""

# ╔═╡ 1fbf9f02-7ad2-4ddd-ac2c-92af79c9ed02
md"""
## KF equation
"""

# ╔═╡ 2e09462e-6030-4e33-bc7a-9d2faed4ca74
md"""
$$0 = - \frac{d}{da}[s_1(a)g_1(a)] - \lambda_1 g_1(a) + \lambda_2 g_2(a)$$
$$0 = - \frac{d}{da}[s_2(a)g_2(a)] - \lambda_2 g_2(a) + \lambda_1 g_1(a)$$

where $s_j(a) + z_j + ra - c_j(a)$

$$1 = \int_{\bar{a}}^\infty g_1(a)da + \int_{\bar{a}}^\infty g_2(a)da$$

**Algorithm**

A finite difference approximation of the KF equation results into the matrix equation $A^T g = 0$ where $A$ is the matrix from implicit algorithm for the HJB equation.
"""

# ╔═╡ 4ebeb5e8-505d-4ed2-8d37-4a64ceb7d796
#using QuantEcon: gth_solve

# ╔═╡ 6c895dce-9c0a-499c-9e95-45c29ee5a459
md"""
```math
\begin{align}
\dot g = (1-m) A g_t + m g_t \\
g_{t+\Delta} - g_t \approx A g_t \Delta \\
g_{t+\Delta} = (1-m)(I + \Delta A) g_t + m g_0 \\
g^* = (1-m)(I + \Delta A) g^* + m g_0 \\
(I - (1-m)(I + \Delta A)) g^* = m g_0 \\
(I - \frac{1-m}{m}\Delta A)) g^* = g_0 \\
\end{align}
```
"""

# ╔═╡ a8279bcc-4cdf-4a89-bb75-e7e1f617643d
function solve_KF_moll(A)
	N = size(A, 1)
	AT = copy(transpose(A))
	b = zeros(N)
	
	i_fix = 1
	b[i_fix] = .1
	AT[i_fix,:] .= 0.0
	AT[i_fix,i_fix] = 1.0

	g = AT\b

	g ./ sum(g)
end

# ╔═╡ eb565042-1823-4d5a-b4d1-ee314dccd4e0
function solve_KF((; N_a, Δa), A)

	#N_a = an
	
	AT = copy(transpose(A))
	b = zeros(2*N_a, 1)

	
	i_fix = 1
	b[i_fix] = .1
	AT[i_fix,:] .= vec(hcat(zeros(1, i_fix-1), 1., zeros(1,2*N_a-i_fix)))

	try
		global g_stacked = AT\b
	catch e
		if e isa SingularException
			@warn "SingularException – added noise"
			
			global g_stacked = (AT + I * √eps())\b
		else
			rethrow(e)
		end
	end
	
	g_sum = sum(g_stacked) * Δa
	g_stacked_norm = g_stacked ./ g_sum

	@assert sum(g_stacked_norm) * Δa ≈ 1
	
	g = reshape(g_stacked_norm, N_a, 2)

end

# ╔═╡ e5cfdbc2-f592-40bb-ba5d-07f455dd5bd4
md"""
## Putting everything together
"""

# ╔═╡ 34779464-624a-446d-b741-81524509aee6
#using QuantEcon: gth_solve

# ╔═╡ a2da91e4-cedd-446f-81d4-adb008c01d5b
function stationary_distribution(A; δ = 0.0, ψ = InfinitesimalGenerators.Zeros(size(A, 1)))
    δ >= 0 ||  throw(ArgumentError("δ needs to be positive"))
    if δ > 0
        g = abs.((δ * I - A') \ (δ * ψ))
    else
        η, g = InfinitesimalGenerators.principal_eigenvalue(A')
        abs(η) <= 1e-5 || @warn "Principal Eigenvalue does not seem to be zero"
    end
    g ./ sum(g)
end

# ╔═╡ dd89d622-6066-4d52-a4ca-366822c63661
function solve_KF_death(A)
	N = size(A, 1)
	g₀ = fill(1/N, N)
	
	stationary_distribution(A, δ = 1e-14, ψ = g₀)
end

# ╔═╡ 33a6780e-7f5a-4677-a51c-c002b7b86dc3
function solve_KF_eigs(A)	
	stationary_distribution(A)
end

# ╔═╡ eade1b48-7d1c-4d17-aaf5-14d111a13556
function solve_KF_iterate(A, Δ, g₀=fill(1/size(A,1), size(A,1)))
	g = copy(g₀)
	B = (I + Δ * A')
	for i ∈ 1:50000
		g_new = B * g

		crit = maximum(abs, g_new - g)
		i % 1000 == 0 && @info crit
		if crit < 1e-12
			@info "converged after $i iterations"
			return g
		end
		g .= g_new
	end
	g
end

# ╔═╡ deeff3d5-3856-43cb-8469-2a4d6d7fca4f
md"""
## Equilibrium interest rate
"""

# ╔═╡ b6101102-2054-4932-b6b6-5070cd84f2be
md"""
$$0 = \int_{\bar{a}}^\infty ag_1(a)da + \int_{\bar{a}}^\infty ag_2(a)da = S(r)$$
"""

# ╔═╡ a8f1dc13-b73c-44b0-8e63-82431f904313
initial_bracket = (0.01, 0.03)

# ╔═╡ e4035abe-11a1-4c2c-8312-4d1c71e2f9ab
md"""
# Appendix 
"""

# ╔═╡ ab69ab43-99bf-4a92-9698-70e170761e82
md"""
## HJB equation (explicit method)
"""

# ╔═╡ 5c7c3dc0-7848-4431-829e-508664d9c5af
md"""
**Basic dea**

Start with some initial guess $v_{i,j}^0$ and update $v_{i,j}^n$ as follows:

$$\frac{v_{i,j}^{n+1} - v_{i,j}^n}{\Delta} + \rho v_{i,j}^n = u(c_{i,j}^n) + (v_{i,j}^n)'(z_j + ra_i - c_{i,j}^n) + \lambda_i (v_{i,-j}^n - v_{i,j}^n)$$

In contrast to the implicit method, we can rearrange for $v_{i,j}^{n+1}$ in the equation above.

The disadvantage of the explicit method is that it converges only if $\Delta$ is not too large.
"""

# ╔═╡ b099dbbf-9648-44c5-984c-fdd80ee81469
function solve_HJB_explicit(m::Moll, r; maxit = 100000, crit = 1e-6)

	(; σ, ρ, z, λ, N_a, aₘᵢₙ, aₘₐₓ, Δa) = m
	da = Δa

	# construct asset grid
	a = construct_a(m)
	Λ = generator(m)
	
	# initialize arrays for forward and backward difference
	dvf = zeros(N_a, 2)
	dvb = zeros(N_a, 2)

	# initial guess for value function
	v₀ = zeros(N_a, 2)
	for (i, zᵢ) in enumerate(z)
		v₀[:,i] = (zᵢ .+ r * a).^(1-σ) / (1-σ) / ρ
	end
	v = v₀

	# initialize vector that keeps track of convergence
	dist = - ones(maxit)

	# step size for updating the value function
	Δ = .9 * da / (z[2] .+ r.*aₘₐₓ)

	for it in range(1, maxit)

		# forward difference
		dvf[1:N_a-1,:] = (v[2:N_a,:] - v[1:N_a-1,:]) / da
		dvf[N_a,:] = (z .+ r * aₘₐₓ) .^ (-σ) # boundary condition a  <= a_max

		# backward difference
		dvb[2:N_a,:] = (v[2:N_a,:] - v[1:N_a-1,:]) / da
		dvb[1,:] = (z .+ r * aₘᵢₙ) .^ (-σ) # boundary condition a >= a_min
	
		I_concave = dvb .> dvf # problems if value function not concave

		# consumption and savings with forward difference
		cf = dvf .^ (-1/σ)
		ȧf = z .+ r .* a - cf

		# consumption and savings with backward difference
		cb = dvb .^ (-1/σ)
		ȧb = z .+ r .* a - cb

		# consumption and derivate of value function at steady state
		c0 = z .+ r .* a
		dv0 = c0 .^ (-σ)

		If = ȧf .> 0 # positive drift => forward difference
		Ib = ȧb .< 0 # negative drift => backward difference
		Ib[N_a,:] .= 1. # make sure backward difference is used at last grid point
		If[N_a,:] .= 0.
		I0 = (1 .- If .- Ib) # steady state
	
		dv_upwind = dvf.*If + dvb.*Ib + dv0.*I0
	
		c = dv_upwind .^ (-1/σ)
		ȧ = z .+ r.*a - c
		u = c.^(1-σ)/(1-σ)
		
		# HJB equation
		v_change = u + dv_upwind .* ȧ + v * Λ' - ρ*v

		#v_switch = zeros(N_a, 2)
		#v_switch[:,2] = v[:,1]
		#v_switch[:,1] = v[:,2]
		#v_change = u + dv_upwind .* ȧ + ones(N_a,1)*λ.*(v_switch - v) - ρ*v

		# updating the value function
		v .= v + Δ * v_change

		dist[it] = maximum(abs.(v_change))

		if dist[it] < crit
			
			return (; v, c, ȧ, it_last=it, dist)

		end

	end

	error("Algorithm did not converge")

end

# ╔═╡ 4bfa5d92-b177-4442-b045-e05cc48b6cc4
t_expl = @elapsed solve_HJB_explicit(m, 0.03; crit=1e-6)

# ╔═╡ bae3209a-1ad6-4e1f-833d-f1074e8239dc


# ╔═╡ fd0fb774-a805-4739-b570-0d2e191a3294
md"""
## Implicit vs. explicit method
"""

# ╔═╡ 8987f1ce-e290-4930-909a-3e09a9113a7a
md"""
On my computer, it takes 0.4 seconds to reach convergence with the implicit method (assuming a tolerance of $10^{-6}$), while it takes approximately 6 seconds with the explicit method.
"""

# ╔═╡ eb7bbd16-04d3-4d7d-b4af-e77f89e4180e
md"""
**Implicit method**
"""

# ╔═╡ bd2e4d28-c68b-45c9-b68e-b477b44fcd75
md"""
**Explicit method**
"""

# ╔═╡ 6e607792-e297-483f-8917-c871fa0c26d0
md"""
## Helper functions
"""

# ╔═╡ 532d0b24-6ace-4579-bf2a-d12c07ee9436
function results_to_df(m; v, c, ȧ, g=nothing)

	N_z = 2
	
	(; N_a, z) = m

	a = construct_a(m)

	df = DataFrame()
	df.a = a * ones(1, N_z) |> vec
	df.z = ones(N_a, 1) * z |> vec
	df.c = c |> vec
	df.ȧ = ȧ |> vec
	df.v = v |> vec
	
	if ! isnothing(g)
		df.g = g |> vec
	end

	df
	
end

# ╔═╡ f3e0b42f-370d-4887-b6a1-d9ecd83c6275
md"""
## Julification of the code
"""

# ╔═╡ b9efca1a-b978-4794-b421-ee9df889a6a8
function statespace(m)
	(; N_a, aₘᵢₙ, aₘₐₓ, z) = m
	a_grid = range(aₘᵢₙ, aₘₐₓ, N_a)
	z_grid = z |> vec

	[(; a, z) for a ∈ a_grid, z ∈ z_grid]
end

# ╔═╡ 01234bcc-a428-4922-969f-9a1ba49c8d62
function statespace_inds(m)
	(; N_a, z) = m
	N_z = length(z)

	[(; i_a, i_z) for i_a ∈ 1:N_a, i_z ∈ 1:N_z]
end

# ╔═╡ 1971dce8-ca39-4ced-88af-c65105977ac4
v₀((; z, a), (; r, σ, ρ)) = u(z + r*a, (; σ))/ρ

# ╔═╡ e12a5f08-7e5b-418c-a08f-8d179e912774
begin
	
Base.@kwdef struct EconPDEsFast
    # income process parameters
	zgrid::Vector{Float64}=[0.5, 1.5]
	Λ::Matrix{Float64}= [-0.2 0.2; 0.2 -0.2]
	
    # utility parameters
    σ::Float64=2.0
	ρ::Float64=0.04
	r::Float64=0.03
	w::Float64=1.0
	
    aₘᵢₙ::Float64=-0.1
    aₘₐₓ::Float64=1.0
	an::Int=500
end

function EconPDEsFast(m::Moll, r; σ=m.σ, ρ=m.ρ, aₘᵢₙ=m.aₘᵢₙ, aₘₐₓ=m.aₘₐₓ, an = m.N_a)
	Λ = generator(m)
	EconPDEsFast(; σ, ρ, r, w=1, aₘᵢₙ, aₘₐₓ, an, zgrid=vec(m.z), Λ)
end
	
function (m::EconPDEsFast)(state::NamedTuple, value::NamedTuple)
    (; zgrid, Λ, ρ) = m
#	zs=zgrid
	nz = length(zgrid)
    (; a) = state
    (; v1, v1a_up, v1a_down, v2, v2a_up, v2a_down) = value

	T = eltype(value)
	
	vs = [v1, v2]
	dvf = [v1a_up, v2a_up]
	dvb = [v1a_down, v2a_down]
	states = [(; a, z=zgrid[1]), (; a, z=zgrid[2])]

	cs = Vector{T}(undef, nz)
	ȧs = Vector{T}(undef, nz)
	dvs = Vector{T}(undef, nz)

	for i ∈ 1:2
		(; c, ȧ, dv) = consumption_and_drift_upwind(states[i], dvf[i], dvb[i], m)
		cs[i] = c
		ȧs[i] = ȧ
		dvs[i] = dv
	end
		
	vts  = ρ .* vs .- (u.(cs, Ref(m)) .+ ȧs .* dvs .+ Λ * vs)

    return (; v1t=vts[1], v2t=vts[2]), (; s1=ȧs[1], s2=ȧs[2], c1=cs[1], c2=cs[2])
end

	pkgtest0 = let
		m = EconPDEsFast()

		#agrid = m.amin .+ range(0, (m.amax - m.amin)^0.8, length = m.an).^(1/0.8)
		agrid = range(m.aₘᵢₙ, m.aₘₐₓ, length=m.an)
		stategrid = OrderedDict(:a => agrid)
		
		yend = OrderedDict(
			(Symbol(:v, i) => [v₀((; z, a), m) for a ∈ agrid]) for (i, z) ∈ enumerate(m.zgrid)
		)

		(; m, stategrid, yend)
	end
end

# ╔═╡ 6213dcf6-cd7a-47cd-95d9-316ba538f4c6
let
	(; m, stategrid, yend) = pkgtest0
	pdesolve(m, stategrid, yend)
	t = @elapsed result = pdesolve(m, stategrid, yend)
	@info t
	@assert result.residual_norm <= 1e-5
end

# ╔═╡ 6b6e49cb-d4b0-49c7-bc82-62bd3cbfabc2
function _solve_HJB_econpdes(M, r; crit = √eps(), v₀=v₀)
	#agrid = M.aₘᵢₙ .+ range(0, (M.aₘₐₓ - M.aₘᵢₙ)^0.8, length = M.an).^(1/0.8)
	agrid = range(M.aₘᵢₙ, M.aₘₐₓ, M.an)
	stategrid = OrderedDict(:a => agrid)

	solend = OrderedDict(
		(Symbol(:v, i) => [v₀((; z, a), M) for a ∈ agrid]) for (i, z) ∈ enumerate(M.zgrid)
	)

	(; residual_norm, optional) = pdesolve(M, stategrid, solend; maxdist = crit)

	(; residual_norm, optional, agrid, M.zgrid)
end

# ╔═╡ 66c76cdb-499d-4141-a32e-879468d91291
function solve_HJB_econpdes(m, r; crit = √eps(), v₀=v₀)
	M = EconPDEsFast(m, r)
	
	(; residual_norm, optional, agrid, zgrid) = _solve_HJB_econpdes(M, r; crit, v₀)

	df = optional_to_df(optional, agrid, zgrid)

	(; df.v, df.c, df.ȧ, df, dist=residual_norm)
end

# ╔═╡ b7024377-6ffe-4f8a-a58a-0fa809683fd2
abstract type Scheme end

# ╔═╡ 2517dcff-c838-4b33-9f0e-e6d4bd8a259e
Base.@kwdef struct Implicit <: Scheme
	Δ=1000
	maxit=100
end

# ╔═╡ 0a1c908e-6432-4e8a-8c8c-7cc903fbca6c
begin
	Base.@kwdef struct Explicit <: Scheme
		Δ
		maxit
	end
	function Explicit((; Δa, z, aₘₐₓ); maxit=100_000)
		r = 0.02
		Δ = .9 * Δa / (maximum(z) + r * aₘₐₓ)
		Explicit(Δ, maxit)
	end
end

# ╔═╡ 3e4b82aa-a6ec-4097-afb5-927b7e16262a
function consumption_and_drift_upwind_vec(ss, dvf, dvb, par)
	consumption_and_drift_upwind.(ss, dvf, dvb, Ref(par)) |> StructArray
end

# ╔═╡ 12cbc2a7-ddf3-42c5-b592-a438ecdc5165
begin
	Base.@kwdef struct EconPDEsSlow
		σ::Float64 = 2.0
		ρ::Float64 = 0.05
		r::Float64 = 0.03
		w::Float64 = 1.0
		aₘᵢₙ::Float64 = -0.1
		aₘₐₓ::Float64 = 1.0
		an::Int = 500
		zgrid::Vector{Float64} =  [0.1, 0.2]
		Λ::Matrix{Float64} = [-0.02 0.02;
			0.03 -0.03]
	end

	function EconPDEsSlow(m::Moll, r; σ=m.σ, ρ=m.ρ, aₘᵢₙ=m.aₘᵢₙ, aₘₐₓ=m.aₘₐₓ, an = m.N_a)
		Λ = generator(m)
		EconPDEsSlow(; σ, ρ, r, w=1, aₘᵢₙ, aₘₐₓ, an, zgrid=vec(m.z), Λ)
	end
	
	function (M::EconPDEsSlow)(state::NamedTuple, sol::NamedTuple)
	 	(; zgrid, r, ρ, σ, Λ) = M
		(; dvf, dvb, v) = dvs_v(sol)

		(; a) = state

		(; dv, c, ȧ) = consumption_and_drift_upwind_vec(states.(a, zgrid), dvf, dvb, M)

		endo = u.(c, Ref(M)) .+ dv .* ȧ	

		vt = ρ * v - (endo + Λ * v)

		nz = length(zgrid)
		return NamedTuple{Tuple(Symbol.(:v, 1:nz, :t))}(vt),
			NamedTuple{Tuple([Symbol.(:c, 1:nz); Symbol.(:s, 1:nz)])}([c; ȧ])
	end	

	pkgtest = let
		m = EconPDEsSlow()

		#agrid = m.amin .+ range(0, (m.amax - m.amin)^0.8, length = m.an).^(1/0.8)
		agrid = range(m.aₘᵢₙ, m.aₘₐₓ, length=m.an)
		stategrid = OrderedDict(:a => agrid)
		
		yend = OrderedDict(
			(Symbol(:v, i) => [v₀((; z, a), m) for a ∈ agrid]) for (i, z) ∈ enumerate(m.zgrid)
		)

		(; m, stategrid, yend)
	end
end

# ╔═╡ 03a5eeac-ebb0-426b-a17d-3ef3e69627f6
let
	(; m, stategrid, yend) = pkgtest
	pdesolve(m, stategrid, yend)
	t = @elapsed result = pdesolve(m, stategrid, yend)
	@info t
	@assert result.residual_norm <= 1e-5
end

# ╔═╡ 90d6eba4-fa47-4717-b742-8b38f0b330ac
function inner_function(ss, v, (; aₘₐₓ, aₘᵢₙ, r, σ, N_a, Δa, z))
	# initialize arrays for forward and backward difference
	dvf = zeros(N_a, 2)
	dvb = zeros(N_a, 2)

	# forward difference
	dvf[1:N_a-1,:] .= (v[2:N_a,:] - v[1:N_a-1,:]) / Δa
	#dvf[N_a,:] .= (vec(z) .+ r * aₘₐₓ) .^ (-σ) # boundary condition a  <= a_max

	# backward difference
	dvb[2:N_a,:] .= (v[2:N_a,:] - v[1:N_a-1,:]) / Δa
	#dvb[1,:] .= (vec(z) .+ r * aₘᵢₙ) .^ (-σ) # boundary condition a >= a_min
	
	#I_concave = dvb .> dvf # problems if value function not concave

	out = consumption_and_drift_upwind_vec(ss, dvf, dvb, (; σ, r, aₘᵢₙ, aₘₐₓ))	
	u = out.c .^ (1-σ)/(1-σ)

	(; u, out.c, out.dv, out.ȧ, out.ȧf, out.ȧb)
end

# ╔═╡ 005988d8-a327-410e-9e53-09ba99dcee9a
function update_v(ss, v, par, _, Λ, (; Δ)::Explicit)
	(; ρ) = par
	(; u, c, ȧ, ȧf, ȧb, dv) = inner_function(ss, v, par)
		
	# HJB equation
	v_change = u + dv .* ȧ + v * Λ' - ρ*v

	v_new = v + Δ * v_change
	(; v_new, v_change, u, c, ȧ, ȧf, ȧb)
end

# ╔═╡ b57ea0a2-bee4-450f-99fa-cc79054fc963
function construct_A_alt(ȧfs, ȧbs, da, N_a, N_z)
	T = typeof((; I_from=1, I_to=1, λ=0.0))
	list = T[]
	
	car_inds = CartesianIndices((N_a, N_z)) |> collect .|> Tuple
	lin_inds = LinearIndices((N_a, N_z))

	for (I_from, (ȧf, ȧb)) in enumerate(zip(ȧfs, ȧbs))
		i_a, i_z = car_inds[I_from]

		if ȧf > 0 && i_a < N_a
			I_to = lin_inds[i_a + 1, i_z]
			λ = ȧf / da
			push!(list, (; I_from, I_to, λ))
		end
		if ȧb < 0 && i_a > 1
			I_to = lin_inds[i_a - 1, i_z]
			λ = -ȧb / da
			push!(list, (; I_from, I_to, λ))
		end
	end
	
	sa = StructArray(list)
	NN = N_a * 2
	A₀ = sparse(sa.I_from, sa.I_to, sa.λ, NN, NN)
	for i ∈ 1:NN
		A₀[i, i] = - sum(A₀[i,:])
	end

	A₀
end

# ╔═╡ 8777505c-6db4-4e75-a57b-59f4620b0051
function construct_A(ȧs, da, N_a, N_z)


	size_ss = (N_a, N_z)
	N = N_a * N_z
	car_inds = CartesianIndices(size_ss) |> collect .|> Tuple
	lin_inds = LinearIndices(size_ss)

	# initialize list of entries A
	T = typeof((; I_from=1, I_to=1, λ=0.0))
	list = T[]

	# create list of entries of A
	for (I_from, ȧ) in enumerate(ȧs)
		i_a, i_z = car_inds[I_from]

		i_a_next = nothing
		λ = nothing

		# in which direction to move?
		if ȧ > 0 && i_a < N_a
			i_a_next = i_a + 1
			λ = ȧ / da
		elseif ȧ < 0 && i_a > 1
			i_a_next = i_a - 1
			λ = -ȧ / da
		end

		if !isnothing(i_a_next)
			I_to = lin_inds[i_a_next, i_z]
			push!(list, (; I_from, I_to, λ))
		end
	end

	# construct sparse matrix from list of entries
	list_sa = StructArray(list)
	A = sparse(list_sa.I_from, list_sa.I_to, list_sa.λ, N, N)
	# fill diagonal
	A[diagind(A)] .= - vec(sum(A, dims=2))

	A
end

# ╔═╡ 434ecb8d-eae2-4913-97f8-3bdbefbf78ff
function construct_A_diag(ȧf, ȧb, da, N_a)
	X = - min.(ȧb,0)/da
	Z =   max.(ȧf,0)/da

	A11 = spdiagm(-1 => X[2:N_a,1], 1 => Z[1:N_a-1,1])
	A22 = spdiagm(-1 => X[2:N_a,2], 1 => Z[1:N_a-1,2])
	A = cat(A11, A22, dims=(1,2))
	A[diagind(A)] .= - vec(sum(A, dims=2))
	
	A
end

# ╔═╡ b77b52ac-73ad-456c-94e5-5f72e310268f
function construct_A_moll(ȧf, ȧb, da, N_a)
	X = - min.(ȧb,0)/da
	Y = - max.(ȧf,0)/da + min.(ȧb,0)/da
	Z =   max.(ȧf,0)/da

	A11 = spdiagm(-1 => X[2:N_a,1], 0 => Y[:,1], 1 => Z[1:N_a-1,1])
	A22 = spdiagm(-1 => X[2:N_a,2], 0 => Y[:,2], 1 => Z[1:N_a-1,2])
	A = cat(A11, A22, dims=(1,2))

	A
end

# ╔═╡ 4d7ee33f-ff78-4ea9-838b-a8320df4651f
function solve_HJB_implicit(m::Moll, r; maxit = 100, crit = √eps(), Δ = 1000)

	(; σ, ρ, z, λ, N_a, aₘᵢₙ, aₘₐₓ, Δa) = m
	da = Δa
	N_z = length(z)
	
	# construct asset grid
	a = range(aₘᵢₙ, aₘₐₓ, N_a)

	# initialize arrays for forward and backward difference
	dvf = zeros(N_a, N_z)
	dvb = zeros(N_a, N_z)

	# precompute A_switch matrix
	id = sparse(I, N_a, N_a)
	A_switch_1 = hcat(-λ[1] * id,  λ[1] * id)
	A_switch_2 = hcat( λ[2] * id, -λ[2] * id)
	A_switch   = vcat(A_switch_1, A_switch_2)

	# initial guess for value function
	v₀ = zeros(N_a, N_z)
	for (i, zᵢ) in enumerate(z)
		v₀[:,i] = (zᵢ .+ r * a).^(1-σ) / (1-σ) / ρ
	end
	v = v₀

	# initialize vector that keeps track of convergence
	dist = - ones(maxit)

	for it in range(1, maxit)

		# STEP 1

		# forward difference
		dvf[1:N_a-1,:] = (v[2:N_a,:] - v[1:N_a-1,:]) / da
		dvf[N_a,:] = (z .+ r * aₘₐₓ) .^ (-σ) # boundary condition a  <= a_max

		# backward difference
		dvb[2:N_a,:] = (v[2:N_a,:] - v[1:N_a-1,:]) / da
		dvb[1,:] = (z .+ r * aₘᵢₙ) .^ (-σ) # boundary condition a >= a_min
	
		I_concave = dvb .> dvf # problems if value function not concave

		# consumption and savings with forward difference
		cf = dvf .^ (-1/σ)
		ȧf = z .+ r .* a - cf

		# consumption and savings with backward difference
		cb = dvb .^ (-1/σ)
		ȧb = z .+ r .* a - cb

		# consumption and derivate of value function at steady state
		c0 = z .+ r .* a
		dv0 = c0 .^ (-σ)

		If = ȧf .> 0 # positive drift => forward difference
		Ib = ȧb .< 0 # negative drift => backward difference
		Ib[N_a,:] .= 1. # make sure backward difference is used at last grid point
		If[N_a,:] .= 0.
		I0 = (1 .- If .- Ib) # steady state
	
		dv_upwind = dvf.*If + dvb.*Ib + dv0.*I0

		# STEP 2
	
		c = dv_upwind .^ (-1/σ)
		u = c.^(1-σ)/(1-σ)

		# STEP 3
		A = construct_A_moll(ȧf, ȧb, Δa, N_a) + A_switch

		B = (ρ + 1/Δ) * I - A
		b = vec(u) + vec(v)/Δ
		v_new_stacked = B \ b
		v_new = reshape(v_new_stacked, N_a, N_z)

		# STEP 4

		v_change = v_new - v
		dist[it] = maximum(abs.(v_change))
		
		v = v_new

		if dist[it] < crit
			
			ȧ = z .+ r.*a - c

			ss = tuple.(a, z)
			return (; v, c, ȧ, a=first.(ss), z=last.(ss), A, it_last=it, dist)
			
		end

	end

	@error("Algorithm did not converge")

end

# ╔═╡ e1376e99-a636-4d28-b711-2dd4be66374f
function solve_df(m::Moll, r; maxit = 100, crit = 1e-6, Δ = 1000)
	(; v, c, ȧ, A, it_last, dist) = solve_HJB_implicit(m, r; crit, Δ)
	# scheme = Implicit(; maxit, Δ)
	#(; v, c, ȧ, A, it_last, dist) = solve_HJB_julian(m, r, scheme; crit)
	g = solve_KF(m, A)
	df = results_to_df(m; v, c, ȧ, g)
	
	return df, it_last, dist
	
end

# ╔═╡ 4de350bb-d2be-46b1-8558-a49f54a29dd1
(df, it_last, dist) = solve_df(m, 0.03; maxit = 100);

# ╔═╡ 1503b0d4-6d0e-4b5d-896f-e13c093ad3d4
let 
	
	df.g_max = min.(df.g,df.g[2]);
	
	@chain df begin
		stack(Not([:a, :z, :g]))
		data(_) * mapping(:a, :value, layout = :variable, color = :z => nonnumeric) * visual(Lines)
		draw(; facet = (linkyaxes = false, ), legend = (position = :top, titleposition = :left))
	end

end

# ╔═╡ 19e28c66-a389-4903-82a5-c963cf0b90b9
function excess_demand(m::Moll, r; maxit = 100, crit = 1e-6, Δ = 1000)
	(; Δa) = m
	(df, it_last, dist) = solve_df(m, r; maxit, crit, Δ)
	A = dot(df.a, df.g) * Δa

end

# ╔═╡ bd353706-be2d-480d-8ebb-cf20cab0dbec
r_eq = find_zero(r -> excess_demand(m, r), initial_bracket, Brent())

# ╔═╡ 3a92eac0-643f-47ec-b6dc-6d541873ac9a
solve_HJB_implicit(m, 0.03)

# ╔═╡ bf01f8b6-b618-4ce8-9ff0-3bdf37901ce1
function update_v(ss, v, par, A_switch, _, (; Δ)::Implicit)
	(; ρ, N_a, Δa) = par
	(; u, c, ȧ, ȧf, ȧb, dv) = inner_function(ss, v, par)
		
	#A = construct_A_alt(ȧf, ȧb, da, N_a, 2) + A_switch
	#A = construct_A_diag(ȧf, ȧb, da, N_a) + A_switch
	A = construct_A_moll(ȧf, ȧb, Δa, N_a) + A_switch
	# A = construct_A(ȧ, da, N_a, 2) + A_switch

	B = (ρ + 1/Δ) * I - A
	b = vec(u) + vec(v)/Δ
	v_new_stacked = B \ b
	v_new = reshape(v_new_stacked, N_a, 2)
	v_change = v_new - v
	(; v_new, v_change, u, c, ȧ, ȧf, ȧb)
end

# ╔═╡ e204ae15-fefd-4f01-8d5e-3772aefe9b0f
function solve_HJB_julian(m::Moll, r, scheme::Scheme; v₀=v₀, crit = 1e-6)

	(; maxit, Δ) = scheme
	
	(; σ, ρ, z, λ, N_a, aₘᵢₙ, aₘₐₓ, Δa) = m
	par = (; σ, ρ, z, N_a, aₘₐₓ, m.aₘᵢₙ, Δa, r)
	Λ = generator(m)
	
	# construct asset grid
	ss = statespace(m)
	
	# precompute A_switch matrix
	A_switch = construct_A_switch(m)

	# initial guess for value function
	v = v₀.(ss, Ref((; r, σ, ρ)))

	# initialize vector that keeps track of convergence
	dists = []

	for it in range(1, maxit)

		# updating the value function
		(; v_new, v_change, c, ȧ, ȧf, ȧb) = update_v(ss, v, par, A_switch, Λ, scheme)

		dist = maximum(abs.(v_change))
		v .= v_new
		push!(dists, dist)
		
		if dist < crit
			A = construct_A_moll(ȧf, ȧb, Δa, N_a) + construct_A_switch(m)
			return (; A, v, c, ȧ, it_last=it, dist=dists)
		end
	end

	@error("Algorithm did not converge")

end

# ╔═╡ 7a98402f-8b61-4228-bc11-33da201e6d82
compare = let
	r = 0.02
	crit = 1e-11
	maxit = 100_000
	scheme = Implicit()
	
	t1 = @elapsed HJB_moll     = solve_HJB_implicit(m, r; crit)
	@info t1
	#t2 = @elapsed HJB_explicit = solve_HJB_explicit(m, r; crit)
	#@info t2
	t3 = @elapsed HJB_julian   = solve_HJB_julian(m, r, scheme; crit)
	@info t3
	t4 = @elapsed HJB_econpdes = solve_HJB_econpdes(m, r; crit)
	@info t4

	@test vec(HJB_moll.v) ≈ vec(HJB_julian.v) ≈ HJB_econpdes.v
	@test vec(HJB_moll.ȧ) ≈ vec(HJB_julian.ȧ) ≈ HJB_econpdes.ȧ
	@test vec(HJB_moll.c) ≈ vec(HJB_julian.c) ≈ HJB_econpdes.c
	
	#@info HJB_moll.it_last, HJB_julian.it_last

	(; HJB_moll, #=HJB_explicit,=# HJB_julian, HJB_econpdes)
end

# ╔═╡ 708ba777-847b-460f-85e0-0615866c18eb
let
	(; HJB_moll, HJB_julian, HJB_econpdes) = compare

	i_z = CartesianIndices(HJB_moll.ȧ) .|> Tuple .|> last |> vec
	#i_a = CartesianIndices(HJB_moll.ȧ) .|> Tuple .|> first |> vec
	#a = 
	df_moll = DataFrame(; 
		ȧ=vec(HJB_moll.ȧ), 
		c=vec(HJB_moll.c), 
		v=vec(HJB_moll.v),
		a=vec(HJB_moll.a),
		z=vec(HJB_moll.z),
	)
	
	df_gomez= HJB_econpdes.df

	df = vcat(df_moll, df_gomez, source = :method => ["Moll", "Gomez"])

	@chain df begin
		stack(Not([:a, :z, :method]))
		data(_) * mapping(
			:a, :value,
			layout = :variable, color=:z => nonnumeric,
			linestyle = :method => nonnumeric
		) * visual(Lines)
		draw(_, facet = (linkyaxes=false, ))
		as_svg
	end
end

# ╔═╡ 9e1600df-9b6a-4b84-a152-1fa636f25fe5
let
	r = 0.03
	crit = 1e-6
	scheme = Implicit()#; maxit, Δ)
	(; v, c, ȧ, A, it_last, dist) = solve_HJB_julian(m, r, scheme; crit)

	N = size(A, 1)
	g₀ = fill(1/N, N)
	
	t0 = @elapsed g0 = solve_KF_iterate(A, 0.05)
	
	
	t1 = @elapsed g1 = solve_KF_moll(A)
	t2 = @elapsed g2 = solve_KF_death(A)
	t3 = @elapsed g3 = solve_KF_eigs(A)
	t4 = @elapsed g4 = gth_solve(Matrix(A))
	

	@info t0, t1, t2, t3, t4
	
	maximum(abs, g0 - g1), maximum(abs, g1 - g2), maximum(abs, g2 - g3),
	maximum(abs, g3 - g4), maximum(abs, g4 - g0)
	
	#norm(g1 - g_stacked)/norm(g1)
	#g1 ./ (sum(g1) * m.Δa)

#	fig = Figure()
#	ax = Axis(fig[1,1])
	
#	lines!(ax, g1)
#	lines!(ax, g_stacked)

#	fig
end

# ╔═╡ 98a7d8e0-1741-4447-b612-e40491fa8673
t_impl = @elapsed solve_HJB_julian(m, 0.03, Implicit())

# ╔═╡ ae8fd43c-f658-47b1-8781-6f699ade6bdb
let 
	df_dist = DataFrame(
		iteration = range(1, it_last), 
		time = range(0, t_impl, it_last),
		log10_dist = log10.(dist[1:it_last])
		)
	figure = (; resolution = (500, 250))
	@chain df_dist begin
		data(_) * mapping(:time, :log10_dist) * visual(Lines)
		draw(; figure)
		as_svg
	end
end

# ╔═╡ f8fbff0d-15d4-43ca-9f9c-29788ff793ec
function solve_explicit_df(m::Moll, r; maxit = 100000, crit = √eps())
	scheme = Explicit(m; maxit)
	(; v, c, ȧ, it_last, dist) = solve_HJB_julian(m, r, scheme; crit)
	df = results_to_df(m; v, c, ȧ)
	
	return df, it_last, dist
	
end

# ╔═╡ 4894a323-c8c8-4f34-b311-20c64176b89d
(df2, it_last2, dist2) = solve_explicit_df(m, 0.03; crit=√eps());

# ╔═╡ 264fc65e-d09a-4c67-94de-f845d42d18a3
let 
	df_dist = DataFrame(
		iteration = range(1, it_last2), 
		time = range(0, t_expl, it_last2),
		log10_dist = log10.(dist2[1:it_last2])
		)
	figure = (; resolution = (500, 250))
	@chain df_dist begin
		data(_) * mapping(:time, :log10_dist) * visual(Lines)
		draw(; figure) 
		as_svg
	end
end

# ╔═╡ 1cc42384-d262-496a-801a-033e7b40e6e9
md"""
### Differentiate
"""

# ╔═╡ 9f6051fa-2bd7-40b1-b2f8-5e7b1c477387
begin
	Δy_up(y, i, Δx) = (y[i+1] - y[i]) / Δx
	Δy_down(y, i, Δx) = (y[i] - y[i-1]) / Δx
	Δy_central(y, i, Δx) = (y[i+1] - y[i-1]) / Δx
	
	function Δgrid(grid, i)
	    last = length(grid)
	    @inbounds down = grid[max(i, 2)]      - grid[max(i-1, 1)]
	    @inbounds up   = grid[min(i+1, last)] - grid[min(i, last-1)]
	    central = (up + down)
	    avg = central / 2
	
	    (; up, down, avg, central)
	end
	
	function Δy(y, bc, i, Δx, fun_name, state_name)
	    up       = i != length(y) ? Δy_up(y, i, Δx.up)     : bc[i]
	    down     = i != 1         ? Δy_down(y, i, Δx.down) : bc[i]
	    second = (up - down) / Δx.avg
	    NamedTuple{deriv_names(fun_name, state_name)}((up, down, second))
	end
	
	deriv_names(fun_name, state_name) = (Symbol(fun_name, state_name, "_", :up), Symbol(fun_name, state_name, "_", :down), Symbol(fun_name, state_name, state_name))
end

# ╔═╡ 518eec58-eb4b-4294-9090-b6645f51a337
md"""
## Imported packages
"""

# ╔═╡ 5c1387ed-973c-4109-9ace-c473a4efe9ee
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
Arpack = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
EconPDEs = "a3315474-fad9-5060-8696-cee5f38a87b7"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
InfinitesimalGenerators = "2fce0c6f-5f0b-5c85-85c9-2ffe1d5ee30d"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuantEcon = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
AlgebraOfGraphics = "~0.6.18"
Arpack = "~0.5.4"
CairoMakie = "~0.11.8"
Chain = "~0.5.0"
DataFrameMacros = "~0.4.1"
DataFrames = "~1.6.1"
EconPDEs = "~1.0.3"
HypertextLiteral = "~0.9.5"
InfinitesimalGenerators = "~0.5.1"
MarkdownLiteral = "~0.1.1"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.55"
QuantEcon = "~0.16.6"
Roots = "~2.1.2"
StructArrays = "~0.6.17"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "a864c2ae73c2fedfe491b1a6e2f7f198a9faf1c0"

[[deps.ADTypes]]
git-tree-sha1 = "41c37aa88889c171f1300ceac1313c06e891d245"
uuid = "47edcb42-4c32-4615-8424-f2b9edc5f35b"
version = "0.2.6"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractLattices]]
git-tree-sha1 = "222ee9e50b98f51b5d78feb93dd928880df35f06"
uuid = "398f06c4-4d28-53ec-89ca-5b2656b7603d"
version = "0.3.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "c278dfab760520b8bb7e9511b968bf4ba38b7acc"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.2.3"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Accessors]]
deps = ["CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Test"]
git-tree-sha1 = "cb96992f1bec110ad211b7e410e57ddf7944c16f"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.35"

    [deps.Accessors.extensions]
    AccessorsAxisKeysExt = "AxisKeys"
    AccessorsIntervalSetsExt = "IntervalSets"
    AccessorsStaticArraysExt = "StaticArrays"
    AccessorsStructArraysExt = "StructArrays"

    [deps.Accessors.weakdeps]
    AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
    IntervalSets = "8197267c-284f-5f27-9208-e0e47529a953"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "0fb305e0253fd4e833d486914367a2ee2c2e78d0"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.1"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "PrecompileTools", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "3fbdee81b0cdc2b106b681dd2b9d4bdc60ca35a2"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.18"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArnoldiMethod]]
deps = ["LinearAlgebra", "Random", "StaticArrays"]
git-tree-sha1 = "62e51b39331de8911e4a7ff6f5aaf38a5f4cc0ae"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.2.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterface]]
deps = ["Adapt", "LinearAlgebra", "Requires", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "bbec08a37f8722786d87bedf84eae19c020c4efa"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "7.7.0"

    [deps.ArrayInterface.extensions]
    ArrayInterfaceBandedMatricesExt = "BandedMatrices"
    ArrayInterfaceBlockBandedMatricesExt = "BlockBandedMatrices"
    ArrayInterfaceCUDAExt = "CUDA"
    ArrayInterfaceGPUArraysCoreExt = "GPUArraysCore"
    ArrayInterfaceStaticArraysCoreExt = "StaticArraysCore"
    ArrayInterfaceTrackerExt = "Tracker"

    [deps.ArrayInterface.weakdeps]
    BandedMatrices = "aae01518-5342-5314-be14-df237901396f"
    BlockBandedMatrices = "ffab5731-97b5-5995-9138-79e8c1846df0"
    CUDA = "052768ef-5323-5732-b1bb-66c8b64840ba"
    GPUArraysCore = "46192b85-c4d5-4398-a991-12ede77f4527"
    StaticArraysCore = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.ArrayLayouts]]
deps = ["FillArrays", "LinearAlgebra"]
git-tree-sha1 = "64d582bcb9c93ac741234789eeb4f16812413efb"
uuid = "4c555306-a7a7-4459-81d9-ec55ddd5c99a"
version = "1.6.0"
weakdeps = ["SparseArrays"]

    [deps.ArrayLayouts.extensions]
    ArrayLayoutsSparseArraysExt = "SparseArrays"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Automa]]
deps = ["PrecompileTools", "TranscodingStreams"]
git-tree-sha1 = "588e0d680ad1d7201d4c6a804dcb1cd9cba79fbb"
uuid = "67c07d97-cdcb-5c2c-af73-a7f9c32a568b"
version = "1.0.3"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BandedMatrices]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra", "PrecompileTools"]
git-tree-sha1 = "931f3f49902e9b6b527fd7cd02d1cd7b4a84264c"
uuid = "aae01518-5342-5314-be14-df237901396f"
version = "1.5.0"
weakdeps = ["SparseArrays"]

    [deps.BandedMatrices.extensions]
    BandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "f1f03a9fa24271160ed7e73051fba3c1a759b53f"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.4.0"

[[deps.BlockArrays]]
deps = ["ArrayLayouts", "FillArrays", "LinearAlgebra"]
git-tree-sha1 = "fc69cbdb4277042f72c6e59cbc7024fbe3034b89"
uuid = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
version = "0.16.39"

[[deps.BlockBandedMatrices]]
deps = ["ArrayLayouts", "BandedMatrices", "BlockArrays", "FillArrays", "LinearAlgebra", "MatrixFactorizations"]
git-tree-sha1 = "b75b1edc92654ceb2bc3f3a4622d68e12fb2e32b"
uuid = "ffab5731-97b5-5995-9138-79e8c1846df0"
version = "0.12.9"
weakdeps = ["SparseArrays"]

    [deps.BlockBandedMatrices.extensions]
    BlockBandedMatricesSparseArraysExt = "SparseArrays"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.CEnum]]
git-tree-sha1 = "389ad5c84de1ae7cf0e28e381131c98ea87d54fc"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.5.0"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CRlibm]]
deps = ["CRlibm_jll"]
git-tree-sha1 = "32abd86e3c2025db5172aa182b982debed519834"
uuid = "96374032-68de-5a5b-8d9e-752f78720389"
version = "1.0.1"

[[deps.CRlibm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e329286945d0cfc04456972ea732551869af1cfc"
uuid = "4e9b3aee-d8a1-5a3d-ad8b-7d824db253f0"
version = "1.0.1+0"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["CRC32c", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "PrecompileTools"]
git-tree-sha1 = "a80d49ed3333f5f78df8ffe76d07e88cc35e9172"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.11.8"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "ad25e7d21ce10e01de973cdc68ad0f850a953c52"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.21.1"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.CodecBzip2]]
deps = ["Bzip2_jll", "Libdl", "TranscodingStreams"]
git-tree-sha1 = "9b1ca1aa6ce3f71b3d1840c538a8210a043625eb"
uuid = "523fee87-0ab8-5b00-afb7-3ecf72e48cfd"
version = "0.8.2"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "PrecompileTools", "URIs"]
git-tree-sha1 = "532c4185d3c9037c0237546d817858b23cf9e071"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.12"

[[deps.CommonSolve]]
git-tree-sha1 = "0eee5eb66b1cf62cd6ad1b460238e60e4b09400c"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.4"

[[deps.CommonSubexpressions]]
deps = ["MacroTools", "Test"]
git-tree-sha1 = "7b8a93dba8af7e3b42fecabf646260105ac373f7"
uuid = "bbf7d656-a473-5ed7-a52c-81e309532950"
version = "0.3.0"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "75bd5b6fc5089df449b5d35fa501c846c9b6549b"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.12.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "c53fc348ca4d40d7b371e71fd52251839080cbc9"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.4"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "f7f4319567fe769debfcf7f8c03d8da1dd4e2fb0"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.9"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataFrameMacros]]
deps = ["DataFrames", "MacroTools"]
git-tree-sha1 = "5275530d05af21f7778e3ef8f167fb493999eea1"
uuid = "75880514-38bc-4a95-a458-c2aea5a3a702"
version = "0.4.1"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "DataStructures", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrecompileTools", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "04c738083f29f86e62c8afc341f0967d8717bdb8"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.6.1"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "ac67408d9ddf207de5cfa9a97e114352430f01ed"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.16"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelaunayTriangulation]]
deps = ["DataStructures", "EnumX", "ExactPredicates", "Random", "SimpleGraphs"]
git-tree-sha1 = "d4e9dc4c6106b8d44e40cd4faf8261a678552c7c"
uuid = "927a84f5-c5f4-47a5-9785-b46e178433df"
version = "0.8.12"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "1f3b7b0d321641c1f2e519f7aed77f8e1f6cb133"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.29"

[[deps.DiffResults]]
deps = ["StaticArraysCore"]
git-tree-sha1 = "782dd5f4561f5d267313f23853baaaa4c52ea621"
uuid = "163ba53b-c6d8-5494-b064-1a9d43ac40c5"
version = "1.1.0"

[[deps.DiffRules]]
deps = ["IrrationalConstants", "LogExpFunctions", "NaNMath", "Random", "SpecialFunctions"]
git-tree-sha1 = "23163d55f885173722d1e4cf0f6110cdbaf7e272"
uuid = "b552c78f-8df3-52c6-915a-8e097449b14b"
version = "1.15.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "7c302d7a5fec5214eb8a5a4c466dcf7a51fcf169"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.107"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e3290f2d49e661fbd94046d7e3726ffcb2d41053"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.4+0"

[[deps.EconPDEs]]
deps = ["BlockBandedMatrices", "FiniteDiff", "LinearAlgebra", "NLsolve", "OrderedCollections", "Printf", "SparseArrays", "SparseDiffTools"]
git-tree-sha1 = "17d709798933f040b3258ba06bbec45c761348f7"
uuid = "a3315474-fad9-5060-8696-cee5f38a87b7"
version = "1.0.3"

[[deps.EnumX]]
git-tree-sha1 = "bdb1942cd4c45e3c678fd11569d5cccd80976237"
uuid = "4e289a0a-7415-4d19-859d-a7e5c4648b56"
version = "1.0.4"

[[deps.ExactPredicates]]
deps = ["IntervalArithmetic", "Random", "StaticArrays"]
git-tree-sha1 = "b3f2ff58735b5f024c392fde763f29b057e4b025"
uuid = "429591f6-91af-11e9-00e2-59fbe8cec110"
version = "2.2.8"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.Extents]]
git-tree-sha1 = "2140cd04483da90b2da7f99b2add0750504fc39c"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.2"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "c5c28c245101bd59154f649e19b038d15901b5dc"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.2"

[[deps.FilePaths]]
deps = ["FilePathsBase", "MacroTools", "Reexport", "Requires"]
git-tree-sha1 = "919d9412dbf53a2e6fe74af62a73ceed0bce0629"
uuid = "8fc22ac5-c921-52a6-82fd-178b2807b824"
version = "0.8.3"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "9f00e42f8d99fdde64d40c8ea5d14269a2e2c1aa"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.21"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random"]
git-tree-sha1 = "5b93957f6dcd33fc343044af3d48c215be2562f1"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.9.3"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FiniteDiff]]
deps = ["ArrayInterface", "LinearAlgebra", "Requires", "Setfield", "SparseArrays"]
git-tree-sha1 = "73d1214fec245096717847c62d389a5d2ac86504"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.22.0"
weakdeps = ["BandedMatrices", "BlockBandedMatrices", "StaticArrays"]

    [deps.FiniteDiff.extensions]
    FiniteDiffBandedMatricesExt = "BandedMatrices"
    FiniteDiffBlockBandedMatricesExt = "BlockBandedMatrices"
    FiniteDiffStaticArraysExt = "StaticArrays"

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
deps = ["CommonSubexpressions", "DiffResults", "DiffRules", "LinearAlgebra", "LogExpFunctions", "NaNMath", "Preferences", "Printf", "Random", "SpecialFunctions"]
git-tree-sha1 = "cf0fe81336da9fb90944683b8c41984b08793dad"
uuid = "f6369f11-7733-5829-9624-2563aa707210"
version = "0.10.36"
weakdeps = ["StaticArrays"]

    [deps.ForwardDiff.extensions]
    ForwardDiffStaticArraysExt = "StaticArrays"

[[deps.FreeType]]
deps = ["CEnum", "FreeType2_jll"]
git-tree-sha1 = "907369da0f8e80728ab49c1c7e09327bf0d6d999"
uuid = "b38be410-82b0-50bf-ab77-7b57e271db43"
version = "4.1.1"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FreeTypeAbstraction]]
deps = ["ColorVectorSpace", "Colors", "FreeType", "GeometryBasics"]
git-tree-sha1 = "055626e1a35f6771fe99060e835b72ca61a52621"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.1"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "273bd1cd30768a2fddfa3fd63bbc746ed7249e5f"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.9.0"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "ec632f177c0d990e64d955ccc1b8c04c485a0950"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.6"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "d4f85701f569584f2cff7ba67a137d03f0cfb7d0"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.3"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "Extents", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "5694b56ccf9d15addedc35e9a4ba9c317721b788"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.10"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "e94c92c7bf4819685eb80186d51c43e71d4afa17"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.76.5+0"

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

[[deps.Graphs]]
deps = ["ArnoldiMethod", "Compat", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "899050ace26649433ef1af25bc17a815b3db52b7"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.9.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "af13a277efd8a6e716d79ef635d5342ccb75be61"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.10.0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions"]
git-tree-sha1 = "f218fe3736ddf977e0e772bc9a586b2383da2685"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.23"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "2e4520d67b0cef90865b3ef727594d2a58e0e1f8"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.11"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "eb49b82c172811fd2c86759fa0553a2221feb909"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.7"

[[deps.ImageCore]]
deps = ["ColorVectorSpace", "Colors", "FixedPointNumbers", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "PrecompileTools", "Reexport"]
git-tree-sha1 = "b2a7eaa169c13f5bcae8131a83bc30eff8f71be0"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.10.2"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "bca20b2f5d00c4fbc192c3212da8fa79f4688009"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.7"

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "355e2b974f2e3212a75dfb60519de21361ad3cb7"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.9"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3d09a9f60edf77f8a4d99f9e015e8fbf9989605d"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.7+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.InfinitesimalGenerators]]
deps = ["Arpack", "Distributions", "FillArrays", "KrylovKit", "LinearAlgebra", "Roots"]
git-tree-sha1 = "c93b89437ba4a89a91691c4ec28ba3d5266614c8"
uuid = "2fce0c6f-5f0b-5c85-85c9-2ffe1d5ee30d"
version = "0.5.1"

[[deps.Inflate]]
git-tree-sha1 = "ea8031dea4aff6bd41f1df8f2fdfb25b33626381"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.4"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

[[deps.IntegerMathUtils]]
git-tree-sha1 = "b8ffb903da9f7b8cf695a8bead8e01814aa24b30"
uuid = "18e54dd8-cb9d-406c-a71d-865a43cbb235"
version = "0.1.2"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5fdf2fe6724d8caabf43b557b84ce53f3b7e2f6b"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.0.2+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

    [deps.Interpolations.weakdeps]
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.IntervalArithmetic]]
deps = ["CRlibm", "RoundingEmulator"]
git-tree-sha1 = "c274ec586ea58eb7b42afd0c5d67e50ff50229b5"
uuid = "d1acc4aa-44c8-5952-acd4-ba5d80a2a253"
version = "0.22.5"
weakdeps = ["DiffRules", "RecipesBase"]

    [deps.IntervalArithmetic.extensions]
    IntervalArithmeticDiffRulesExt = "DiffRules"
    IntervalArithmeticRecipesBaseExt = "RecipesBase"

[[deps.IntervalSets]]
git-tree-sha1 = "581191b15bcb56a2aa257e9c160085d0f128a380"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.9"
weakdeps = ["Random", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "68772f49f54b479fa88ace904f6127f0a3bb2e46"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.12"

[[deps.InvertedIndices]]
git-tree-sha1 = "0dc7b50b8d436461be01300fd8cd45aa0274b038"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.3.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.Isoband]]
deps = ["isoband_jll"]
git-tree-sha1 = "f9b6d97355599074dc867318950adaa6f9946137"
uuid = "f1662d9f-8043-43de-a69a-05efc1cc6ff4"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "fa6d0bcff8583bac20f1ffa708c3913ca605c611"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.5"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60b1194df0a3298f460063de985eae7b01bc011a"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.1+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "fee018a29b60733876eb557804b5b109dd3dd8a7"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.8"

[[deps.KrylovKit]]
deps = ["ChainRulesCore", "GPUArraysCore", "LinearAlgebra", "Printf"]
git-tree-sha1 = "5cebb47f472f086f7dd31fb8e738a8db728f1f84"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.6.1"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Lazy]]
deps = ["MacroTools"]
git-tree-sha1 = "1370f8202dac30758f3c345f9909b97f53d87d3f"
uuid = "50d2b5c4-7a5e-59d5-8109-a42b560f39c0"
version = "0.15.1"

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
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

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

[[deps.LightXML]]
deps = ["Libdl", "XML2_jll"]
git-tree-sha1 = "3a994404d3f6709610701c7dabfc03fed87a81f8"
uuid = "9c8b4983-aa76-5018-a973-4c85ecc9e179"
version = "0.9.1"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "7bbea35cec17305fc70a0e5b4641477dc0789d9d"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.2.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LinearAlgebraX]]
deps = ["LinearAlgebra", "Mods", "Primes", "SimplePolynomials"]
git-tree-sha1 = "d76cec8007ec123c2b681269d40f94b053473fcf"
uuid = "9b3f67b0-2d00-526e-9884-9e4938f8fb88"
version = "0.2.7"

[[deps.Loess]]
deps = ["Distances", "LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "a113a8be4c6d0c64e217b472fb6e61c760eb4022"
uuid = "4345ca2d-374a-55d4-8d30-97f9976e7612"
version = "0.6.3"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "7d6dd4e9212aebaeed356de34ccf262a3cd415aa"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.26"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl"]
git-tree-sha1 = "72dc3cf284559eb8f53aa593fe62cb33f83ed0c0"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.0.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Makie]]
deps = ["Animations", "Base64", "CRC32c", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "DelaunayTriangulation", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG_jll", "FileIO", "FilePaths", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalArithmetic", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MacroTools", "MakieCore", "Markdown", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "PrecompileTools", "Printf", "REPL", "Random", "RelocatableFolders", "Scratch", "ShaderAbstractions", "Showoff", "SignedDistanceFields", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "40c5dfbb99c91835171536cd571fe6f1ba18ff97"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.20.7"

[[deps.MakieCore]]
deps = ["Observables", "REPL"]
git-tree-sha1 = "248b7a4be0f92b497f7a331aed02c1e9a878f46b"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.7.3"

[[deps.MappedArrays]]
git-tree-sha1 = "2dab0221fe2b0f2cb6754eaa743cc266339f527e"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.2"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "0d3fa2dd374934b62ee16a4721fe68c418b92899"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.1"

[[deps.MathOptInterface]]
deps = ["BenchmarkTools", "CodecBzip2", "CodecZlib", "DataStructures", "ForwardDiff", "JSON", "LinearAlgebra", "MutableArithmetics", "NaNMath", "OrderedCollections", "PrecompileTools", "Printf", "SparseArrays", "SpecialFunctions", "Test", "Unicode"]
git-tree-sha1 = "8b40681684df46785a0012d352982e22ac3be59e"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.25.2"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "UnicodeFun"]
git-tree-sha1 = "96ca8a313eb6437db5ffe946c457a401bbb8ce1d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.7"

[[deps.MatrixFactorizations]]
deps = ["ArrayLayouts", "LinearAlgebra", "Printf", "Random"]
git-tree-sha1 = "78f6e33434939b0ac9ba1df81e6d005ee85a7396"
uuid = "a3b82374-2e81-5b9e-98ce-41277c0e4c87"
version = "2.1.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mods]]
git-tree-sha1 = "924f962b524a71eef7a21dae1e6853817f9b658f"
uuid = "7475f97c-0381-53b1-977b-4c60186c8d62"
version = "2.2.4"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.Multisets]]
git-tree-sha1 = "8d852646862c96e226367ad10c8af56099b4047e"
uuid = "3b2b4ff1-bcff-5658-a3ee-dbcf1ce5ac09"
version = "0.4.4"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "806eea990fb41f9b36f1253e5697aa645bf6a9f8"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.4.0"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "a0b464d183da839699f4c79e7606d9d186ec172c"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.3"

[[deps.NLopt]]
deps = ["NLopt_jll"]
git-tree-sha1 = "d1d09c342c3dd9b3bae985b088bd928632e4d79e"
uuid = "76087f3c-5699-56af-9a33-bf431cd00edd"
version = "1.0.1"
weakdeps = ["MathOptInterface"]

    [deps.NLopt.extensions]
    NLoptMathOptInterfaceExt = ["MathOptInterface"]

[[deps.NLopt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9b1f15a08f9d00cdb2761dcfa6f453f5d0d6f973"
uuid = "079eb43e-fd8e-5478-9966-2cf3e3edb778"
version = "2.7.1+0"

[[deps.NLsolve]]
deps = ["Distances", "LineSearches", "LinearAlgebra", "NLSolversBase", "Printf", "Reexport"]
git-tree-sha1 = "019f12e9a1a7880459d0173c182e6a99365d7ac1"
uuid = "2774e3e8-f4cf-5e23-947b-6d7e65073b56"
version = "4.5.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "d92b107dbb887293622df7697a2223f9f8176fcd"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "6a731f2b5c03157418a20c12195eb4b74c8f8621"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.13.0"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "60e3045590bd104a16fefb12836c00c0ef8c7f8c"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "MathOptInterface", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "d024bfb56144d947d4fafcd9cb5cafbe3410b133"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.9.2"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "67186a2bc9a90f9f85ff3cc8277868961fb57cbd"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.4.3"

[[deps.PackageExtensionCompat]]
git-tree-sha1 = "fb28e33b8a95c4cee25ce296c817d89cc2e53518"
uuid = "65ce6f38-6b18-4e1d-a461-8949797d7930"
version = "1.0.2"
weakdeps = ["Requires", "TOML"]

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "0fac6313486baae819364c52b4f483450a9d793f"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.12"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4745216e94f71cb768d58330b059c9b76f32cb66"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.14+0"

[[deps.Parameters]]
deps = ["OrderedCollections", "UnPack"]
git-tree-sha1 = "34c0e9ad262e5f7fc75b10a9952ca7692cfc5fbe"
uuid = "d96e819e-fc66-5662-9728-84c9c7592b0a"
version = "0.12.3"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Permutations]]
deps = ["Combinatorics", "LinearAlgebra", "Random"]
git-tree-sha1 = "eb3f9df2457819bf0a9019bd93cc451697a0751e"
uuid = "2ae35dd2-176d-5d53-8349-f30d82d94d4f"
version = "0.4.20"

[[deps.PikaParser]]
deps = ["DocStringExtensions"]
git-tree-sha1 = "d6ff87de27ff3082131f31a714d25ab6d0a88abf"
uuid = "3bbf5609-3e7b-44cd-8549-7c69f321e792"
version = "0.6.1"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f9501cc0430a26bc3d156ae1b5b0c1b47af4d6da"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.3"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "862942baf5663da528f66d24996eb6da85218e76"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.0"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "68723afdb616445c6caaef6255067a8339f91325"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.55"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "RecipesBase", "Setfield", "SparseArrays"]
git-tree-sha1 = "a9c7a523d5ed375be3983db190f6a5874ae9286d"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "4.0.6"
weakdeps = ["ChainRulesCore", "FFTW", "MakieCore", "MutableArithmetics"]

    [deps.Polynomials.extensions]
    PolynomialsChainRulesCoreExt = "ChainRulesCore"
    PolynomialsFFTWExt = "FFTW"
    PolynomialsMakieCoreExt = "MakieCore"
    PolynomialsMutableArithmeticsExt = "MutableArithmetics"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "36d8b4b899628fb92c2749eb488d884a926614d3"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.3"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "03b4c25b43cb84cee5c90aa9b5ea0a78fd848d2f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.0"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "00805cd429dcb4870060ff49ef443486c262e38e"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.1"

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "88b895d13d53b5577fd53379d913b9ab9ac82660"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.1"

[[deps.Primes]]
deps = ["IntegerMathUtils"]
git-tree-sha1 = "1d05623b5952aed1307bf8b43bec8b8d1ef94b6e"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.5"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "00099623ffee15972c16111bcf84c58a0051257c"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.9.0"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "9b23c31e76e333e6fb4c1595ae6afa74966a729e"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.9.4"

[[deps.QuantEcon]]
deps = ["DSP", "DataStructures", "Distributions", "FFTW", "Graphs", "LinearAlgebra", "Markdown", "NLopt", "Optim", "Pkg", "Primes", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "034293b29fdbcae73aeb7ca0b2755e693f04701b"
uuid = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
version = "0.16.6"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.RingLists]]
deps = ["Random"]
git-tree-sha1 = "f39da63aa6d2d88e0c1bd20ed6a3ff9ea7171ada"
uuid = "286e9d63-9694-5540-9e3c-4e6708fa07b2"
version = "0.2.8"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6ed52fdd3382cf21947b15e8870ac0ddbff736da"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.0+0"

[[deps.Roots]]
deps = ["Accessors", "ChainRulesCore", "CommonSolve", "Printf"]
git-tree-sha1 = "754acd3031a9f2eaf6632ba4850b1c01fe4460c1"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.1.2"

    [deps.Roots.extensions]
    RootsForwardDiffExt = "ForwardDiff"
    RootsIntervalRootFindingExt = "IntervalRootFinding"
    RootsSymPyExt = "SymPy"
    RootsSymPyPythonCallExt = "SymPyPythonCall"

    [deps.Roots.weakdeps]
    ForwardDiff = "f6369f11-7733-5829-9624-2563aa707210"
    IntervalRootFinding = "d2bf35a9-74e0-55ec-b149-d360ff49b807"
    SymPy = "24249f21-da20-56a4-8eb1-6a02cf4ae2e6"
    SymPyPythonCall = "bc8888f7-b21e-4b7c-a06a-5d9c9496438c"

[[deps.RoundingEmulator]]
git-tree-sha1 = "40b9edad2e5287e05bd413a38f61a8ff55b9557b"
uuid = "5eaf0fd0-dfba-4ccb-bf02-d820a40db705"
version = "0.2.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SciMLOperators]]
deps = ["ArrayInterface", "DocStringExtensions", "Lazy", "LinearAlgebra", "Setfield", "SparseArrays", "StaticArraysCore", "Tricks"]
git-tree-sha1 = "51ae235ff058a64815e0a2c34b1db7578a06813d"
uuid = "c0aeaf25-5076-4817-a8d5-81caf7dfa961"
version = "0.3.7"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "0e7508ff27ba32f26cd459474ca2ede1bc10991f"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.ShaderAbstractions]]
deps = ["ColorTypes", "FixedPointNumbers", "GeometryBasics", "LinearAlgebra", "Observables", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "79123bc60c5507f035e6d1d9e563bb2971954ec8"
uuid = "65257c39-d410-5151-9873-9b3e5be5013e"
version = "0.4.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.ShiftedArrays]]
git-tree-sha1 = "503688b59397b3307443af35cd953a13e8005c16"
uuid = "1277b4bf-5013-50f5-be3d-901d8477a67a"
version = "2.0.0"

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

[[deps.SimpleGraphs]]
deps = ["AbstractLattices", "Combinatorics", "DataStructures", "IterTools", "LightXML", "LinearAlgebra", "LinearAlgebraX", "Optim", "Primes", "Random", "RingLists", "SimplePartitions", "SimplePolynomials", "SimpleRandom", "SparseArrays", "Statistics"]
git-tree-sha1 = "f65caa24a622f985cc341de81d3f9744435d0d0f"
uuid = "55797a34-41de-5266-9ec1-32ac4eb504d3"
version = "0.8.6"

[[deps.SimplePartitions]]
deps = ["AbstractLattices", "DataStructures", "Permutations"]
git-tree-sha1 = "e9330391d04241eafdc358713b48396619c83bcb"
uuid = "ec83eff0-a5b5-5643-ae32-5cbf6eedec9d"
version = "0.3.1"

[[deps.SimplePolynomials]]
deps = ["Mods", "Multisets", "Polynomials", "Primes"]
git-tree-sha1 = "7063828369cafa93f3187b3d0159f05582011405"
uuid = "cc47b68c-3164-5771-a705-2bc0097375a0"
version = "0.2.17"

[[deps.SimpleRandom]]
deps = ["Distributions", "LinearAlgebra", "Random"]
git-tree-sha1 = "3a6fb395e37afab81aeea85bae48a4db5cd7244a"
uuid = "a6525b86-64cd-54fa-8f65-62fc48bdc0e8"
version = "0.3.1"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "2da10356e31327c7096832eb9cd86307a50b1eb6"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SparseDiffTools]]
deps = ["ADTypes", "Adapt", "ArrayInterface", "Compat", "DataStructures", "FiniteDiff", "ForwardDiff", "Graphs", "LinearAlgebra", "PackageExtensionCompat", "Random", "Reexport", "SciMLOperators", "Setfield", "SparseArrays", "StaticArrayInterface", "StaticArrays", "Tricks", "UnPack", "VertexSafeGraphs"]
git-tree-sha1 = "3b38ae7a1cbe9b8b1344359599753957644b03d4"
uuid = "47a9eef4-7e08-11e9-0b38-333d64bd3804"
version = "2.16.0"

    [deps.SparseDiffTools.extensions]
    SparseDiffToolsEnzymeExt = "Enzyme"
    SparseDiffToolsPolyesterForwardDiffExt = "PolyesterForwardDiff"
    SparseDiffToolsSymbolicsExt = "Symbolics"
    SparseDiffToolsZygoteExt = "Zygote"

    [deps.SparseDiffTools.weakdeps]
    Enzyme = "7da242da-08ed-463a-9acd-ee780be4f1d9"
    PolyesterForwardDiff = "98d1487c-24ca-40b6-b7ab-df2af84e126b"
    Symbolics = "0c5d862f-8b57-4792-8d23-62f2024744c7"
    Zygote = "e88e6eb3-aa80-5325-afca-941959d7151f"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "e2cfc4012a19088254b3950b85c3c1d8882d864d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.3.1"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.StableHashTraits]]
deps = ["Compat", "PikaParser", "SHA", "Tables", "TupleTools"]
git-tree-sha1 = "662f56ffe22b3985f3be7474f0aecbaf214ecf0f"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "1.1.6"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "b366eb1eb68075745777d80861c6706c33f588ae"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.8.9"

[[deps.StaticArrayInterface]]
deps = ["ArrayInterface", "Compat", "IfElse", "LinearAlgebra", "PrecompileTools", "Requires", "SparseArrays", "Static", "SuiteSparse"]
git-tree-sha1 = "5d66818a39bb04bf328e92bc933ec5b4ee88e436"
uuid = "0d7ed370-da01-4f52-bd93-41d350b8b718"
version = "1.5.0"
weakdeps = ["OffsetArrays", "StaticArrays"]

    [deps.StaticArrayInterface.extensions]
    StaticArrayInterfaceOffsetArraysExt = "OffsetArrays"
    StaticArrayInterfaceStaticArraysExt = "StaticArrays"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "7b0e9c14c624e435076d19aea1e5cbdec2b9ca37"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.2"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "36b3d696ce6366023a0ea192b4cd442268995a0d"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.2"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "1d77abd07f617c4868c33d4f5b9e1dbb2643c9cf"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.2"

[[deps.StatsFuns]]
deps = ["HypergeometricFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"
weakdeps = ["ChainRulesCore", "InverseFunctions"]

    [deps.StatsFuns.extensions]
    StatsFunsChainRulesCoreExt = "ChainRulesCore"
    StatsFunsInverseFunctionsExt = "InverseFunctions"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsAPI", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "5cf6c4583533ee38639f73b880f35fc85f2941e0"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.3"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.StructArrays]]
deps = ["Adapt", "ConstructionBase", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "1b0b1205a56dc288b71b1961d48e351520702e24"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.17"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "cb76cf677714c095e535e3501ac7954732aeea2d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.11.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "34cc045dd0aaa59b8bbe86c644679bc57f1d5bd0"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.8"

[[deps.TranscodingStreams]]
git-tree-sha1 = "54194d92959d8ebaa8e26227dbe3cdefcdcd594f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.3"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "155515ed4c4236db30049ac1495e2969cc06be9d"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.4.3"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

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

[[deps.VertexSafeGraphs]]
deps = ["Graphs"]
git-tree-sha1 = "8351f8d73d7e880bfc042a8b6922684ebeafb35c"
uuid = "19fa3120-7c27-5ec5-8db8-b0b0aa330d6f"
version = "0.2.0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "801cbe47eae69adc50f36c3caec4758d2650741b"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.2+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

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
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

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
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "93284c28274d9e75218a416c65ec49d0e0fcdf3d"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.40+0"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "libpng_jll"]
git-tree-sha1 = "d4f63314c8aa1e48cd22aa0c17ed76cd1ae48c3c"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.10.3+0"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

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
# ╟─cdd24518-1cd3-43af-9a88-f0a8cd3cd6ed
# ╟─0edf9214-7213-4b1d-8fa6-54a365525d29
# ╟─1c23af35-64f5-4b96-9eab-98ee9b1ecb4c
# ╠═64d44910-988c-4e24-baf6-6408da65bd21
# ╠═95764a2a-0a92-447c-98a8-df22167abda6
# ╠═a5233ee4-ba60-4428-8de5-4b012ed43408
# ╟─ac27dd69-1216-4a02-ba62-ebb098b33fa1
# ╠═3dd0c186-04f0-425c-93dc-825c7d4b606e
# ╠═23a991be-7c8b-45e2-bd75-af5e146fc6b0
# ╟─2289c7a7-3493-4bfb-8ffe-31b074d17b14
# ╟─0caee0bf-77c0-4542-a68e-9052d230ca74
# ╟─38555a5e-5d39-4c61-9381-94c8fb67a257
# ╠═4d7ee33f-ff78-4ea9-838b-a8320df4651f
# ╟─fb88f286-0bce-4cb0-9d52-83b373f6fbcf
# ╠═5cbd313e-b357-4b00-bb7d-4c87f8522cdd
# ╟─e507dfb0-08a4-45a8-b342-09ca989e05f5
# ╠═6213dcf6-cd7a-47cd-95d9-316ba538f4c6
# ╠═e12a5f08-7e5b-418c-a08f-8d179e912774
# ╠═12cbc2a7-ddf3-42c5-b592-a438ecdc5165
# ╠═69583d36-6190-4483-87fb-291d9cf82c77
# ╠═db13a6b9-2de2-47ce-8af0-99d868183354
# ╠═03a5eeac-ebb0-426b-a17d-3ef3e69627f6
# ╠═56234834-6de9-4e44-95af-c5e3fb828a9d
# ╠═a10ef170-045c-439e-a23b-122e812853aa
# ╠═68aeb6cf-cff8-49e9-8c39-d94aadd5a444
# ╠═e197b24a-7e51-45e4-8186-11f88bf48de6
# ╠═7734c75e-5f2b-4807-8b9a-82e7e010edac
# ╠═e57a2dfa-0937-49cf-9161-fa1e39fb5e80
# ╠═91c8dce8-07b8-46d2-8c61-e6bccced64e4
# ╠═66c76cdb-499d-4141-a32e-879468d91291
# ╠═ad4b46b0-b6fb-4711-9501-0912eda3d6c3
# ╠═6b6e49cb-d4b0-49c7-bc82-62bd3cbfabc2
# ╠═aeefc434-3375-400e-a360-4dc263b96252
# ╟─16b48db6-ef8d-4d7b-9e9f-b281143d0a90
# ╠═7a98402f-8b61-4228-bc11-33da201e6d82
# ╠═708ba777-847b-460f-85e0-0615866c18eb
# ╠═02c84fbb-f2cb-4051-8bda-62f5e2a4c6d1
# ╠═d704e278-6ae8-4166-bade-599e71e54336
# ╟─1fbf9f02-7ad2-4ddd-ac2c-92af79c9ed02
# ╟─2e09462e-6030-4e33-bc7a-9d2faed4ca74
# ╠═4ebeb5e8-505d-4ed2-8d37-4a64ceb7d796
# ╟─6c895dce-9c0a-499c-9e95-45c29ee5a459
# ╠═a8279bcc-4cdf-4a89-bb75-e7e1f617643d
# ╠═dd89d622-6066-4d52-a4ca-366822c63661
# ╠═33a6780e-7f5a-4677-a51c-c002b7b86dc3
# ╠═eb565042-1823-4d5a-b4d1-ee314dccd4e0
# ╟─e5cfdbc2-f592-40bb-ba5d-07f455dd5bd4
# ╠═34779464-624a-446d-b741-81524509aee6
# ╠═a2da91e4-cedd-446f-81d4-adb008c01d5b
# ╠═d3612723-d4eb-4887-b902-fa47e6a26226
# ╠═5bc621c2-1d87-49a9-8dc6-d81d5f998640
# ╠═9e1600df-9b6a-4b84-a152-1fa636f25fe5
# ╠═abff0506-1792-420d-80d8-cb3645483443
# ╠═eade1b48-7d1c-4d17-aaf5-14d111a13556
# ╠═e1376e99-a636-4d28-b711-2dd4be66374f
# ╠═4de350bb-d2be-46b1-8558-a49f54a29dd1
# ╠═98a7d8e0-1741-4447-b612-e40491fa8673
# ╠═3a92eac0-643f-47ec-b6dc-6d541873ac9a
# ╠═1503b0d4-6d0e-4b5d-896f-e13c093ad3d4
# ╟─deeff3d5-3856-43cb-8469-2a4d6d7fca4f
# ╟─b6101102-2054-4932-b6b6-5070cd84f2be
# ╠═19e28c66-a389-4903-82a5-c963cf0b90b9
# ╠═a8f1dc13-b73c-44b0-8e63-82431f904313
# ╠═bd353706-be2d-480d-8ebb-cf20cab0dbec
# ╟─e4035abe-11a1-4c2c-8312-4d1c71e2f9ab
# ╟─ab69ab43-99bf-4a92-9698-70e170761e82
# ╟─5c7c3dc0-7848-4431-829e-508664d9c5af
# ╠═b099dbbf-9648-44c5-984c-fdd80ee81469
# ╠═f8fbff0d-15d4-43ca-9f9c-29788ff793ec
# ╠═4894a323-c8c8-4f34-b311-20c64176b89d
# ╠═4bfa5d92-b177-4442-b045-e05cc48b6cc4
# ╠═bae3209a-1ad6-4e1f-833d-f1074e8239dc
# ╟─fd0fb774-a805-4739-b570-0d2e191a3294
# ╟─8987f1ce-e290-4930-909a-3e09a9113a7a
# ╟─eb7bbd16-04d3-4d7d-b4af-e77f89e4180e
# ╟─ae8fd43c-f658-47b1-8781-6f699ade6bdb
# ╟─bd2e4d28-c68b-45c9-b68e-b477b44fcd75
# ╟─264fc65e-d09a-4c67-94de-f845d42d18a3
# ╟─6e607792-e297-483f-8917-c871fa0c26d0
# ╠═532d0b24-6ace-4579-bf2a-d12c07ee9436
# ╟─f3e0b42f-370d-4887-b6a1-d9ecd83c6275
# ╠═70930f84-fe69-4fa3-a37c-eec9a377526c
# ╠═b9efca1a-b978-4794-b421-ee9df889a6a8
# ╠═01234bcc-a428-4922-969f-9a1ba49c8d62
# ╠═1971dce8-ca39-4ced-88af-c65105977ac4
# ╠═b7024377-6ffe-4f8a-a58a-0fa809683fd2
# ╠═2517dcff-c838-4b33-9f0e-e6d4bd8a259e
# ╠═0a1c908e-6432-4e8a-8c8c-7cc903fbca6c
# ╠═e204ae15-fefd-4f01-8d5e-3772aefe9b0f
# ╠═005988d8-a327-410e-9e53-09ba99dcee9a
# ╠═bf01f8b6-b618-4ce8-9ff0-3bdf37901ce1
# ╠═90d6eba4-fa47-4717-b742-8b38f0b330ac
# ╠═3e4b82aa-a6ec-4097-afb5-927b7e16262a
# ╠═b57ea0a2-bee4-450f-99fa-cc79054fc963
# ╠═8777505c-6db4-4e75-a57b-59f4620b0051
# ╠═434ecb8d-eae2-4913-97f8-3bdbefbf78ff
# ╠═b77b52ac-73ad-456c-94e5-5f72e310268f
# ╟─1cc42384-d262-496a-801a-033e7b40e6e9
# ╠═9f6051fa-2bd7-40b1-b2f8-5e7b1c477387
# ╟─518eec58-eb4b-4294-9090-b6645f51a337
# ╠═5c1387ed-973c-4109-9ace-c473a4efe9ee
# ╠═6f92c837-1760-423e-9777-9db9ad758475
# ╠═563c267e-0c6a-4b90-81a7-fc8ff4c73c75
# ╠═9fa33ac0-bb2f-4055-b124-8c7c19621226
# ╠═276e171b-271d-4610-b9bb-01b212193123
# ╠═aecea3fe-f0ee-4953-857b-29d6a5640530
# ╠═9353566b-a70e-4dbe-8f02-b16d2c0570d2
# ╠═0b9c3978-0eb9-4f80-b51d-540d4ed88c3e
# ╠═639eeead-5c9e-445b-8e40-c811634c998a
# ╠═674afcf7-925e-4cc1-ac4e-9dba4cdc4693
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

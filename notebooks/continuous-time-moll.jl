### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

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

# ╔═╡ 6d5523f2-a517-48df-9f31-bbd516e1208e
using Parameters

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

# ╔═╡ cdd24518-1cd3-43af-9a88-f0a8cd3cd6ed
md"""
`continuous-time-moll.jl` | **Version 1.0** | *last updated: June 21, 2022* | *created by [Daniel Schmidt](https://github.com/danieljschmidt)*
"""

# ╔═╡ 0edf9214-7213-4b1d-8fa6-54a365525d29
md"""
# Huggett model in continuous time

In this notebook, we consider a Huggett economy in continuous time. Income is a Poisson process with two states. 

More information on the algorithms can be found in the [Online Appendix of Achdou et al. (2021)](https://benjaminmoll.com/wp-content/uploads/2020/02/HACT_Numerical_Appendix.pdf). The Julia code in this notebook follows the closely the [Matlab code snippets](https://benjaminmoll.com/codes/) that Ben Moll uploaded on his website."""

# ╔═╡ 1c23af35-64f5-4b96-9eab-98ee9b1ecb4c
md"""
## Model
"""

# ╔═╡ 64d44910-988c-4e24-baf6-6408da65bd21
@with_kw struct HuggettPoisson
	
	σ::Float64 = 2.   # risk aversion coefficient u'(c) = c^(-σ)
	ρ::Float64 = 0.05 # rate of time preference

	z::Matrix{Float64} = reshape([0.1, 0.2], 1, :)   # income state (row vector)
	λ::Matrix{Float64} = reshape([0.02, 0.03], 1, :) # intensities (row vector)

	# asset grid parameters
	N_a::Int64 = 500
	aₘᵢₙ::Float64 = - 0.1
	aₘₐₓ::Float64 = 1.
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
	[aᵢ for aᵢ in range(aₘᵢₙ, aₘₐₓ, N_a)] |> hcat # asset grid (column vector)
end

# ╔═╡ 23a991be-7c8b-45e2-bd75-af5e146fc6b0
m = HuggettPoisson();

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

# ╔═╡ eb565042-1823-4d5a-b4d1-ee314dccd4e0
function solve_KF(m::HuggettPoisson, A)

	(; N_a, Δa) = m
	
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

	g ./ sum(g_stacked)
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
function solve_HJB_explicit(m::HuggettPoisson, r; maxit = 100000, crit = 1e-6)

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

# ╔═╡ e57a2dfa-0937-49cf-9161-fa1e39fb5e80
function consumption_and_drift((; a, z), dv, (; r, σ))
	#dv = max(dv, eps(0.0))
	c = dv ^ (-1/σ)
	ȧ = z + r * a - c

	(; c, ȧ, dv)
end

# ╔═╡ 91c8dce8-07b8-46d2-8c61-e6bccced64e4
function consumption_and_drift₀((; a, z), (; r, σ))
	ȧ = 0.0
	c = z + r * a

	dv = c^(-σ)

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

# ╔═╡ 3e4b82aa-a6ec-4097-afb5-927b7e16262a
function consumption_and_drift_upwind_vec(ss, dvf, dvb, par)
	consumption_and_drift_upwind.(ss, dvf, dvb, Ref(par)) |> StructArray
end

# ╔═╡ 1971dce8-ca39-4ced-88af-c65105977ac4
v₀((; z, a), (; r, σ, ρ)) = (z + r*a)^(1-σ)/(1-σ)/ρ

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

# ╔═╡ 90d6eba4-fa47-4717-b742-8b38f0b330ac
function inner_function(ss, v, (; aₘₐₓ, aₘᵢₙ, r, σ, N_a, Δa, z))
	# initialize arrays for forward and backward difference
	dvf = zeros(N_a, 2)
	dvb = zeros(N_a, 2)

	# forward difference
	dvf[1:N_a-1,:] .= (v[2:N_a,:] - v[1:N_a-1,:]) / Δa
	dvf[N_a,:] .= (vec(z) .+ r * aₘₐₓ) .^ (-σ) # boundary condition a  <= a_max

	# backward difference
	dvb[2:N_a,:] .= (v[2:N_a,:] - v[1:N_a-1,:]) / Δa
	dvb[1,:] .= (vec(z) .+ r * aₘᵢₙ) .^ (-σ) # boundary condition a >= a_min
	
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
function solve_HJB_implicit(m::HuggettPoisson, r; maxit = 100, crit = 1e-6, Δ = 1000)

	(; σ, ρ, z, λ, N_a, aₘᵢₙ, aₘₐₓ, Δa) = m
	da = Δa

	# construct asset grid
	a = construct_a(m)

	# initialize arrays for forward and backward difference
	dvf = zeros(N_a, 2)
	dvb = zeros(N_a, 2)

	# precompute A_switch matrix
	id = sparse(I, N_a, N_a)
	A_switch_1 = hcat(-λ[1] * id,  λ[1] * id)
	A_switch_2 = hcat( λ[2] * id, -λ[2] * id)
	A_switch   = vcat(A_switch_1, A_switch_2)

	# initial guess for value function
	v₀ = zeros(N_a, 2)
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
		v_new = reshape(v_new_stacked, N_a, 2)

		# STEP 4

		v_change = v_new - v
		dist[it] = maximum(abs.(v_change))
		
		v = v_new

		if dist[it] < crit

			ȧ = z .+ r.*a - c

			return v, c, ȧ, A, it, dist

		end

	end

	@error("Algorithm did not converge")

end

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
function solve_HJB_julian(m::HuggettPoisson, r, scheme::Scheme; v₀=v₀, crit = 1e-6)

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
	t4 = @elapsed g4 = gth_solve(A)
	

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

# ╔═╡ e1376e99-a636-4d28-b711-2dd4be66374f
function solve_df(m::HuggettPoisson, r; maxit = 100, crit = 1e-6, Δ = 1000)
	scheme = Implicit(; maxit, Δ)
	(; v, c, ȧ, A, it_last, dist) = solve_HJB_julian(m, r, scheme; crit)
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
function excess_demand(m::HuggettPoisson, r; maxit = 100, crit = 1e-6, Δ = 1000)
	(; Δa) = m
	(df, it_last, dist) = solve_df(m, r; maxit, crit, Δ)
	A = dot(df.a, df.g) * Δa

end

# ╔═╡ bd353706-be2d-480d-8ebb-cf20cab0dbec
# ╠═╡ disabled = true
#=╠═╡
r_eq = find_zero(r -> excess_demand(m, r), initial_bracket, Brent())
  ╠═╡ =#

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
	end
end

# ╔═╡ f8fbff0d-15d4-43ca-9f9c-29788ff793ec
function solve_explicit_df(m::HuggettPoisson, r; maxit = 100000, crit = 1e-6)
	scheme = Explicit(m; maxit)
	(; v, c, ȧ, it_last, dist) = solve_HJB_julian(m, r, scheme; crit)
	df = results_to_df(m; v, c, ȧ)
	
	return df, it_last, dist
	
end

# ╔═╡ 4894a323-c8c8-4f34-b311-20c64176b89d
(df2, it_last2, dist2) = solve_explicit_df(m, 0.03; crit=1e-6);

# ╔═╡ ebd88b43-e277-4c79-b443-1661a3c438b8
(df2[:,[:c, :ȧ, :v]] .- df[:,[:c, :ȧ, :v]]) ./ df[:,[:c, :ȧ, :v]]

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
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
InfinitesimalGenerators = "2fce0c6f-5f0b-5c85-85c9-2ffe1d5ee30d"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Parameters = "d96e819e-fc66-5662-9728-84c9c7592b0a"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
QuantEcon = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
AlgebraOfGraphics = "~0.6.8"
Arpack = "~0.5.3"
CairoMakie = "~0.8.7"
Chain = "~0.4.10"
DataFrames = "~1.3.4"
InfinitesimalGenerators = "~0.4.0"
Parameters = "~0.12.3"
PlutoUI = "~0.7.39"
QuantEcon = "~0.16.3"
Roots = "~1.4.1"
StructArrays = "~0.6.8"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0-rc1"
manifest_format = "2.0"
project_hash = "36a5a0077662c9a75f41d5ffd35e39a948eff97e"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "6f1d9bc1c08f9f4a8fa92e3ea3cb50153a1b40d4"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.1.0"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "896c6cd357f08890d407e36c1301087c42a13e61"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.0"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "d8eeef6f37c3c96f5d3d36595090222c45003bfe"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.8"

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
git-tree-sha1 = "f87e559f87a45bece9c9ed97458d3afe98b1ebb9"
uuid = "ec485272-7323-5ecc-a04f-4719b315124d"
version = "0.1.0"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "91ca22c4b8437da89b030f08d71db55a379ce958"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.3"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "d618d3cf75e8ed5064670e939289698ecf426c7f"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.12"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

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

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BenchmarkTools]]
deps = ["JSON", "Logging", "Printf", "Profile", "Statistics", "UUIDs"]
git-tree-sha1 = "4c10eee4af024676200bc7752e536f858c6b8f93"
uuid = "6e4b80f9-dd63-53aa-95a3-0cdb28fa8baf"
version = "1.3.1"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA"]
git-tree-sha1 = "9510066f04dd9363adbc70432523edc1550363a6"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.7"

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

[[deps.Chain]]
git-tree-sha1 = "339237319ef4712e6e5df7758d0bccddf5c237d9"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.4.10"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

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
git-tree-sha1 = "7297381ccb5df764549818d9a7d57e45f1057d30"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.18.0"

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
version = "0.5.2+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DSP]]
deps = ["Compat", "FFTW", "IterTools", "LinearAlgebra", "Polynomials", "Random", "Reexport", "SpecialFunctions", "Statistics"]
git-tree-sha1 = "3fb5d9183b38fdee997151f723da42fb83d1c6f2"
uuid = "717857b8-e6f2-59f4-9121-6e50c889abd2"
version = "0.7.6"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

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

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

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
git-tree-sha1 = "7669d53b75e9f9e2fa32d5215cb2af348b2c13e2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.21"

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
git-tree-sha1 = "0ec161f87bf4ab164ff96dfacf4be8ffff2375fd"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.62"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

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
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "a087a23129ac079d43ba6b534c6350325fcd41c9"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "d8a578692e3077ac998b50c0217dfd67f21d1e5f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.0+0"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "505876577b5481e50d089c1c68899dfb6faebc62"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.4.6"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[deps.FiniteDiff]]
deps = ["ArrayInterfaceCore", "LinearAlgebra", "Requires", "SparseArrays", "StaticArrays"]
git-tree-sha1 = "ee13c773ce60d9e95a6c6ea134f25605dce2eda3"
uuid = "6a86dc24-6348-571c-b903-95158fe2bd41"
version = "2.13.0"

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

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "de6980e052d67c0da1872dfdb2c49fb7d3f56b07"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.0"

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
git-tree-sha1 = "d778af2dcb083169807d43aa9d15f9f7e3909d4c"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.7.7"

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
deps = ["DualNumbers", "LinearAlgebra", "SpecialFunctions", "Test"]
git-tree-sha1 = "cb7099a0109939f16a4d3b572ba8396b1f6c7c31"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.10"

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

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "d9a03ffc2f6650bd4c831b285637929d99a4efb5"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.5"

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

[[deps.InfinitesimalGenerators]]
deps = ["Arpack", "Distributions", "FillArrays", "FiniteDiff", "KrylovKit", "LinearAlgebra", "Roots"]
git-tree-sha1 = "3f450b067ee224ebd0e97ffbfac3e1615cc832f5"
uuid = "2fce0c6f-5f0b-5c85-85c9-2ffe1d5ee30d"
version = "0.4.0"

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

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
deps = ["AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "b7bc05649af456efc75d178846f47006c2c4c3c7"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.13.6"

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
git-tree-sha1 = "591e8dc09ad18386189610acafb970032c519707"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.3"

[[deps.KrylovKit]]
deps = ["LinearAlgebra", "Printf"]
git-tree-sha1 = "49b0c1dd5c292870577b8f58c51072bd558febb9"
uuid = "0b1a1467-8014-51b9-945f-bf0ae24f4b77"
version = "0.5.4"

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
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.81.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

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

[[deps.LightGraphs]]
deps = ["ArnoldiMethod", "DataStructures", "Distributed", "Inflate", "LinearAlgebra", "Random", "SharedArrays", "SimpleTraits", "SparseArrays", "Statistics"]
git-tree-sha1 = "432428df5f360964040ed60418dd5601ecd240b6"
uuid = "093fc24a-ae57-5d10-9952-331d41423f4d"
version = "1.3.5"

[[deps.LineSearches]]
deps = ["LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "Printf"]
git-tree-sha1 = "f27132e551e959b3667d8c93eae90973225032dd"
uuid = "d3d80556-e9d4-5f37-9878-2ab0fcc64255"
version = "7.1.1"

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
git-tree-sha1 = "09e4b894ce6a976c354a69041a04748180d43637"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.15"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

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
git-tree-sha1 = "390df26c6a7a90c7c3b94fee62bdb4492f791273"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.7"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "c8fdab2a70daafaa81aaff28659eba4831e22ef9"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.3.3"

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
git-tree-sha1 = "c167b0d6d165ce49f35fbe2ee1aea8844e7c7cea"
uuid = "b8f27783-ece8-5eb3-8dc8-9495eed66fee"
version = "1.4.0"

[[deps.MathProgBase]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "9abbe463a1e9fc507f12a69e7f29346c2cdc472c"
uuid = "fdba3010-5040-5b88-9595-932c9decdf73"
version = "0.7.8"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "e8610e4631e395c42cffeb4937683a37bf8ffd53"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.2"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

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
version = "2022.2.1"

[[deps.MutableArithmetics]]
deps = ["LinearAlgebra", "SparseArrays", "Test"]
git-tree-sha1 = "4e675d6e9ec02061800d6cfb695812becbd03cdf"
uuid = "d8a4904e-b15c-11e9-3269-09a3773c0cb0"
version = "1.0.4"

[[deps.NLSolversBase]]
deps = ["DiffResults", "Distributed", "FiniteDiff", "ForwardDiff"]
git-tree-sha1 = "50310f934e55e5ca3912fb941dec199b49ca9b68"
uuid = "d41bc354-129a-5804-8e4c-c37616107c6c"
version = "7.8.2"

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
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "dfd8d34871bc3ad08cd16026c1828e271d554db9"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.1"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "ec2e30596282d722f018ae784b7f44f3b88065e4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.6"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

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
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ab05aa4cc89736e95915b01e7279e61b1bfe33b8"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Optim]]
deps = ["Compat", "FillArrays", "ForwardDiff", "LineSearches", "LinearAlgebra", "NLSolversBase", "NaNMath", "Parameters", "PositiveFactorizations", "Printf", "SparseArrays", "StatsBase"]
git-tree-sha1 = "7a28efc8e34d5df89fc87343318b0a8add2c4021"
uuid = "429524aa-4258-5aef-a3af-852621145aeb"
version = "1.7.0"

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
git-tree-sha1 = "7f4869861f8dac4990d6808b66b57e5a425cfd99"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.13"

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
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "a7a7e1a88853564e551e4eba8650f8c38df79b37"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.1.1"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "bb16469fd5224100e422f0b027d26c5a25de1200"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.2.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.Polynomials]]
deps = ["LinearAlgebra", "MutableArithmetics", "RecipesBase"]
git-tree-sha1 = "a8d37fbaba422166e9f5354b6d8f6197e1f74fe5"
uuid = "f27b6e38-b328-58d1-80ce-0feddd5e7a45"
version = "3.1.3"

[[deps.PooledArrays]]
deps = ["DataAPI", "Future"]
git-tree-sha1 = "a6062fe4063cdafe78f4a0a81cfffb89721b30e7"
uuid = "2dfb63ee-cc39-5dd5-95bd-886bf059d720"
version = "1.4.2"

[[deps.PositiveFactorizations]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "17275485f373e6673f7e7f97051f703ed5b15b20"
uuid = "85a6dd25-e78a-55b7-8502-1745935b8125"
version = "0.2.4"

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
git-tree-sha1 = "747f4261ebe38a2bc6abf0850ea8c6d9027ccd07"
uuid = "27ebfcd6-29c5-5fa9-bf4b-fb8fc14df3ae"
version = "0.5.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Profile]]
deps = ["Printf"]
uuid = "9abbd945-dff8-562f-b5e8-e1ebf5ef1b79"

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

[[deps.QuantEcon]]
deps = ["DSP", "DataStructures", "Distributions", "FFTW", "LightGraphs", "LinearAlgebra", "Markdown", "NLopt", "Optim", "Pkg", "Primes", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "Test"]
git-tree-sha1 = "d777434be1b3536821caea3fc5c4d9fd9d350c4f"
uuid = "fcd29c91-0bd7-5a09-975d-7ac3f643a60c"
version = "0.16.3"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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
git-tree-sha1 = "838b60ee62bebc794864c880a47e331e00c47505"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "1.4.1"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "c49318f1b9ca3d927ae576d323fa6f724d01ba53"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.1"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

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

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

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
git-tree-sha1 = "69fa1bef454c483646e8a250f384e589fd76562b"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "1.8.6"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "2bbd9f2e40afd197a1379aef05e0d85dba649951"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.4.7"

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
git-tree-sha1 = "8977b17906b0a1cc74ab2e3a05faa16cf08a8291"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.16"

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
git-tree-sha1 = "9abba8f8fb8458e9adf07c8a2377a070674a24f1"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.8"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

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
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f90022b44b7bf97952756a6b6737d1a0024a3233"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.5"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

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
version = "1.2.12+3"

[[deps.isoband_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51b5eeb3f98367157a7a12a1fb0aa5328946c03c"
uuid = "9a68df92-36a6-505f-a73e-abb412b6bfb4"
version = "0.2.3+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.0+0"

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
version = "1.41.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

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
# ╟─1fbf9f02-7ad2-4ddd-ac2c-92af79c9ed02
# ╟─2e09462e-6030-4e33-bc7a-9d2faed4ca74
# ╠═4ebeb5e8-505d-4ed2-8d37-4a64ceb7d796
# ╠═6c895dce-9c0a-499c-9e95-45c29ee5a459
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
# ╠═ebd88b43-e277-4c79-b443-1661a3c438b8
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
# ╠═e57a2dfa-0937-49cf-9161-fa1e39fb5e80
# ╠═91c8dce8-07b8-46d2-8c61-e6bccced64e4
# ╠═7734c75e-5f2b-4807-8b9a-82e7e010edac
# ╠═3e4b82aa-a6ec-4097-afb5-927b7e16262a
# ╠═1971dce8-ca39-4ced-88af-c65105977ac4
# ╠═b7024377-6ffe-4f8a-a58a-0fa809683fd2
# ╠═2517dcff-c838-4b33-9f0e-e6d4bd8a259e
# ╠═0a1c908e-6432-4e8a-8c8c-7cc903fbca6c
# ╠═e204ae15-fefd-4f01-8d5e-3772aefe9b0f
# ╠═005988d8-a327-410e-9e53-09ba99dcee9a
# ╠═bf01f8b6-b618-4ce8-9ff0-3bdf37901ce1
# ╠═90d6eba4-fa47-4717-b742-8b38f0b330ac
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
# ╠═6d5523f2-a517-48df-9f31-bbd516e1208e
# ╠═9fa33ac0-bb2f-4055-b124-8c7c19621226
# ╠═276e171b-271d-4610-b9bb-01b212193123
# ╠═aecea3fe-f0ee-4953-857b-29d6a5640530
# ╠═9353566b-a70e-4dbe-8f02-b16d2c0570d2
# ╠═0b9c3978-0eb9-4f80-b51d-540d4ed88c3e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

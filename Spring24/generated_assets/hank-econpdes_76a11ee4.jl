### A Pluto.jl notebook ###
# v0.19.38

#> [frontmatter]
#> chapter = 8
#> section = 2
#> order = 2
#> title = "HANK with EconPDEs.jl"
#> layout = "layout.jlhtml"
#> tags = ["unfinished"]
#> description = ""

using Markdown
using InteractiveUtils

# ╔═╡ cdbd5837-456a-4ab8-99e7-a291ecaedb68
using Distributions

# ╔═╡ eb9ebd9c-d1aa-4b9e-b211-31737b11aeef
using CairoMakie

# ╔═╡ 8cce558b-91bd-40cc-909f-53f0fd5716c9
using PlutoUI

# ╔═╡ ed10146f-1397-43e4-9426-ed15a3130067
using EconPDEs

# ╔═╡ 5eb30630-253e-421c-906b-2345266d53f5
using EconPDEs: StateGrid

# ╔═╡ 48b79f5f-0ca9-478f-b117-baeb3b3bc8a4
md"""
`hank-econpdes.jl` | **Version 0.1** | *last updated: May 10 2022*
"""

# ╔═╡ 5aa0904d-5ab5-4c3d-8234-7053053b3a59
md"""
# HANK
"""

# ╔═╡ 8501d4f2-e9bf-423e-a7ae-8337577cae6d
u(c, (; γ)) = c > 0 ? γ == 1 ? log(c) : c^(1-γ)/(1-γ) : 100 * c - 1000

# ╔═╡ 8e3792bb-a977-4a22-bbbb-8fbc839a15b1
u_prime(c, (; γ)) = c^(-γ)

# ╔═╡ 1c12771b-51ce-43f2-b397-6e9eab92dfe2
u_prime_inv(x, (; γ)) = x^(-1/γ)

# ╔═╡ b2450bcf-e13c-4948-b068-c4f6d69d61fb
rb(b, (; rb_pos, rb_neg)) = b < 0 ? rb_neg : rb_pos

# ╔═╡ c59cbd45-aedd-4bdf-a1e8-af6321c9b0fb
ra(a, (; ra, τ, amax)) = ra * (1 - (1.33*amax/a)^(1-τ))

# ╔═╡ e9fc109b-fc8f-45c5-a2ea-89e59e8990d8
χ(d, a, (; χ₀, χ₁)) = χ₀ * abs(d) + χ₁ * d^2/2 * (max(a,10^(-5)))^(-1)

# ╔═╡ 8a3aca3a-54bf-46a4-bbd1-a5bd5edd6a6f
function two_asset_kinked_foc(va, vb, a, (; χ₀, χ₁))
	d = min(va/vb - 1 + χ₀, 0) * a/ χ₁ +  max(va / vb - 1 - χ₀,0) * a/ χ₁
end

# ╔═╡ 95de908f-0023-4515-9fc3-381e766ecf7f
function c_bdot(vb, (; b, z), param)
	(; ξ, w) = param
	vb = max(vb, √√eps())

	c_star = u_prime_inv(vb, param)
	b_dot_c = (1-ξ) * w * z + b * rb(b, param) - c_star # bdrift1

	Hc = u(c_star, param) + vb * b_dot_c
	
	(; vb_c = vb, c_star, b_dot_c, Hc)
end

# ╔═╡ 49f4011c-018d-430c-a69b-0bd01d23c022
function c_bdot₀((; b, z), param)
	(; ξ, w) = param

	c_star = (1-ξ) * w * z + rb(b, param) * b # bdrift1(0)
	vb_c = u_prime(c_star, param)
	b_dot_c = 0.0

	Hc = u(c_star, param)
	
	(; vb_c, c_star, b_dot_c, Hc)
end

# ╔═╡ e848cc15-3cb1-449b-afb5-3fb5cd6026c1
function d_bdot(va, vb, (; a, b, z), param)
	va = max(va, √eps())
	vb = max(vb, √eps())

	d_star = two_asset_kinked_foc(va, vb, a, param)

	b_dot_d = - d_star - χ(d_star, a, param) # bdrift2
	a_dot_d = d_star # adrift2 = d

	Hd = va * a_dot_d + vb * b_dot_d
	(; va_d = va, vb_d = vb, d_star, b_dot_d, a_dot_d, Hd)
end

# ╔═╡ 8c79d03b-a4d9-41fc-b432-a807f532b7af
function d_bdot₀((; a, b, z), param)
	a_dot_d = 0.0
	d_star = 0.0
	b_dot_d = 0.0
	(; va_d = 0.0, vb_d = 0.0, d_star, b_dot_d, a_dot_d)
end

# ╔═╡ a94b4061-cba6-4d6d-9805-d903fe8f23dd
function policies_drifts₀((; a, b, z), param)
	(; ξ, w) = param
	
	c_star = (1-ξ) * w * z + rb(b, param) * b
#	vk = u_prime(c_star; γ)
	b_dot_c = 0.0
	
	(; c_star, b_dot_c)
end

# ╔═╡ eae9c1a1-ab55-407a-80e3-d41881521125
begin
	Base.@kwdef struct SimpleHANK
		γ = 2.0 # ga = 2; %CRRA utility with parameter gamma
		ra = 0.05
		rb_pos = 0.03
		rb_neg = 0.12
		ρ = 0.06 #; %discount rate
		χ₀ = 0.03
		χ₁ = 2.0
		τ = 10
		ξ = 0.1 #; %fraction of income that is automatically deposited
	#	if ra - 1/cχ₁ > 0
	#	    @warn "Warning: ra - 1/χ₁ > 0"
	#	end
		# Income process (two-state Poisson process):
		w = 4
		Nz = 2
		z  = [.8, 1.3]
		la_mat = [-1/3 1/3; 1/3 -1/3]
		bmin = -2.0 # 0
		bmax = 40
		amin = 0
		amax = 70

		κz = 0.1
		zbar = 1.0
		σz = 0.07
	end

	function (m::SimpleHANK)(state::NamedTuple, y::NamedTuple)
		(; ρ, ξ, w, amax, amin, bmin) = m  
		
		# Only relevant for RBC
		if haskey(state, :z)
			(; z) = state
			(; vz_up, vz_down, vzz) = y
			(; κz, zbar, σz) = m
			μz = κz * (zbar - z)
    		vz = (μz >= 0) ? vz_up : vz_down
			exo = μz * vz + 0.5 * vzz * σz^2
		else
			exo = 0
			state = (; z = 2, state...)
		end
		
		(; a, b, z) = state
		(; v, va_up, va_down, vb_up, vb_down, vaa, vbb, vab) = y

		# consumption
		U = c_bdot(vb_up, state, m)
		D = c_bdot(vb_down, state, m)
		O = c_bdot₀(state, m)
		if U.b_dot_c > 0.0 #&& (D.b_dot_c ≥ 0.0 || U.Hc ≥ D.Hc) && U.Hc ≥ O.Hc
			out_c = U
		elseif D.b_dot_c < 0.0 #&& (U.b_dot_c < 0.0 || D.Hc ≥ U.Hc) && D.Hc ≥ O.Hc
			out_c = D
		else
			out_c = O
		end
		if (b ≈ bmin) && (out_c.b_dot_c ≤ 0.0)
			out_c = O
    	end
		(; vb_c, c_star, b_dot_c) = out_c
		
		# deposits - try all combinations of va and vb
		
		UU = d_bdot(va_up, vb_up, state, m)
		DD = d_bdot(va_down, vb_down, state, m)
		UD = d_bdot(va_up, vb_down, state, m)
		DU = d_bdot(va_down, vb_up, state, m)
		OO = d_bdot₀(state, m) #vb_d = max((vb_up + vb_down)/2, eps())

		validUU = UU.a_dot_d > 0 && UU.b_dot_d > 0 && UU.Hd > 0 && a != amax
		validDD = DD.a_dot_d < 0 && DD.b_dot_d < 0 && DD.Hd > 0 && a != amin
		validUD = UD.a_dot_d > 0 && UD.b_dot_d < 0 && UD.Hd > 0 && a != amax
		validDU = DU.a_dot_d < 0 && DU.b_dot_d > 0 && DU.Hd > 0 && a != amin

		valid_d = sum([validUU, validDD, validUD, validDD])
		
		if 		validUD &&
				(!validDU || UD.Hd ≥ DU.Hd) &&
				(!validDD || UD.Hd ≥ DD.Hd)
			out = UD
			dir = 1
		elseif 	validDU &&
				(!validUD || DU.Hd ≥ UD.Hd) &&
				(!validDD || DU.Hd ≥ DD.Hd)
			out = DU
			dir = 2
		elseif  validDD &&
				(!validUD || DD.Hd ≥ UD.Hd) && 
				(!validDU || DD.Hd ≥ DU.Hd)
			out = DD
			dir = 3
		elseif  validUU && 
				(!validUD || DD.Hd ≥ UD.Hd) && 
				(!validDU || DD.Hd ≥ DU.Hd)
			out = UU
			dir = 4
		#elseif DU.a_dot < 0 && DU.b_dot_d > 0
		#	out = DU
		#	dir = 4
		elseif !validDD && !validUD && !validDU
			out = OO
			dir = 5
		else
			@error "Unreachable reached"
		end

		(; va_d, vb_d, d_star, a_dot_d, b_dot_d) = out

		va_a = a < 0 ? va_down : va_up
	
		endo = u(c_star, m) + vb_c * b_dot_c +
					vb_d * b_dot_d + 
					#va_d * a_dot_d +
					va_a * ra(a, m) * a

		vt = ρ * v - (endo + exo)

		(; vt), (; c_star, d_star, dir, valid_d)
	end
end

# ╔═╡ e6f9702a-f816-42a7-908a-8b62e841ec5b
grids((; bmin, bmax, amin, amax, κz, zbar, σz)) = let
	I = 50 #100
	#bmax = 40
	b_grid = range(bmin,bmax,I)
	
	J = 25 #50
	#amin = 0
	#amax = 70
	a_grid = range(amin, amax, J)

	#distribution = Gamma(2 * κz * zbar / σz^2, σz^2 / (2 * κz))
	#z_grid = range(quantile(distribution, 0.01), quantile(distribution, 0.99), length = 5)
	(; a_grid, b_grid)#, z_grid)
end

# ╔═╡ 6d7b011f-b9d0-493e-a1ac-62b8432a18f2
v0(a, b, z, (; ρ, ξ, ra, rb_neg, γ, w)) =
	u((1-ξ)*w*z + ra * a + rb_neg * b, (; γ))/ρ

# ╔═╡ 6dbdd3ce-d690-44f8-9861-263154c063b3
begin
	m = SimpleHANK(amax = 100, bmax = 50)
	g = grids(m)
	yend = OrderedDict(:v => [v0(a, b, 1, m) for a ∈ g.a_grid, b ∈ g.b_grid])#, z ∈ g.z_grid])
	stategrid = OrderedDict(:a => g.a_grid, :b => g.b_grid)#, :z => g.z_grid)
	result = pdesolve(m, stategrid, yend) #, 0.0001:0.0001:0.0010)
	#@assert result.residual_norm <= 1e-5
end

# ╔═╡ 4dc4d556-32e1-4678-91a4-f64bf6beb7c5
surface(result.zero[:v], axis = (type = Axis3, )) #.zero[:v] #|> allunique

# ╔═╡ 6ac298ff-fdd3-49a0-a295-8c28bc8f93ce
surface(result.optional[:c_star], axis = (type = Axis3, ))

# ╔═╡ 8a852666-b506-444d-a57b-534f9d7733d4
surface(g.a_grid, g.b_grid, result.optional[:d_star], axis = (type = Axis3, ))

# ╔═╡ b3113a37-f9b7-4e3c-b253-91d6241b6d4c
hist(vec(result.optional[:dir]))

# ╔═╡ 6edc36e7-8404-4efc-828d-1c9f9cff1fa7
hist(vec(result.optional[:valid_d]))

# ╔═╡ b99aeb8d-313e-43c4-8682-e2179421f64a
md"""
# Appendix
"""

# ╔═╡ 64696f48-2209-4550-bc4a-c86e9deb9c73


# ╔═╡ e3ffb03b-d9bb-42dc-b9f5-1797f7b1d705
TableOfContents()

# ╔═╡ 9fe0fb70-486b-4b24-9805-e21f2cae0d38
begin
	admonition(kind, title, text) = Markdown.MD(Markdown.Admonition(kind, title, [text]))
	hint(   text, title="Hint")    = admonition("hint",    title, text)
	warning(text, title="Warning") = admonition("warning", title, text)
	danger( text, title="Danger")  = admonition("danger",  title, text)
	correct(text, title="Correct") = admonition("correct", title, text)

	almost(text) = warning(text, "Almost there!")
	keep_working(text=md"The answer is not quite right.") = danger(text, "Keep working on it!")
	yays = [md"Great!", md"Yay ❤", md"Great! 🎉", md"Well done!", md"Keep it up!", md"Good job!", md"Awesome!", md"You got the right answer!", md"Let's move on to the next section."]
	got_it(text=rand(yays)) = correct(text, "Got it!")
end

# ╔═╡ d0c842f3-3ed9-45f5-b8dc-d6d1d91b041c
danger(
	md"**This notebook is not ready for public consumption.** Use at your own risk.", "Under construction!"
)

# ╔═╡ 25312bea-4adc-478c-b51e-e92ecc3a0705
md"""
## Patching `EconPDEs.jl`

EconPDEs doesn't allow three state variables
"""

# ╔═╡ a38744fc-f87e-4eba-8aaa-8c8a561cbced
Δy_up(y, i, Δx) = (y[i+1] - y[i]) / Δx

# ╔═╡ 62ebb755-4f34-4243-92ba-cd09e90a3ef9
Δy_down(y, i, Δx) = (y[i] - y[i-1]) / Δx

# ╔═╡ 84d5eb68-42a6-4342-abb8-b003afd0f32f
Δy_central(y, i, Δx) = (y[i+1] - y[i-1]) / Δx

# ╔═╡ 5f5ac3cb-03ce-4d64-a186-1da4e09046ee
function Δgrid(grid, i)
    last = length(grid)
    @inbounds down = grid[max(i, 2)]      - grid[max(i-1, 1)]
    @inbounds up   = grid[min(i+1, last)] - grid[min(i, last-1)]
    central = (up + down)
    avg = central / 2

    (; up, down, avg, central)
end

# ╔═╡ 5479a846-43a0-4d65-84a2-a45c1aa21145
deriv_names(fun_name, state_name) = (Symbol(fun_name, state_name, "_", :up), Symbol(fun_name, state_name, "_", :down), Symbol(fun_name, state_name, state_name))

# ╔═╡ 70532bb4-eee5-47f7-8b23-3d5fa387fa0f
function Δy(y, bc, i, Δx, fun_name, state_name)
    up       = i != length(y) ? Δy_up(y, i, Δx.up)     : bc[i]
    down     = i != 1         ? Δy_down(y, i, Δx.down) : bc[i]
    second = (up - down) / Δx.avg
    NamedTuple{deriv_names(fun_name, state_name)}((up, down, second))
end

# ╔═╡ fecb3988-fbe6-4b88-b234-a7fbfded0989
function cross_difference(y, grids, inds)
    @assert length(grids) == length(inds) == length(size(y)) == 2
    
    i1, i2 = Tuple(inds)
    grid1, grid2 = grids

    Δx1 = Δgrid(grid1, i1)
    Δx2 = Δgrid(grid2, i2)
    
    i1_lo = max(i1-1, 1)
    i1_hi = min(i1+1, length(grid1))
    
    if i2 == 1
        a = Δy_up(view(y, i1_hi, :), i2, Δx2.central) # use Δx2.up?
        b = Δy_up(view(y, i1_lo, :), i2, Δx2.central) # use Δx2.up?
    elseif i2 == length(grid2)
        a = Δy_down(view(y, i1_hi, :), i2, Δx2.central) # use Δx2.down?
        b = Δy_down(view(y, i1_lo, :), i2, Δx2.central) # use Δx2.down?
    else
        a = Δy_central(view(y, i1_hi, :), i2, Δx2.central)
        b = Δy_central(view(y, i1_lo, :), i2, Δx2.central)
    end

    vab = (a - b) / Δx1.central # adjust Δx1.central when i1 is adjusted?
end

# ╔═╡ 58e350ff-70cd-4022-9c6f-f540a69b97e9
function select_all_but_one_dim(y0, dim_inds_drop)
    y = reshape(view(y0, :), size(y0))

    for (dim_drop, i_drop) ∈ reverse(dim_inds_drop)
        y = selectdim(y, dim_drop, i_drop)
    end
    y
end

# ╔═╡ 70650555-1750-4ad0-9506-ef35813d11e1
# three state variables
function differentiate(Tsolution, grid::StateGrid{<: Any, 3, <: NamedTuple}, y, inds, bc)
    solnames = Tsolution.parameters[1]
    grids = grid.x
    statenames = collect(keys(grids))
    dim_inds = [(; dim, i) for (dim, i) ∈ enumerate(Tuple(inds))]

    n_states = length(grids)
    
    nts = map(enumerate(solnames)) do (k, solname)

        yk = selectdim(y, n_states+1, k)
        bck = selectdim(bc, n_states+1, k)

        # upwind differences for each state        
        nts1 = map(enumerate(statenames)) do (dim, statename)
            i = inds[dim]
            grid = grids[statename]
            Δx = Δgrid(grid, i)

            dim_inds_drop = dim_inds[Not(dim)]

            y_sub = select_all_but_one_dim(yk, dim_inds_drop)
            bc_sub = select_all_but_one_dim(bck, dim_inds_drop)
    
            va = Δy(y_sub, bc_sub, i, Δx, solname, statename)    
        end

        # upwind cross-differences for each combination of state
        nts2 = map(dim_inds) do (dim_drop, i_drop)
            state_drop = statenames[dim_drop]
            
            sub_grids = delete(grids, state_drop)
            sub_inds = [Tuple(inds)...][Not(dim_drop)]
            sub_y  = selectdim(yk,  dim_drop, i_drop)
            sub_bc = selectdim(bck, dim_drop, i_drop)
            sub_statenames = filter(!=(state_drop), statenames)
            
            vab = cross_difference(sub_y, sub_grids, sub_inds)
            cross_name = Symbol(solname, sub_statenames...)

            (; Symbol(solname, sub_statenames...) => vab)
        end    

        (; solname => yk[Tuple(inds)...], merge(nts1...)..., merge(nts2...)...)
    end
    merge(nts...)
end


# ╔═╡ Cell order:
# ╟─d0c842f3-3ed9-45f5-b8dc-d6d1d91b041c
# ╟─48b79f5f-0ca9-478f-b117-baeb3b3bc8a4
# ╟─5aa0904d-5ab5-4c3d-8234-7053053b3a59
# ╠═8501d4f2-e9bf-423e-a7ae-8337577cae6d
# ╠═8e3792bb-a977-4a22-bbbb-8fbc839a15b1
# ╠═1c12771b-51ce-43f2-b397-6e9eab92dfe2
# ╠═b2450bcf-e13c-4948-b068-c4f6d69d61fb
# ╠═c59cbd45-aedd-4bdf-a1e8-af6321c9b0fb
# ╠═e9fc109b-fc8f-45c5-a2ea-89e59e8990d8
# ╠═8a3aca3a-54bf-46a4-bbd1-a5bd5edd6a6f
# ╠═95de908f-0023-4515-9fc3-381e766ecf7f
# ╠═49f4011c-018d-430c-a69b-0bd01d23c022
# ╠═e848cc15-3cb1-449b-afb5-3fb5cd6026c1
# ╠═8c79d03b-a4d9-41fc-b432-a807f532b7af
# ╠═a94b4061-cba6-4d6d-9805-d903fe8f23dd
# ╠═eae9c1a1-ab55-407a-80e3-d41881521125
# ╠═6dbdd3ce-d690-44f8-9861-263154c063b3
# ╠═4dc4d556-32e1-4678-91a4-f64bf6beb7c5
# ╠═6ac298ff-fdd3-49a0-a295-8c28bc8f93ce
# ╠═8a852666-b506-444d-a57b-534f9d7733d4
# ╠═b3113a37-f9b7-4e3c-b253-91d6241b6d4c
# ╠═6edc36e7-8404-4efc-828d-1c9f9cff1fa7
# ╠═e6f9702a-f816-42a7-908a-8b62e841ec5b
# ╠═6d7b011f-b9d0-493e-a1ac-62b8432a18f2
# ╟─b99aeb8d-313e-43c4-8682-e2179421f64a
# ╠═64696f48-2209-4550-bc4a-c86e9deb9c73
# ╠═cdbd5837-456a-4ab8-99e7-a291ecaedb68
# ╠═eb9ebd9c-d1aa-4b9e-b211-31737b11aeef
# ╠═8cce558b-91bd-40cc-909f-53f0fd5716c9
# ╠═e3ffb03b-d9bb-42dc-b9f5-1797f7b1d705
# ╠═ed10146f-1397-43e4-9426-ed15a3130067
# ╠═9fe0fb70-486b-4b24-9805-e21f2cae0d38
# ╟─25312bea-4adc-478c-b51e-e92ecc3a0705
# ╠═5eb30630-253e-421c-906b-2345266d53f5
# ╠═70650555-1750-4ad0-9506-ef35813d11e1
# ╠═a38744fc-f87e-4eba-8aaa-8c8a561cbced
# ╠═62ebb755-4f34-4243-92ba-cd09e90a3ef9
# ╠═84d5eb68-42a6-4342-abb8-b003afd0f32f
# ╠═5f5ac3cb-03ce-4d64-a186-1da4e09046ee
# ╠═70532bb4-eee5-47f7-8b23-3d5fa387fa0f
# ╠═5479a846-43a0-4d65-84a2-a45c1aa21145
# ╠═fecb3988-fbe6-4b88-b234-a7fbfded0989
# ╠═58e350ff-70cd-4022-9c6f-f540a69b97e9

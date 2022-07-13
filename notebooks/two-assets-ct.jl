### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ b48b0674-1bcb-48e5-9b05-57dea5877715
using LinearAlgebra

# ╔═╡ 1ee8ccc3-d498-48c6-b299-1032165e4ab9
using StructArrays

# ╔═╡ 9aa61364-51a3-45d0-b1c2-757b864de132
using PlutoTest

# ╔═╡ 026cfe16-ff0f-4f68-b412-b1f6c1902824
using CairoMakie, AlgebraOfGraphics

# ╔═╡ 9d99416e-8119-4a40-b577-0135050a0e4e
using SparseArrays

# ╔═╡ 6188aab9-86bf-4ec4-bb10-43b59f71e3e2
using DataFrames, DataFrameMacros, Chain

# ╔═╡ 9696e6ca-6953-43e2-8d47-fbfe24ba4250
using PlutoUI

# ╔═╡ c89f1918-01ee-11ed-22fa-edd66e0f6c59
md"""
# Solving the HJB: Implicit scheme
"""

# ╔═╡ 629b8291-0f13-419e-b1c0-d10d5e708720
md"""
## Setup
"""

# ╔═╡ 828dee22-1ee7-41c4-b68b-f88facea86d9
md"""
## HJB Moll
"""

# ╔═╡ 3ab0c985-5317-4ea4-bddc-6289ab90bcad
function two_asset_kinked_cost_new(d,a, (; χ₀, χ₁))
	χ₀ * abs(d) + χ₁ * d^2/2 *(max(a,10^(-5)))^(-1)
end

# ╔═╡ f2ce6352-450e-4cde-a2fe-3586461c3bdf
function two_asset_kinked_FOC_new(pa, pb, a, (; χ₀, χ₁))
	min(pa / pb - 1 +  χ₀, 0.0) * a / χ₁ + max(pa/pb - 1 -  χ₀, 0.0) * a / χ₁
end

# ╔═╡ 9c1ef8d1-57bc-4da2-83ec-fb8f1a8ce296
function two_asset_kinked_cost(d,a, (; chi0, chi1))
	chi0 * abs(d) + chi1 * d^2/2 *(max(a,10^(-5)))^(-1)
end

# ╔═╡ 91b63bfc-f4a4-41c1-a472-7d13df27b93c
function two_asset_kinked_FOC(pa,pb,a, (; chi0, chi1))
	min(pa / pb - 1 + chi0, 0.0) * a / chi1 + max(pa/pb - 1 - chi0, 0.0) * a / chi1
end

# ╔═╡ 68a97aab-7924-4472-aa47-7903add8aea4
function solve_HJB_base(maxit = 35)
	ga = 2 #CRRA utility with parameter gamma
	ra = 0.05
	rb_pos = 0.03
	rb_neg = 0.12
	rho = 0.06 #discount rate
	chi0 = 0.03
	chi1 = 2;


	xi = 0.1 #fraction of income that is automatically deposited

	if ra - 1/chi1 > 0
    	@warn("Warning: ra - 1/chi1 > 0")
	end


	#Income process (two-state Poisson process):
	w = 4;
	Nz = 2;
	z      = [.8, 1.3]
	la_mat = [-1/3 1/3; 1/3 -1/3];

	crit = 10^(-5);
	Delta = 100;
	maxit = 35

	#grids
	I = 100;
	bmin = -2;
	#bmin = 0;
	bmax = 40;
	b = range(bmin,bmax,length=I)
	db = (bmax-bmin)/(I-1);

	J= 50;
	amin = 0;
	amax = 70;
	a = range(amin,amax,length=J);
	da = (amax-amin)/(J-1);

	bb = b * ones(1,J)
	aa = ones(I,1) * a'
	zz = ones(J,1) * z'

	dist = zeros(maxit)


	bbb = zeros(I,J,Nz)
	aaa = zeros(I,J,Nz)
	zzz = zeros(I,J,Nz)
	for nz = 1:Nz
	    bbb[:,:,nz] .= bb
	    aaa[:,:,nz] .= aa
	    zzz[:,:,nz] .= z[nz]
	end

	
	Bswitch = [
	    LinearAlgebra.I(I*J)*la_mat[1,1] LinearAlgebra.I(I*J)*la_mat[1,2];
	    LinearAlgebra.I(I*J)*la_mat[2,1] LinearAlgebra.I(I*J)*la_mat[2,2]
	]
	
	#Preallocation
	VbF = zeros(I,J,Nz);
	VbB = zeros(I,J,Nz);
	VaF = zeros(I,J,Nz);
	VaB = zeros(I,J,Nz);
	c = zeros(I,J,Nz);
	updiag = zeros(I*J,Nz);
	lowdiag = zeros(I*J,Nz);
	centdiag = zeros(I*J,Nz);
	AAi = Array{AbstractArray}(undef, Nz)
	BBi = Array{AbstractArray}(undef, Nz)

	d_B = zeros(I,J,Nz)
	d_F = zeros(I,J,Nz)
	Id_B = zeros(I,J,Nz)
	Id_F = zeros(I,J,Nz)
	c = zeros(I,J,Nz)
	u = zeros(I,J,Nz)

	
	#INITIAL GUESS
	v0 = (((1-xi)*w*zzz + ra.*aaa + rb_neg.*bbb).^(1-ga))/(1-ga)/rho
	v = copy(v0)


	#return at different points in state space
	#matrix of liquid returns
	Rb = rb_pos .* (bbb .> 0) .+ rb_neg .* (bbb .< 0)
	raa = ra .* ones(1,J)
	#if ra>>rb, impose tax on ra*a at high a, otherwise some households
	#accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
	tau = 10
	raa = ra .* (1 .- (1.33 .* amax ./ a) .^ (1-tau))#; plot(a,raa.*a)
	#matrix of illiquid returns

	Ra = zeros(I,J,Nz)
	Ra[:,:,1] .= raa'
	Ra[:,:,2] .= raa'

	for n=1:maxit
	    V = v;   
	    #DERIVATIVES W.R.T. b
	    # forward difference
	    VbF[1:I-1,:,:] .= (V[2:I,:,:] .- V[1:I-1,:,:]) ./ db;
	    VbF[I,:,:] = ((1-xi)*w*zzz[I,:,:] + Rb[I,:,:] .* bmax).^(-ga); #state constraint boundary condition
			
	    # backward difference
	    VbB[2:I,:,:] = (V[2:I,:,:]-V[1:I-1,:,:])/db;
	    VbB[1,:,:] = ((1-xi)*w*zzz[1,:,:] + Rb[1,:,:].*bmin).^(-ga); #state constraint boundary condition
	
	    #DERIVATIVES W.R.T. a
	    # forward difference
	    VaF[:,1:J-1,:] = (V[:,2:J,:]-V[:,1:J-1,:])/da;
	    # backward difference
	    VaB[:,2:J,:] = (V[:,2:J,:]-V[:,1:J-1,:])/da;
	 
	    #useful quantities
	    c_B = max.(VbB,10^(-6)).^(-1/ga);
	    c_F = max.(VbF,10^(-6)).^(-1/ga); 
		dBB = two_asset_kinked_FOC.(VaB,VbB,aaa, Ref((; chi0, chi1)))
	    dFB = two_asset_kinked_FOC.(VaB,VbF,aaa, Ref((; chi0, chi1)))
	    #VaF(:,J,:) = VbB(:,J,:).*(1-ra.*chi1 - chi1*w*zzz(:,J,:)./a(:,J,:));
	    dBF = two_asset_kinked_FOC.(VaF,VbB,aaa, Ref((; chi0, chi1)))
	    #VaF(:,J,:) = VbF(:,J,:).*(1-ra.*chi1 - chi1*w*zzz(:,J,:)./a(:,J,:));
	    dFF = two_asset_kinked_FOC.(VaF,VbF,aaa, Ref((; chi0, chi1)))
	    
	    #UPWIND SCHEME
	    d_B .= (dBF .> 0) .* dBF .+ (dBB .< 0) .* dBB;
	   
		#state constraints at amin and amax
	    d_B[:,1,:] = (dBF[:,1,:] .> 10^(-12)) .* dBF[:,1,:] #make sure d>=0 at amax, don't use VaB(:,1,:)
	    d_B[:,J,:] = (dBB[:,J,:] .< -10^(-12)) .* dBB[:,J,:] #make sure d<=0 at amax, don't use VaF(:,J,:)
	    d_B[1,1,:] = max.(d_B[1,1,:],0)
	    #split drift of b and upwind separately
	    sc_B = (1-xi) .* w .* zzz .+ Rb .* bbb .- c_B;
	    sd_B = (-d_B - two_asset_kinked_cost.(d_B, aaa, Ref((; chi0, chi1))))
	    
	    d_F .= (dFF .> 0) .* dFF + (dFB .< 0) .* dFB
	    #state constraints at amin and amax
	    d_F[:,1,:] = (dFF[:,1,:] .> 10^(-12)) .* dFF[:,1,:] #make sure d>=0 at amin, don't use VaB(:,1,:)
	    d_F[:,J,:] = (dFB[:,J,:] .< -10^(-12)) .* dFB[:,J,:] #make sure d<=0 at amax, don't use VaF(:,J,:)
	
	    #split drift of b and upwind separately
	    sc_F = (1-xi)*w*zzz .+ Rb.*bbb .- c_F;
	    sd_F = (-d_F .- two_asset_kinked_cost.(d_F,aaa, Ref((; chi0, chi1))));
	    sd_F[I,:,:] = min.(sd_F[I,:,:],0.0)
	    
	    Ic_B = (sc_B .< -10^(-12))
	    Ic_F = (sc_F .> 10^(-12)) .* (1 .- Ic_B)
	    Ic_0 = 1 .- Ic_F .- Ic_B
	    
	    Id_F .= (sd_F .> 10^(-12))
	    Id_B .= (sd_B .< -10^(-12)) .* (1 .- Id_F)
	    Id_B[1,:,:] .= 0
	    Id_F[I,:,:] .= 0
		Id_B[I,:,:] .= 1 #don't use VbF at bmax so as not to pick up articial state constraint
	    Id_0 = 1 .- Id_F .- Id_B
	    
	    c_0 = (1-xi) * w * zzz + Rb .* bbb
	  
	    c .= c_F .* Ic_F + c_B .* Ic_B + c_0 .* Ic_0
	    u .= c .^ (1-ga) ./(1-ga)
	    
	    #CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
	    X = -Ic_B .* sc_B ./db .- Id_B .* sd_B ./ db
	    Y = (Ic_B .* sc_B .- Ic_F .* sc_F) ./db .+ (Id_B .* sd_B .- Id_F .* sd_F) ./db;
	    Z = Ic_F.*sc_F/db + Id_F.*sd_F/db;
	    
	    for i = 1:Nz
	        centdiag[:,i] = reshape(Y[:,:,i],I*J,1)
	    end
	
	    lowdiag[1:I-1,:] = X[2:I,1,:]
	    updiag[2:I,:] = Z[1:I-1,1,:]
	    for j = 2:J
	        lowdiag[1:j*I,:] = [lowdiag[1:(j-1)*I,:]; X[2:I,j,:]; zeros(1,Nz)]
	        updiag[1:j*I,:] = [updiag[1:(j-1)*I,:]; zeros(1,Nz); Z[1:I-1,j,:]];
	    end
	
	    for nz = 1:Nz
	    	BBi[nz] = spdiagm(
				I*J, I*J, 
				0 => centdiag[:,nz],
				1 => updiag[2:end,nz],
				-1 => lowdiag[1:end-1,nz]
						 )
	    end
	
	    BB = cat(BBi..., dims = (1,2))
	
	
	    #CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
	    dB = Id_B .* dBB .+ Id_F .* dFB
	    dF = Id_B .* dBF .+ Id_F .* dFF
	    MB = min.(dB,0.0)
	    MF = max.(dF,0.0) .+ xi .* w .* zzz .+ Ra .* aaa
	    MB[:,J,:] = xi .* w .* zzz[:,J,:] .+ dB[:,J,:] .+ Ra[:,J,:] .* amax #this is hopefully negative
	    MF[:,J,:] .= 0.0
	    chi = -MB ./ da
	    yy =  (MB - MF) ./da
	    zeta = MF ./ da
	
	    # MATRIX AAi
	    for nz=1:Nz
	        #This will be the upperdiagonal of the matrix AAi
	        AAupdiag = zeros(I,1); #This is necessary because of the peculiar way spdiags is defined.
	        for j=1:J
	            AAupdiag=[AAupdiag; zeta[:,j,nz]]
	        end
	        
	        #This will be the center diagonal of the matrix AAi
	        AAcentdiag = yy[:,1,nz]
	        for j=2:J-1
	            AAcentdiag = [AAcentdiag; yy[:,j,nz]];
	        end
	        AAcentdiag = [AAcentdiag; yy[:,J,nz]];
	        
	        #This will be the lower diagonal of the matrix AAi
	        AAlowdiag = chi[:,2,nz]
	        for j=3:J
	            AAlowdiag = [AAlowdiag; chi[:,j,nz]]
	        end
	
			#@info AAcentdiag
			#@info AAlowdiag
			#@info AAupdiag
	
		    #Add up the upper, center, and lower diagonal into a sparse matrix
	        AAi[nz] = spdiagm(
				I*J, I*J,
				0 => AAcentdiag,
				-I => AAlowdiag,
				I => AAupdiag[begin+I:end-I]
			)
	
	    end
	
		AA = cat(AAi..., dims = (1,2))
	    
	    A = AA + BB + Bswitch
	
		
	    if maximum(abs, sum(A,dims=2)) > 10^(-12)
	        @warn("Improper Transition Matrix")
	        break
	    end
	    
	#    if maximum(abs, sum(A, dims=2)) > 10^(-9)
	#       @warn("Improper Transition Matrix")
	#       break
	#    end
	    
	    B = (1/Delta + rho)*LinearAlgebra.I(I*J*Nz) - A
	    
	    u_stacked = reshape(u,I*J*Nz,1)
	    V_stacked = reshape(V,I*J*Nz,1)
	    
	    vec = u_stacked + V_stacked/Delta;
	    
	    V_stacked = B\vec #SOLVE SYSTEM OF EQUATIONS
	        
	    V = reshape(V_stacked,I,J,Nz)   
	    
	    
	    Vchange = V - v
	    v = V
	    	   
	    dist[n] = maximum(abs, Vchange)
	    @info "Value Function, Iteration $n | max Vchange = $(dist[n])"
	    if dist[n]<crit
	        @info("Value Function Converged, Iteration = $n")
	        break
	    end 
	end

	d = Id_B .* d_B + Id_F .* d_F
	m = d + xi*w*zzz + Ra.*aaa;
	s = (1-xi)*w*zzz + Rb.*bbb - d - two_asset_kinked_cost.(d, aaa, Ref((; chi0, chi1))) - c

	sc = (1-xi)*w*zzz + Rb.*bbb - c;
	sd = - d - two_asset_kinked_cost.(d,aaa, Ref((; chi0, chi1)))

	df = DataFrame(
		a = vec(aaa),
		b = vec(bbb),
		z = vec(zzz),
		c = vec(c), 
		d = vec(d),
		s = vec(s),
		m = vec(m),
		u = vec(u),
		sc = vec(sc),
		sd = vec(sd),
		v = vec(v)
	)
end

# ╔═╡ 75c6bed0-86a2-4393-83a7-fbd7862a3975
df_base = solve_HJB_base()

# ╔═╡ b8a7269f-27ae-4c37-b73e-7decb8333ea9
# ╠═╡ disabled = true
#=╠═╡
@chain df begin
	stack([:s, :d])
	data(_) * mapping(:b, :a, :value, col = :z => nonnumeric, row = :variable) * visual(Surface)
	draw(axis = (type = Axis3, ))
end
  ╠═╡ =#

# ╔═╡ 830448b7-1700-4312-91ce-55f86aaa33a4
md"""
## HJB Greimel
"""

# ╔═╡ ed6045c0-b76c-4691-9f05-c943c542d13f
Base.@kwdef struct TwoAssets
	γ = 2 #CRRA utility with parameter gamma
	ra = 0.05
	rb_pos = 0.03
	rb_neg = 0.12
	rho = 0.06 #discount rate
	χ₀ = 0.03
	χ₁ = 2
	xi = 0.1 #fraction of income that is automatically deposited
	#Income process (two-state Poisson process):
	w = 4
	Nz = 2
	z      = [.8, 1.3]
	la_mat = [-1/3 1/3; 1/3 -1/3]
	crit = 10^(-5)
	Delta = 100
	#grids
	I = 100
	bmin = -2
	bmax = 40
	b = range(bmin,bmax,length=I)
	db = (bmax-bmin)/(I-1)
	J = 50
	amin = 0
	amax = 70
	a = range(amin,amax,length=J)
	da = (amax-amin)/(J-1)
	τ = 10
end

# ╔═╡ 03ec6276-09a4-4f66-a864-19e2e0d825eb
md"""
if ``r_a \gg r_b``, impose tax on ``ra \cdot a`` at high ``a``, otherwise some households accumulate infinite illiquid wealth (not needed if ``r_a`` is close to or less than ``r_b``)
"""

# ╔═╡ 98f11cbd-8b05-464b-be44-3b76277c6d0d
R_a(a, (; ra, amax, τ)) = ra * (1 - (1.33 * amax / a) ^ (1-τ)) * a

# ╔═╡ 8a6cd6dd-49f6-496a-bb50-e43e06c4a1db
R_b(b, (; rb_pos, rb_neg)) = (b ≥ 0 ? rb_pos : rb_neg) * b

# ╔═╡ a6356900-5530-494f-9d01-03041805ebe6
function check((; ra, χ₁))
	if ra - 1/χ₁ > 0
    	@warn("Warning: ra - 1/χ₁ > 0")
	end
end

# ╔═╡ f6c0329e-1a02-4292-96b7-11c8cd9c3b54
util(c, (; γ)) = c^(1-γ)/(1-γ)

# ╔═╡ 5a6c37d8-af66-46a1-9c93-581057c41f94
u_prime(c, (; γ)) = c^(-γ)

# ╔═╡ cdcca65f-59df-4990-97b7-2b511bdb61e6
u_prime_inv(x, (; γ)) = x^(-1/γ)

# ╔═╡ 920e3393-fd38-4154-8d90-ce9dc712ed1a
function get_d(VaB, VaF, Vb, a, b, model)
	(; amax, amin, bmin, bmax) = model
	dxB = two_asset_kinked_FOC_new(VaB,Vb,a, model)
	dxF = two_asset_kinked_FOC_new(VaF,Vb,a, model)

	if dxF > 0 && dxB < 0
		d = dxF + dxB
	elseif dxF > 0
		d = dxF
	elseif dxB < 0
		d = dxB
	else
		d = 0.0
	end

	if a == amin
		d = dxF * (dxF > 10^(-12))
	end
	if a == amax
		d = dxB * (dxB < -10^(-12))
	end
	if a == amin && b == bmin
		d = max(d, 0.0)
	end

	sd = -d - two_asset_kinked_cost_new(d, a, model)

	(; d, dxB, dxF, sd)
end

# ╔═╡ 9c1544ef-fdc9-4c9e-a36c-fe5b3ea89728
function get_d_upwind(VaB, VaF, VbB, VbF, a, b, model)
	(; amax, amin, bmin, bmax) = model

	outB = get_d(VaB, VaF, VbB, a, b, model)
	outF = get_d(VaB, VaF, VbF, a, b, model)

	d_B = outB.d
	d_F = outF.d
	sd_B = outB.sd
	sd_F = outF.sd
	
	if b == bmax
		sd_F = min(sd_F, 0.0)
	end

	Id_F = false
	Id_B = false
	Id_0 = false

	if b == bmax
		Id_B = true
		dB = outB.dxB
		dF = outB.dxF
		d = d_B
	elseif sd_F > 10^(-12) #&& b < bmax
		Id_F = true
		dB = outF.dxB
		dF = outF.dxF
		d = d_F
	elseif sd_B < -10^(-12) && b > bmin
		Id_B = true
		dB = outB.dxB
		dF = outB.dxF
		d = d_B
	else
		Id_0 = true
		dB = 0.0
		dF = 0.0
		d = 0.0
	end

	(; d_B, d_F, dB, dF, sd_B, sd_F, Id_B, Id_F, Id_0, d)
end	

# ╔═╡ 4023a748-5e4f-4137-811b-0e93567021dd
function get_c(Vb, b, z, model)
	(; γ, xi, w) = model
	c = u_prime_inv.(max.(Vb,10^(-6)), Ref((; γ)))
	sc = (1-xi) * w * z + R_b(b, model) - c

	(; c, sc)
end

# ╔═╡ 541a5b4f-0d49-4c36-85c8-799e3260c77d
function get_c₀(b, z, model)
	(; γ, xi, w) = model
	sc = 0.0
	c = (1-xi) * w * z + R_b(b, model)
	(; c, sc)
end

# ╔═╡ 09a8d045-6867-43a9-a021-5a39c22171a5
function get_c_upwind(VbB, VbF, b, z, model)
	out_F = get_c(VbF, b, z, model)
	if out_F.sc > 10^(-12)
		return out_F
	end
	out_B = get_c(VbB, b, z, model)
	if out_B.sc < -10^(-12)
		return out_B
	end
	return get_c₀(b, z, model)
end

# ╔═╡ 2fb709ca-5327-41e4-916b-4a0098859c3e
function solve_HJB_new(model, maxit = 35)
	(; rho, ra, rb_pos, rb_neg, xi, w) = model
	(; Delta, crit) = model
	(; a, b, z, I, J, Nz, la_mat, amin, amax, bmin, bmax, db, da) = model

	bb = b * ones(1,J)
	aa = ones(I,1) * a'
	zz = ones(J,1) * z'

	dist = zeros(maxit)


	bbb = zeros(I,J,Nz)
	aaa = zeros(I,J,Nz)
	zzz = zeros(I,J,Nz)
	for nz = 1:Nz
	    bbb[:,:,nz] .= bb
	    aaa[:,:,nz] .= aa
	    zzz[:,:,nz] .= z[nz]
	end

	
	Bswitch = [
	    LinearAlgebra.I(I*J)*la_mat[1,1] LinearAlgebra.I(I*J)*la_mat[1,2];
	    LinearAlgebra.I(I*J)*la_mat[2,1] LinearAlgebra.I(I*J)*la_mat[2,2]
	]
	
	#Preallocation
	VbF = zeros(I,J,Nz);
	VbB = zeros(I,J,Nz);
	VaF = zeros(I,J,Nz);
	VaB = zeros(I,J,Nz);
	c = zeros(I,J,Nz);
	updiag = zeros(I*J,Nz);
	lowdiag = zeros(I*J,Nz);
	centdiag = zeros(I*J,Nz);
	AAi = Array{AbstractArray}(undef, Nz)
	BBi = Array{AbstractArray}(undef, Nz)

	d = zeros(I,J,Nz)
	#d_F = zeros(I,J,Nz)
	#Id_B = zeros(I,J,Nz)
	#Id_F = zeros(I,J,Nz)
	#Id_0 = zeros(I,J,Nz)
	c = zeros(I,J,Nz)
	u = zeros(I,J,Nz)
	
	#INITIAL GUESS
	v0 = util.((1-xi)*w*zzz + ra.*aaa + rb_neg.*bbb, Ref(model)) ./ rho
	v = copy(v0)


	#return at different points in state space
	#matrix of liquid returns
	Rb = rb_pos .* (bbb .> 0) .+ rb_neg .* (bbb .< 0)
	raa = ra .* ones(1,J)
	#if ra>>rb, impose tax on ra*a at high a, otherwise some households
	#accumulate infinite illiquid wealth (not needed if ra is close to or less than rb)
	tau = 10
	raa = ra .* (1 .- (1.33 .* amax ./ a) .^ (1-tau))#; plot(a,raa.*a)
	#matrix of illiquid returns

	Ra = zeros(I,J,Nz)
	Ra[:,:,1] .= raa'
	Ra[:,:,2] .= raa'

	for n=1:maxit
	    V = v;   
	    #DERIVATIVES W.R.T. b
	    # forward difference
	    VbF[1:I-1,:,:] .= (V[2:I,:,:] .- V[1:I-1,:,:]) ./ db;
	    VbF[I,:,:] = u_prime.((1-xi)*w*zzz[I,:,:] + Rb[I,:,:] .* bmax, Ref(model)) #state constraint boundary condition
			
	    # backward difference
	    VbB[2:I,:,:] = (V[2:I,:,:]-V[1:I-1,:,:])/db;
	    VbB[1,:,:] = u_prime.((1-xi)*w*zzz[1,:,:] + Rb[1,:,:].*bmin, Ref(model)) #state constraint boundary condition
	
	    #DERIVATIVES W.R.T. a
	    # forward difference
	    VaF[:,1:J-1,:] = (V[:,2:J,:]-V[:,1:J-1,:])/da;
	    # backward difference
	    VaB[:,2:J,:] = (V[:,2:J,:]-V[:,1:J-1,:])/da;
	 
		out_d = get_d_upwind.(VaB, VaF, VbB, VbF, aaa, bbb, Ref(model)) |> StructArray
		(; d_B, d_F, Id_B, Id_F, sd_B, sd_F) = out_d
		d .= out_d.d

		out_c = get_c_upwind.(VbB, VbF, bbb, zzz, Ref(model)) |> StructArray
		(; sc) = out_c

	    c .= out_c.c
	    u .= util.(c, Ref(model))
	    
	    #CONSTRUCT MATRIX BB SUMMARING EVOLUTION OF b
	    X = - min.(sc, 0.0) ./db .- Id_B .* sd_B ./ db
	    Y = (min.(sc, 0.0) - max.(sc, 0.0)) ./db .+ (Id_B .* sd_B .- Id_F .* sd_F) ./db;
	    Z = max.(sc, 0.0)/db + Id_F.*sd_F/db;
	    
	    for i = 1:Nz
	        centdiag[:,i] = reshape(Y[:,:,i],I*J,1)
	    end
	
	    lowdiag[1:I-1,:] = X[2:I,1,:]
	    updiag[2:I,:] = Z[1:I-1,1,:]
	    for j = 2:J
	        lowdiag[1:j*I,:] = [lowdiag[1:(j-1)*I,:]; X[2:I,j,:]; zeros(1,Nz)]
	        updiag[1:j*I,:] = [updiag[1:(j-1)*I,:]; zeros(1,Nz); Z[1:I-1,j,:]];
	    end
	
	    for nz = 1:Nz
	    	BBi[nz] = spdiagm(
				I*J, I*J, 
				0 => centdiag[:,nz],
				1 => updiag[2:end,nz],
				-1 => lowdiag[1:end-1,nz]
						 )
	    end
	
	    BB = cat(BBi..., dims = (1,2))
	
	
	    #CONSTRUCT MATRIX AA SUMMARIZING EVOLUTION OF a
		(; dB, dF) = out_d
	    MB = min.(dB,0.0)
	    MF = max.(dF,0.0) .+ xi .* w .* zzz .+ Ra .* aaa
	    MB[:,J,:] = xi .* w .* zzz[:,J,:] .+ dB[:,J,:] .+ Ra[:,J,:] .* amax #this is hopefully negative
	    MF[:,J,:] .= 0.0
	    chi = -MB ./ da
	    yy =  (MB - MF) ./da
	    zeta = MF ./ da
	
	    # MATRIX AAi
	    for nz=1:Nz
	        #This will be the upperdiagonal of the matrix AAi
	        AAupdiag = zeros(I,1); #This is necessary because of the peculiar way spdiags is defined.
	        for j=1:J
	            AAupdiag=[AAupdiag; zeta[:,j,nz]]
	        end
	        
	        #This will be the center diagonal of the matrix AAi
	        AAcentdiag = yy[:,1,nz]
	        for j=2:J-1
	            AAcentdiag = [AAcentdiag; yy[:,j,nz]];
	        end
	        AAcentdiag = [AAcentdiag; yy[:,J,nz]];
	        
	        #This will be the lower diagonal of the matrix AAi
	        AAlowdiag = chi[:,2,nz]
	        for j=3:J
	            AAlowdiag = [AAlowdiag; chi[:,j,nz]]
	        end
		
		    #Add up the upper, center, and lower diagonal into a sparse matrix
	        AAi[nz] = spdiagm(
				I*J, I*J,
				0 => AAcentdiag,
				-I => AAlowdiag,
				I => AAupdiag[begin+I:end-I]
			)
	
	    end
	
		AA = cat(AAi..., dims = (1,2))
	    
	    A = AA + BB + Bswitch
	
		
	    if maximum(abs, sum(A,dims=2)) > 10^(-12)
	        @warn("Improper Transition Matrix")
	        break
	    end
	    
	#    if maximum(abs, sum(A, dims=2)) > 10^(-9)
	#       @warn("Improper Transition Matrix")
	#       break
	#    end
	    
	    B = (1/Delta + rho)*LinearAlgebra.I(I*J*Nz) - A
	    
	    u_stacked = reshape(u,I*J*Nz,1)
	    V_stacked = reshape(V,I*J*Nz,1)
	    
	    vec = u_stacked + V_stacked/Delta;
	    
	    V_stacked = B\vec #SOLVE SYSTEM OF EQUATIONS
	        
	    V = reshape(V_stacked,I,J,Nz)   
	    
	    
	    Vchange = V - v
	    v .= V
	    	   
	    dist[n] = maximum(abs, Vchange)
	    @info "Value Function, Iteration $n | max Vchange = $(dist[n])"
	    if dist[n]<crit
	        @info("Value Function Converged, Iteration = $n")
	        break
	    end 
	end

	#d = Id_B .* d_B + Id_F .* d_F
	m = d + xi*w*zzz + Ra.*aaa;
	s = (1-xi)*w*zzz + Rb.*bbb - d - two_asset_kinked_cost_new.(d, aaa, Ref(model)) - c

	sc = (1-xi)*w*zzz + Rb.*bbb - c;
	sd = - d - two_asset_kinked_cost_new.(d,aaa, Ref(model))

	df = DataFrame(
		a = vec(aaa),
		b = vec(bbb),
		z = vec(zzz),
		c = vec(c), 
		d = vec(d),
		s = vec(s),
		m = vec(m),
		u = vec(u),
		sc = vec(sc),
		sd = vec(sd),
		v = vec(v)
	)
end

# ╔═╡ 48ecc7ee-b943-4497-bc15-1b62d78e9271
begin
	m = TwoAssets()
	df_new = solve_HJB_new(m)
end

# ╔═╡ 8f3da06b-2887-4564-87c7-12a798580f53
@test df_new.v ≈ df_base.v

# ╔═╡ 0cf5bd0e-51e8-437b-bea6-2027b898a579
@test df_new.c ≈ df_base.c

# ╔═╡ 87074572-87c8-4f01-8895-a8ebcbeef9a0
@test df_new.d ≈ df_base.d

# ╔═╡ 4fe6acdd-521e-44b1-a7b6-97f78f72986d
@test df_new.s ≈ df_base.s

# ╔═╡ 3a43bab0-c058-4665-8939-a3920c9986d1
md"""
# Appendix
"""

# ╔═╡ 311f99f6-baeb-402e-8512-999d91829ec9
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
AlgebraOfGraphics = "~0.6.9"
CairoMakie = "~0.8.8"
Chain = "~0.5.0"
DataFrameMacros = "~0.2.1"
DataFrames = "~1.3.4"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.39"
StructArrays = "~0.6.11"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"

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
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "8a8f4d8eddc2e8c4ab71c1855b91b7d762ef05fe"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.9"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"

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
git-tree-sha1 = "76d499235febafad126b229e8d5027800a91874d"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.8"

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
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "2dd813e5f2f7eec2d1268c57cf2373d3ee91fcea"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.1"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

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

[[deps.Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "9be8be1d8a6f44b96482c8af52238ea7987da3e3"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.45.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

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
git-tree-sha1 = "d530092b57aef8b96b27694e51c575b09c7f0b2e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.64"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "b19534d1895d702889b219c382a6e18010797f0b"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.8.6"

[[deps.Downloads]]
deps = ["ArgTools", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"

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
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

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
git-tree-sha1 = "9267e5f50b0e12fdfd5a2455534345c4cf2c7f7a"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.14.0"

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
git-tree-sha1 = "a88992eaa3073e65c970ce73f4636c080b68e21e"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.7.9"

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

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Serialization", "Showoff", "SignedDistanceFields", "SparseArrays", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "UnicodeFun"]
git-tree-sha1 = "b0946fd8f4f981210980bef0a7ed63ab5fb4206f"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.8"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "469221640e5e798b52877fd12c596204cee05df1"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.3.4"

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

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "114ef48a73aea632b8aebcb84f796afcc510ac7c"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

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

[[deps.NaNMath]]
git-tree-sha1 = "737a5957f387b17e74d4ad2f440eb330b39a62c5"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.0"

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
git-tree-sha1 = "9a36165cf84cff35851809a40a928e1103702013"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.16+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

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
git-tree-sha1 = "ca433b9e2f5ca3a0ce6702a032fce95a3b6e1e48"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.14"

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

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

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

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

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

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.SIMD]]
git-tree-sha1 = "7dbc15af7ed5f751a82bf3ed37757adf76c32402"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.1"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "8c3e2c64dac132efa8828b1b045a47cbf0881def"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.2"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "0b4b7f1393cff97c33891da2a0bf69c6ed241fda"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "9f8a5dc5944dc7fbbe6eb4180660935653b0a9d9"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.0"

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
git-tree-sha1 = "48598584bacbebf7d30e20880438ed1d24b7c7d6"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.18"

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

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"

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
git-tree-sha1 = "fcf41697256f2b759de9380a7e8196d6516f0310"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.0"

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
# ╠═c89f1918-01ee-11ed-22fa-edd66e0f6c59
# ╟─629b8291-0f13-419e-b1c0-d10d5e708720
# ╠═b48b0674-1bcb-48e5-9b05-57dea5877715
# ╟─828dee22-1ee7-41c4-b68b-f88facea86d9
# ╠═68a97aab-7924-4472-aa47-7903add8aea4
# ╠═3ab0c985-5317-4ea4-bddc-6289ab90bcad
# ╠═f2ce6352-450e-4cde-a2fe-3586461c3bdf
# ╠═9c1ef8d1-57bc-4da2-83ec-fb8f1a8ce296
# ╠═91b63bfc-f4a4-41c1-a472-7d13df27b93c
# ╠═75c6bed0-86a2-4393-83a7-fbd7862a3975
# ╠═b8a7269f-27ae-4c37-b73e-7decb8333ea9
# ╟─830448b7-1700-4312-91ce-55f86aaa33a4
# ╠═ed6045c0-b76c-4691-9f05-c943c542d13f
# ╟─03ec6276-09a4-4f66-a864-19e2e0d825eb
# ╠═98f11cbd-8b05-464b-be44-3b76277c6d0d
# ╠═8a6cd6dd-49f6-496a-bb50-e43e06c4a1db
# ╠═a6356900-5530-494f-9d01-03041805ebe6
# ╠═f6c0329e-1a02-4292-96b7-11c8cd9c3b54
# ╠═5a6c37d8-af66-46a1-9c93-581057c41f94
# ╠═cdcca65f-59df-4990-97b7-2b511bdb61e6
# ╠═920e3393-fd38-4154-8d90-ce9dc712ed1a
# ╠═9c1544ef-fdc9-4c9e-a36c-fe5b3ea89728
# ╠═4023a748-5e4f-4137-811b-0e93567021dd
# ╠═541a5b4f-0d49-4c36-85c8-799e3260c77d
# ╠═09a8d045-6867-43a9-a021-5a39c22171a5
# ╠═1ee8ccc3-d498-48c6-b299-1032165e4ab9
# ╠═2fb709ca-5327-41e4-916b-4a0098859c3e
# ╠═48ecc7ee-b943-4497-bc15-1b62d78e9271
# ╠═8f3da06b-2887-4564-87c7-12a798580f53
# ╠═0cf5bd0e-51e8-437b-bea6-2027b898a579
# ╠═87074572-87c8-4f01-8895-a8ebcbeef9a0
# ╠═4fe6acdd-521e-44b1-a7b6-97f78f72986d
# ╟─3a43bab0-c058-4665-8939-a3920c9986d1
# ╠═9aa61364-51a3-45d0-b1c2-757b864de132
# ╠═026cfe16-ff0f-4f68-b412-b1f6c1902824
# ╠═9d99416e-8119-4a40-b577-0135050a0e4e
# ╠═6188aab9-86bf-4ec4-bb10-43b59f71e3e2
# ╠═9696e6ca-6953-43e2-8d47-fbfe24ba4250
# ╠═311f99f6-baeb-402e-8512-999d91829ec9
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

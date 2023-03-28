### A Pluto.jl notebook ###
# v0.18.4

using Markdown
using InteractiveUtils

# ╔═╡ 8d268ee9-6a3b-4bf3-a0d1-59d86fe8f263
md"""
`social-connectedness.jl` | **Version 1.2** | _last updated on June 20, 2022_
"""

# ╔═╡ 662cad72-af63-41c9-a895-5633be486f3f
md"""
# Assignment: Social networks and the housing market

The goal is to replicate some of the analysis in *[__The Economic Effects of Social Networks: Evidence from the Housing Market__](https://www.journals.uchicago.edu/doi/abs/10.1086/700073) (Bailey, Cao, Kuchler & Stroebel, 2018, JPE)* using publicly available (county-level) data.

We will use 

* Zillow house price data
* HMDA mortgage data
* the Facebook Social Connectedness Index
"""

# ╔═╡ cdfa2f05-45cc-4050-baf4-15d66910c4e9
md"""
## Part 1: Friends' house price experiences (3 points)

For a pair of regions ``(A, B)`` the Facebook Social Connectedness is defined as follows

```math
\text{SCI}_{AB} = \frac{\text{\# FB friendships}_{AB}}{\text{\# FB users}_A \cdot \text{\#FB users}_B}
```

For ``A\neq B`` the index measures how likely two random users ``(a, b) \in A \times B`` are friends with each other. Under the plausible assumption that Facebook friendships are representative for all friendships, we can use this measure to approxmiate how socially connected two regions are.

We use the SCI to calculate the average house price change in the home counties of county ``c``'s friends. If ``p_{\tilde ct}`` is the change in log-house prices in period ``t``, the average friend of an inhabitant of county ``c`` will have experienced
```math
p^\text{friends}_{ct} =\frac{1}{\sum_c \text{\# friendships}_{c\tilde c}}\sum_c \text{\# friendships}_{c\tilde c} \cdot p_{\tilde ct}.
```
``p^\text{friends}_{ct}`` is a weighted average of house price changes, using friendship links as weights.

The figure below shows that there is a strong correlation between ``p_{ct}`` and ``p_{ct}^{\text{friends}}``.
"""

# ╔═╡ 2bc4cd24-a12e-46c2-97d4-03d0bdaca3a3
@chain zillow_and_friends_df begin
	@transform(:pop = log(:population2010/1000))
	@subset(!ismissing(:Δ_log_hpi))
	disallowmissing!
	@subset(:year ∈ 2001:5:2021)
	data(_) * mapping(:Δ_log_hpi,
		:friends_exp, # all friends, incl own county
		#:friends_exp_000, # distant friends (all other counties)
		#:friends_exp_100, # distant friends (> 100 miles away)
		#:friends_exp_300, # distant friends (> 300 miles away)
		layout = :year=> nonnumeric, markersize = :pop) * visual(Scatter)
	draw
end

# ╔═╡ c3b70e50-420b-4956-9c36-d4596fe07634
@mdx """
 Look at lines 7–10 of the code that generates this plot. The variable 
* `friends_exp_000` omits friendship links within the county,
* `friends_exp_x00` omits friendship links when the other county is ``≤ x00`` miles away.
* <details> <summary> Click here if you are wondering why you should care about these variables. </summary> 

	These variables can be used as instruments to mimic the IV strategy in the paper:
	
	> To exploit only variation in friends’ house price experiences that is orthogonal to
	> a person’s own experiences, we instrument for the house price experiences of all
	> friends with the experiences of only her friends in geographically distant housing
	> markets.
	
	</details>


👉 Look at the variations of the plot using `friends_exp_x00`, describe what you discover and try to make sense of it.
"""

# ╔═╡ 6959807f-03c5-4a9d-8c4d-c04c06c6e0a5
md"""
Your answer goes here ...
"""

# ╔═╡ 63d1da26-9ac4-411f-8f47-02a02b37f142
md"""
## Part 2: Social connectedness and the housing market (7 points)

We want to see if we can replicate the main results of _Bailey et al. (2018)._ Do friends' house price experiences drive mortgage choices?

Let counties be indexed by ``c`` and time ``t`` run from 2001 to 2016. We run the following regressions.

```math
\begin{align}
\text{originations}_{ct} &= \alpha \cdot \Delta\log(\text{hpi}_{ct}) + \beta \cdot \Delta\log(\text{hpi}^{\text{friends}}_{ct}) + \delta_c + \delta_t + \varepsilon_{ct} \\
\log(\text{amount}_{ct}) &= \alpha \cdot \Delta\log(\text{hpi}_{ct}) + \beta \cdot \Delta\log(\text{hpi}^{\text{friends}}_{ct}) + \delta_c + \delta_t + \varepsilon_{ct}
\end{align}
```

* ``\text{originations}_{ct}`` is the number (the count) of mortgages originated (normalized by population) in county ``c`` and year ``t``
* ``\text{amount}_{ct}`` is the average size of originated mortgages in county ``c`` and year ``t``
* ``\delta_c`` and ``\delta_t`` are region and time fixed effects.
"""

# ╔═╡ dbba5ead-cd2a-4714-9e56-dd34c7d72039
md"""
👉 Look at the output of the regressions below. Interpret the results. Are you convinced? You might want to play around with the specification to check robustness.

(It might be helpful to read the introduction of [_Bailey et al. (2018)_](https://www.journals.uchicago.edu/doi/abs/10.1086/700073) and refer to it in your discussion.)
"""

# ╔═╡ 8001cd65-2542-4a1a-93b8-2e7243226956
md"""
Your answer goes here ...
"""

# ╔═╡ 2357bd9e-93e2-468e-8bd4-7720f07eed5e
using FixedEffectModels

# ╔═╡ 2eb3a5db-f5a3-489f-b712-146934229772
reg(joined_df, 
	@formula(amount ~ Δ_log_hpi + (friends_exp_000 ~ friends_exp_100) + fe(year)),
	#Vcov.cluster(:state),
	weights = :population2010
)

# ╔═╡ 7297dc74-44b4-4d2d-98f5-7120b74cb196
reg(joined_df, 
	@formula(count_per_pop ~ Δ_log_hpi + (friends_exp_000 ~ friends_exp_100) + fe(year)),
	#Vcov.cluster(:state),
	weights = :population2010
)

# ╔═╡ ca4b0557-9360-4089-ab64-a62e07d23de9
joined_df = @chain hmda_df begin
	@select(:income, :amount, :count, :year, :fips, :state)
	innerjoin(zillow_and_friends_df, on = [:year, :fips])
	@transform!(:count_per_pop = :count / :population2010)
end

# ╔═╡ 24013463-f8ea-463f-81e9-06542880dc8b
md"""
# Appendix: Getting all the data
"""

# ╔═╡ cc13875f-0407-4749-a30e-39737c28c8b3
md"""
## Getting the shapes for plotting maps
"""

# ╔═╡ aa8593ef-3f6f-4b0c-8046-f88aef680918
using GeoTables

# ╔═╡ 1af8b2a8-f7fb-4afa-a746-cedbe1fea1a6
using GeoTables.Meshes

# ╔═╡ 301662e7-2d65-4b3b-9128-c415918e0058
using AlgebraOfGraphics.Makie.GeometryBasics

# ╔═╡ 94686899-25fd-4d28-9894-5dcdc2650efe
shapes0 = GeoTables.gadm("USA", depth=2);

# ╔═╡ e8183f41-5921-4bfc-a262-6ac6b165870c
shapes_df = @chain shapes0 begin
	DataFrame
	@subset(:TYPE_2 ∉ ["Water body"])
	@subset!(:NAME_1 ∉ ["Alaska", "Hawaii"])
	@select!(
		:center = point_to_point(Meshes.centroid(:geometry)),		
		#:area = Meshes.area(:geometry),
		:geometry = multi_to_multi(:geometry),
		:state_name = :NAME_1, #:county_name = :NAME_2,
		#:GID_2 => ByRow(GID2_to_fips) => AsTable
		:county_matching = clean_county_name_for_matching(:NAME_2)
	)
end

# ╔═╡ 85768c95-3f0c-43ad-9e37-c479c9bb6a5b
shapes_pop_df = @chain pop_df begin
	@select(:fips, :state, :county, :population2010, :state_name, :county_name,
		:county_matching = clean_county_name_for_matching(:county_name)
	)
	innerjoin(shapes_df, on = [:state_name, :county_matching], makeunique=true)
end

# ╔═╡ edf2ff0d-ab91-47ff-ae03-cee7507601a1
function clean_county_name_for_matching(county)
	county_match = county
	
	repl = [
		" County"      => "",
		" Parish"      => "",
		" City"        => "", " city"    => "",
		" and Borough" => "", " Borough" => "",
		" Census Area" => "",
		" Municipality"=> "",
		"\xf1"         => "n", # ñ => n
		"St."          => "Saint",
		"Ste."         => "Sainte",
		"Oglala Lakota"=> "Shannon",
		" " => ""
	]
		
	for r in repl
		county_match = replace(county_match, r)
	end
		
	county_match = lowercase(county_match)
end

# ╔═╡ c9f03768-f6bd-4b52-b18f-a4b6698b23ca
function GID2_to_fips(code)
	_, state, county_x = split(code, ".")
	county, _ = split(county_x, "_")

	fips = state * lpad(county, 3, '0')
	(; state = parse(Int, state), county = parse(Int, county), fips = parse(Int, fips))
end

# ╔═╡ ac0699ba-d578-45e9-996a-a7c88111aa4a
begin
	point_to_point = GeometryBasics.Point ∘ Meshes.coordinates 
	
	function chain_to_polygon(chain)
		chain |> Meshes.vertices .|>
		point_to_point |>
		#Meshes.coordinates .|> GeometryBasics.Point |>
		GeometryBasics.Polygon
	end

	function multi_to_multi(multi_polygon)
		multi_polygon |> Meshes.chains .|> chain_to_polygon |> GeometryBasics.MultiPolygon
	end
end

# ╔═╡ 58df19fa-9f1b-49db-9302-b8003ff8244e
@chain shapes_pop_df begin
	data(_) * (mapping(:geometry, color = :population2010 => log) * visual(Poly) + mapping(:center) * visual(Scatter, markersize = 1))
	draw
end

# ╔═╡ e889d127-4796-4a89-ab36-b39bae718b6d
md"""
## Zillow real estate prices
"""

# ╔═╡ 9b4cc5c3-ce13-4f90-b41e-d7a7398f3e65
@chain zillow_df begin
	@subset(!ismissing(:hpi), :month == 3, :year ∈ 2000:5:2020)
	disallowmissing
	innerjoin(shapes_pop_df, on = [:fips, :state, :county])
	data(_) * mapping(:geometry, color = :hpi => log, layout = :year => nonnumeric) * visual(Poly)
	draw
end

# ╔═╡ decc7aec-daa8-42e2-a62d-e09c2bf559cd
@chain zillow_df begin
	leftjoin(pop_df, on = [:state, :county, :fips])
	@subset(!ismissing(:hpi))
	@groupby(:date = Date(:year, :month))
	@aside total = sum(pop_df.population2010)
	@combine(:coverage = sum(:population2010) / total)
	data(_) * mapping(:date => "", :coverage) * visual(Lines)
	draw(axis=(title="What fraction of the US population do Zillow data cover?",))
end

# ╔═╡ 288cafba-1854-49e0-89ac-350282c3e249
@chain zillow_df begin
	leftjoin(pop_df, on = [:fips, :state, :county])
	@subset(!ismissing(:hpi))
	disallowmissing(:population2010)
	@groupby(:date = Date(:year, :month), :state)
	@combine(:hpi = mean(:hpi, weights(:population2010)), :pop = sum(:population2010))
	data(_) * mapping(:date => "", :hpi => log, group=:state => nonnumeric) * visual(Lines)
	draw(axis=(; title = "(log) House prices across states over time"))
end

# ╔═╡ e33a1edc-d414-4e87-91aa-ed3ef9fa1ca2
zillow_growth_df = @chain zillow_pop_df begin
	sort([:fips, :year])
	@groupby(:fips)
	@transform(
		:Δ_hpi = @c([missing; diff(:hpi)]),
		:Δ_log_hpi = @c([missing; diff(log.(:hpi))])
	)
	@groupby(:fips)
	@transform(
		:l_hpi = @c(lag(:hpi)),
		:l_Δ_hpi = @c(lag(:Δ_hpi)),
		:l_Δ_log_hpi = @c(lag(:Δ_log_hpi))
	)
end

# ╔═╡ f3df7592-bc70-4cc4-8ad0-778d48d75dc5
# ╠═╡ disabled = true
#=╠═╡
@chain zillow_growth_df begin
	@transform(:pop = log(:population2010 / 1000))
	@subset(!ismissing(:Δ_log_hpi), !ismissing(:l_Δ_log_hpi))
	data(_) * mapping(:Δ_log_hpi, :l_Δ_log_hpi) * (visual(Scatter, color=(:blue, 0.1)) * mapping(markersize = :pop) + AlgebraOfGraphics.density() * mapping(weights = :population2010) * visual(Contour))
	draw
end
  ╠═╡ =#

# ╔═╡ 9d975887-a482-4d30-9600-1d9ac7e3821e
zillow_pop_df = @chain zillow_df begin
	@groupby(:fips, :year, :state, :county)
	@combine(:hpi = mean(:hpi))
	leftjoin(pop_df, on = [:fips, :state, :county])
	@select(:year, :hpi, :fips, :population2010)
end

# ╔═╡ c1bbfe1e-27c0-4cd3-9679-e95dd7e8be57
zillow_df = @chain zillow_df0 begin
	rename(:StateCodeFIPS => :state, :MunicipalCodeFIPS => :county)
	stack(r"^20", [:state, :county], value_name = :hpi, variable_name = :date )
	@transform(:date = Date(:date))
	@select(:year = year(:date), :month = month(:date), Not(:date), :fips = parse(Int, string(:state)*lpad(:county, 3, '0')))
end

# ╔═╡ d6e02e8c-5783-4be7-aedc-7605f843d011
zillow_df0 = joinpath(datadep"zillow", "County_zhvi_uc_sfrcondo_tier_0.33_0.67_sm_sa_month.csv") |> CSV.File |> DataFrame

# ╔═╡ 19bc4c50-66ed-49ea-ac30-283137a1aedf
md"""
## The Social Connectedness Index
"""

# ╔═╡ f08a1afd-0ad6-4c85-98c1-28754356e73f
md"""
### SCI with population and distances
"""

# ╔═╡ dd918f80-ac1f-46e9-ba6a-416e4bb2e39f
sci_pop_distance_df = @chain county_df begin
	@aside centers_df = @select(shapes_pop_df, :fips, :center, :pop = :population2010)
	innerjoin(centers_df, on = :user_loc => :fips)
	rename!(:center => :user_center, :pop => :user_pop)
	innerjoin(centers_df, on = :fr_loc => :fips)
	rename!(:center => :fr_center, :pop => :fr_pop)
	@transform!(:distance = norm(:fr_center - :user_center))
	@select!(:user_loc, :fr_loc, :scaled_sci, 
		:user_pop, :fr_pop,
		:distance, :distance_mi = :distance * 55
	)
	dropmissing
end

# ╔═╡ 397ec8b4-5be1-4d13-a1ec-a9d1ddf8a00f
# ╠═╡ disabled = true
#=╠═╡
@chain sci_pop_distance_df begin
	@transform(
		:pop_pop = :fr_pop * :user_pop
	)
	data(_) * mapping(:distance_mi, weights = :pop_pop) * AlgebraOfGraphics.density()
	draw
end
  ╠═╡ =#

# ╔═╡ 6afa713a-291c-41bc-94fc-466d64d73eef
county_df = SCI_data(:US_counties)

# ╔═╡ bb4d5193-6b4a-4443-b3e4-a50a0e370eb2
md"""
## Price experiences of friends

How much did house prices grow where friends live?
"""

# ╔═╡ f8c180f5-6976-4ce4-8d7e-4a162c400288
@chain zillow_and_friends_df begin
	@transform(:pop = log(:population2010/1000))
	@subset(!ismissing(:Δ_log_hpi))
	disallowmissing!
	@subset(:year ∈ 2001:5:2021)
	data(_) * mapping(:Δ_log_hpi, :friends_exp_100, layout = :year=> nonnumeric, markersize = :pop) * visual(Scatter)
	draw
end

# ╔═╡ 6cbf8d8c-29e1-4a8d-9359-5a6145ea8bd1
zillow_and_friends_df = @chain zillow_growth_df begin
	@select(:year, :fips, :population2010, :hpi, :Δ_log_hpi)
	@subset(!ismissing(:Δ_log_hpi))
	innerjoin(friends_experience_df, on = [:year, :fips])
end

# ╔═╡ c485a1ce-54b5-4f9b-8ab0-f6739af6344a
friends_experience_df = @chain sci_pop_distance_df begin
	@transform!(:sci_pop = :scaled_sci * :fr_pop)
	sort!(:user_loc)
	#@subset(:user_loc < 20000)
	@groupby(:fips = :user_loc)
	combine([:fips, :fr_loc, :sci_pop, :distance_mi] => friends_exp(zillow_growth_df, :Δ_log_hpi) => AsTable)
end

# ╔═╡ 3c0bbf62-3bc2-429c-a8d6-47854cc55da1
function friends_exp(x_df0, x)
	x_df = @chain x_df0 begin
		@select(:fips, :year, :x = $(x))
		@subset(!ismissing(:x))
	end
		
	function (user_loc, fr_loc, sci_pop, distance_mi)
		user_loc = only(unique(user_loc))
	
		@chain begin
			DataFrame(; fr_loc, sci_pop, distance_mi)
			#@subset!(:distance_mi > 500)
			sort!(:fr_loc)
			rightjoin(_, x_df, on = :fr_loc => :fips, makeunique=true)
			#@subset!(:fr_loc != user_loc)
			@subset!(!ismissing(:sci_pop))
			disallowmissing!([:sci_pop, :x])
			@groupby(:year)
			@combine(
				:friends_exp_500 = mean(:x, weights(:sci_pop .* (:distance_mi .> 500))),
				:friends_exp_300 = mean(:x, weights(:sci_pop .* (:distance_mi .> 300))),
				:friends_exp_100 = mean(:x, weights(:sci_pop .* (:distance_mi .> 100))),
				:friends_exp_000 = mean(:x, weights(:sci_pop .* (:distance_mi .> 0))),
				:friends_exp = mean(:x, weights(:sci_pop))
			)
		end
	end
end

# ╔═╡ f2329bd8-b585-41a5-b286-8e3a8233a4cd
md"""
## The HMDA mortgage database
"""

# ╔═╡ 550295bd-ae15-404c-810a-f2f2561fcf5f
hmda_df0 = RData.load(datadep"hmda-panel/hmda_big.RData")["hmda_big"]

# ╔═╡ cdc99479-37cd-4ada-9ce0-08e136281a42
hmda_panel_url = "https://gitlab.com/drechsel-grau-greimel/hmda/-/raw/master/data-processed/hmda_big.RData"

# ╔═╡ fefa0c22-b554-49d5-a024-2f5307f13b31
hmda_df = @chain hmda_df0 begin
	@subset(!ismissing(:owner_occupied), !ismissing(:purpose_loan), !ismissing(:county))
	@subset!(:owner_occupied, :purpose_loan == "purchase")
	@transform!(:state = Int(:state), :county = tryparse(Int, replace(get(:county), "O" => "0")))
	@subset(!isnothing(:county))
	@transform!(:fips = parse(Int, string(:state) * lpad(string(:county), 3, '0')))
end

# ╔═╡ a6e7fdcb-0ad4-4b04-a032-ff14d3f63926
@chain hmda_df begin
	@select(:amount, :count, :year, :fips)
	innerjoin(shapes_pop_df, on = :fips)
	@transform!(:count_per_pop = :count / :population2010)
	@subset!(:year ∈ 2000:5:2015)
	data(_) * mapping(:geometry, color = :count_per_pop, layout = :year => nonnumeric) * visual(Poly)
	draw
end

# ╔═╡ d5dc551e-c120-4a20-8c94-15e2e8300b7a
@chain hmda_df begin
	@select(:amount = log(:amount / 1000), :count, :year, :fips)
	innerjoin(shapes_pop_df, on = :fips)
	@transform!(:count_per_pop = :count / :population2010)
	@subset!(:year ∈ 2000:5:2015)
	data(_) * mapping(:geometry, color = :amount, layout = :year => nonnumeric) * visual(Poly)
	draw
end

# ╔═╡ 844c432e-c804-4a47-adad-bef2887f7dc0
md"""
# Appendix
"""

# ╔═╡ 6871566b-fb9a-4d7a-b6f0-6b8d38fc7127
using MarkdownLiteral: @mdx

# ╔═╡ e85d31e6-0068-44df-a88d-f245d90e21d1
using Dates

# ╔═╡ c1c01e0e-d155-11ec-14a0-1173f9a37f8f
using RData: RData

# ╔═╡ 67697ae7-049e-42e2-b476-49c4fd807823
using LinearAlgebra: norm

# ╔═╡ 54c7ee04-19f9-4083-8869-f6a9d07ba51f
using DataFrames, DataFrameMacros, Chain

# ╔═╡ e3a83114-db7d-4488-9db7-14ed919b0f33
using CairoMakie, AlgebraOfGraphics

# ╔═╡ cb6557ad-1cc6-417b-8668-38e2399fd6af
using StatsBase: weights

# ╔═╡ a5455cd9-fe6c-4bbe-a0c6-4537ddf585b5
using Statistics: mean

# ╔═╡ 93c33060-d39f-4ece-ab2c-c929cc249a5b
using PlutoUI

# ╔═╡ 4953ba71-27ae-4197-9d2d-997c7f513159
TableOfContents()

# ╔═╡ af98fe40-d533-46fc-a3da-2cee27a40b9e
using ShiftedArrays

# ╔═╡ 1f40fbb1-feb7-41ea-8cfe-8d00c0c30a54
using HTTP: HTTP

# ╔═╡ e6848c06-8895-444b-864f-6e4c86e488c1
md"""
## DataDeps
"""

# ╔═╡ d43c8098-2c1e-48bc-8d82-64b110451ac4
dist_urls = [
	"500mi" => "https://nber.org/distance/2010/sf1/county/sf12010countydistance500miles.csv.zip",
	"100mi" => "https://data.nber.org/distance/2010/sf1/county/sf12010countydistance100miles.csv.zip"]

# ╔═╡ 48cc091a-afab-446a-a0dc-22b4f5572162
function county_dist_data(id)
	urls = dist_urls
	id = string(id)
	files = urls_to_files(urls)
	valid_ids = first.(urls)
	if id ∉ valid_ids
		ArgumentError("provide one of $valid_ids") |> throw
	end
	path = joinpath(@datadep_str("county_dist_$id"), files[id])
	
	CSV.File(path) |> DataFrame
end

# ╔═╡ b0449f77-f2e5-47da-ae96-976b467388ba
sci_url_pre = "https://data.humdata.org/dataset/e9988552-74e4-4ff4-943f-c782ac8bca87/resource/"

# ╔═╡ 4b7380ac-d9d4-44c4-ac82-55d4859d117b
sci_urls = [
	:countries =>
		"35ca6ade-a5bd-4782-b266-797169dca74b/download/countries-countries-fb-social-connectedness-index-october-2021.tsv",
	:US_counties =>
		"c59fd5ac-0458-4e83-b6be-5334f0ea9a69/download/us-counties-us-counties-fb-social-connectedness-index-october-2021.zip",
	:US_counties__countries => 
		"868a2fdb-f5c8-4a98-af7c-cfc8bf0daeb3/download/us-counties-countries-fb-social-connectedness-index-october-2021.tsv",
	:GADM_NUTS2 =>
		"cc5b6046-c417-4e25-930a-3d31538dffc5/download/gadm1_nuts2-gadm1_nuts2-fb-social-connectedness-index-october-2021.zip",
	:GADM_NUTS3_counties =>
		"18bf46fe-7f84-47b7-9d7e-b79a9c491f52/download/gadm1_nuts3_counties-gadm1_nuts3_counties-fb-social-connectedness-index-october-2021.zip"
	]

# ╔═╡ 2bd5e792-3e38-423d-a684-f716a646294c
sci_checksums = Dict(
	:US_counties =>
		"023cf8a522a4e15ca321113adf9dcda85b7796fee3a5688f825ffc71a0eeaa1f",
	:countries =>
		"bc6269eb10da945e2327beb725e3ef4fa955fd430535b0357cc118a0b3c7cfd6",
	:US_counties__countries => nothing,
	:GADM_NUTS2 => nothing,
	:GADM_NUTS3_counties => nothing
	)

# ╔═╡ 3160a516-005c-4126-a750-bbf98c54d722
function urls_to_files(urls)
	map(urls) do (id, url)
		file = split(url, "/") |> last |> string
		id => replace(file, ".zip" => "")
	end |> Dict
end

# ╔═╡ 92bea217-2d43-40e4-ad04-bc35c55ab52d
#sci_files = urls_to_files(sci_urls)

# ╔═╡ f3fbc391-0dea-4637-8911-126e9a52b913
function SCI_data(id)
	urls = sci_urls
	id = Symbol(id)
	files = urls_to_files(urls)
	valid_ids = first.(urls)
	if id ∉ valid_ids
		ArgumentError("provide one of $valid_ids") |> throw
	end
	if id == :US_counties
		path = joinpath(@datadep_str("SCI_$id"), "county_county.tsv")
	else
		path = joinpath(@datadep_str("SCI_$id"), files[id])
	end
	
	CSV.File(path) |> DataFrame
end

# ╔═╡ cc203ba7-ab85-452f-b9c2-34dc0c0d937e
import CSV

# ╔═╡ fe876989-2cb8-4846-a27f-67ed07e93335
begin
	using DataDeps: DataDeps, DataDep, @datadep_str, register, unpack
	ENV["DATADEPS_ALWAYS_ACCEPT"] = "true"
	
	for (id, url) in sci_urls
		register(DataDep(
    		"SCI_$id",
    		"""
	
			""",
    		sci_url_pre * url,
	    	sci_checksums[id],
		 	post_fetch_method = id ∉ [:US_counties__countries, :countries] ? unpack : identity 
		))
	end
	
	for (id, url) in dist_urls
		register(DataDep(
    		"county_dist_$id",
    		"""
	
			""",
    		url,
	    	#sci_checksums[id],
		 	post_fetch_method = unpack
		))
	end

	register(DataDep(
		"county_pop",
		" ",
		"https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv"
	))

	register(DataDep(
		"zillow",
		" ",
		[prices_zip_url, prices_county_url, rentals_zip_url]#,
		#["2a0792d21d0134f95f029ecf8b01a272bfdac79bb231ab20c8176502dc917048"]
	))

	register(DataDep(
		"hmda-panel",
		" ",
		hmda_panel_url,
		"95e6e5d87915cf324912eb9a40b0ecfeeec3b8346ceb77d4b3cae16394db762d"
	))
end

# ╔═╡ 89e50510-636b-47ae-8f4a-7ef18317c284
md"""
### Population data
"""

# ╔═╡ fa730aed-a6f6-4a4a-94b2-76633eb8551b
function get_pop(args...; kwargs...)
	file = joinpath(datadep"county_pop", "co-est2019-alldata.csv")
	df = CSV.File(file, args...; kwargs...) |> DataFrame
end

# ╔═╡ 8fb93eef-88c0-43b9-b193-4a8911ec8f15
pop_url = "https://www2.census.gov/programs-surveys/popest/datasets/2010-2019/counties/totals/co-est2019-alldata.csv"

# ╔═╡ 0176205f-b123-4b1d-b38a-d871b1c6bbe4
pop_df = @chain begin
	get_pop(select = [:REGION, :DIVISION, :STATE, :COUNTY, :STNAME, :CTYNAME, :CENSUS2010POP])
	@transform(:fips = string(:STATE) * lpad(string(:COUNTY), 3, '0'))
	@select(
		:fips = Meta.parse(:fips),
		:state_name = :STNAME, :county_name = :CTYNAME,
		:state = :STATE, :county = :COUNTY,
		:population2010 = :CENSUS2010POP,
	#	#:divisor = :STNAME * " " * :CTYNAME

	)
	#@transform(:county_match = clean_county_name_for_matching(:county))
	@subset(:county != 0)
end

# ╔═╡ afb0db63-970a-469e-a7f4-f109a289d06f
md"""
## House price data
"""

# ╔═╡ bd17661f-6d02-47e3-b2df-8f6bb3d99876
prices_county_url = "https://files.zillowstatic.com/research/public_csvs/zhvi/County_zhvi_uc_sfrcondo_tier_0.33_0.67_sm_sa_month.csv"

# ╔═╡ 914522a3-ace0-49f1-8b87-50fa3af2dbb3
rentals_zip_url = "https://files.zillowstatic.com/research/public_csvs/zori/Zip_ZORI_AllHomesPlusMultifamily_Smoothed.csv"

# ╔═╡ 872cb78f-27d3-40f5-b9e9-e116132a473c
prices_zip_url = "https://files.zillowstatic.com/research/public_csvs/zhvi/Zip_zhvi_uc_sfrcondo_tier_0.33_0.67_sm_sa_month.csv"

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CSV = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataDeps = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
Dates = "ade2ca70-3891-5945-98fb-dc099432e06a"
FixedEffectModels = "9d5cd8c9-2029-5cab-9928-427838db53e3"
GeoTables = "e502b557-6362-48c1-8219-d30d308dcdb0"
HTTP = "cd3eb016-35fb-5094-929b-558a96fad6f3"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
RData = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
ShiftedArrays = "1277b4bf-5013-50f5-be3d-901d8477a67a"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"

[compat]
AlgebraOfGraphics = "~0.6.7"
CSV = "~0.10.4"
CairoMakie = "~0.8.5"
Chain = "~0.4.10"
DataDeps = "~0.7.8"
DataFrameMacros = "~0.2.1"
DataFrames = "~1.3.4"
FixedEffectModels = "~1.6.6"
GeoTables = "~0.4.0"
HTTP = "~0.9.17"
MarkdownLiteral = "~0.1.1"
PlutoUI = "~0.7.39"
RData = "~0.8.3"
ShiftedArrays = "~1.0.0"
StatsBase = "~0.33.16"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.2"
manifest_format = "2.0"
project_hash = "bbb4519dbc3ae745081d2f5b4840dabe2f75673a"

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
git-tree-sha1 = "03e0550477d86222521d254b741d470ba17ea0b5"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.3.4"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "af92965fb30777147966f58acb05da51c5616b5f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.3.3"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "4c245508a6c1f05c42b507c8fc85f3a17bf88206"
repo-rev = "master"
repo-url = "https://github.com/JuliaPlots/AlgebraOfGraphics.jl.git"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.7"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArchGDAL]]
deps = ["CEnum", "ColorTypes", "Dates", "DiskArrays", "Extents", "GDAL", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "ImageCore", "Tables"]
git-tree-sha1 = "f5592638035b4709a42a9b27e31c13309be80eb9"
uuid = "c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3"
version = "0.9.0"

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
git-tree-sha1 = "e65431fdcee9883cfac2ea2fd388fe6e4372fa56"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.8.5"

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
git-tree-sha1 = "9489214b993cd42d17f44c36e359bf6a7c919abf"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.0"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "1e315e3f4b0b7ce40feded39c73049692126cf53"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.3"

[[deps.CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "3587fdbecba8c44f7e7285a1957182711b95f580"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.3.1"

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

[[deps.Combinatorics]]
git-tree-sha1 = "08c8b6831dc00bfea825826be0bc8336fc369860"
uuid = "861a8166-3701-5b0c-9a16-15d98fcdc6aa"
version = "1.0.2"

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "URIs"]
git-tree-sha1 = "4cd7063c9bdebdbd55ede1af70f3c2f48fab4215"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.6"

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

[[deps.DBFTables]]
deps = ["Printf", "Tables", "WeakRefStrings"]
git-tree-sha1 = "f5b78d021b90307fb7170c4b013f350e6abe8fed"
uuid = "75c7ada1-017a-5fb6-b8c7-2125ff2d6c93"
version = "1.0.0"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataDeps]]
deps = ["BinaryProvider", "HTTP", "Libdl", "Reexport", "SHA", "p7zip_jll"]
git-tree-sha1 = "e299d8267135ef2f9c941a764006697082c1e7e8"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.8"

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

[[deps.DiskArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "230d999fc78652ea070312373ed1bfe2489e4fe5"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.3.6"

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

[[deps.ExprTools]]
git-tree-sha1 = "56559bbef6ca5ea0c0818fa5c90320398a6fbf8d"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.8"

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

[[deps.FixedEffectModels]]
deps = ["DataFrames", "FixedEffects", "LinearAlgebra", "Printf", "Reexport", "Statistics", "StatsBase", "StatsFuns", "StatsModels", "Tables", "Vcov"]
git-tree-sha1 = "f6bd6f55724d239abcb833678183da20ff4f2c54"
uuid = "9d5cd8c9-2029-5cab-9928-427838db53e3"
version = "1.6.6"

[[deps.FixedEffects]]
deps = ["GroupedArrays", "LinearAlgebra", "Printf", "Requires", "StatsBase"]
git-tree-sha1 = "06c114eaad4566df6287c5d303c194309f923efb"
uuid = "c8885935-8500-56a7-9867-7708b20db0eb"
version = "2.1.1"

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

[[deps.GADM]]
deps = ["ArchGDAL", "DataDeps", "GeoInterface", "Logging", "Tables"]
git-tree-sha1 = "5c7b6420512082cbfad419641ea68c151bda0526"
uuid = "a8dd9ffe-31dc-4cf5-a379-ea69100a8233"
version = "0.4.1"

[[deps.GDAL]]
deps = ["CEnum", "GDAL_jll", "NetworkOptions", "PROJ_jll"]
git-tree-sha1 = "9ce70502472a9f23f8889f0f9e2be8451413fe7b"
uuid = "add2ef01-049f-52c4-9ee2-e494f65e021a"
version = "1.4.0"

[[deps.GDAL_jll]]
deps = ["Artifacts", "Expat_jll", "GEOS_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "OpenJpeg_jll", "PROJ_jll", "Pkg", "SQLite_jll", "Zlib_jll", "Zstd_jll", "libgeotiff_jll"]
git-tree-sha1 = "756a15a73ded80cf194e7458abeb6f559d5070e2"
uuid = "a7073274-a066-55f0-b90d-d619367d196c"
version = "300.500.0+1"

[[deps.GEOS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4ceb4cdae127931b852ced4d3782bb51ab5e2632"
uuid = "d604d12d-fa86-5845-992e-78dc15976526"
version = "3.10.2+0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "039118892476c2bf045a43b88fcb75ed566000ff"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.0"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "434166198434a5c2fcc0a1a59d22c3b0ad460889"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.1"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "de6980e052d67c0da1872dfdb2c49fb7d3f56b07"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.0"

[[deps.GeoInterfaceRecipes]]
deps = ["GeoInterface", "RecipesBase"]
git-tree-sha1 = "29e1ec25cfb6762f503a19495aec347acf867a9e"
uuid = "0329782f-3d07-4b52-b9f6-d3137cf03c7a"
version = "1.0.0"

[[deps.GeoTables]]
deps = ["ArchGDAL", "GADM", "GeoInterface", "Meshes", "Shapefile", "Tables"]
git-tree-sha1 = "47b68dd4c3d1214b8c489846547512817dcde9a7"
uuid = "e502b557-6362-48c1-8219-d30d308dcdb0"
version = "0.4.0"

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

[[deps.GroupedArrays]]
deps = ["DataAPI", "Missings"]
git-tree-sha1 = "44c812379b629eea08b6d25a196010f1f4b001e3"
uuid = "6407cd72-fade-4a84-8a1e-56e431fc1533"
version = "0.3.3"

[[deps.HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

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
git-tree-sha1 = "9a5c62f231e5bba35695a20988fc7cd6de7eeb5a"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.3"

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

[[deps.Inflate]]
git-tree-sha1 = "f5fc07d4e706b84f72d54eedcc1c13d92fb0871c"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.2"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "61feba885fac3a407465726d0c330b3055df897f"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.1.2"

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

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

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

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LittleCMS_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pkg"]
git-tree-sha1 = "110897e7db2d6836be22c18bffd9422218ee6284"
uuid = "d3a379c0-f9a3-5b72-a4c0-6bf4d2e8af0f"
version = "2.12.0+0"

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
git-tree-sha1 = "85b7bbc5ba64ee9c68fe32d908aa29f0a97a6145"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.17.5"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "76b079514ddae2cb3bf839c499678e1ef6c5df78"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.3.2"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MarkdownLiteral]]
deps = ["CommonMark", "HypertextLiteral"]
git-tree-sha1 = "0d3fa2dd374934b62ee16a4721fe68c418b92899"
uuid = "736d6165-7244-6769-4267-6b50796e6954"
version = "0.1.1"

[[deps.Match]]
git-tree-sha1 = "1d9bc5c1a6e7ee24effb93f175c9342f9154d97f"
uuid = "7eb4fadd-790c-5f42-8a69-bfa0b872bfbf"
version = "1.2.0"

[[deps.MathTeXEngine]]
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test"]
git-tree-sha1 = "5c1e3d66b3a36029de4e5ac07ab8bafd5a8041e5"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.4.1"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"

[[deps.Meshes]]
deps = ["CategoricalArrays", "CircularArrays", "Distances", "IterTools", "IteratorInterfaceExtensions", "LinearAlgebra", "NearestNeighbors", "Random", "RecipesBase", "ReferenceFrameRotations", "SimpleTraits", "SparseArrays", "SpecialFunctions", "StaticArrays", "StatsBase", "TableTraits", "Tables"]
git-tree-sha1 = "9beca5fb86e768ebab0ee9b7f0da68c2b39cfdc8"
uuid = "eacbb407-ea5a-433e-ab97-5258b1ca43fa"
version = "0.22.8"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "29714d0a7a8083bba8427a4fbfb00a540c681ce7"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.3"

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

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "ded92de95031d4a8c61dfb6ba9adb6f1d8016ddd"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.10"

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

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"

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
git-tree-sha1 = "3e32c8dbbbe1159a5057c80b8a463369a78dd8d8"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.12"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "Pkg", "SQLite_jll"]
git-tree-sha1 = "12bd68665a0c3cb4635c4359d3fa9e2769ed59e5"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "900.0.0+0"

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

[[deps.RData]]
deps = ["CategoricalArrays", "CodecZlib", "DataFrames", "Dates", "FileIO", "Requires", "TimeZones", "Unicode"]
git-tree-sha1 = "19e47a495dfb7240eb44dc6971d660f7e4244a72"
uuid = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
version = "0.8.3"

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

[[deps.ReferenceFrameRotations]]
deps = ["Crayons", "LinearAlgebra", "Printf", "Random", "StaticArrays"]
git-tree-sha1 = "ec9bde2e30bc221e05e20fcec9a36a9c315e04a6"
uuid = "74f56ac7-18b3-5285-802d-d4bd4f104033"
version = "3.0.0"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "307761d71804208c0c62abdbd0ea6822aa5bbefd"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.2.0"

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

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "b6d006c4c57278d532de38912e16adf626c949c7"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.38.4+0"

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

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "db8481cf5d6278a121184809e9eb1628943c7704"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.13"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Shapefile]]
deps = ["DBFTables", "Extents", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "RecipesBase", "Tables"]
git-tree-sha1 = "2f400236c85ba357dfdc2a56af80c939dc118f02"
uuid = "8e980c4a-a4fe-5da2-b3a7-4b4b0353a2f4"
version = "0.8.0"

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
git-tree-sha1 = "a9e798cae4867e3a41cae2dd9eb60c047f1212db"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.6"

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
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "f90022b44b7bf97952756a6b6737d1a0024a3233"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.5.5"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Serialization", "Unicode"]
git-tree-sha1 = "0a359b0ee27e4fbc90d9b3da1f48ddc6f98a0c9e"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.7.3"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "216b95ea110b5972db65aa90f88d8d89dcb8851c"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.6"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

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

[[deps.Vcov]]
deps = ["Combinatorics", "GroupedArrays", "LinearAlgebra", "StatsBase", "Tables"]
git-tree-sha1 = "1233e304ac41897767c137f6ee281391d6ebfb0e"
uuid = "ec2bfdc2-55df-4fc9-b9ae-4958c2cf2486"
version = "0.5.2"

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

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

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

[[deps.libgeotiff_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "PROJ_jll", "Pkg"]
git-tree-sha1 = "e51bca193c8a4774dc1d2e5d40d5c4491c1b4fd4"
uuid = "06c338fa-64ff-565b-ac2f-249532af990e"
version = "1.7.1+0"

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
# ╟─f5450eab-0f9f-4b7f-9b80-992d3c553ba9
# ╟─8d268ee9-6a3b-4bf3-a0d1-59d86fe8f263
# ╟─662cad72-af63-41c9-a895-5633be486f3f
# ╟─cdfa2f05-45cc-4050-baf4-15d66910c4e9
# ╠═2bc4cd24-a12e-46c2-97d4-03d0bdaca3a3
# ╟─c3b70e50-420b-4956-9c36-d4596fe07634
# ╠═6959807f-03c5-4a9d-8c4d-c04c06c6e0a5
# ╟─63d1da26-9ac4-411f-8f47-02a02b37f142
# ╟─dbba5ead-cd2a-4714-9e56-dd34c7d72039
# ╠═8001cd65-2542-4a1a-93b8-2e7243226956
# ╠═2357bd9e-93e2-468e-8bd4-7720f07eed5e
# ╠═2eb3a5db-f5a3-489f-b712-146934229772
# ╠═7297dc74-44b4-4d2d-98f5-7120b74cb196
# ╠═ca4b0557-9360-4089-ab64-a62e07d23de9
# ╟─24013463-f8ea-463f-81e9-06542880dc8b
# ╟─cc13875f-0407-4749-a30e-39737c28c8b3
# ╠═aa8593ef-3f6f-4b0c-8046-f88aef680918
# ╠═1af8b2a8-f7fb-4afa-a746-cedbe1fea1a6
# ╠═301662e7-2d65-4b3b-9128-c415918e0058
# ╠═94686899-25fd-4d28-9894-5dcdc2650efe
# ╠═e8183f41-5921-4bfc-a262-6ac6b165870c
# ╠═85768c95-3f0c-43ad-9e37-c479c9bb6a5b
# ╠═edf2ff0d-ab91-47ff-ae03-cee7507601a1
# ╠═c9f03768-f6bd-4b52-b18f-a4b6698b23ca
# ╠═ac0699ba-d578-45e9-996a-a7c88111aa4a
# ╠═58df19fa-9f1b-49db-9302-b8003ff8244e
# ╟─e889d127-4796-4a89-ab36-b39bae718b6d
# ╠═9b4cc5c3-ce13-4f90-b41e-d7a7398f3e65
# ╟─decc7aec-daa8-42e2-a62d-e09c2bf559cd
# ╟─288cafba-1854-49e0-89ac-350282c3e249
# ╠═e33a1edc-d414-4e87-91aa-ed3ef9fa1ca2
# ╠═f3df7592-bc70-4cc4-8ad0-778d48d75dc5
# ╠═9d975887-a482-4d30-9600-1d9ac7e3821e
# ╠═c1bbfe1e-27c0-4cd3-9679-e95dd7e8be57
# ╠═d6e02e8c-5783-4be7-aedc-7605f843d011
# ╟─19bc4c50-66ed-49ea-ac30-283137a1aedf
# ╟─f08a1afd-0ad6-4c85-98c1-28754356e73f
# ╠═dd918f80-ac1f-46e9-ba6a-416e4bb2e39f
# ╠═397ec8b4-5be1-4d13-a1ec-a9d1ddf8a00f
# ╠═6afa713a-291c-41bc-94fc-466d64d73eef
# ╟─bb4d5193-6b4a-4443-b3e4-a50a0e370eb2
# ╠═f8c180f5-6976-4ce4-8d7e-4a162c400288
# ╠═6cbf8d8c-29e1-4a8d-9359-5a6145ea8bd1
# ╠═c485a1ce-54b5-4f9b-8ab0-f6739af6344a
# ╠═3c0bbf62-3bc2-429c-a8d6-47854cc55da1
# ╟─f2329bd8-b585-41a5-b286-8e3a8233a4cd
# ╠═550295bd-ae15-404c-810a-f2f2561fcf5f
# ╠═cdc99479-37cd-4ada-9ce0-08e136281a42
# ╠═fefa0c22-b554-49d5-a024-2f5307f13b31
# ╠═a6e7fdcb-0ad4-4b04-a032-ff14d3f63926
# ╠═d5dc551e-c120-4a20-8c94-15e2e8300b7a
# ╟─844c432e-c804-4a47-adad-bef2887f7dc0
# ╠═6871566b-fb9a-4d7a-b6f0-6b8d38fc7127
# ╠═e85d31e6-0068-44df-a88d-f245d90e21d1
# ╠═c1c01e0e-d155-11ec-14a0-1173f9a37f8f
# ╠═67697ae7-049e-42e2-b476-49c4fd807823
# ╠═54c7ee04-19f9-4083-8869-f6a9d07ba51f
# ╠═e3a83114-db7d-4488-9db7-14ed919b0f33
# ╠═cb6557ad-1cc6-417b-8668-38e2399fd6af
# ╠═a5455cd9-fe6c-4bbe-a0c6-4537ddf585b5
# ╠═93c33060-d39f-4ece-ab2c-c929cc249a5b
# ╠═4953ba71-27ae-4197-9d2d-997c7f513159
# ╠═af98fe40-d533-46fc-a3da-2cee27a40b9e
# ╠═1f40fbb1-feb7-41ea-8cfe-8d00c0c30a54
# ╟─e6848c06-8895-444b-864f-6e4c86e488c1
# ╟─d43c8098-2c1e-48bc-8d82-64b110451ac4
# ╠═48cc091a-afab-446a-a0dc-22b4f5572162
# ╠═b0449f77-f2e5-47da-ae96-976b467388ba
# ╠═4b7380ac-d9d4-44c4-ac82-55d4859d117b
# ╠═2bd5e792-3e38-423d-a684-f716a646294c
# ╠═3160a516-005c-4126-a750-bbf98c54d722
# ╠═92bea217-2d43-40e4-ad04-bc35c55ab52d
# ╠═f3fbc391-0dea-4637-8911-126e9a52b913
# ╠═cc203ba7-ab85-452f-b9c2-34dc0c0d937e
# ╠═fe876989-2cb8-4846-a27f-67ed07e93335
# ╟─89e50510-636b-47ae-8f4a-7ef18317c284
# ╠═fa730aed-a6f6-4a4a-94b2-76633eb8551b
# ╠═8fb93eef-88c0-43b9-b193-4a8911ec8f15
# ╠═0176205f-b123-4b1d-b38a-d871b1c6bbe4
# ╟─afb0db63-970a-469e-a7f4-f109a289d06f
# ╠═bd17661f-6d02-47e3-b2df-8f6bb3d99876
# ╠═914522a3-ace0-49f1-8b87-50fa3af2dbb3
# ╠═872cb78f-27d3-40f5-b9e9-e116132a473c
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

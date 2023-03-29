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
AlgebraOfGraphics = "~0.6.14"
CSV = "~0.10.9"
CairoMakie = "~0.10.3"
Chain = "~0.5.0"
DataDeps = "~0.7.10"
DataFrameMacros = "~0.4.1"
DataFrames = "~1.5.0"
FixedEffectModels = "~1.9.0"
GeoTables = "~1.2.2"
HTTP = "~1.7.4"
MarkdownLiteral = "~0.1.1"
PlutoUI = "~0.7.50"
RData = "~1.0.0"
ShiftedArrays = "~2.0.0"
StatsBase = "~0.33.21"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "273841b421da5425338bf1c23e0b9f9bfbc969a0"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "16b6dbc4cf7caee4e1e75c49485ec67b667098a0"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.3.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.AbstractTrees]]
git-tree-sha1 = "faa260e4cb5aba097a73fab382dd4b5819d8ec8c"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.4"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "cc37d689f599e8df4f464b2fa3870ff7db7492ef"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.6.1"

[[deps.AlgebraOfGraphics]]
deps = ["Colors", "Dates", "Dictionaries", "FileIO", "GLM", "GeoInterface", "GeometryBasics", "GridLayoutBase", "KernelDensity", "Loess", "Makie", "PlotUtils", "PooledArrays", "RelocatableFolders", "SnoopPrecompile", "StatsBase", "StructArrays", "Tables"]
git-tree-sha1 = "43c2ef89ca0cdaf77373401a989abae4410c7b8a"
repo-rev = "master"
repo-url = "https://github.com/JuliaPlots/AlgebraOfGraphics.jl.git"
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.14"

[[deps.Animations]]
deps = ["Colors"]
git-tree-sha1 = "e81c509d2c8e49592413bfb0bb3b08150056c79d"
uuid = "27a7e980-b3e6-11e9-2bcd-0b925532e340"
version = "0.4.1"

[[deps.ArchGDAL]]
deps = ["CEnum", "ColorTypes", "Dates", "DiskArrays", "Extents", "GDAL", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "ImageCore", "Tables"]
git-tree-sha1 = "70908bb727c9a0ba863c5145aa48ee838cc29b84"
uuid = "c9ce4bd3-c3d5-55b8-8973-c0e20141b8c3"
version = "0.9.3"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arrow_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Lz4_jll", "Pkg", "Thrift_jll", "Zlib_jll", "boost_jll", "snappy_jll"]
git-tree-sha1 = "d64cb60c0e6a138fbe5ea65bcbeea47813a9a700"
uuid = "8ce61222-c28f-5041-a97a-c2198fb817bf"
version = "10.0.0+1"

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

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bessels]]
git-tree-sha1 = "4435559dc39793d53a9e3d278e185e920b4619ef"
uuid = "0e736298-9ec6-45e8-9647-e4fc86a2fe38"
version = "0.2.8"

[[deps.BitFlags]]
git-tree-sha1 = "43b1a4a8f797c1cddadf60499a8a077d4af2cd2d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.7"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CRC32c]]
uuid = "8bf52ea8-c179-5cab-976a-9e18b702a9bc"

[[deps.CSV]]
deps = ["CodecZlib", "Dates", "FilePathsBase", "InlineStrings", "Mmap", "Parsers", "PooledArrays", "SentinelArrays", "SnoopPrecompile", "Tables", "Unicode", "WeakRefStrings", "WorkerUtilities"]
git-tree-sha1 = "c700cce799b51c9045473de751e9319bdd1c6e94"
uuid = "336ed68f-0bac-5ca0-87d4-7b16caf5d00b"
version = "0.10.9"

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "7a6a830076a6eb2a8289e751e7237c04a1ce0ddd"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.3"

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

[[deps.CategoricalArrays]]
deps = ["DataAPI", "Future", "Missings", "Printf", "Requires", "Statistics", "Unicode"]
git-tree-sha1 = "5084cc1a28976dd1642c9f337b28a3cb03e0f7d2"
uuid = "324d7699-5711-5eae-9e2f-1d82baa6b597"
version = "0.10.7"

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c6d890a52d2c4d55d326439580c3b8d0875a77d9"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.7"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "485193efd2176b88e6622a39a246f8c5b600e74e"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.6"

[[deps.CircularArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "3587fdbecba8c44f7e7285a1957182711b95f580"
uuid = "7a955b69-7140-5f4e-a0ed-f168c5e2e749"
version = "1.3.1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "9c209fb7536406834aa938fb149964b985de6c83"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.1"

[[deps.ColorBrewer]]
deps = ["Colors", "JSON", "Test"]
git-tree-sha1 = "61c5334f33d91e570e1d0c3eb5465835242582c4"
uuid = "a2cac450-b92f-5266-8821-25eda20663c8"
version = "0.4.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random", "SnoopPrecompile"]
git-tree-sha1 = "aa3edc8f8dea6cbfa176ee12f7c2fc82f0608ed3"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.20.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "600cc5508d66b78aae350f7accdb58763ac18589"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.10"

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
deps = ["Crayons", "JSON", "SnoopPrecompile", "URIs"]
git-tree-sha1 = "e2f4627b0d3f2c1876360e0b242a7c23923b469d"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.10"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "7a60c856b9fa189eb34f5f8a6f6b5529b7942957"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "89a9db8d28102b094992472d333674bd1a83ce2a"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DBFTables]]
deps = ["Dates", "Tables", "WeakRefStrings"]
git-tree-sha1 = "6c6cb6614e5ff0769662a144d0a62d443b80be43"
uuid = "75c7ada1-017a-5fb6-b8c7-2125ff2d6c93"
version = "1.2.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataDeps]]
deps = ["HTTP", "Libdl", "Reexport", "SHA", "p7zip_jll"]
git-tree-sha1 = "bc0a264d3e7b3eeb0b6fc9f6481f970697f29805"
uuid = "124859b0-ceae-595e-8997-d05f6a7a8dfe"
version = "0.7.10"

[[deps.DataFrameMacros]]
deps = ["DataFrames", "MacroTools"]
git-tree-sha1 = "5275530d05af21f7778e3ef8f167fb493999eea1"
uuid = "75880514-38bc-4a95-a458-c2aea5a3a702"
version = "0.4.1"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InlineStrings", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SentinelArrays", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "aa51303df86f8626a962fccb878430cdb0a97eee"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.5.0"

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

[[deps.DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "e82c3c97b5b4ec111f3c1b55228cebc7510525a2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.25"

[[deps.DiskArrays]]
deps = ["OffsetArrays"]
git-tree-sha1 = "3f87990e0882e44c0f4e5c9699d09a0edbfa25c8"
uuid = "3c3547ce-8d99-4f5e-a174-61eb10b00ae3"
version = "0.3.9"

[[deps.Distances]]
deps = ["LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "49eba9ad9f7ead780bfb7ee319f962c811c6d3b2"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.8"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "da9e1a9058f8d3eec3a8c9fe4faacfb89180066b"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.86"

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

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.ExprTools]]
git-tree-sha1 = "c1d06d129da9f55715c6c212866f5b1bddc5fa00"
uuid = "e2ba6199-217a-4e67-a87a-7c52f15ade04"
version = "0.1.9"

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
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "74faea50c1d007c85837327f6775bea60b5492dd"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+2"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "f9818144ce7c8c41edf5c4c179c684d92aa4d9fe"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.6.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FilePathsBase]]
deps = ["Compat", "Dates", "Mmap", "Printf", "Test", "UUIDs"]
git-tree-sha1 = "e27c4ebe80e8699540f2d6c805cc12203b614f12"
uuid = "48062228-2e41-5def-b9a4-89aafe57970f"
version = "0.9.20"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "7072f1e3e5a8be51d525d64f63d3ec1287ff2790"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.11"

[[deps.FixedEffectModels]]
deps = ["DataFrames", "FixedEffects", "LinearAlgebra", "Printf", "Reexport", "SnoopPrecompile", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels", "Tables", "Vcov"]
git-tree-sha1 = "510c54fe919707cc4d754f90726cfa4d8c2a8182"
uuid = "9d5cd8c9-2029-5cab-9928-427838db53e3"
version = "1.9.0"

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
git-tree-sha1 = "38a92e40157100e796690421e34a11c107205c86"
uuid = "663a7486-cb36-511b-a19d-713bb74d65c9"
version = "0.10.0"

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
git-tree-sha1 = "744de52d59bd1f7ae18b81b72a7de9ecc884df19"
uuid = "a8dd9ffe-31dc-4cf5-a379-ea69100a8233"
version = "1.0.1"

[[deps.GDAL]]
deps = ["CEnum", "GDAL_jll", "NetworkOptions", "PROJ_jll"]
git-tree-sha1 = "aa6f8ca2f7a0eb46f4d8353eb725c717de40da6e"
uuid = "add2ef01-049f-52c4-9ee2-e494f65e021a"
version = "1.5.1"

[[deps.GDAL_jll]]
deps = ["Arrow_jll", "Artifacts", "Expat_jll", "GEOS_jll", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "LibPQ_jll", "Libdl", "Libtiff_jll", "NetCDF_jll", "OpenJpeg_jll", "PROJ_jll", "Pkg", "SQLite_jll", "Zlib_jll", "Zstd_jll", "libgeotiff_jll"]
git-tree-sha1 = "aa913bff49c25482fe3db2c357cb5f8127a6d2ba"
uuid = "a7073274-a066-55f0-b90d-d619367d196c"
version = "301.600.200+0"

[[deps.GEOS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "818ab247de98d8848a022c7be084b1283d912326"
uuid = "d604d12d-fa86-5845-992e-78dc15976526"
version = "3.11.2+0"

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "cd3e314957dc11c4c905d54d1f5a65c979e4748a"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.2"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoFormatTypes]]
git-tree-sha1 = "434166198434a5c2fcc0a1a59d22c3b0ad460889"
uuid = "68eda718-8dee-11e9-39e7-89f7f65f511f"
version = "0.4.1"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "0eb6de0b312688f852f347171aba888658e29f20"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.3.0"

[[deps.GeoInterfaceRecipes]]
deps = ["GeoInterface", "RecipesBase"]
git-tree-sha1 = "29e1ec25cfb6762f503a19495aec347acf867a9e"
uuid = "0329782f-3d07-4b52-b9f6-d3137cf03c7a"
version = "1.0.0"

[[deps.GeoJSON]]
deps = ["Extents", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "JSON3", "Tables"]
git-tree-sha1 = "624eb2bb45428b1ca3f6b5428aa105028912c71c"
uuid = "61d90e0f-e114-555e-ac52-39dfb47a3ef9"
version = "0.6.4"

[[deps.GeoTables]]
deps = ["ArchGDAL", "GADM", "GeoInterface", "GeoJSON", "Meshes", "Shapefile", "Tables"]
git-tree-sha1 = "358acbc723191444ceaf968c5bbf4370cc0d1828"
uuid = "e502b557-6362-48c1-8219-d30d308dcdb0"
version = "1.2.2"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "303202358e38d2b01ba46844b92e48a3c238fd9e"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.6"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "d3b3624125c1474292d0d8ed0f65554ac37ddb23"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.74.0+2"

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
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.GroupedArrays]]
deps = ["DataAPI", "Missings"]
git-tree-sha1 = "44c812379b629eea08b6d25a196010f1f4b001e3"
uuid = "6407cd72-fade-4a84-8a1e-56e431fc1533"
version = "0.3.3"

[[deps.HDF5_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "OpenSSL_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "4cc2bb72df6ff40b055295fdef6d92955f9dede8"
uuid = "0234f1f7-429e-5d53-9886-15a909be8d59"
version = "1.12.2+2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "37e4657cd56b11abe3d10cd4a1ec5fbdb4180263"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.7.4"

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

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

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

[[deps.ImageMetadata]]
deps = ["AxisArrays", "ImageAxes", "ImageBase", "ImageCore"]
git-tree-sha1 = "36cbaebed194b292590cba2593da27b34763804a"
uuid = "bc367c6b-8a6b-528e-b4bd-a4b897500b49"
version = "0.9.8"

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

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InlineStrings]]
deps = ["Parsers"]
git-tree-sha1 = "9cc2baf75c6d09f9da536ddf58eb2f29dedaf461"
uuid = "842dd82b-1e85-43dc-bf29-5d0ee9dffc48"
version = "1.4.0"

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
git-tree-sha1 = "721ec2cf720536ad005cb38f50dbba7b02419a15"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.7"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "16c0cc91853084cb5f58a78bd209513900206ce6"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.4"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "49510dfcb407e572524ba94aeae2fced1f3feb0f"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.8"

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

[[deps.JSON3]]
deps = ["Dates", "Mmap", "Parsers", "SnoopPrecompile", "StructTypes", "UUIDs"]
git-tree-sha1 = "84b10656a41ef564c39d2d477d7236966d2b5683"
uuid = "0f8b85d8-7281-11e9-16c2-39a750bddbf1"
version = "1.12.0"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "106b6aa272f294ba47e96bd3acbabdc0407b5c60"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.2"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6f2675ef130a300a112286de91973805fcc5ffbc"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.91+0"

[[deps.Kerberos_krb5_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "60274b4ab38e8d1248216fe6b6ace75ae09b0502"
uuid = "b39eb1a6-c29a-53d7-8c32-632cd16f18da"
version = "1.19.3+0"

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
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibPQ_jll]]
deps = ["Artifacts", "JLLWrappers", "Kerberos_krb5_jll", "Libdl", "OpenSSL_jll", "Pkg"]
git-tree-sha1 = "a299629703a93d8efcefccfc16b18ad9a073d131"
uuid = "08be9ffa-1c94-5ee5-a977-46a84ec9b350"
version = "14.3.0+1"

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
git-tree-sha1 = "c7cb1f5d892775ba13767a87c7ada0b980ea0a71"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+2"

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
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "cedb76b37bc5a6c702ade66be44f831fa23c681e"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.0"

[[deps.Lz4_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "5d494bc6e85c4c9b626ee0cab05daa4085486ab1"
uuid = "5ced341a-0733-55b8-9ab6-a4889d929147"
version = "1.9.3+0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "Pkg"]
git-tree-sha1 = "2ce8695e1e699b68702c03402672a69f54b8aca9"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2022.2.0+0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Makie]]
deps = ["Animations", "Base64", "ColorBrewer", "ColorSchemes", "ColorTypes", "Colors", "Contour", "Distributions", "DocStringExtensions", "Downloads", "FFMPEG", "FileIO", "FixedPointNumbers", "Formatting", "FreeType", "FreeTypeAbstraction", "GeometryBasics", "GridLayoutBase", "ImageIO", "InteractiveUtils", "IntervalSets", "Isoband", "KernelDensity", "LaTeXStrings", "LinearAlgebra", "MakieCore", "Markdown", "Match", "MathTeXEngine", "MiniQhull", "Observables", "OffsetArrays", "Packing", "PlotUtils", "PolygonOps", "Printf", "Random", "RelocatableFolders", "Setfield", "Showoff", "SignedDistanceFields", "SnoopPrecompile", "SparseArrays", "StableHashTraits", "Statistics", "StatsBase", "StatsFuns", "StructArrays", "TriplotBase", "UnicodeFun"]
git-tree-sha1 = "e7b6e3eebbadcdfd9f40ad99be84044968a562ee"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.3"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "9926529455a331ed73c19ff06d16906737a876ed"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.3"

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
deps = ["AbstractTrees", "Automa", "DataStructures", "FreeTypeAbstraction", "GeometryBasics", "LaTeXStrings", "REPL", "RelocatableFolders", "Test", "UnicodeFun"]
git-tree-sha1 = "64890e1e8087b71c03bd6b8af99b49c805b2a78d"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.5"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "03a9b9718f5682ecb107ac9f7308991db4ce395b"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.7"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Meshes]]
deps = ["Bessels", "CircularArrays", "Distances", "IterTools", "LinearAlgebra", "NearestNeighbors", "Random", "ReferenceFrameRotations", "SparseArrays", "StaticArrays", "StatsBase", "Tables", "TransformsBase"]
git-tree-sha1 = "8350fc1b783af43c4fb7491acca5dbe2bc5dd8ec"
uuid = "eacbb407-ea5a-433e-ab97-5258b1ca43fa"
version = "0.28.0"

[[deps.MiniQhull]]
deps = ["QhullMiniWrapper_jll"]
git-tree-sha1 = "9dc837d180ee49eeb7c8b77bb1c860452634b0d1"
uuid = "978d7f02-9e05-4691-894f-ae31a51d76ca"
version = "0.4.0"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.Mocking]]
deps = ["Compat", "ExprTools"]
git-tree-sha1 = "782e258e80d68a73d8c916e55f8ced1de00c2cea"
uuid = "78c3b35d-d492-501b-9361-3d52fe80e533"
version = "0.7.6"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "7b86a5d4d70a9f5cdf2dacb3cbe6d251d1a61dbe"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.4"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "2c3726ceb3388917602169bed973dbc97f1b51a8"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.13"

[[deps.NetCDF_jll]]
deps = ["Artifacts", "HDF5_jll", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "072f8371f74c3b9e1b26679de7fbf059d45ea221"
uuid = "7243133f-43d8-5620-bbf4-c2c921802cf3"
version = "400.902.5+1"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "6862738f9796b3edc1c09d0890afce4eca9e7e93"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.4"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "82d7c9e310fe55aa54996e6f7f94674e2a38fcb4"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.9"

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
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "a4ca623df1ae99d09bc9868b008262d0c0ac1e4f"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.4+0"

[[deps.OpenJpeg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libtiff_jll", "LittleCMS_jll", "Pkg", "libpng_jll"]
git-tree-sha1 = "76374b6e7f632c130e78100b166e5a48464256f8"
uuid = "643b3616-a352-519d-856d-80112ee9badc"
version = "2.4.0+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "6503b77492fd7fcb9379bf73cd31035670e3c509"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.3.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9ff31d101d987eb9d66bd8b176ac7c277beccd09"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.20+0"

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
git-tree-sha1 = "d78db6df34313deaca15c5c0b9ff562c704fe1ab"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.5.0"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.40.0+0"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "67eae2738d63117a196f497d7db789821bce61d1"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.17"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "f809158b27eba0c18c269cf2a2be6ed751d3e81d"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.17"

[[deps.PROJ_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "Pkg", "SQLite_jll"]
git-tree-sha1 = "fcb3f39ae1184a056ecc415863d46d2109aa6947"
uuid = "58948b4f-47e0-5654-a9ad-f609743f8632"
version = "900.100.0+0"

[[deps.Packing]]
deps = ["GeometryBasics"]
git-tree-sha1 = "ec3edfe723df33528e085e632414499f26650501"
uuid = "19eb6ba3-879d-56ad-ad62-d5c202156566"
version = "0.5.0"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "84a314e3926ba9ec66ac097e3635e270986b0f10"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.50.9+0"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "478ac6c952fddd4399e71d4779797c538d0ff2bf"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.8"

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
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "SnoopPrecompile", "Statistics"]
git-tree-sha1 = "c95373e73290cf50a8a22c3375e4625ded5c5280"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.4"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

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
deps = ["Crayons", "Formatting", "LaTeXStrings", "Markdown", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "548793c7859e28ef026dba514752275ee871169f"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.3"

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

[[deps.QhullMiniWrapper_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Qhull_jll"]
git-tree-sha1 = "607cf73c03f8a9f83b36db0b86a3a9c14179621f"
uuid = "460c41e3-6112-5d7f-b78c-b6823adb3f2d"
version = "1.0.0+1"

[[deps.Qhull_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "238dd7e2cc577281976b9681702174850f8d4cbc"
uuid = "784f63db-0788-585a-bace-daefebcd302b"
version = "8.0.1001+0"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "6ec7ac8412e83d57e313393220879ede1740f9ee"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.2"

[[deps.RData]]
deps = ["CategoricalArrays", "CodecZlib", "DataAPI", "DataFrames", "Dates", "FileIO", "Requires", "TimeZones", "Unicode"]
git-tree-sha1 = "9a6220c8f59c38ddf6217638042ae6788973f617"
uuid = "df47a6cb-8c03-5eed-afd8-b6050d6c41da"
version = "1.0.0"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "261dddd3b862bd2c940cf6ca4d1c8fe593e457c8"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.3"

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
git-tree-sha1 = "90bc7a7c96410424509e4263e277e43250c05691"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.SQLite_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "54d66b0f69f4578f4988fc08d579783fcdcd764f"
uuid = "76ed43ae-9a5d-5a62-8c75-30186b810ce8"
version = "3.41.0+0"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "30449ee12237627992a99d5e30ae63e4d78cd24a"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.0"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "77d3c4726515dca71f6d80fbb5e251088defe305"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.3.18"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.Shapefile]]
deps = ["DBFTables", "Extents", "GeoFormatTypes", "GeoInterface", "GeoInterfaceRecipes", "OrderedCollections", "RecipesBase", "Tables"]
git-tree-sha1 = "806a1cc22939a77d58e028a098780ee5fc2f130e"
uuid = "8e980c4a-a4fe-5da2-b3a7-4b4b0353a2f4"
version = "0.8.1"

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

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

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

[[deps.SnoopPrecompile]]
deps = ["Preferences"]
git-tree-sha1 = "e760a70afdcd461cf01a575947738d359234665c"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.3"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "a4ada03f999bd01b3a25dcaa30b2d929fe537e00"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.1.0"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "ef28127915f4229c971eb43f3fc075dd3fe91880"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.2.0"

[[deps.StableHashTraits]]
deps = ["CRC32c", "Compat", "Dates", "SHA", "Tables", "TupleTools", "UUIDs"]
git-tree-sha1 = "0b8b801b8f03a329a4e86b44c5e8a7d7f4fe10a3"
uuid = "c5dd0088-6c3f-4803-b00e-f31a60c170fa"
version = "0.3.1"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "b8d897fe7fa688e93aef573711cb207c08c9e11e"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.19"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "45a7769a04a3cf80da1c1c7c60caf932e6f4c9f7"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.6.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "f625d686d5a88bcd2b15cd81f18f98186fdc0c9a"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.3.0"

[[deps.StatsModels]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "Printf", "REPL", "ShiftedArrays", "SparseArrays", "StatsBase", "StatsFuns", "Tables"]
git-tree-sha1 = "06a230063087c11910e9bbd17ccbf5af792a27a4"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.7.0"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "521a0e828e98bb69042fec1809c1b5a680eb7389"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.15"

[[deps.StructTypes]]
deps = ["Dates", "UUIDs"]
git-tree-sha1 = "ca4bccb03acf9faaf4137a9abc1881ed1841aa70"
uuid = "856f2bd8-1eba-4b0a-8007-ebc267875bd4"
version = "1.10.0"

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
git-tree-sha1 = "1544b926975372da01227b382066ab70e574a3ec"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.1"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Thrift_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "boost_jll"]
git-tree-sha1 = "fd7da49fae680c18aa59f421f0ba468e658a2d7a"
uuid = "e0b8ae26-5307-5830-91fd-398402328850"
version = "0.16.0+0"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "8621f5c499a8aa4aa970b1ae381aae0ef1576966"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.4"

[[deps.TimeZones]]
deps = ["Dates", "Downloads", "InlineStrings", "LazyArtifacts", "Mocking", "Printf", "RecipesBase", "Scratch", "Unicode"]
git-tree-sha1 = "a92ec4466fc6e3dd704e2668b5e7f24add36d242"
uuid = "f269a46b-ccf7-5d73-abea-4c690281aa53"
version = "1.9.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

[[deps.TransformsBase]]
deps = ["AbstractTrees"]
git-tree-sha1 = "2412fb54902b0063c69c2bcfbec6b571120cc856"
uuid = "28dd2a49-a57a-4bfb-84ca-1a49db9b96b8"
version = "0.1.2"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.TriplotBase]]
git-tree-sha1 = "4d4ed7f294cda19382ff7de4c137d24d16adc89b"
uuid = "981d1d27-644d-49a2-9326-4793e63143c3"
version = "0.1.0"

[[deps.TupleTools]]
git-tree-sha1 = "3c712976c47707ff893cf6ba4354aa14db1d8938"
uuid = "9d95972d-f1c8-5527-a6e0-b4b365fa01f6"
version = "1.3.0"

[[deps.URIs]]
git-tree-sha1 = "074f993b0ca030848b897beff716d93aca60f06a"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.2"

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
deps = ["Combinatorics", "GroupedArrays", "LinearAlgebra", "StatsAPI", "StatsBase", "Tables"]
git-tree-sha1 = "5ef9c8f67948b2b5e9d93a3da052bab0c1515e7c"
uuid = "ec2bfdc2-55df-4fc9-b9ae-4958c2cf2486"
version = "0.7.0"

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

[[deps.WorkerUtilities]]
git-tree-sha1 = "cd1659ba0d57b71a464a29e64dbc67cfe83d54e7"
uuid = "76eceee3-57b5-4d4a-8e66-0e911cebbf60"
version = "1.6.1"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "93c41695bc1c08c46c5899f4fe06d6ead504bb73"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.10.3+0"

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

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c6edfe154ad7b313c01aceca188c05c835c67360"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.4+0"

[[deps.boost_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "7a89efe0137720ca82f99e8daa526d23120d0d37"
uuid = "28df3c45-c428-5900-9ff8-a3135698ca75"
version = "1.76.0+1"

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
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libgeotiff_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Libtiff_jll", "PROJ_jll", "Pkg"]
git-tree-sha1 = "13dfba87a1fe301c4b40f991d0ec990bbee59bbe"
uuid = "06c338fa-64ff-565b-ac2f-249532af990e"
version = "100.700.100+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

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
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.snappy_jll]]
deps = ["Artifacts", "JLLWrappers", "LZO_jll", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "985c1da710b0e43f7c52f037441021dfd0e3be14"
uuid = "fe1e1685-f7be-5f59-ac9f-4ca204017dfd"
version = "1.1.9+1"

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

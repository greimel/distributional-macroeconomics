### A Pluto.jl notebook ###
# v0.19.14

using Markdown
using InteractiveUtils

# ╔═╡ 597b5e3e-c281-4315-8251-bcf6525a6eaf
using PlutoUI

# ╔═╡ 15beb495-81ee-4d80-89db-2b6cb4500f3b
using PlutoTeachingTools

# ╔═╡ 30be491a-662c-4baf-9242-47bf724f7c94
present_button()

# ╔═╡ 07590a45-f524-4f9c-bb98-718ea373e418
md"""
# The Income and Wealth Distribution in Macroeconomics
"""

# ╔═╡ de8afbd6-6e5f-4155-be79-e55fcd1bc307
md"""
# Measuring Inequality
"""

# ╔═╡ e4134887-c3a9-4360-b2f1-8faf9e5fc42d
md"""
* [Distributional National Accounts](https://gabriel-zucman.eu/usdina/): detailed _micro-data_ on income and wealth for the US, 1946-now (Piketty, Zucman & Saez, 2018, QJE)
* [World Inequality Database](https://wid.world): comparable data on inequality for many countries in the World
* even better: admin data (CBS or similar)

#### Main insights

* income and wealth distributions have Pareto tails
* Pareto tail coefficient is strongly increasing over time in many countries
"""

# ╔═╡ cb32ed19-f6d7-441d-935d-d4f343a6d2cf
md"""
# Modeling Inequality
"""

# ╔═╡ 3be29490-6e7f-48fb-a4f9-b442efac6bb7
md"""
## Modeling Inequality
* regard the economy as a distribution
  * e.g. joint distribution of income ``y`` and wealth ``w``
* how does this distribution look?
* how does this distribution evolve over time?
"""

# ╔═╡ 1ea84570-6433-11ed-3978-955b4b70b598
md"""
## General Approach I
  * distribution is determined by aggregating individual choices
  * individual choices depend on _exogenous_ stochastic process (e.g. exogenous incomes ``y``)
  * critical assumption: _ex-ante identical_ agents (but _ex-post_ heterogeneous)
  * simulate many individuals and look at cross-sectional distribution

## General Approach II
  * determine _state variables_ of interest (e.g. income ``y`` and wealth ``w``)
  * determine _equilibrium objects_ (e.g. interest rate ``r``)
  * Solution strategy
    * derive transition function of an _individual_ ``(y_t, w_t) \stackrel{A(r)}{\mapsto} (y_{t+1}, w_{t+1})``
    * simulate _stationary distribution_ ``F_r(y, w)``
    * iterate on _equilibrium objects_ until _equilibrium conditions_ hold \
      e.g. market clearing ``\int w \; d F_r(y, w) \stackrel{!}{=} K``
"""

# ╔═╡ 9b24fd2b-776f-4c69-b89c-c2f464e73d38
md"""
## The Benchmark Model
* individual agents solve _Income Fluctuations Problem_
```math
\begin{align}
&\max_{(c_t, w_{t+1})_t} E_0\sum_{t=0}^\infty \beta^t u(c_t) \\
&\begin{aligned}\text{s.t. }& c_t + w_{t+1} ≤ y_t + (1+r) w_t \\
& y_t \sim \text{Markov chain} \\
& w_{t+1} \geq \bar w \\
& w_{-1} \text{ given}
\end{aligned}
\end{align}
```
* model solution: policy functions ``c_t^*(w_t, y_t), w_{t+1}^*(w_t, y_t)``
* transition function ``\Pr((w_t, y_t), (w_{t+1}, y_{t+1})) = \Pr(y_t, y_{t+1}) \cdot 1_{w^*_{t+1}(w_t, y_t)}(w_{t+1})``
"""

# ╔═╡ bf6ab362-16e4-42a0-89ef-c6543e84e4b1
md"""
## Modeling Inequality: Outlook

* _Income and Wealth Distribution in Macroeconomics: A Continuous-Time Approach_ (Achdou, Han, Lasry, Lions, Moll, REStud, 2021)
* Solving models with aggregate uncertainty [(Auclert, Bardoczy, Rognlie, Straub, 2021, Ecma)](https://onlinelibrary.wiley.com/doi/abs/10.3982/ECTA17434)
* [EconPDEs.jl](https://github.com/matthieugomez/EconPDEs.jl)
"""

# ╔═╡ 44062795-2591-4cda-a48f-cfd5512b8ca9
md"""
# Recent Research 1: Causes of Inequality
"""

# ╔═╡ 1fae9911-c120-42f0-bd76-31a5b80c8c85
md"""

> How to model the exogenous income process?

* [Gabaix, Lasry, Lions & Moll (2016, Ecma)](https://scholar.harvard.edu/files/xgabaix/files/dynamics_of_inequality.pdf) Superstar entrepreneurs (type dependence) necessary to generate fast increase of Pareto Tail
* [Kim & Jones (2018, JPE)](https://web.stanford.edu/~chadj/papers.html#inequality) Schumpeterian Creative Destruction leads to Pareto Distribution

> How to model returns on wealth?

* [Moll & Restrepo (2022, Ecma)](https://www.econometricsociety.org/publications/econometrica/2022/07/01/uneven-growth-automation’s-impact-income-and-wealth-inequality): Automation's impact on income inequality (return to wealth increasing in wealth)

"""

# ╔═╡ c60d372e-cca2-4c26-8ac3-6f5a75415f73
md"""
# Recent Research 2: Consequences of Rising Inequality
"""

# ╔═╡ 8e57dcb5-d202-4ee5-94d5-b99e1cadd9d5
md"""
## Household Debt and Financial Crises

> What if inequality changes exogenously?

* Kumhof, Ranciere & Winant (2015, AER): financial crises become for likely
* Mian, Sufi & Straub (2021, QJE): Savings of the increases, Debt of non-rich increases, interest rates go up
* Drechsel-Grau & Greimel (2021): Rising house prices and mortgage debt

"""

# ╔═╡ 61d761ca-3c7f-471a-b0a3-182a1071ac09
md"""
# Recent Research 3: Distributional Macroeconomics

* Consumer Credit and Default
* Macroeconomics of Housing Markets
* Environmental Macroeconomics
* Family Macroeconomics (intra-household decision-making)
* Distributional Effects of Monetary Policy
* Optimal Taxation
"""

# ╔═╡ f98217db-d290-4e4f-9e83-bce6584b3f79
md"""
# Appendix
"""

# ╔═╡ 438d5210-374b-41df-887d-bc3f3c625fa8
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoTeachingTools = "~0.2.5"
PlutoUI = "~0.7.48"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.2"
manifest_format = "2.0"
project_hash = "a084e953ecc8d9d6d85465607b1a1d49f1e47fec"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "cc4bd91eba9cdbbb4df4746124c22c0832a460d6"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.1"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

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

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0f960b1404abb0b244c1ece579a0ec78d056a5d1"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.15"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "ab9aa169d2160129beb241cb2750ca499b4e90e9"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.17"

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

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "42324d08725e200c23d4dfb549e0d5d89dede2d2"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.10"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.Parsers]]
deps = ["Dates", "SnoopPrecompile"]
git-tree-sha1 = "cceb0257b662528ecdf0b4b4302eb00e767b38e7"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoHooks]]
deps = ["InteractiveUtils", "Markdown", "UUIDs"]
git-tree-sha1 = "072cdf20c9b0507fdd977d7d246d90030609674b"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0774"
version = "0.0.5"

[[deps.PlutoLinks]]
deps = ["FileWatching", "InteractiveUtils", "Markdown", "PlutoHooks", "Revise", "UUIDs"]
git-tree-sha1 = "8f5fa7056e6dcfb23ac5211de38e6c03f6367794"
uuid = "0ff47ea0-7a50-410d-8455-4348d5de0420"
version = "0.1.6"

[[deps.PlutoTeachingTools]]
deps = ["Downloads", "HypertextLiteral", "LaTeXStrings", "Latexify", "Markdown", "PlutoLinks", "PlutoUI", "Random"]
git-tree-sha1 = "ea3e4ac2e49e3438815f8946fa7673b658e35bdb"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.5"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "efc140104e6d0ae3e7e30d56c98c4a927154d684"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.48"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SnoopPrecompile]]
git-tree-sha1 = "f604441450a3c0569830946e5b33b78c928e1a85"
uuid = "66db9d55-30c0-4569-8b51-7e840670fc0c"
version = "1.0.1"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

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

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"
"""

# ╔═╡ Cell order:
# ╠═30be491a-662c-4baf-9242-47bf724f7c94
# ╟─07590a45-f524-4f9c-bb98-718ea373e418
# ╟─de8afbd6-6e5f-4155-be79-e55fcd1bc307
# ╟─e4134887-c3a9-4360-b2f1-8faf9e5fc42d
# ╟─cb32ed19-f6d7-441d-935d-d4f343a6d2cf
# ╟─3be29490-6e7f-48fb-a4f9-b442efac6bb7
# ╟─1ea84570-6433-11ed-3978-955b4b70b598
# ╟─9b24fd2b-776f-4c69-b89c-c2f464e73d38
# ╟─bf6ab362-16e4-42a0-89ef-c6543e84e4b1
# ╟─44062795-2591-4cda-a48f-cfd5512b8ca9
# ╟─1fae9911-c120-42f0-bd76-31a5b80c8c85
# ╟─c60d372e-cca2-4c26-8ac3-6f5a75415f73
# ╟─8e57dcb5-d202-4ee5-94d5-b99e1cadd9d5
# ╟─61d761ca-3c7f-471a-b0a3-182a1071ac09
# ╟─f98217db-d290-4e4f-9e83-bce6584b3f79
# ╠═597b5e3e-c281-4315-8251-bcf6525a6eaf
# ╠═15beb495-81ee-4d80-89db-2b6cb4500f3b
# ╠═438d5210-374b-41df-887d-bc3f3c625fa8
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

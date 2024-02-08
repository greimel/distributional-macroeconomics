### A Pluto.jl notebook ###
# v0.19.22

#> [frontmatter]
#> chapter = 1
#> section = 2
#> order = 1
#> title = "More Julia"
#> layout = "layout.jlhtml"
#> description = ""
#> tags = ["julia-basics"]

using Markdown
using InteractiveUtils

# ╔═╡ be621538-c0dd-4d71-8aa5-1a3c7d6bd5f0
using DataFrames

# ╔═╡ 56f2b729-df63-457e-899c-180e4be46bf0
using DataFrameMacros

# ╔═╡ 3580cde2-f96d-4229-a6ed-fcb046eb9234
using Chain: @chain

# ╔═╡ b9aea038-b0f5-4101-b8d8-840caf80e9b2
using PlutoUI

# ╔═╡ ba649465-3d71-4854-9151-bf6147225865
md"""
`more-julia.jl` | **Version 1.1** | *last updated: May 10 2022* | *created by [Daniel Schmidt](https://github.com/danieljschmidt)*
"""

# ╔═╡ c704103c-3093-4f50-971e-c37eb14546de
md"""
# More Julia
"""

# ╔═╡ e44404db-3327-4c52-8841-b5c10287a08e
md"""
The purpose of this notebook is to make you familiar with Julia syntax that is frequently used in the course material. It also points out the peculiarities of Pluto notebooks.
"""

# ╔═╡ ae686dfa-a6e1-4d30-89bf-03d1ea809849
md"""
## Named tuples and keyword arguments

### Named tuples

Like many other programming languages, Julia allows you to create tuples. The elements of tuples can be accessed by using the corresponding index:
"""

# ╔═╡ 39f1cb09-3c61-464b-8258-2de5d446217b
standard_tuple = (1, 2)

# ╔═╡ 29f29203-054d-4353-9431-59e0e5bb575a
standard_tuple[1]

# ╔═╡ 0498aa7b-14bf-441f-9760-c75c8c39382d
md"""
Julia also has named tuples. Just like standard tuples, named tuples allow you to access an element using its index, but alternatively you can also access an element by its name.
"""

# ╔═╡ 85088a9a-18d0-431a-b56f-31a8e3863fd8
named_tuple = (a = 1, b = 2)

# ╔═╡ b9912c59-59db-4b06-9603-68b14b2e0efd
named_tuple.a

# ╔═╡ b85abb17-1e62-481e-9134-16c1c6d7ff04
md"""
This is convenient because it means that you do not need to remember if some parameter is the first element in a tuple or the second one.

Below you can see an alternative way of creating named tuples that is frequently used in the course material:
"""

# ╔═╡ 63e5adb0-529f-4396-866d-85c14828e56a
let 
	a = 1
	named_tuple2 = (; a, b=2) # (; a) is equivalent to (; a=a)
end

# ╔═╡ 3f80dd01-7203-4404-b38c-32757b2f4223
md"""
### Keyword arguments

A similar syntax with a semicolon is used for keyword arguments that are identified by name (and not by their position as normal function arguments):
"""

# ╔═╡ aef759a7-a0ed-41e5-b91b-f276853bd30b
function some_function(; a, b = 1)
	return a + b
end

# ╔═╡ 746beef7-2a9d-43c1-9570-82f2fdf79f0c
some_function(; a = 3, b = 5)

# ╔═╡ 654400a3-598b-41d1-823e-3824f04a6542
let
	a = 3
	some_function(; a, b = 5) # (; a) is equivalent to (; a=a)
end

# ╔═╡ d705a9d5-97aa-4e7f-adcc-24ad53d9079b
md"""
## Vectorization with dot syntax

You can apply a function to all elements of a vector by using the dot syntax:
"""

# ╔═╡ 8e6f94b0-d7ce-4461-b40b-8f8f45de6c80
[1,2,3] .^ 2

# ╔═╡ 77244481-f942-4d09-a859-e3db7185895a
log.([1,2,3])

# ╔═╡ 76d105ee-7d73-43a3-946d-5c4b7181a2a3
md"""
Note how Julia usually does things in a way that's mathematically consistent. Look at the following code.
"""

# ╔═╡ 594f37a2-654b-4eb4-8ee9-2b8609eda73f
A = ones(3, 3)

# ╔═╡ 56402ce6-964c-4ca0-b3a9-917987c976ea
exp(A)

# ╔═╡ 9172b664-e56b-4c26-a950-dc901fa68482
exp.(A)

# ╔═╡ a2039109-41ab-4358-8b28-3e1be577d8ec
A^2

# ╔═╡ d4cca863-6b61-4abc-93df-8d9553760cf6
A .^ 2

# ╔═╡ 88b83d9e-7f2a-4539-a133-577f98d32889
md"""
What would Python, R or Matlab do?
"""

# ╔═╡ 4c983a8b-abe9-4cc5-a7ab-102a0d05eedd
md"""
## The pipe operator

The pipe operator |> makes it possible to write down nested function calls in a more readable way. For example, the two expressions below do the same thing:
"""

# ╔═╡ ae60f715-86b7-45ff-8080-4b961cdfb465
round(log(named_tuple.b))

# ╔═╡ 19638c7f-1cef-41af-bdf1-b48e81503904
named_tuple.b |> log |> round

# ╔═╡ 5d0b2c96-0fea-4e56-be58-0beba790d88a
md"""
(In case you use R: The |> operator in Julia is similar to the %>% operator in R.)
"""

# ╔═╡ c6abdf8a-9f2a-469f-992a-ed64589bd4a0
md"""
# Working with data

Here are a few simple examples of working with DataFrames in Julia. If you are familiar with pandas, dplyr or Stata, have a look at [this cheatsheet](https://dataframes.juliadata.org/stable/man/comparisons/). Also, have a look at the documentation of [DataFrameMacros.jl](https://jkrumbiegel.github.io/DataFrameMacros.jl/stable) which makes working with DataFrames much more convenient.
"""

# ╔═╡ dbcb3e8e-2261-4861-9f0d-cd208e084718
md"""
## Creating a `DataFrame`
"""

# ╔═╡ 1c8ead52-c270-4a8a-9ad4-6b25af1603e4
df1 = DataFrame(x = 1:10)

# ╔═╡ 62e73de9-e6f1-4e30-b739-394e59d54879
md"""
## Transforming (mutating) a column
"""

# ╔═╡ c7fa1dbc-3133-4032-8bc5-f8bfefe61dd2
transform(df1, :x => ByRow(x -> sqrt(x)))

# ╔═╡ e521ab22-cc5a-4575-b20d-852ac239e644
@transform(df1, sqrt(:x)) # uses @transform from DataFrameMacros.jl

# ╔═╡ 7f6fc125-2001-4396-bcb2-8a2e1e0cf036
md"""
## The ```@chain``` macro
"""

# ╔═╡ b4c63d80-7e17-40be-bb0d-678c6a1a6a0a
md"""
The ```@chain``` macro works similar to the pipe operator. In the code for this course, the ```@chain``` macro is often applied to data frames together with macros from the ```DataFrameMacros``` package.
"""

# ╔═╡ 9751ec14-ef27-4d67-89a1-4f1be103a81d
df2 = DataFrame(A=[1,2,2,1], B=randn(4))

# ╔═╡ 4682cb1e-77e9-4a23-ae4b-cd18f4c58872
md"""
Consider the data frame above. Let's say you would like to
- add up the values in the B column separately for each value of A
- take the absolute value of the resulting sums of B values.

Using the ```@chain``` macro, we can perform this task with relatively concise code:
"""

# ╔═╡ cecea5f2-5446-4bf9-b186-11834ac25d37
@chain df2 begin
	@groupby(:A)
	@combine(
		sum(:B), # automatically named
		:C = sum(:B) # specify name
	)
	@transform(abs(:C))
end

# ╔═╡ 925436fb-e535-4735-888a-b047693c9132
md"""
Without the ```@chain``` macro, the code would look like this:
"""

# ╔═╡ f0d4415c-4f2b-415e-a600-cc5784811ee9
begin
	df_groups = groupby(df2, :A)
	df_sum = @combine(df_groups, sum(:B), :C = sum(:B))
	@transform(df_sum, abs(:C))
end

# ╔═╡ 5cc63149-7edf-45b0-a198-d0f308a52618
md"""
## Unicode characters

You can use Greek letters and other Unicode characters in your Julia code. For example, type "\alpha" in the cell below (without the quotation marks) and press Tab on your keyboard. This should create an $\alpha$ symbol.
"""

# ╔═╡ aa02b393-9b84-4b9f-9922-c79d880113e6


# ╔═╡ bd064547-0ad0-4cf6-91a1-6c7e82e86b24
md"""
See the [Julia documention](https://docs.julialang.org/en/v1/manual/unicode-input/) for a list of supported Unicode characters. For Greek letters, the abbreviations are the same as in LaTeX.
"""

# ╔═╡ 9e327f16-5320-4ca4-9a5c-f2049afa9fac
md"""
### The $\in$ symbol
An elegant way of writing loops is to use the $\in$ symbol instead of writing "in". The $\in$ symbol can be created by typing "\in" and pressing Tab.
"""

# ╔═╡ 1ed23fc8-090b-476b-a6bf-d79f3602133b
[i^2 for i ∈ 1:5]

# ╔═╡ d5ddeb01-d83e-4a94-98c9-bb9459668b14
md"""
# Pluto notebooks

## General advice

* Press F1 to see shortcuts for Pluto notebooks.

* Ctrl + Click on an underlined variable or function to jump to its definition.

* Use the Live docs in the bottom right corner to get more information about any Julia function or object.

* Check the [Github wiki](https://github.com/fonsp/Pluto.jl/wiki) for more information on Pluto notebooks.

## Automated updating of cells

When changing a function or variable, Pluto automatically updates all affected cells. For example, change the value of $d$ to some other number and see how the following cell updates automatically:
"""

# ╔═╡ 542a9921-e9a2-448f-877b-e27c418eee35
d = 5

# ╔═╡ 7c2cf481-6754-4ddf-9bfa-9a04818cea4a
d + 10

# ╔═╡ de8f50c0-544c-4430-8548-76bff0b622cb
md"""
This is different from jupyter notebooks in which you have to update related cells manually. 

The automated updating the advantage that you do not have to keep in mind the order in which to evaluate cells. However, it can also be annoying if some of the affected cells take a long time to run. At the time of writing, there is no way to turn off the automated updating but you can always manually disable cells with long run times by clicking on the three dots in the top right corner of a cell.
"""

# ╔═╡ 70046f23-5792-47ba-bfdd-52549fb6e68a
md"""
## Only one expression per cell

Pluto notebooks only allow one expression per cell. If you nevertheless want to place several expressions into the same cell, you have to use a begin ... end block:
"""

# ╔═╡ 43760b9f-d48c-44e3-be34-2f728b770260
begin
	e = 10
	f = e + 5
end

# ╔═╡ da466543-2824-4af3-81ce-98e55259802a
md"""
## Cannot reuse variable names
"""

# ╔═╡ 45468076-5bcc-442b-8f84-5c152a679129
md"""
Use `let` blocks to specify variable names locally.
"""

# ╔═╡ 553beaca-357b-4ef3-9d66-6b8146eb9a2e
let
	g = 3
end

# ╔═╡ 4e6ffc08-02ee-41e9-ba47-59e465138b2c
g

# ╔═╡ 9f32f1e8-a0b3-4a57-9147-30a9c610a4ee
md"""
## Deactivate cells
"""

# ╔═╡ e2e7d66b-18fc-4b29-919e-bb81b1e30b4d
# ╠═╡ disabled = true
#=╠═╡
sleep(5)
  ╠═╡ =#

# ╔═╡ f0fe4fcd-bbd9-40ad-987c-62a3eb299f21
md"""
# Imported Packages
"""

# ╔═╡ 07374ec1-2f35-431a-9001-4f9e47cf40c5
TableOfContents()

# ╔═╡ 54416ff1-ef91-4ff3-968f-22acae580f4f
g = 1

# ╔═╡ 7b78d185-12c6-4baa-8ffc-b5a806cda32f
g = 2

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Chain = "~0.5.0"
DataFrameMacros = "~0.4.0"
DataFrames = "~1.4.4"
PlutoUI = "~0.7.49"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "26781b0a3b356446156365442bf774b44842e35e"

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

[[deps.Chain]]
git-tree-sha1 = "8c4920235f6c561e401dfe569beb8b924adad003"
uuid = "8be319e6-bccf-4806-a6f7-6fae938471bc"
version = "0.5.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.1+0"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

[[deps.DataFrameMacros]]
deps = ["DataFrames", "MacroTools"]
git-tree-sha1 = "92ae44e8d08667be722ca197c97e60bcff1db968"
uuid = "75880514-38bc-4a95-a458-c2aea5a3a702"
version = "0.4.0"

[[deps.DataFrames]]
deps = ["Compat", "DataAPI", "Future", "InvertedIndices", "IteratorInterfaceExtensions", "LinearAlgebra", "Markdown", "Missings", "PooledArrays", "PrettyTables", "Printf", "REPL", "Random", "Reexport", "SnoopPrecompile", "SortingAlgorithms", "Statistics", "TableTraits", "Tables", "Unicode"]
git-tree-sha1 = "d4f69885afa5e6149d0cab3818491565cf41446d"
uuid = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
version = "1.4.4"

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

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

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

[[deps.InvertedIndices]]
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

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

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "f66bdc5de519e8f8ae43bdc598782d35a25b1272"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.1.0"

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
git-tree-sha1 = "151d91d63d8d6c1a5789ecb7de51547e00480f1b"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.4"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eadad7b14cf046de6eb41f13c9275e5aa2711ab6"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.49"

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
git-tree-sha1 = "96f6db03ab535bdb901300f88335257b0018689d"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.2.2"

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

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

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

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

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
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

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
git-tree-sha1 = "ac00576f90d8a259f2c9d823e91d1de3fd44d348"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.1"

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
# ╟─ba649465-3d71-4854-9151-bf6147225865
# ╟─c704103c-3093-4f50-971e-c37eb14546de
# ╟─e44404db-3327-4c52-8841-b5c10287a08e
# ╟─ae686dfa-a6e1-4d30-89bf-03d1ea809849
# ╠═39f1cb09-3c61-464b-8258-2de5d446217b
# ╠═29f29203-054d-4353-9431-59e0e5bb575a
# ╟─0498aa7b-14bf-441f-9760-c75c8c39382d
# ╠═85088a9a-18d0-431a-b56f-31a8e3863fd8
# ╠═b9912c59-59db-4b06-9603-68b14b2e0efd
# ╟─b85abb17-1e62-481e-9134-16c1c6d7ff04
# ╠═63e5adb0-529f-4396-866d-85c14828e56a
# ╟─3f80dd01-7203-4404-b38c-32757b2f4223
# ╠═aef759a7-a0ed-41e5-b91b-f276853bd30b
# ╠═746beef7-2a9d-43c1-9570-82f2fdf79f0c
# ╠═654400a3-598b-41d1-823e-3824f04a6542
# ╟─d705a9d5-97aa-4e7f-adcc-24ad53d9079b
# ╠═8e6f94b0-d7ce-4461-b40b-8f8f45de6c80
# ╠═77244481-f942-4d09-a859-e3db7185895a
# ╟─76d105ee-7d73-43a3-946d-5c4b7181a2a3
# ╠═594f37a2-654b-4eb4-8ee9-2b8609eda73f
# ╠═56402ce6-964c-4ca0-b3a9-917987c976ea
# ╠═9172b664-e56b-4c26-a950-dc901fa68482
# ╠═a2039109-41ab-4358-8b28-3e1be577d8ec
# ╠═d4cca863-6b61-4abc-93df-8d9553760cf6
# ╟─88b83d9e-7f2a-4539-a133-577f98d32889
# ╟─4c983a8b-abe9-4cc5-a7ab-102a0d05eedd
# ╠═ae60f715-86b7-45ff-8080-4b961cdfb465
# ╠═19638c7f-1cef-41af-bdf1-b48e81503904
# ╟─5d0b2c96-0fea-4e56-be58-0beba790d88a
# ╟─5cc63149-7edf-45b0-a198-d0f308a52618
# ╠═aa02b393-9b84-4b9f-9922-c79d880113e6
# ╟─bd064547-0ad0-4cf6-91a1-6c7e82e86b24
# ╟─9e327f16-5320-4ca4-9a5c-f2049afa9fac
# ╠═1ed23fc8-090b-476b-a6bf-d79f3602133b
# ╟─c6abdf8a-9f2a-469f-992a-ed64589bd4a0
# ╠═be621538-c0dd-4d71-8aa5-1a3c7d6bd5f0
# ╟─dbcb3e8e-2261-4861-9f0d-cd208e084718
# ╠═1c8ead52-c270-4a8a-9ad4-6b25af1603e4
# ╟─62e73de9-e6f1-4e30-b739-394e59d54879
# ╠═c7fa1dbc-3133-4032-8bc5-f8bfefe61dd2
# ╠═56f2b729-df63-457e-899c-180e4be46bf0
# ╠═e521ab22-cc5a-4575-b20d-852ac239e644
# ╟─7f6fc125-2001-4396-bcb2-8a2e1e0cf036
# ╟─b4c63d80-7e17-40be-bb0d-678c6a1a6a0a
# ╠═3580cde2-f96d-4229-a6ed-fcb046eb9234
# ╠═9751ec14-ef27-4d67-89a1-4f1be103a81d
# ╟─4682cb1e-77e9-4a23-ae4b-cd18f4c58872
# ╠═cecea5f2-5446-4bf9-b186-11834ac25d37
# ╟─925436fb-e535-4735-888a-b047693c9132
# ╠═f0d4415c-4f2b-415e-a600-cc5784811ee9
# ╟─d5ddeb01-d83e-4a94-98c9-bb9459668b14
# ╠═542a9921-e9a2-448f-877b-e27c418eee35
# ╠═7c2cf481-6754-4ddf-9bfa-9a04818cea4a
# ╟─de8f50c0-544c-4430-8548-76bff0b622cb
# ╟─70046f23-5792-47ba-bfdd-52549fb6e68a
# ╠═43760b9f-d48c-44e3-be34-2f728b770260
# ╟─da466543-2824-4af3-81ce-98e55259802a
# ╠═54416ff1-ef91-4ff3-968f-22acae580f4f
# ╠═7b78d185-12c6-4baa-8ffc-b5a806cda32f
# ╟─45468076-5bcc-442b-8f84-5c152a679129
# ╠═553beaca-357b-4ef3-9d66-6b8146eb9a2e
# ╠═4e6ffc08-02ee-41e9-ba47-59e465138b2c
# ╟─9f32f1e8-a0b3-4a57-9147-30a9c610a4ee
# ╠═e2e7d66b-18fc-4b29-919e-bb81b1e30b4d
# ╟─f0fe4fcd-bbd9-40ad-987c-62a3eb299f21
# ╠═b9aea038-b0f5-4101-b8d8-840caf80e9b2
# ╠═07374ec1-2f35-431a-9001-4f9e47cf40c5
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

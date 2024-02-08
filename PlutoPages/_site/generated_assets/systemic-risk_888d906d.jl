### A Pluto.jl notebook ###
# v0.19.22

#> [frontmatter]
#> chapter = 5
#> section = 2
#> order = 1
#> title = "Systemic risk in financial networks"
#> layout = "layout.jlhtml"
#> description = ""
#> tags = ["financial-networks"]

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

# ╔═╡ 36ac25b9-27e5-4728-bbcf-920f231ff6ab
using Revise

# ╔═╡ 4878d013-053d-4490-b7c6-db0de4bf0a82
using PlutoUI: NumberField

# ╔═╡ 74262581-3c64-4e5b-9328-416c4e1efc91
using MarkdownLiteral: @markdown

# ╔═╡ 2b405f2b-3256-4c47-8334-c2d93713d409
using PlutoUI: Select

# ╔═╡ 103babf9-bbc3-4c77-9185-72f792a09706
using LinearAlgebra: dot

# ╔═╡ bb0e41cf-66fb-4cae-8cd9-6ad771d1acf4
using HypertextLiteral

# ╔═╡ 969aa5eb-c56f-4115-83df-bb0ec911b6aa
using MarkdownLiteral: @mdx

# ╔═╡ 935da683-a4e2-4ddc-999f-55cb61390f39
using StructArrays

# ╔═╡ 2d1da571-0561-4c28-bb63-35ca5f9538d5
using PlutoUI: PlutoUI, TableOfContents, Slider, CheckBox, as_svg

# ╔═╡ 7c068d65-f472-44c2-af56-581bf9309bd5
using StatsBase: weights

# ╔═╡ 631ce85b-db76-4f3b-bda5-0c51ffb3bc43
using CategoricalArrays

# ╔═╡ fa01788e-7585-4f3f-93b2-a44541065820
using DataAPI: refarray

# ╔═╡ 6a89953d-310b-49c9-89e1-8d51c8b75be0
using Graphs, SimpleWeightedGraphs

# ╔═╡ 50ac0182-32fe-4e21-8c6d-895ffc67ec27
using NetworkLayout

# ╔═╡ 8b75bc15-e07b-43e5-adb3-c5d6481ee9d8
using LinearAlgebra: I

# ╔═╡ cf61ebba-6800-4bfa-bb7c-cb9a6c846b65
using PlutoLinks

# ╔═╡ 954d2dde-8088-4e91-bed3-f8339090b77b
using PlutoTeachingTools

# ╔═╡ 2872c686-8e4f-4230-a07a-5d988aba39b7
using Chain: @chain

# ╔═╡ 3198bb79-f0c7-4b01-8dce-ef7629e8d7e6
using GeometryBasics: Rect, Vec2f0

# ╔═╡ e9365f8d-6f58-458c-855d-c0444c6f409f
using DataFrames, DataFrameMacros

# ╔═╡ 955c190a-1d00-4c78-bdfb-daac31edf76f
using Roots: find_zero

# ╔═╡ fc01674e-73de-4cbf-9c80-fa1ea47fcb21
using Colors: @colorant_str

# ╔═╡ 00a4054a-dd10-4c5c-ae20-c0dc176e8e18
using Statistics: mean

# ╔═╡ 7d2d3618-7810-444d-90c0-0d592f8eba8c
using GraphMakie

# ╔═╡ 12968f5e-6e75-4429-bff8-0d1c644330a7
using NamedTupleTools

# ╔═╡ 8f15ccf5-860a-47ca-b5c2-842c4e6f861a
using CairoMakie: @key_str, @lift, @L_str

# ╔═╡ 3758f272-297c-4325-9a7c-042b4f41d615
using CairoMakie, AlgebraOfGraphics, Colors

# ╔═╡ a349099a-64aa-4e86-9bf2-b6157feff394
using LaTeXStrings

# ╔═╡ c9775558-94f2-4c8c-9d0e-ab948fa5ead4
using PlutoTest: @test

# ╔═╡ 52052d98-0c41-45ec-95bf-d936b1c43e81
md"""
`systemic-risk.jl` | **Version 2.3** | *last updated: March 10, 2023*
"""

# ╔═╡ 72e25b9c-89e3-441b-bf89-c1122535318a
present_button()

# ╔═╡ a4dba7da-9d63-4d0d-a736-565d6ccf7d09
space = 0

# ╔═╡ cdfc35e1-5082-4959-9d30-3abb2dc8a7ac
md"""
# Systemic Risk in Financial Networks

based on _[Acemoglu, Ozdaglar & Tahbaz-Salehi, 2015](https://www.aeaweb.org/articles?id=10.1257/aer.20130456), American Economic Review_
"""

# ╔═╡ 944b07ed-20e5-4068-ac5e-8a48f866fdd2
md"""
## Overview: Contagion on Financial Networks

> _Financial network_:  banks linked by interbank lending
"""

# ╔═╡ 144cfaf2-b78c-4f87-8b35-559914abf532
aside_figure(fig) = aside(md"$(as_svg(fig))")

# ╔═╡ eba94b09-f061-432a-80ce-a68be83e6a99
md"""
# I. Model Setup
"""

# ╔═╡ ffc743af-b97c-4083-a02e-ea5725829f2a
md"""
## Model I: Banks
"""

# ╔═╡ acfdf348-e8bf-41ee-aa19-fe7ec29087cc
md"""
## Model II: Firms

* fixed initial size
* dividend ``\delta = a - ε`` paid to bank (``ε`` is a shock)
* fraction ``\tilde ℓ`` can be liquidated (recover ``\tilde ℓ ζ A``)
* final pay-off ``(1-\tilde ℓ) A``

"""

# ╔═╡ 2f69a836-7fed-42e6-a509-096bc8cabbc2
dividend((; a), ε) = max(a - ε, 0.0)

# ╔═╡ 26a5c4f4-75c1-4a4c-a813-5fad0c571da1
recovery((; ζ, A), ℓ) = ℓ * ζ   * A

# ╔═╡ 5d283f89-4419-4b14-81ba-9826b4f1689e
function cashflow₁(firm, ε, ℓ)
	div = dividend(firm, ε)
	rec = recovery(firm, ℓ)
	cf = div + rec
	(; cf, div, rec)
end

# ╔═╡ a0ac4092-cd27-4523-85f9-f4a4d81456b3
md"""
# II. Payment equilibrium
"""

# ╔═╡ 01eb57fd-a815-41f2-9e25-7730bff7917d
md"""
### Banks' choice

1. payables > receivables: repay all debt (**_solvent_**)
2. not enough? liquidate (part) of the firm ``ℓ > 0`` (**_insolvent_**)
3. still not enough? (partially) default on interbank debt (**_bankrupt_**)
4. still not enough? (partially) default on external liabilities ``ν``
"""

# ╔═╡ 6af1e245-aaf8-4b7d-977a-1d76998f7cf8


# ╔═╡ 5b6f26db-f63f-4db5-ae3e-d41ae476948f


# ╔═╡ 0e4a53b5-1751-4723-a05d-a5504e427e3c


# ╔═╡ 523e177c-f74d-4805-8a2b-9c27a4b0bc63
md"""
* size of shock $(@bind ε_cash Slider(0.0:0.05:1.0, show_value = true, default = 0.0))
* show illiquid firm value ``A`` $(@bind show_illiquid CheckBox(default = false))
* recovery rate ``ζ`` $(@bind recovery_rate Slider(0:0.1:0.5, show_value = true, default = 0.0))
""" |> aside

# ╔═╡ 7f058e4a-9b12-41c9-9fd7-4ad023058a14
md"""
# III. Financial Contagion

* default on interbank market 
* ``\implies`` lower interbank receivables for other banks
* ``\implies`` reduces distance to default
"""

# ╔═╡ 96878ebb-fbc0-4d53-998a-210e13a42492
md"""
| variable            | value                                                        |
|:------------------- | :----------------------------------------------------------- |
| ring vs complete ``γ`` | $(@bind γ2 Slider(0:0.01:1, default=0.5, show_value=true))|
| dividend ``a``      | $(@bind a Slider(0:0.01:1, default=0.7, show_value=true))    |
| recovery rate ``ζ`` | $(@bind ζ Slider(0:0.01:1, default=0.3, show_value=true))    |
| **external**        | |
| cash ``c``          | $(@bind c Slider(0:0.1:2, default=1.2, show_value=true))     |
| shock ``ε``         | $(@bind ε Slider(0:0.1:2, default=0.0, show_value=true))     |
| deposits ``ν``      | $(@bind ν Slider(0.0:0.1:2, default=1.8, show_value=true))   |
| **internal**        | |
| interbank ``y``     | $(@bind y Slider(0:0.1:2, default=1.0, show_value=true))     |
""" |> aside

# ╔═╡ 9446b4db-4d93-4153-84e2-73d01fb31254


# ╔═╡ da7c558d-a2b5-41a8-9c78-3e39a00dfd31
md"""
# IV. Stability and Resilience
"""

# ╔═╡ bf719f30-72c1-488e-ba77-c183effb7c60
md"""
### 1. More interbank lending leads to less stability and less resilience

Formally, for a given regular financial network ``(y_{ij})`` let ``\tilde y_{ij} = \beta y_{ij}`` for ``\beta > 1``. Financial network ``\tilde y`` is less resilient and less stable *(see __Proposition 3__)*.

* _Proposition 3_: More interbank lending leads to a more fragile system. (For a given network topology, the defaults are increasing in ``\bar y``.)
"""

# ╔═╡ a0767d80-0857-47ef-90a1-72bc34064716
md"""
### 2. Densely connected networks are _robust, yet fragile_

> Our results thus confirm a conjecture of Andrew Haldane (2009, then Executive Director for Financial Stability at the Bank of England), who suggested that highly interconnected financial networks may be “robust-yet-fragile” in the sense that “within a certain range, connections serve as shock-absorbers [and] connectivity engenders robustness.” However, beyond that range, interconnections start to serve as a mechanism for the propagation of shocks, “the system [flips to] the wrong side of the knife-edge,” and fragility prevails.

* Upto some ``\varepsilon^*``, the shock does not spread
* Above this value, all banks default
* This is an example for a _phase transition_: it flips from being the most to the least stable and resilient network

__Compare *Propositions 4 and 6*__:
* _Proposition 4_: For small shocks and big ``\bar y``, the complete network ist the most resilitient and most stable financial network and the ring is the least resilient and least stable.

* _Proposition 6_: For big shocks, the ring and complete networks are the least resilient and least stable.

"""

# ╔═╡ c920d82c-cfe9-462a-bacd-436f01c314cf
A = 0.5

# ╔═╡ 72e5f4d3-c3e4-464d-b35c-2cf19fa9d4b5
md"""
# Assignment
"""

# ╔═╡ d6345a8b-6d8f-4fd2-972b-199412cbdc26
md"""
## Task 1 (6 Points)
"""

# ╔═╡ 51bfc57f-7b06-4e27-af32-51df38af30a1
md"""
👉 (1.1 | 2 points) Consider financial network ``\tilde y``. What is the minimal number of shocks ``\tilde p`` that have to occur to wipe out the whole financial system. Which banks have to be hit? How big would the shocks ``\tilde \varepsilon`` have to be?
"""

# ╔═╡ b817cdf6-dfe6-4607-98da-2299a33d4906
answer11 = md"""
Your answer goes here ... You can type math like this: ``\tilde p = 17``, ``\tilde \varepsilon = 1.1``
"""

# ╔═╡ 07e9c77f-e05e-452c-8c47-cdd9dfc8e2fc
md"""
👉 (1.2 | 2 points) Consider financial network ``\hat y``. What is the minimal number of shocks ``\hat p`` that have to occur to wipe out the whole financial system. Which banks have to be hit? How big would the shocks ``\hat \varepsilon`` have to be?
"""

# ╔═╡ f5293bee-9413-4507-8568-54836eb6d4a2
answer12 = md"""
Your answer goes here ... You can type math like this: ``\hat p = 17``, ``\hat \varepsilon = 1.1``
"""

# ╔═╡ 49c2fb2d-de6e-4ab2-a558-59fb153cf703
md"""
👉 (1.3 | 2 points) Now consider ``\varepsilon > \max\{\tilde \varepsilon, \hat \varepsilon \}`` and ``p = 1``. Which network is more stable? Which network is more resilient?
"""

# ╔═╡ c0c711de-6916-4ab9-ab73-c476654332c4
answer13 = md"""
Your answer goes here ... You can type math like this: ``\hat p = 17``, ``\tilde \varepsilon = 1.1``
"""

# ╔═╡ a95431eb-14a0-4dc3-bbe6-9c409f6cc596
md"""
## Task 2 (4 Points)

Consider the model of systemic risk by _Acemoglu, Ozdaglar & Tahbaz-Salehi (2015)_ with $n = 5$ banks. There are two scenarios that differ by the structure of the financial network ($\tilde y$ vs. $\hat y$).  Both $\tilde y$ and $\hat y$ are regular financial networks with parameter $y$.
"""

# ╔═╡ f00d9e1a-b111-4b6a-95f5-b9736329befe
md"""
Each bank owns a project that pays $z_i = a - \varepsilon_i$ in the intermediate period. $\varepsilon_i > 0$ is a shock. Let $\tilde p$ be the minimal number of shocks, and $\tilde \varepsilon$  be the minimal shock size that can wipe out the whole system (that is, a simultaneous default of all banks) under the financial network $\tilde y$.
"""

# ╔═╡ f8eb242f-a974-48aa-9173-b0bc7ff697d5
md"""
👉 (2.1 | 2 points) What is $\tilde p$? Explain which banks would have to be hit to wipe out the whole system.
"""

# ╔═╡ c2633df1-2e30-4387-8749-de3280b0602d
answer21 = md"""
Your answer goes here ... You can type math like this: ``p = 17``, ``\varepsilon = 1.1``
"""

# ╔═╡ 253ab06f-6284-4cbf-b2a2-232ff99548c9
md"""
👉 (2.2 | 2 points) Let $\delta > 0$ be very small and shock size $\varepsilon' \in (\tilde \varepsilon - \delta, \tilde \varepsilon)$. Conditional on $p' = 1$ and $\varepsilon'$, which network is more resilient and which network is more stable? Explain.
"""

# ╔═╡ 1d058f8b-16f5-4744-8425-452876006c47
answer22 = md"""
Your answer goes here ... You can type math like this: ``p = 17``, ``\varepsilon = 1.1``
"""

# ╔═╡ 27fadf93-0b17-446e-8001-d8394b7befaa
md"""
### Before you submit ...

👉 Make sure you have added **your names** and **your group number** in the cells below.

👉 Make sure that that **all group members proofread** your submission (especially your little essays).

👉 Go to the very top of the notebook and click on the symbol in the very top-right corner. **Export a static html file** of this notebook for submission. (The source code is embedded in the html file.)
"""

# ╔═╡ aed99485-cec3-4bf3-b05d-4d20572ec907
group_members = ([
	(firstname = "Ella-Louise", lastname = "Flores"),
	(firstname = "Padraig", 	lastname = "Cope"),
	(firstname = "Christy",  	lastname = "Denton")
	]);

# ╔═╡ db841316-9106-40bb-9ca3-ae6f8b975404
group_number = 99

# ╔═╡ 639f04a2-0979-41e2-b465-a4f52049166d
if group_number == 99 || (group_members[1].firstname == "Ella-Louise" && group_members[1].lastname == "Flores")
	md"""
!!! danger "Note!"
    **Before you submit**, please replace the randomly generated names above by the names of your group and put the right group number in the appropriate cell.
	"""
end

# ╔═╡ 900a4b24-04dc-4a1b-9829-a166cf9eb7fb
md"""
## Task 3 (not graded): Avoiding a bank run
"""

# ╔═╡ 871f33f0-4882-4ff0-bbde-eb954059e907
md"""
Consider the setup of Allen & Gale with banks ``i \in \{1, 2, 3, 4\}``. Banks know that a fraction ``\gamma`` of the population are _early types_. In the social optimum, banks offer deposit contracts ``(c_1, c_2)``. The fraction of early types ``\omega_i`` in each bank is random. There are three possible states ``S_j = (\omega_{1j}, \omega_{2j}, \omega_{3j}, \omega_{4j})``

```math
\begin{align}
S_1 &= (\gamma, \gamma, \gamma, \gamma) \\
S_2 &= (\gamma + \varepsilon, \gamma + \varepsilon, \gamma - \varepsilon, \gamma - \varepsilon) \\
S_3 &= (\gamma - \varepsilon, \gamma - \varepsilon, \gamma + \varepsilon, \gamma + \varepsilon) \\
\end{align}
```

The states are shown in the figure below. The red dots mean "more early customers" (``\gamma + \varepsilon``), the green dots mean "more late customers" (``\gamma - \varepsilon``) and the gray dots mean "no shock" (``\gamma``).
"""

# ╔═╡ 4242531a-e74f-4618-939f-2adf9d6e1db2
states = [
	[:none, :none, :none, :none],
	[:early, :early, :late, :late],
	[:late, :late, :early, :early],
]

# ╔═╡ cb435758-d2a8-4203-9915-971e041d4319
md"""
Select state (the ``i`` in ``S_i``): $(@bind i_state NumberField(1:length(states))).
"""

# ╔═╡ 6f28bcfb-2b56-4548-bbd1-9528876525dd
md"""
👉 **(a)** What is the minimal number of edges that will prevent a bank run in period ``t=1`` in state ``S_1``? Explain briefly.
"""

# ╔═╡ fa50a9b9-cdc4-4d84-ae3b-db039f1609e4
answer_a = md"""
Your answer goes here ...
"""

# ╔═╡ cfa98250-1d4a-43a6-a99f-cf106001f3cb
md"""
👉 **(b)** What is the minimal number of edges that will prevent a bank run in period ``t=1`` in all possible states? Explain and adjust the adjacency matrix `G_minimal` accordingly.
"""

# ╔═╡ 9b0da913-fc5e-42ea-bc5f-e37bd59f2cd2
G_minimal = [
	0 1 1 1;
    0 0 1 1;
	1 0 0 1;
	1 1 0 0
]

# ╔═╡ a7e36f2c-d588-4fc5-a247-4323d646a51b
answer_b = md"""
Your answer goes here ...
"""

# ╔═╡ 77b5b0ea-bfae-406f-8fcf-472165bdcd1d
md"""
👉 **(c)** Assume that your minimal network from **(a)** has _uniform weights_. What is the lower bound ``y_\text{min}`` for that weight that will allow the socially optimal allocation in all states?
"""

# ╔═╡ dbfd2f13-89f4-4932-a778-b2d375d45ac6
answer_c = md"""
Your answer goes here ...
"""

# ╔═╡ 63a8c85f-9cd1-4cf5-9f58-e482494f8d24
md"""
👉 **(d)** What will happen if ``y < y_\text{min}``?
"""

# ╔═╡ 55904c61-4531-4984-b73c-1065a7114772
answer_d = md"""
Your answer goes here ...
"""

# ╔═╡ 264e3358-babf-4bf4-9b57-f436676aa02a
md"""
👉 **(e)** Assume that there is a complete interbank network with a uniform weights to ensure the socially optimal allocation in all states. What would be an alternative state ``S_4`` in which the complete interbank network has a better outcome?
"""

# ╔═╡ 36e610ff-1f42-4d58-b0a7-1bc33bd0d4af
answer_e = md"""
Your answer goes here ...
"""

# ╔═╡ e8637286-ea8b-49c8-b49f-1ab556b83f0c
md"""
## Functions for Task 3
"""

# ╔═╡ a047aeaa-fa54-4cbf-90f4-42d0537b7d06
exX = let
	S = states[i_state]
	n = length(S)

	node_styles = (title = "shock", df = DataFrame(label = string.([:early, :late, :none]), color = ["lightgreen", "tomato", "lightgray"]))

	df = @chain begin
		DataFrame(bank = 1:n, label = string.(S))
		leftjoin(_, node_styles.df, on = :label)
	end

	(; n, color = df.color, node_styles)
end;

# ╔═╡ b1c0d43c-f483-4290-998f-177ce79f41fa
function node_legend(figpos, node_styles, title = "")
	
	elems = [MarkerElement(; color, markersize = 15, marker = :circle) for color ∈ node_styles.color]

	if length(title) == 0
		title_tuple = ()
	else
		title_tuple = (title, )
	end
	
	Legend(
		figpos,
    	elems, node_styles.label, title_tuple...;
		orientation=:horizontal, titleposition=:left, framevisible=false
	)
end

# ╔═╡ 4e9b785f-ad74-4aa8-ad48-89fa8b236939
md"""
# Appendix
"""

# ╔═╡ 1a997e44-f29c-4c55-a953-a9039f096d47
TableOfContents()

# ╔═╡ 78bedfcc-3671-4852-985b-3e1b5aaade5a
md"""
## Helpers
"""

# ╔═╡ 25e84f19-9cd8-43ad-ae6a-d500b8ac74b6
md"""
## Packages
"""

# ╔═╡ 11d6ac4f-e910-4a9f-9ee4-bdd270e9400b
md"""
## HTML Helpers
"""

# ╔═╡ c1b1a22b-8e18-4d19-bb9b-14d2853f0b72
@htl("
<style>
.blurry-text {
	color: transparent;
	text-shadow: 0 0 5px rgba(0,0,0,0.5);
	white-space: nowrap;
}

.blurry-text:hover {
	color: transparent;
	text-shadow: 0 0 0px rgba(0,0,0.5);
	white-space: nowrap;
}

</style>
")

# ╔═╡ cc1ff1e6-2968-4baf-b513-e963ab2ff1b4
blur(text) = @htl("""<span class="blurry-text">$text</span>""")

# ╔═╡ 1d5d8c8a-8d86-426f-bb17-bd2279d91ff1
md"""
* ring vs complete ``(γ)``: $(@bind _γ_ Slider(0:0.1:1.0, show_value = true, default = 0.5))
* shock to bank ``1`` ``(ε)``: $(@bind _ε_ Slider(0.0:0.05:1, show_value = true, default = 0.0))
*  $(blur("number of islands:")) $(@bind n_islands Slider(1:3, show_value = true, default = 1))
"""

# ╔═╡ c65323f5-211d-4a95-aed3-d6129bdd083e
strike(text) = @htl("<s>$text</s>")

# ╔═╡ bafc7db1-ac1d-4314-be83-0b9f6c63b5fc
md"""
### Bank

* external debt (deposits) ``ν``
* external assets (cash) ``c``
* borrowing and lending on _**financial network**_ (see below)
* investment in $(strike("loan to")) _**firm**_ (see below)

### Financial network

* ``n`` banks
* ``y_{ij} > 0 ⟺`` bank ``i`` has debt with bank ``j``
* the network is **_regular_** \
  ``\forall i \sum_{j} y_{ij} = \sum_{j} y_{ji} = \bar y``
"""

# ╔═╡ 5d263432-856b-4e9b-a303-a476222e8963
vertical_space(space=space) = @htl("""
<p style="line-height:$(space)cm;"><br></p>
""")

# ╔═╡ e11a99df-d0f2-4838-b325-473d3043be98
vertical_space()

# ╔═╡ ba39c037-c4e1-465e-a232-e5194caa14ba
vertical_space()

# ╔═╡ 834ff5a0-3bae-44ca-a413-f45d3394b511
vertical_space()

# ╔═╡ 4ae1b6c2-bb63-4ca8-b8ec-057c8d2a371f
vertical_space()

# ╔═╡ aaffd333-3aa1-48ee-b5db-d577bd7da830
vertical_space()

# ╔═╡ 25a01039-25e9-44b0-afd0-c3df37c7598f
md"""
## Theorems, etc
"""

# ╔═╡ 71231141-c2f5-4695-ade0-548a0039f511
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

# ╔═╡ c993c1d8-8823-4db4-bf6e-bf2c21ea3d39
begin
	admonition(kind, title, text) = Markdown.MD(Markdown.Admonition(kind, title, [text]))
	proposition(text; number = nothing, title = nothing) = 
		admonition("correct", theorem_header("Proposition", number, title), text)
	corollary(text; number = nothing, title = nothing) =
		admonition("note", theorem_header("Corollary", number, title), text)
	definition(text; number = nothing, title = nothing) =
		admonition("warning", theorem_header("Definition", number, title), text)
	assumption(text; number = nothing, title = nothing) =
		admonition("danger", theorem_header("Assumption", number, title), text)
#	danger(text, title="Danger")   = admonition("danger",  title, text)
#	correct(text, title="Correct") = admonition("hint", title, text)

end

# ╔═╡ 7b876239-8ddc-4929-ad52-752edb72c0eb
Foldable("Take-away", 
	proposition(md"""
**Small Shocks:** more complete ``\implies`` fewer defaults

**Large shocks:** split up network ``\implies`` fewer defaults $(Foldable("in their words", md"> Not all financial systems, however, are as fragile in the presence of large shocks. In fact, as part ``(ii)`` shows, for small enough values of ``δ``, any ``δ``-connected financial network is strictly more stable and resilient than both the complete and the ring networks. The presence of such _weakly connected_ components in the network guarantees that the losses—rather than being transmitted to all other banks—are borne in part by the distressed bank’s senior creditors."
	))
""",
title = "Acemoglu, Ozdaglar & Tahbaz-Salehi (2015)", )
)

# ╔═╡ f24387ee-c1cf-4ec0-a34e-4b4f33ee9010
md"""
## Packages
"""

# ╔═╡ fbf40f06-f63a-40ca-bbe3-78104d39ee71
md"""
## Line patterns in AlgebraOfGraphics
"""

# ╔═╡ f863066f-270f-4972-8046-7f51cb697ac5
hatch_dict = Dict(
	:/ => Vec2f0(1),
	:\ => Vec2f0(1, -1),
	:- => Vec2f0(1, 0),
	:| => Vec2f0(0, 1),
	:x => [Vec2f0(1), Vec2f0(1, -1)],
	:+ => [Vec2f0(1, 0), Vec2f0(0, 1)]
)

# ╔═╡ d285879a-bdfd-4efa-aa5d-9dacf08a2dc6
hatches = [:\, :/, :-, :|, :x, :+]

# ╔═╡ 70bbb90e-06f1-4b60-a952-275866945c58
line_pattern(hatch; width=1.5, tilesize=(10,10), linecolor=:black, color=:white) = Makie.LinePattern(; direction=hatch_dict[Symbol(hatch)], width, tilesize, linecolor, background_color=color)

# ╔═╡ 3040803d-e95a-40bc-aa72-92c7d158e226
md"""
# FinancialNetworks.jl
"""

# ╔═╡ f13e010a-d394-4e40-9535-c5e2e3e226aa
md"""
## Payments and Equilibrium
"""

# ╔═╡ 7871b479-2aaf-42c1-ad84-42ac17cfc6e1
function repay((; c, ν), (; x, y), firm, ε)
	
	assets(ℓ) = x + c + cashflow₁(firm, ε, ℓ).cf

	if ν + y ≤ assets(0.0)
		out = (; ℓ=0, y_pc=1.0, ν_pc=1.0)
	elseif assets(0.0) < ν + y ≤ assets(1.0)
		ℓ = find_zero(ℓ_ -> ν + y - assets(ℓ_), (0.0, 1.0))
		out = (; ℓ,   y_pc=1.0, ν_pc=1.0)
	else
		ℓ = 1.0
		ass = assets(ℓ)
		if ν ≤ ass < ν + y
			out = (; ℓ, y_pc = (ass - ν)/y, ν_pc = 1.0)
		else # assets < ν
			out = (; ℓ, y_pc = 0.0, ν_pc = ass/ν)
		end
	end		

	(; out..., assets = assets(out.ℓ), cashflow₁(firm, ε, out.ℓ)..., firm_pc = 1 - out.ℓ, x, y = y * out.y_pc)
	
end

# ╔═╡ 7c026a42-1c05-4968-b068-c8561ca5a2db
md"""
## Visualize bank-firm-network
"""

# ╔═╡ 3998e822-d7e0-40f2-a866-71a3d81c4ca3
function add_legend!(figpos; kwargs...)
	dict = [:lightgray => "solvent", :orange => "insolvent", :red => "bankrupt"]
	
	elements = [
		MarkerElement(marker = :circle, strokewidth=1, markersize = 20, color = c) for c in first.(dict)
	]

	Legend(figpos, elements, last.(dict); kwargs...)

	figpos
end

# ╔═╡ 0d207cbe-5b97-4bd7-accf-fc49ce4522e9
md"""
## Visualize bank balance sheet
"""

# ╔═╡ 5104f0a5-a551-4f4c-8f89-2b8f834b3587
function balance_sheet_palette()
	palette_df = DataFrame(
		label = ["external", "interbank", "firm", "liquidated", "shortfall"],
		color = [Makie.wong_colors()[[5, 2, 4, 3, ]]..., Pattern("/")]
	)
	
	palette_df.label .=> palette_df.color
end

# ╔═╡ f63aee78-1209-40dd-9c9d-2699194807d8
function _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))
	receivable = DataFrame([
		(color = "external",   val=c,   lab=L"cash $c$"),
		(color = "interbank",  val=x,   lab=L"IB deposit $x$"),
		(color = "firm",       val=div, lab=L"dividend $δ$"),
		(color = "firm",       val=ill, lab="illiquid"),
		(color = "liquidated", val=rec, lab=L"recovered $ρ$"),
	])

	payable = DataFrame([
		(color = "external",  val=ν_paid,    lab=L"deposits $ν$"),
		(color = "interbank", val=y_paid,    lab=L"IB debt $y$"),
		(color = "shortfall", val=shortfall, lab=""),
	])

	nt = (; receivable, payable)

	vcat(nt..., source = :side => collect(string.(keys(nt))))
end

# ╔═╡ 6a2b3b9c-df52-41c3-908b-5dbf052ad107
function balance_sheet_df_new(bank, firm, (; x, y); ε_cash=get(bank, :ε, 0.0))
	c₀ = bank.c
	c = max(c₀ - ε_cash, 0.0)
	bank = (; c, bank.ν)
	
	(; y_pc, ν_pc, ℓ, rec, div) = 
		repay(bank, (; x, y), firm, ε)

	(; a, A) = firm
	(; ν, c) = bank
	y_paid = y_pc * y
	ν_paid = ν_pc * ν
	shortfall = (1-ν_pc) * ν + (1-y_pc) * y
	ill = (1-ℓ) * A
	
	df = _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))

	bs_max = max((A + a) + x + c₀, y + ν)
	
	(; df, bs_max, ℓ, shortfall)
end

# ╔═╡ 1b9c8f2b-8011-46f1-b46d-479031eb9ac3
function _visualize_simple_balance_sheet_(bank, firm, (; x, y); show_illiquid=true, fontsize=15)

	(; df, bs_max, ℓ, shortfall) = balance_sheet_df_new(bank, firm, (; x, y))

	@transform!(df, :stack = :lab == "illiquid" ? :lab : :color)
	
	if !show_illiquid || ℓ ≈ 1
		@subset!(df, :lab ≠ "illiquid")
	end
	if ℓ * firm.ζ ≈ 0
		@subset!(df, :color ≠ "liquidated")
	end
	@subset!(df, :val > 0)
	@transform!(df, @subset(:val < 0.3), :lab = "")

	df.color = categorical(df.color)
	
	ordered = @chain begin
		["external", "interbank", "firm", "liquidated", "shortfall"]
		filter(∈(unique(df.color)), _)
	end
	levels!(df.color, ordered)

	df.stack = categorical(df.stack)
	ordered = @chain begin
		["external", "interbank", "firm", "liquidated", "illiquid", "shortfall"]
		filter(∈(unique(df.stack)), _)
	end
	levels!(df.stack, ordered)
	
	plt = data(df) * mapping(
		:side => sorter("receivable", "payable") => "", 
		:val => "",
		stack=:stack, color=:color => "",
		bar_labels=:lab => verbatim
	) * visual(BarPlot, flip_labels_at = 0, label_size=fontsize, strokewidth=0.5)

	(; plt, bs_max)
end

# ╔═╡ ebf41ae2-3985-439b-b193-eabfab701d16
function visualize_simple_balance_sheet(args...; figure=(;), legend=(;), kwargs...)
	(; plt, bs_max) = _visualize_simple_balance_sheet_(args...; kwargs...)

	if isempty(legend)
		legend = (; position = :top, titleposition=:left, framevisible=false)
	end
	
	draw(
		plt;
		axis = (limits = (nothing, nothing, 0, 1.05 * bs_max),),
		figure,
		palettes = (; color=balance_sheet_palette()),
		legend,
	)
end

# ╔═╡ b7d74679-904c-44c4-bedd-89f1b68a5e42
function visualize_simple_balance_sheet!(figpos, legpos, args...; legend=(;), kwargs...)
	(; plt, bs_max) = _visualize_simple_balance_sheet_(args...; kwargs...)

	fg = draw!(
		figpos,
		plt,
		palettes = (; color=balance_sheet_palette())
	)
	
	legend!(legpos, fg; legend...)
end

# ╔═╡ ab1cf7ab-80fb-4423-a924-1d6e24e9c9bc
function balance_sheet_df((; ν, c), (; a, A, ζ), (; ℓ, y_pc, ν_pc, x, y, div, rec))

	y_paid = y
	y = y / y_pc
	ν_paid = ν_pc * ν
	shortfall = (1-ν_pc) * ν + (1-y_pc) * y
	ill = (1-ℓ) * A
	
	df = _balance_sheet_df_((; c, x, div, ill, rec, ν_paid, y_paid, shortfall))

	@subset!(df, :lab ≠ "illiquid")
	
	bs_max = max(a + x + c, y + ν) * 1.05

	(; df, bs_max = isfinite(bs_max) ? bs_max : zero(bs_max))
end

# ╔═╡ 0ff81eb8-d843-49a4-af93-ec2414797e87
function visualize_balance_sheets!(figpos, bank_df, banks, firm)
	
	bank_dfs = [
		balance_sheet_df(banks[i], firm, bank_df[i,:]) for i ∈ 1:length(banks)
	] |> StructArray

	combined_bank_df = vcat(
		bank_dfs.df..., source = :bank
	)

	@transform!(combined_bank_df, :bank = "bank $(:bank)")
	
	bs_max = maximum(bank_dfs.bs_max)

	plt = @chain combined_bank_df begin
		data(_) * mapping(
			:side => sorter("receivable", "payable") => "",
			:val => "",
			color = :color => "", stack = :color,
			layout = :bank
		) * visual(BarPlot, strokewidth=0.5)
	end

	fg = draw!(figpos[2,1], plt, 
		axis = (limits = (nothing, nothing, 0, 1.05 * bs_max),),
		palettes=(color=balance_sheet_palette(),)
	)
	legend!(figpos[1,1], fg, orientation=:horizontal, titleposition=:left, framevisible=false)

	nothing
end

# ╔═╡ 0f03f537-f589-4abd-9587-0bb18835d9b9
function visualize_balance_sheets(out, banks, firms, shares)
	fig = Figure()
	visualize_balance_sheets!(fig[1,1], out, banks, firms, shares)
	fig
end

# ╔═╡ ffe50130-a9cb-4ba9-a861-c247bf688873
md"""
# NetworkUtils.jl
"""

# ╔═╡ d3e5b8f2-51b1-4ba4-97b4-2be156a74643
md"""
## Regular interbank market
"""

# ╔═╡ d5c128e0-6371-4ad2-8bfc-c17faadc520b
weighted_adjacency_matrix(graph) = Graphs.weights(graph) .* (adjacency_matrix(graph) .> 0)

# ╔═╡ 38a1f960-7178-4968-89e4-6b659b64baa2
function payments_received((; promises, pc_served))
	Y = weighted_adjacency_matrix(promises)
	vec(pc_served' * Y)
end

# ╔═╡ f689810f-8f43-4368-a822-0ee4f3015271
function payments_promised((; promises))
	Y = weighted_adjacency_matrix(promises)
	dropdims(sum(Y, dims=2), dims=2)
end

# ╔═╡ 601b1453-20bd-4e93-abb7-4fad4f5adf8c
function iterate_payments(banks, promises, pc_served, firm, εs_firm, ℓs = zeros(length(banks)))
	ys = payments_promised((; promises))
	
	x_new = copy(y)
	out = map(enumerate(banks)) do (i, bank)
		# compute repayment
		xs = payments_received((; promises, pc_served))
		internal = (; x=xs[i], y=ys[i])
		
		rpy = repay(bank, internal, firm, εs_firm[i])
		(; y_pc, ℓ) = rpy
	
		ℓs[i] = ℓ
		pc_served[i] = y_pc
		# return bank's choices
		(; bank = i, rpy...)
	end

	(; pc_served, out, ℓs)
end

# ╔═╡ 9cdc97e5-67f8-4e71-97d8-9cda5e9d7bd8
function equilibrium(banks, promises, firm, ε_firms; maxit = 100)
	n_banks = length(banks)
	
	pc_served = ones(n_banks)
	ℓs = zeros(n_banks)

	for it ∈ 1:maxit
		ℓs_old = copy(ℓs)
		pc_served_old = copy(pc_served)
		
		(; pc_served, out, ℓs) = iterate_payments(banks, promises, pc_served, firm, ε_firms, ℓs)
		converged = pc_served_old ≈ pc_served && ℓs_old ≈ ℓs
		
		if converged || it == maxit
			bank_df = DataFrame(out)
			
			return (; bank_df, it, success = it != maxit, banks)
		end
	end
	
end

# ╔═╡ 5bde1113-5297-4e93-b052-7a0f93f4ea84
function is_regular(Y)
	vec(sum(Y, dims=1)) ≈ vec(sum(Y, dims=2))
end

# ╔═╡ 73456a8f-8534-4b55-8af6-ff7952bd3a3a
function is_regular(Y, ȳ)
	is_regular(Y) && mean(sum(Y, dims=1)) ≈ ȳ
end

# ╔═╡ e6def365-4438-4591-acbd-60d9de466b0a
initial_network(interbank_market) = adjacency_matrix(interbank_market.network)

# ╔═╡ 055a811c-bc4a-4959-88f1-bc71e8749313
updated_network(interbank_market) = interbank_market.payments

# ╔═╡ 3c211dad-73b0-4715-b066-e10005f8f3a7
function complete_network(n, ȳ)
	Y = fill(ȳ/(n-1), n, n)
		
	for i in 1:n
		Y[i,i] = 0.0
	end
	
	Y
end

# ╔═╡ e2045a53-21d3-4c7c-ad95-3fc8d2444821
function CompleteNetwork(n, ȳ)
	SimpleWeightedDiGraph(complete_network(n, ȳ))
end

# ╔═╡ 2a5f02e7-fec1-41b3-b8b4-b33b1cc1232c
function ring_network(n, ȳ)
	Y = zeros(n, n)
	
	for i in 1:(n-1)
		Y[i, i+1] = ȳ
	end
	Y[n, 1] = ȳ
	
	Y
end

# ╔═╡ ee10e134-84be-41f8-9e7e-c1cd68227977
function RingNetwork(n, ȳ)
	SimpleWeightedDiGraph(ring_network(n, ȳ))
end

# ╔═╡ 47118d1b-8471-4288-85fe-d3e94667dc96
function γ_network(n, ȳ, γ)
	Y = γ * ring_network(n, ȳ) + (1-γ) * complete_network(n, ȳ)
end

# ╔═╡ dffce961-6844-4a44-beea-ab085c2f9f3f
function γNetwork(n, ȳ, γ)
	SimpleWeightedDiGraph(γ_network(n, ȳ, γ))
end

# ╔═╡ f708e6c0-cfac-4b4d-a3ed-69b98883294a
function island_network(n_islands, n_banks_per_island, ȳ, γ)
	blocks = (γ_network(n_banks_per_island, ȳ, γ) for _ in 1:n_islands)
	
	cat(blocks...,dims=(1,2))
end

# ╔═╡ d6cb95c1-a075-4544-9031-58aef65c7577
function island_network(n_banks::AbstractVector, ȳ, γ)
	blocks = (γ_network(n, ȳ, γ) for n ∈ n_banks)
	
	cat(blocks...,dims=(1,2))
end

# ╔═╡ 1bb841e0-ddd2-4571-83a5-d929e0a8a69c
function IslandNetwork(n_islands, n_banks_per_island, ȳ; γ=0.0)
	Y = island_network(n_islands, n_banks_per_island, ȳ, γ)	
	SimpleWeightedDiGraph(Y)
end

# ╔═╡ 7fbcfbde-0b5e-4bf2-9eda-9b15a4dd6bec
function IslandNetwork(n_banks::AbstractVector, ȳ; γ=0.0)
	Y = island_network(n_banks, ȳ, γ)	
	SimpleWeightedDiGraph(Y)
end

# ╔═╡ 8a2f3e4d-9c61-4a21-a229-58f731964181
initial_analysis = let
	n_banks = 6
	n_firms = n_banks

	ν = 2.3
	c = 2.4

	banks = [(; ν, c = max(c - (i==1)*_ε_, 0)) for i ∈ 1:n_banks]

	firm = (; ζ=0.0, a=0.0, A=0.5)
	εs = zeros(n_banks)

	ȳ = 1.2
	promises = IslandNetwork(n_islands, n_banks ÷ n_islands, ȳ; γ=_γ_)
	
	(; bank_df) = equilibrium(banks, promises, firm, εs)

	(; promises, banks, firm, bank_df)
end;

# ╔═╡ a8a8a197-d54d-4c67-b7ce-19bdc8b64401
transmission_analysis = let
	n_islands = 1
	n_banks = 4
	
	banks = [(; ν, c = max(c - (i==1)*ε, 0)) for i ∈ 1:n_banks]

	firm = (; ζ, a, A)
	εs = zeros(n_banks)

	ȳ = y
	IM = IslandNetwork(n_islands, n_banks ÷ n_islands, ȳ; γ=γ2)

	(; bank_df) = equilibrium(banks, IM, firm, εs)

	(; banks, firm, IM, bank_df)
end

# ╔═╡ ff77f7b4-52d1-4fc8-abfc-e623f7bcd423
md"""
## Layout
"""

# ╔═╡ c7c8581c-e82c-428e-a599-3c003cc0151c
begin
	function points_on_circle(n; c = Point(0, 0), r = 1, start = 0)
		if n == 0
			Point2[]
		else
			x = range(0, 2, n + 1) .+ start # .- 0.5
			x = x[begin:end-1]
		
			Point2.(r * sinpi.(x), r * cospi.(x)) .+ Ref(c)
		end
	end

	function nested_circles(inner; r=1.5, start = 0)
		function (g)
			n_nodes = nv(g)
			outer = n_nodes - inner
			[
				points_on_circle(inner);
				points_on_circle(outer; r, start)
			]
		end
	end
	function componentwise_circle(g::AbstractGraph; kwargs...)
		components = connected_components(g)
		componentwise_circle(components; kwargs...)
	end
			
	function componentwise_circle(components; scale=2.5)
		nodes = Int[]
		node_locs = Point2f[]

		N = length(components)
		ncol = ceil(Int, sqrt(N))

		for (i, component) in enumerate(components)
			(row, col) = fldmod1(i, ncol)
			n = length(component)
			append!(nodes, component)
			append!(node_locs, points_on_circle(n, c = Point(scale * (col-1), scale * (1 - row))))
		end

		node_locs[sortperm(nodes)]
	end
end

# ╔═╡ 69698f6a-110d-49ea-bd57-16af0de90844
md"""
# NamedGraphPlot
"""

# ╔═╡ 2c22e86e-1139-4f66-8992-093d38f4c6cb
md"""
## Plotting helpers
"""

# ╔═╡ 60293fac-4d20-409b-bfd2-d5283f189320
attr(; kwargs...) = (; arrow_size = 10, elabels_distance = 15, edge_width=1, nlabels_distance = 10, kwargs...)

# ╔═╡ bef59a5b-6952-42aa-b700-4ad81d848fe4
minimal(; extend_limits=0.05, hidespines=true, kwargs...) = (; 
	xgridvisible=false, xticksvisible=false, xticklabelsvisible=false,
	ygridvisible=false, yticksvisible=false, yticklabelsvisible=false, 
	leftspinevisible=!hidespines, rightspinevisible=!hidespines, topspinevisible=!hidespines, bottomspinevisible=!hidespines,
	xautolimitmargin = (extend_limits, extend_limits),
	yautolimitmargin = (extend_limits, extend_limits),
	kwargs...
)

# ╔═╡ ce563ab5-6324-4a86-be61-7a107ff0e3b3
figure(xscale=1, yscale=xscale) = (; resolution = (xscale * 300, yscale * 300))

# ╔═╡ 4ccb52c2-a177-4b13-8ce5-2af6f286c737
md"""
## Named GraphPlot
"""

# ╔═╡ 4a1fd4a3-71a0-4947-8b75-cc2100fb636b
function text_bbox(textstring::AbstractString, fontsize::Union{AbstractVector, Number}, font, align, rotation, justification, lineheight)
    glyph_collection = Makie.layout_text(
            textstring, fontsize,
            font, nothing, align, rotation, justification, lineheight,
            RGBAf(0,0,0,0), RGBAf(0,0,0,0), 0f0, 100
        )

    return Rect2f(Makie.boundingbox(glyph_collection, Point3f(0), Makie.to_rotation(rotation)))
end

# ╔═╡ 3496c181-c4c3-4b1b-a5e5-83df27182c99
md"""
## `LaTeXStrings` as `bar_labels`
"""

# ╔═╡ e8d198fa-11db-43b2-846e-311abe9f59aa
function get_xshift(lb, ub, align; default=0.5f0)
    if align isa Symbol
        align = align == :left   ? 0.0f0 :
                align == :center ? 0.5f0 :
                align == :right  ? 1.0f0 : default
    end
    lb * (1-align) + ub * align |> Float32
end

# ╔═╡ 2711b045-b879-4225-afae-13328a4e14b3
function get_yshift(lb, ub, align; default=0.5f0)
    if align isa Symbol
        align = align == :bottom ? 0.0f0 :
                align == :center ? 0.5f0 :
                align == :top    ? 1.0f0 : default
    end
    lb * (1-align) + ub * align |> Float32
end

# ╔═╡ c4032c7a-73f4-4342-a5a7-19fd14017402
begin
fonts = (; regular = CairoMakie.Makie.MathTeXEngine.texfont(), bold = CairoMakie.Makie.MathTeXEngine.texfont())
set_theme!(; fonts)
	
using Makie: automatic, compute_x_and_width, barplot_labels, bar_rectangle, generate_tex_elements, TeXChar, FreeTypeAbstraction, GlyphExtent, height_insensitive_boundingbox_with_advance, origin, GlyphCollection, stack_grouped_from_to, ATTRIBUTES, theme, bar_label_formatter
	
function Makie.plot!(p::BarPlot)

    labels = Observable(Tuple{Union{String,LaTeXStrings.LaTeXString}, Point2f}[])
    label_aligns = Observable(Vec2f[])
    label_offsets = Observable(Vec2f[])
    label_colors = Observable(RGBAf[])
    function calculate_bars(xy, fillto, offset, width, dodge, n_dodge, gap, dodge_gap, stack,
                            dir, bar_labels, flip_labels_at, label_color, color_over_background,
                            color_over_bar, label_formatter, label_offset)

        in_y_direction = get((y=true, x=false), dir) do
            error("Invalid direction $dir. Options are :x and :y.")
        end

        x = first.(xy)
        y = last.(xy)

        # by default, `width` is `minimum(diff(sort(unique(x)))`
        if width === automatic
            x_unique = unique(filter(isfinite, x))
            x_diffs = diff(sort(x_unique))
            width = isempty(x_diffs) ? 1.0 : minimum(x_diffs)
        end

        # compute width of bars and x̂ (horizontal position after dodging)
        x̂, barwidth = compute_x_and_width(x, width, gap, dodge, n_dodge, dodge_gap)

        # --------------------------------
        # ----------- Stacking -----------
        # --------------------------------

        if stack === automatic
            if fillto === automatic
                fillto = offset
            end
        elseif eltype(stack) <: Integer
            fillto === automatic || @warn "Ignore keyword fillto when keyword stack is provided"
            if !iszero(offset)
                @warn "Ignore keyword offset when keyword stack is provided"
                offset = 0.0
            end
            i_stack = stack

            from, to = stack_grouped_from_to(i_stack, y, (x = x̂,))
            y, fillto = to, from
        else
            ArgumentError("The keyword argument `stack` currently supports only `AbstractVector{<: Integer}`") |> throw
        end

        # --------------------------------
        # ----------- Labels -------------
        # --------------------------------

        if !isnothing(bar_labels)
            oback = color_over_background === automatic ? label_color : color_over_background
            obar = color_over_bar === automatic ? label_color : color_over_bar
            label_args = barplot_labels(x̂, y, bar_labels, in_y_direction,
                                        flip_labels_at, to_color(oback), to_color(obar),
                                        label_formatter, label_offset)
            labels[], label_aligns[], label_offsets[], label_colors[] = label_args
        end

        return bar_rectangle.(x̂, y .+ offset, barwidth, fillto, in_y_direction)
    end

    bars = lift(calculate_bars, p[1], p.fillto, p.offset, p.width, p.dodge, p.n_dodge, p.gap,
                p.dodge_gap, p.stack, p.direction, p.bar_labels, p.flip_labels_at,
                p.label_color, p.color_over_background, p.color_over_bar, p.label_formatter, p.label_offset)

	kwargs_df = DataFrame(; 
		color=p.color[], colormap = p.colormap[], colorrange = p.colorrange[],
        strokewidth = p.strokewidth[], strokecolor = p.strokecolor[],
		visible = p.visible[],
        inspectable = p.inspectable[], transparency = p.transparency[],
        highclip = p.highclip[], lowclip = p.lowclip[], nan_color = p.nan_color[]
	)
	map(zip(bars[], eachrow(kwargs_df))) do (bar, kwargs)
		poly!(p, bar; kwargs...)
	end

    if !isnothing(p.bar_labels[])
		#@info p.label_size[]
        text!(p, labels; align=label_aligns, offset=label_offsets, color=label_colors, font=p.label_font, fontsize=p.label_size,
		rotation=p.label_rotation)
    end
end

function Makie.texelems_and_glyph_collection(str::LaTeXString, fontscale_px, halign, valign,
        rotation, color, strokecolor, strokewidth, word_wrap_width)

    rot = convert_attribute(rotation, key"rotation"())

    all_els = generate_tex_elements(str)
    els = filter(x -> x[1] isa TeXChar, all_els)

    # hacky, but attr per char needs to be fixed
    fs = Vec2f(first(fontscale_px))

    scales_2d = [Vec2f(x[3] * Vec2f(fs)) for x in els]

    texchars = [x[1] for x in els]
    glyphindices = [FreeTypeAbstraction.glyph_index(texchar) for texchar in texchars]
    fonts = [texchar.font for texchar in texchars]
    extents = GlyphExtent.(texchars)

    bboxes = map(extents, scales_2d) do ext, scale
        unscaled_hi_bb = height_insensitive_boundingbox_with_advance(ext)
        return Rect2f(
            origin(unscaled_hi_bb) * scale,
            widths(unscaled_hi_bb) * scale
        )
    end

    basepositions = [to_ndim(Vec3f, fs, 0) .* to_ndim(Point3f, x[2], 0) for x in els]

    if word_wrap_width > 0
        last_space_idx = 0
        last_newline_idx = 1
        newline_offset = Point3f(basepositions[1][1], 0f0, 0)

        for i in eachindex(texchars)
            basepositions[i] -= newline_offset
            if texchars[i].represented_char == ' ' || i == length(texchars)
                right_pos = basepositions[i][1] + width(bboxes[i])
                if last_space_idx != 0 && right_pos > word_wrap_width
                    section_offset = basepositions[last_space_idx + 1][1]
                    lineheight = maximum((height(bb) for bb in bboxes[last_newline_idx:last_space_idx]))
                    last_newline_idx = last_space_idx+1
                    newline_offset += Point3f(section_offset, lineheight, 0)

                    # TODO: newlines don't really need to represented at all?
                    # chars[last_space_idx] = '\n'
                    for j in last_space_idx+1:i
                        basepositions[j] -= Point3f(section_offset, lineheight, 0)
                    end
                end
                last_space_idx = i
            elseif texchars[i].represented_char == '\n'
                last_space_idx = 0
            end
        end
    end

    bb = isempty(bboxes) ? BBox(0, 0, 0, 0) : begin
        mapreduce(union, zip(bboxes, basepositions)) do (b, pos)
            Rect2f(Rect3f(b) + pos)
        end
    end

    xshift = get_xshift(minimum(bb)[1], maximum(bb)[1], halign)
    yshift = get_yshift(minimum(bb)[2], maximum(bb)[2], valign, default=0f0)
    
    shift = Vec3f(xshift, yshift, 0)
    positions = basepositions .- Ref(shift)
    positions .= Ref(rot) .* positions

    pre_align_gl = GlyphCollection(
        glyphindices,
        fonts,
        Point3f.(positions),
        extents,
        scales_2d,
        rot,
        color,
        strokecolor,
        strokewidth
    )

    all_els, pre_align_gl, Point2f(xshift, yshift)
end

end

# ╔═╡ efbc1a42-ecc2-4907-a981-bd1d29ca0803
balance_sheet = let
	c = 0.7
	ν = 1.2
	x = 1.0
	a = 0.7
	y = 1.0
	A = 1.0
	ζ = recovery_rate
	shortfall = max(y + ν - (c - ε_cash + x + a), 0)

	figure = (; fonts, resolution = (400, 350), fontsize=15)

	firm = (; a, A, ζ)
	visualize_simple_balance_sheet((; c, ν, ε=ε_cash), firm, (; x, y); figure, show_illiquid) |> as_svg
end	

# ╔═╡ 8e6bdcfd-3cb0-4f6b-8dcf-f99f4afe87c2
aside(md"$(balance_sheet)")

# ╔═╡ 6fddaaa4-688d-4f40-a1b3-04bf62024955
aside(md"$balance_sheet")

# ╔═╡ 1f23193e-0bca-4bee-82b1-453100721ab6
begin
	@recipe(NamedGraphPlot, graph) do scene
	    scatter_theme = default_theme(scene, Scatter)
	    lineseg_theme = default_theme(scene, LineSegments)
	    labels_theme = default_theme(scene, Makie.Text)
	    Attributes(
			graphplot_attr = Attributes(),
	        # node attributes (Scatter)
	        node_color = :lightgray,
			node_strokewidth = 0.5,
	        node_size = automatic,
	        node_marker = :circle,
			node_aspect = :regular, # need for a better name
			node_labels = automatic,
			node_font=labels_theme.fonts.regular,
			node_fontsize=labels_theme.fontsize
	    )
	end
	
	function Makie.plot!(ngp::NamedGraphPlot)
		# Extract attributes from plot object
		(; graph, graphplot_attr,
		   node_color, node_strokewidth, node_size, node_aspect,
		   node_marker, node_labels, node_font, node_fontsize) = ngp

		# Compute marker sizes
		out = lift(
				 node_labels, graph, node_size, node_fontsize, node_font, node_aspect
			) do node_labels, graph, node_size, node_fontsize, node_font, node_aspect
			
			if node_labels === automatic
				node_labels = string.(1:nv(graph))
			end
				
			# Extract text sizes and layout accordingly
		    label_sizes = [widths(text_bbox(
				string(label), 
				node_fontsize,
				node_font,
				(:center, :center),
				0f0, 0, 0
			)) for label in node_labels]
		
			# if regular: use circle/square, otherwise elipse/rectangle
			if node_aspect === :regular
				label_sizes = [max(pt...) for pt in label_sizes]
			end
				
			# We can specify markersize as a Vec2f, which is the eltype of label_sizes.
		    # Thus, we can explicitly cause the node drawing to be large enough
			# to accomodate the marker size.
			if node_size === automatic
				node_size = map(x -> 1.7x .+ 0node_fontsize, label_sizes)
			end
	
			(; node_labels, label_sizes, node_size)
		end
	
		@lift (; node_labels, label_sizes, node_size) = $out
		
		node_attr = (; 
	        marker = node_marker,
			color = node_color,
	        strokewidth = node_strokewidth,
			markersize = node_size, 
	        markerspace = :pixel,
	    )
	
		if hasproperty(graphplot_attr, :node_attr)
			@warn "Ignoring :node_attr"
			delete!(graphplot_attr, :node_attr)
		end
	
		gp = graphplot!(ngp, graph; graphplot_attr..., node_attr)
	
		positions = @lift($(gp[:layout])($graph))
		
		text!(ngp, positions; text=node_labels, font=node_font, fontsize=node_fontsize, align = (:center, :center))
		
	    return ngp
	end
end

# ╔═╡ c99e52e2-6711-4fb6-bcc0-8e4f378ed479
out_T1 = let
	y = 2.1
	ỹ = IslandNetwork([3, 3, 5], y)
	ŷ = IslandNetwork([6, 2, 1, 1, 1], y)

	@assert nv(ỹ) == nv(ŷ)
	n = nv(ỹ)

	layout = (_) -> componentwise_circle([1:6, 6 .+ (1:5)])

	fig = Figure(; figure(2.5, 1.0)...)
	ax1 = Axis(fig[1,1]; minimal(hidespines=false, title = L"interbank market $\tilde{y}$")...)
	ax2 = Axis(fig[1,2]; minimal(hidespines=false, title = L"interbank market $\hat{y}$")...)

	namedgraphplot!(ax1, ỹ, graphplot_attr = (; layout))
	namedgraphplot!(ax2, ŷ, graphplot_attr = (; layout))

	(; ỹ, ŷ, fig, n, y)
end; out_T1.fig |> as_svg

# ╔═╡ 15f45669-516b-4f3f-9ec1-f9e2c1d2e71a
@markdown("""
Consider the interbank networks ``\\tilde y`` and ``\\hat y`` of $(out_T1.n) banks as depicted above. For all non-isolated banks the sum of interbank liabilities equal the sum of interbank claims (``y = $(out_T1.y)``).
""")

# ╔═╡ d7111001-f632-4d0d-a2c7-7bbfd67bf87d
md"""
For this exercise you can use the tool below, to simulate the payment equilibrium for a given interbank market, shock size, and the bank that is hit by the shock.

* Which bank is hit? ``i`` $(@bind i_T1 Slider(1:out_T1.n, default = 1, show_value = true))
* Size of the shock ``\varepsilon``  $(@bind ε_T1 Slider(0.0:0.1:3.0, show_value = true, default = 1.0))
* Select ``\tilde y`` or ``\hat y`` $(@bind IM_T1 Select([out_T1.ỹ => "ỹ", out_T1.ŷ => "ŷ"]))
"""

# ╔═╡ c7b99d3c-5d32-45e6-84fa-8a6513e6beb9
out_T2 = let
	ȳ = 2.1
	IM1 = IslandNetwork([3, 2], ȳ; γ=0.0)
	IM2 = IslandNetwork([3, 2], ȳ; γ=1.0)

	n1 = nv(IM1)
	n2 = nv(IM2)
	if n1 == n2
		n = n1
	else
		n = (n1, n2)
	end
	layout = Shell()

	fig = Figure(; figure(2, 1)...)
	ax1 = Axis(fig[1,1]; minimal(hidespines=false, title = L"interbank market $\tilde{y}$")...)
	ax2 = Axis(fig[1,2]; minimal(hidespines=false, title = L"interbank market $\hat{y}$")...)

	namedgraphplot!(ax1, IM1; graphplot_attr=(; layout))
	namedgraphplot!(ax2, IM2; graphplot_attr=(; layout))

	(; IM1, IM2, fig, n)# |> as_svg
end; out_T2.fig |> as_svg

# ╔═╡ 0fb4d187-f03a-435b-b9fc-188925e058f1
md"""
If you have understood the mechanics of the model, you should be able to solve these tasks without simulations. You can use the given tool to verify your answer.

* Which bank is hit? ``i`` $(@bind i_bank_T2 Slider(1:out_T2.n, default = 1, show_value = true))
* Size of the shock ``\varepsilon``  $(@bind ε_T2 Slider(0.0:0.1:3.0, show_value = true, default = 1.0))
* Select ``\tilde y`` or ``\hat y`` $(@bind IM_T2 Select([out_T2.IM1 => "ỹ", out_T2.IM2 => "ŷ"]))
"""

# ╔═╡ 2fe4c931-d4b2-4b4d-8634-73573125cfb5
let
	g = SimpleDiGraph(G_minimal)

	fig, ax, _ = namedgraphplot(g;
		node_color = exX.color,
		graphplot_attr = (; layout = Shell()),
		figure = figure(),
		axis = minimal(title = latexstring("interbank network in state \$ S_$i_state \$"), extend_limits=0.1)
	)

	(; node_styles) = exX
	if !ismissing(node_styles)
		(title, df) = node_styles
		node_legend(fig[end+1,1], df, title)
	end

	fig |> as_svg
end

# ╔═╡ 222665dc-397b-4e58-b3ee-935b115cf13d
function visualize_bank_firm_network!(ax, IM, bank_df; r = 1.4, start = Makie.automatic, layout=Makie.automatic, show_firms=true, kwargs...)

	A = adjacency_matrix(IM)
	n_banks = nv(IM)

	if show_firms
		n_firms = n_banks
		shares = I(n_firms)
		big_A = [A       shares';
	    	     0 * shares 0 * I]
		
	else
		big_A = A
		
	end
	
	g = big_A |> SimpleWeightedDiGraph

	n_banks = size(A, 1)

	edge_attr_df = map(edges(g)) do (; src, dst, weight)
		if src ≤ n_banks && dst ≤ n_banks
			linewidth = weight
			linestyle = Makie.automatic
			arrow_size = 7
		elseif dst > n_banks
			linewidth = 0.5
			arrow_size = 0.0
			linestyle = "--"
		end
		(; linewidth, arrow_size, linestyle)
	end |> DataFrame

	arrow_attr = (; markersize = edge_attr_df.arrow_size)
	edge_attr = (; edge_attr_df.linewidth, edge_attr_df.linestyle)

	node_color = ifelse.(bank_df.y_pc .< 1.0, :red, ifelse.(bank_df.ℓ .> 0.0, :orange, :lightgray))
	nlabels = string.(1:n_banks)
	node_marker = fill(:circle, n_banks)
	
	if show_firms
		nlabels = [nlabels; ["F$i" for i ∈ 1:n_firms]]
		node_marker = [node_marker; fill(:rect, n_firms)]
		
		node_color = [
			node_color;
			fill(colorant"limegreen", n_firms)
		]
	end
	
	start = start === Makie.automatic ? 1/n_banks : start
	layout = layout === Makie.automatic ? nested_circles(n_banks; start, r) : layout

	graphplot_attr = (; 
		layout, edge_attr, arrow_attr,
		edge_plottype = :beziersegments,
		kwargs...
	)
	namedgraphplot!(ax, g;
		node_labels=nlabels,
		node_marker,
		node_color,
		graphplot_attr
	)

	nothing
end

# ╔═╡ e2041c57-0e3d-4dad-9bab-d70434f18509
let
	(; IM, firm, banks, bank_df) = transmission_analysis

	fig = Figure(resolution = (800, 400))
		
	visualize_bank_firm_network!(Axis(fig[1,1]; minimal(extend_limits=0.1)...), IM, bank_df; start = 1/8, show_firms=false)

	visualize_balance_sheets!(fig[1,2:3], bank_df, banks, firm)

	fig |> as_svg
end

# ╔═╡ 7b70a862-faf4-4c42-917c-238718c43708
function visualize_bank_firm_network(IM, bank_df; figure = figure(), add_legend=false, hidespines=true, kwargs...)
	fig = Figure(; figure...)
	ax = Axis(fig[1,1]; minimal(; hidespines, extend_limits=0.1)...)
	visualize_bank_firm_network!(ax, IM, bank_df; kwargs...)

	if add_legend
		add_legend!(fig[0,:], orientation=:horizontal, framevisible=false)
	end
	fig # |> as_svg
end

# ╔═╡ 073982c7-6333-43f6-866a-91a49f8ba7eb
fig = let
	(; promises, bank_df) = initial_analysis
	
	visualize_bank_firm_network(promises, bank_df; figure=figure(), start = 1/6, show_firms=false)
end;

# ╔═╡ 1d13564d-7256-4cf1-919f-6d8298cc476b
md"""
1. _How do shocks propagate_ on a financial network?
2. _How should a network look_ so that shocks don't spread? $(fig |> aside_figure)
"""

# ╔═╡ 73ba4210-d8f6-4a74-bf4d-d7bc0902bb5e
md"""
$(aside_figure(fig))

## Plan for the lecture

> **_Goal_**: Understand model and key results of Acemoglu, Ozdaglar & Tahbaz-Salehi (2015)

* I. Model setup
* II. **insolvency** and **bankruptcy** in the payment equilibrium
* III. **financial contagion**
* IV. **stability** and **resilience** of financial networks
  * more interbank lending leads to higher fragility
  * densely connected networks are **robust, yet fragile**
  * with **big shocks**, we want to have **disconnected components**
* V. Empirical relevance, Outlook

"""

# ╔═╡ 7970090f-eb3a-488e-b225-ce4198494f1d
firm_bank_aside = let
	(; IM, firm, banks, bank_df) = transmission_analysis
		
	visualize_bank_firm_network(IM, bank_df; show_firms=false, start = 1/8) |> as_svg |> x -> md"$(x)" |> aside
end

# ╔═╡ 0d18cdf0-441e-4ca9-98e3-50bc3efa837f
let
	i = i_T1
	IM = IM_T1
	ε = ε_T1
	
	n_banks = nv(IM)

	ν = 3.0
	c = 0.0
	ζ = 0.1
	A = 3.5
	a = 3.0

	
	shares = I(n_banks)
	
	banks = [(; ν, c = c) for i ∈ 1:n_banks]

	firm = (; ζ, a, A)
	εs = zeros(n_banks)
	εs[i] = min(ε, a)
	
	(; bank_df) = equilibrium(banks, IM, firm, εs)

	layout = (_) -> componentwise_circle([1:6, 6 .+ (1:5)])
	
	visualize_bank_firm_network(IM, bank_df; figure=figure(1.25, 1.0), hidespines=false, start = 1/6, layout, add_legend=true, show_firms=false) |> as_svg
end

# ╔═╡ a9d27019-72b7-4257-b72a-12952b516db9
let
	i = i_bank_T2
	#_ε4 = 0.4
	IM = IM_T2
	n_banks = nv(IM)

	ν = 3.0
	c = 0.0
	ζ = 0.0
	A = 0.0
	a = 3.5
	ε = ε_T2
	
	shares = I(n_banks)
	
	banks = [(; ν, c) for i ∈ 1:n_banks]

	firm = (; ζ, a, A)
	εs = zeros(n_banks)
	εs[i] = min(ε, a)
	
	(; bank_df) = equilibrium(banks, IM, firm, εs)

	layout = Shell()
	
	visualize_bank_firm_network(IM, bank_df; figure=figure(1.0, 1.0), hidespines=false, start = 1/6, layout, add_legend=true, show_firms=false) |> aside_figure
end

# ╔═╡ b54dc329-7764-41e6-8716-ef20bef0b29b
md"""
## Assignment infrastructure
"""

# ╔═╡ cc9d4ea3-4d1b-4f6e-8d49-77fec05e2804
begin
	hint(text) = Markdown.MD(Markdown.Admonition("hint", "Hint", [text]))
	almost(text) = Markdown.MD(Markdown.Admonition("warning", "Almost there!", [text]))
	still_missing(text=md"Replace `missing` with your answer.") = Markdown.MD(Markdown.Admonition("warning", "Here we go!", [text]))
	keep_working(text=md"The answer is not quite right.") = Markdown.MD(Markdown.Admonition("danger", "Keep working on it!", [text]))
	yays = [md"Great!", md"Yay ❤", md"Great! 🎉", md"Well done!", md"Keep it up!", md"Good job!", md"Awesome!", md"You got the right answer!", md"Let's move on to the next section."]
	correct(text=rand(yays)) = Markdown.MD(Markdown.Admonition("correct", "Got it!", [text]))
end

# ╔═╡ 1d8be056-aa73-4491-8d6e-57502ccc48be
members = let
	names = map(group_members) do (; firstname, lastname)
		firstname * " " * lastname
	end
	join(names, ", ", " & ")
end

# ╔═╡ 048d7ac7-204c-4bb4-aea0-612af18bb6d2
md"""
*Assignment submitted by* **$members** (*group $(group_number)*)
"""

# ╔═╡ 5ef89a83-8c7c-4fcd-8498-9c2f452a13c8
function wordcount(text)
	stripped_text = strip(replace(string(text), r"\s" => " "))
   	words = split(stripped_text, (' ', '-', '.', ',', ':', '_', '"', ';', '!', '\''))
   	length(filter(!=(""), words))
end

# ╔═╡ 0a0c51ec-9443-4901-a3db-f205c4e94e99
@test wordcount("  Hello,---it's me.  ") == 4

# ╔═╡ 01f0e271-acb6-44f8-85c5-ada44f8d401b
@test wordcount("This;doesn't really matter.") == 5

# ╔═╡ 8ed4dff0-c0b5-4247-a779-59ef7aa500a1
show_words(answer) = md"_approximately $(wordcount(answer)) words_"

# ╔═╡ 55d00dc9-b257-446b-9d60-688a43b79a7f
show_words(answer_a)

# ╔═╡ 6a615560-37dd-4c08-852e-da67e3a6ccf2
show_words(answer_b)

# ╔═╡ c0203246-97cd-4568-93bb-d79898fa7233
show_words(answer_c)

# ╔═╡ 61d2b83f-12a2-46e4-bf41-477053455e4f
show_words(answer_d)

# ╔═╡ 87aae64c-c713-473f-8a8c-d28d5973273f
show_words(answer_e)

# ╔═╡ b5a91cd7-6e0b-4690-9dfa-36a2986ac8db
function show_words_limit(answer, limit)
	count = wordcount(answer)
	if count < 1.02 * limit
		return show_words(answer)
	else
		return almost(md"You are at $count words. Please shorten your text a bit, to get **below $limit words**.")
	end
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AlgebraOfGraphics = "cbdf2221-f076-402e-a563-3d30da359d67"
CairoMakie = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
CategoricalArrays = "324d7699-5711-5eae-9e2f-1d82baa6b597"
Chain = "8be319e6-bccf-4806-a6f7-6fae938471bc"
Colors = "5ae59095-9a9b-59fe-a467-6f913c188581"
DataAPI = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
DataFrameMacros = "75880514-38bc-4a95-a458-c2aea5a3a702"
DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
GraphMakie = "1ecd5474-83a3-4783-bb4f-06765db800d2"
Graphs = "86223c79-3864-5bf0-83f7-82e725a168b6"
HypertextLiteral = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
LaTeXStrings = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Makie = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
MarkdownLiteral = "736d6165-7244-6769-4267-6b50796e6954"
NamedTupleTools = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
NetworkLayout = "46757867-2c16-5918-afeb-47bfcb05e46a"
PlutoLinks = "0ff47ea0-7a50-410d-8455-4348d5de0420"
PlutoTeachingTools = "661c6b06-c737-4d37-b85c-46df65de6f69"
PlutoTest = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
Roots = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
SimpleWeightedGraphs = "47aef6b3-ad0c-573a-a1e2-d07658019622"
Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"

[compat]
AlgebraOfGraphics = "~0.6.14"
CairoMakie = "~0.10.2"
CategoricalArrays = "~0.10.7"
Chain = "~0.5.0"
Colors = "~0.12.10"
DataAPI = "~1.14.0"
DataFrameMacros = "~0.4.1"
DataFrames = "~1.5.0"
GeometryBasics = "~0.4.5"
GraphMakie = "~0.5.3"
Graphs = "~1.8.0"
HypertextLiteral = "~0.9.4"
LaTeXStrings = "~1.3.0"
Makie = "~0.19.2"
MarkdownLiteral = "~0.1.1"
NamedTupleTools = "~0.14.3"
NetworkLayout = "~0.4.4"
PlutoLinks = "~0.1.6"
PlutoTeachingTools = "~0.2.6"
PlutoTest = "~0.2.2"
PlutoUI = "~0.7.50"
Revise = "~3.5.1"
Roots = "~2.0.9"
SimpleWeightedGraphs = "~1.3.0"
StatsBase = "~0.33.21"
StructArrays = "~0.6.14"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.5"
manifest_format = "2.0"
project_hash = "a79ce4a828f0c4da7303c24e89efd794fde31c5e"

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
uuid = "cbdf2221-f076-402e-a563-3d30da359d67"
version = "0.6.14"

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

[[deps.Cairo]]
deps = ["Cairo_jll", "Colors", "Glib_jll", "Graphics", "Libdl", "Pango_jll"]
git-tree-sha1 = "d0b3f8b4ad16cb0a2988c6788646a5e6a17b6b1b"
uuid = "159f3aea-2a34-519c-b102-8c37f9878175"
version = "1.0.5"

[[deps.CairoMakie]]
deps = ["Base64", "Cairo", "Colors", "FFTW", "FileIO", "FreeType", "GeometryBasics", "LinearAlgebra", "Makie", "SHA", "SnoopPrecompile"]
git-tree-sha1 = "abb7df708fe1335367518659989627100a61f3f0"
uuid = "13f3f980-e62b-5c42-98c6-ff1f3baf88f0"
version = "0.10.2"

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

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "4f619d394ac521dc59cb80a2cd8f78578e483a9d"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.2.1"

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

[[deps.CommonMark]]
deps = ["Crayons", "JSON", "SnoopPrecompile", "URIs"]
git-tree-sha1 = "e2f4627b0d3f2c1876360e0b242a7c23923b469d"
uuid = "a80b9123-70ca-4bc0-993e-6e3bcb318db6"
version = "0.8.10"

[[deps.CommonSolve]]
git-tree-sha1 = "9441451ee712d1aec22edad62db1a9af3dc8d852"
uuid = "38540f10-b2f7-11e9-35d8-d573e4eb0ff2"
version = "0.2.3"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "61fdd77467a5c3ad071ef8277ac6bd6af7dd4c04"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.6.0"

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

[[deps.DataAPI]]
git-tree-sha1 = "e8119c1a33d267e16108be441a287a6981ba1630"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.14.0"

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
git-tree-sha1 = "7be5f99f7d15578798f338f5433b6c432ea8037b"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.16.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "d3ba08ab64bdfd27234d3f61956c966266757fe6"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.13.7"

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

[[deps.GLM]]
deps = ["Distributions", "LinearAlgebra", "Printf", "Reexport", "SparseArrays", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns", "StatsModels"]
git-tree-sha1 = "884477b9886a52a84378275737e2823a5c98e349"
uuid = "38e38edf-8417-5370-95a0-9cbb8c7f171a"
version = "1.8.1"

[[deps.GPUArraysCore]]
deps = ["Adapt"]
git-tree-sha1 = "1cd7f0af1aa58abc02ea1d872953a97359cb87fa"
uuid = "46192b85-c4d5-4398-a991-12ede77f4527"
version = "0.1.4"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "e07a1b98ed72e3cdd02c6ceaab94b8a606faca40"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.2.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "fe9aea4ed3ec6afdfbeb5a4f39a2208909b162a6"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.5"

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

[[deps.GraphMakie]]
deps = ["GeometryBasics", "Graphs", "LinearAlgebra", "Makie", "NetworkLayout", "PolynomialRoots", "StaticArrays"]
git-tree-sha1 = "72882a1584f367cfecc83e3e8a232c7720c262cd"
uuid = "1ecd5474-83a3-4783-bb4f-06765db800d2"
version = "0.5.3"

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
git-tree-sha1 = "1cf1d7dcb4bc32d7b4a5add4232db3750c27ecb4"
uuid = "86223c79-3864-5bf0-83f7-82e725a168b6"
version = "1.8.0"

[[deps.GridLayoutBase]]
deps = ["GeometryBasics", "InteractiveUtils", "Observables"]
git-tree-sha1 = "678d136003ed5bceaab05cf64519e3f956ffa4ba"
uuid = "3955a311-db13-416c-9275-1d80ed98e5e9"
version = "0.9.1"

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
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

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
git-tree-sha1 = "82aec7a3dd64f4d9584659dc0b62ef7db2ef3e19"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.2.0"

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

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "d9ae7a9081d9b1a3b2a5c1d3dac5e2fdaafbd538"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.22"

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

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Printf", "Requires"]
git-tree-sha1 = "2422f47b34d4b127720a18f86fa7b1aa2e141f29"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.18"

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
git-tree-sha1 = "0a1b7c2863e44523180fdb3146534e265a91870b"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.23"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "60168780555f3e663c536500aa790b6368adc02a"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.3.0"

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
git-tree-sha1 = "274fa9c60a10b98ab8521886eb4fe22d257dca65"
uuid = "ee78f7c6-11fb-53f2-987a-cfe4a2b5a57a"
version = "0.19.2"

[[deps.MakieCore]]
deps = ["Observables"]
git-tree-sha1 = "2c3fc86d52dfbada1a2e5e150e50f06c30ef149c"
uuid = "20f20a25-4f0e-4fdf-b5d1-57303727442b"
version = "0.6.2"

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
git-tree-sha1 = "f04120d9adf4f49be242db0b905bea0be32198d1"
uuid = "0a4f8689-d25c-4efe-a92b-7142dfc1aa53"
version = "0.5.4"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

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

[[deps.NamedTupleTools]]
git-tree-sha1 = "90914795fc59df44120fe3fff6742bb0d7adb1d0"
uuid = "d9ec5142-1e00-5aa0-9d6a-321866360f50"
version = "0.14.3"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore", "ImageMetadata"]
git-tree-sha1 = "5ae7ca23e13855b3aba94550f26146c01d259267"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.1.0"

[[deps.NetworkLayout]]
deps = ["GeometryBasics", "LinearAlgebra", "Random", "Requires", "SparseArrays"]
git-tree-sha1 = "cac8fc7ba64b699c678094fa630f49b80618f625"
uuid = "46757867-2c16-5918-afeb-47bfcb05e46a"
version = "0.4.4"

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
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

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
git-tree-sha1 = "6f4fbcd1ad45905a5dee3f4256fabb49aa2110c6"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.5.7"

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
git-tree-sha1 = "eb11c2e0586fdf48d5d262ba6e29e438ccc512d9"
uuid = "661c6b06-c737-4d37-b85c-46df65de6f69"
version = "0.2.6"

[[deps.PlutoTest]]
deps = ["HypertextLiteral", "InteractiveUtils", "Markdown", "Test"]
git-tree-sha1 = "17aa9b81106e661cffa1c4c36c17ee1c50a86eda"
uuid = "cb4044da-4d16-4ffa-a6a3-8cad7f73ebdc"
version = "0.2.2"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "5bb5129fdd62a2bbbe17c2756932259acf467386"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.50"

[[deps.PolygonOps]]
git-tree-sha1 = "77b3d3605fc1cd0b42d95eba87dfcd2bf67d5ff6"
uuid = "647866c9-e3ac-4575-94e7-e3d426903924"
version = "0.1.2"

[[deps.PolynomialRoots]]
git-tree-sha1 = "5f807b5345093487f733e520a1b7395ee9324825"
uuid = "3a141323-8675-5d76-9d11-e1df1406c778"
version = "1.0.0"

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
git-tree-sha1 = "786efa36b7eff813723c4849c90456609cf06661"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.8.1"

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

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

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

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "90cb983381a9dc7d3dff5fb2d1ee52cd59877412"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.5.1"

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
deps = ["ChainRulesCore", "CommonSolve", "Printf", "Setfield"]
git-tree-sha1 = "9c2f5d3768804ed465f0c51540c6074ae9f63900"
uuid = "f2b01f46-fcfa-551c-844a-d8ac1e96c665"
version = "2.0.9"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.SIMD]]
deps = ["SnoopPrecompile"]
git-tree-sha1 = "8b20084a97b004588125caebf418d8cab9e393d1"
uuid = "fdea26ae-647d-5447-a871-4b548cad5224"
version = "3.4.4"

[[deps.ScanByte]]
deps = ["Libdl", "SIMD"]
git-tree-sha1 = "2436b15f376005e8790e318329560dcc67188e84"
uuid = "7b38b023-a4d7-4c5e-8d43-3f3097f304eb"
version = "0.3.3"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

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

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.SimpleWeightedGraphs]]
deps = ["Graphs", "LinearAlgebra", "Markdown", "SparseArrays", "Test"]
git-tree-sha1 = "7d0b07df35fccf9b866a94bcab98822a87a3cb6f"
uuid = "47aef6b3-ad0c-573a-a1e2-d07658019622"
version = "1.3.0"

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
git-tree-sha1 = "2d7d9e1ddadc8407ffd460e24218e37ef52dd9a3"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.16"

[[deps.StaticArraysCore]]
git-tree-sha1 = "6b7ba252635a5eff6a0b0664a41ee140a1c9e72a"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

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
git-tree-sha1 = "a5e15f27abd2692ccb61a99e0854dfb7d48017db"
uuid = "3eaba693-59b7-5ba5-a881-562e759f1c8d"
version = "0.6.33"

[[deps.StringManipulation]]
git-tree-sha1 = "46da2434b41f41ac3594ee9816ce5541c6096123"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.0"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "GPUArraysCore", "StaticArraysCore", "Tables"]
git-tree-sha1 = "b03a3b745aa49b566f128977a7dd1be8711c5e71"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.14"

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
git-tree-sha1 = "c79322d36826aa2f4fd8ecfa96ddb47b174ac78d"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.10.0"

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

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "7e6b0e3e571be0b4dd4d2a9a3a83b65c04351ccc"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.3"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "94f38103c984f89cf77c402f2a68dbd870f8165f"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.11"

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

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

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
# ╟─048d7ac7-204c-4bb4-aea0-612af18bb6d2
# ╟─639f04a2-0979-41e2-b465-a4f52049166d
# ╟─52052d98-0c41-45ec-95bf-d936b1c43e81
# ╠═72e25b9c-89e3-441b-bf89-c1122535318a
# ╠═a4dba7da-9d63-4d0d-a736-565d6ccf7d09
# ╟─cdfc35e1-5082-4959-9d30-3abb2dc8a7ac
# ╟─944b07ed-20e5-4068-ac5e-8a48f866fdd2
# ╟─1d13564d-7256-4cf1-919f-6d8298cc476b
# ╟─1d5d8c8a-8d86-426f-bb17-bd2279d91ff1
# ╟─7b876239-8ddc-4929-ad52-752edb72c0eb
# ╠═e11a99df-d0f2-4838-b325-473d3043be98
# ╠═144cfaf2-b78c-4f87-8b35-559914abf532
# ╠═073982c7-6333-43f6-866a-91a49f8ba7eb
# ╠═8a2f3e4d-9c61-4a21-a229-58f731964181
# ╟─73ba4210-d8f6-4a74-bf4d-d7bc0902bb5e
# ╟─eba94b09-f061-432a-80ce-a68be83e6a99
# ╟─ffc743af-b97c-4083-a02e-ea5725829f2a
# ╟─8e6bdcfd-3cb0-4f6b-8dcf-f99f4afe87c2
# ╟─bafc7db1-ac1d-4314-be83-0b9f6c63b5fc
# ╠═ba39c037-c4e1-465e-a232-e5194caa14ba
# ╟─acfdf348-e8bf-41ee-aa19-fe7ec29087cc
# ╠═834ff5a0-3bae-44ca-a413-f45d3394b511
# ╠═2f69a836-7fed-42e6-a509-096bc8cabbc2
# ╠═26a5c4f4-75c1-4a4c-a813-5fad0c571da1
# ╠═5d283f89-4419-4b14-81ba-9826b4f1689e
# ╟─a0ac4092-cd27-4523-85f9-f4a4d81456b3
# ╟─6fddaaa4-688d-4f40-a1b3-04bf62024955
# ╟─01eb57fd-a815-41f2-9e25-7730bff7917d
# ╠═6af1e245-aaf8-4b7d-977a-1d76998f7cf8
# ╠═5b6f26db-f63f-4db5-ae3e-d41ae476948f
# ╠═0e4a53b5-1751-4723-a05d-a5504e427e3c
# ╟─523e177c-f74d-4805-8a2b-9c27a4b0bc63
# ╠═4ae1b6c2-bb63-4ca8-b8ec-057c8d2a371f
# ╟─efbc1a42-ecc2-4907-a981-bd1d29ca0803
# ╟─7f058e4a-9b12-41c9-9fd7-4ad023058a14
# ╟─96878ebb-fbc0-4d53-998a-210e13a42492
# ╟─e2041c57-0e3d-4dad-9bab-d70434f18509
# ╠═9446b4db-4d93-4153-84e2-73d01fb31254
# ╠═aaffd333-3aa1-48ee-b5db-d577bd7da830
# ╠═a8a8a197-d54d-4c67-b7ce-19bdc8b64401
# ╠═7970090f-eb3a-488e-b225-ce4198494f1d
# ╟─da7c558d-a2b5-41a8-9c78-3e39a00dfd31
# ╟─bf719f30-72c1-488e-ba77-c183effb7c60
# ╟─a0767d80-0857-47ef-90a1-72bc34064716
# ╠═c920d82c-cfe9-462a-bacd-436f01c314cf
# ╟─72e5f4d3-c3e4-464d-b35c-2cf19fa9d4b5
# ╟─d6345a8b-6d8f-4fd2-972b-199412cbdc26
# ╟─c99e52e2-6711-4fb6-bcc0-8e4f378ed479
# ╟─15f45669-516b-4f3f-9ec1-f9e2c1d2e71a
# ╟─d7111001-f632-4d0d-a2c7-7bbfd67bf87d
# ╟─0d18cdf0-441e-4ca9-98e3-50bc3efa837f
# ╟─51bfc57f-7b06-4e27-af32-51df38af30a1
# ╠═b817cdf6-dfe6-4607-98da-2299a33d4906
# ╟─07e9c77f-e05e-452c-8c47-cdd9dfc8e2fc
# ╠═f5293bee-9413-4507-8568-54836eb6d4a2
# ╟─49c2fb2d-de6e-4ab2-a558-59fb153cf703
# ╠═c0c711de-6916-4ab9-ab73-c476654332c4
# ╟─a95431eb-14a0-4dc3-bbe6-9c409f6cc596
# ╟─c7b99d3c-5d32-45e6-84fa-8a6513e6beb9
# ╟─f00d9e1a-b111-4b6a-95f5-b9736329befe
# ╟─f8eb242f-a974-48aa-9173-b0bc7ff697d5
# ╠═c2633df1-2e30-4387-8749-de3280b0602d
# ╟─253ab06f-6284-4cbf-b2a2-232ff99548c9
# ╠═1d058f8b-16f5-4744-8425-452876006c47
# ╟─a9d27019-72b7-4257-b72a-12952b516db9
# ╟─0fb4d187-f03a-435b-b9fc-188925e058f1
# ╟─27fadf93-0b17-446e-8001-d8394b7befaa
# ╠═aed99485-cec3-4bf3-b05d-4d20572ec907
# ╠═db841316-9106-40bb-9ca3-ae6f8b975404
# ╟─900a4b24-04dc-4a1b-9829-a166cf9eb7fb
# ╟─871f33f0-4882-4ff0-bbde-eb954059e907
# ╟─cb435758-d2a8-4203-9915-971e041d4319
# ╟─2fe4c931-d4b2-4b4d-8634-73573125cfb5
# ╟─4242531a-e74f-4618-939f-2adf9d6e1db2
# ╟─6f28bcfb-2b56-4548-bbd1-9528876525dd
# ╠═fa50a9b9-cdc4-4d84-ae3b-db039f1609e4
# ╟─55d00dc9-b257-446b-9d60-688a43b79a7f
# ╟─cfa98250-1d4a-43a6-a99f-cf106001f3cb
# ╠═9b0da913-fc5e-42ea-bc5f-e37bd59f2cd2
# ╠═a7e36f2c-d588-4fc5-a247-4323d646a51b
# ╟─6a615560-37dd-4c08-852e-da67e3a6ccf2
# ╟─77b5b0ea-bfae-406f-8fcf-472165bdcd1d
# ╠═dbfd2f13-89f4-4932-a778-b2d375d45ac6
# ╟─c0203246-97cd-4568-93bb-d79898fa7233
# ╟─63a8c85f-9cd1-4cf5-9f58-e482494f8d24
# ╠═55904c61-4531-4984-b73c-1065a7114772
# ╟─61d2b83f-12a2-46e4-bf41-477053455e4f
# ╟─264e3358-babf-4bf4-9b57-f436676aa02a
# ╠═36e610ff-1f42-4d58-b0a7-1bc33bd0d4af
# ╟─87aae64c-c713-473f-8a8c-d28d5973273f
# ╟─e8637286-ea8b-49c8-b49f-1ab556b83f0c
# ╠═4878d013-053d-4490-b7c6-db0de4bf0a82
# ╠═a047aeaa-fa54-4cbf-90f4-42d0537b7d06
# ╠═b1c0d43c-f483-4290-998f-177ce79f41fa
# ╟─4e9b785f-ad74-4aa8-ad48-89fa8b236939
# ╠═1a997e44-f29c-4c55-a953-a9039f096d47
# ╟─78bedfcc-3671-4852-985b-3e1b5aaade5a
# ╠═0f03f537-f589-4abd-9587-0bb18835d9b9
# ╟─25e84f19-9cd8-43ad-ae6a-d500b8ac74b6
# ╠═74262581-3c64-4e5b-9328-416c4e1efc91
# ╠═2b405f2b-3256-4c47-8334-c2d93713d409
# ╠═103babf9-bbc3-4c77-9185-72f792a09706
# ╠═bb0e41cf-66fb-4cae-8cd9-6ad771d1acf4
# ╠═969aa5eb-c56f-4115-83df-bb0ec911b6aa
# ╠═935da683-a4e2-4ddc-999f-55cb61390f39
# ╠═2d1da571-0561-4c28-bb63-35ca5f9538d5
# ╠═7c068d65-f472-44c2-af56-581bf9309bd5
# ╠═631ce85b-db76-4f3b-bda5-0c51ffb3bc43
# ╠═fa01788e-7585-4f3f-93b2-a44541065820
# ╠═6a89953d-310b-49c9-89e1-8d51c8b75be0
# ╠═50ac0182-32fe-4e21-8c6d-895ffc67ec27
# ╠═8b75bc15-e07b-43e5-adb3-c5d6481ee9d8
# ╠═cf61ebba-6800-4bfa-bb7c-cb9a6c846b65
# ╠═36ac25b9-27e5-4728-bbcf-920f231ff6ab
# ╟─11d6ac4f-e910-4a9f-9ee4-bdd270e9400b
# ╠═c1b1a22b-8e18-4d19-bb9b-14d2853f0b72
# ╠═cc1ff1e6-2968-4baf-b513-e963ab2ff1b4
# ╠═c65323f5-211d-4a95-aed3-d6129bdd083e
# ╠═5d263432-856b-4e9b-a303-a476222e8963
# ╟─25a01039-25e9-44b0-afd0-c3df37c7598f
# ╟─71231141-c2f5-4695-ade0-548a0039f511
# ╠═c993c1d8-8823-4db4-bf6e-bf2c21ea3d39
# ╟─f24387ee-c1cf-4ec0-a34e-4b4f33ee9010
# ╠═954d2dde-8088-4e91-bed3-f8339090b77b
# ╠═2872c686-8e4f-4230-a07a-5d988aba39b7
# ╠═3198bb79-f0c7-4b01-8dce-ef7629e8d7e6
# ╠═e9365f8d-6f58-458c-855d-c0444c6f409f
# ╟─fbf40f06-f63a-40ca-bbe3-78104d39ee71
# ╠═f863066f-270f-4972-8046-7f51cb697ac5
# ╠═d285879a-bdfd-4efa-aa5d-9dacf08a2dc6
# ╠═70bbb90e-06f1-4b60-a952-275866945c58
# ╟─3040803d-e95a-40bc-aa72-92c7d158e226
# ╟─f13e010a-d394-4e40-9535-c5e2e3e226aa
# ╠═9cdc97e5-67f8-4e71-97d8-9cda5e9d7bd8
# ╠═601b1453-20bd-4e93-abb7-4fad4f5adf8c
# ╠═955c190a-1d00-4c78-bdfb-daac31edf76f
# ╠═7871b479-2aaf-42c1-ad84-42ac17cfc6e1
# ╟─7c026a42-1c05-4968-b068-c8561ca5a2db
# ╠═7b70a862-faf4-4c42-917c-238718c43708
# ╠═fc01674e-73de-4cbf-9c80-fa1ea47fcb21
# ╠═222665dc-397b-4e58-b3ee-935b115cf13d
# ╠═3998e822-d7e0-40f2-a866-71a3d81c4ca3
# ╟─0d207cbe-5b97-4bd7-accf-fc49ce4522e9
# ╠═ebf41ae2-3985-439b-b193-eabfab701d16
# ╠═b7d74679-904c-44c4-bedd-89f1b68a5e42
# ╠═1b9c8f2b-8011-46f1-b46d-479031eb9ac3
# ╠═5104f0a5-a551-4f4c-8f89-2b8f834b3587
# ╠═f63aee78-1209-40dd-9c9d-2699194807d8
# ╠═6a2b3b9c-df52-41c3-908b-5dbf052ad107
# ╠═ab1cf7ab-80fb-4423-a924-1d6e24e9c9bc
# ╠═0ff81eb8-d843-49a4-af93-ec2414797e87
# ╟─ffe50130-a9cb-4ba9-a861-c247bf688873
# ╟─d3e5b8f2-51b1-4ba4-97b4-2be156a74643
# ╠═00a4054a-dd10-4c5c-ae20-c0dc176e8e18
# ╠═d5c128e0-6371-4ad2-8bfc-c17faadc520b
# ╠═38a1f960-7178-4968-89e4-6b659b64baa2
# ╠═f689810f-8f43-4368-a822-0ee4f3015271
# ╠═5bde1113-5297-4e93-b052-7a0f93f4ea84
# ╠═73456a8f-8534-4b55-8af6-ff7952bd3a3a
# ╠═e6def365-4438-4591-acbd-60d9de466b0a
# ╠═055a811c-bc4a-4959-88f1-bc71e8749313
# ╠═e2045a53-21d3-4c7c-ad95-3fc8d2444821
# ╠═3c211dad-73b0-4715-b066-e10005f8f3a7
# ╠═ee10e134-84be-41f8-9e7e-c1cd68227977
# ╠═2a5f02e7-fec1-41b3-b8b4-b33b1cc1232c
# ╠═47118d1b-8471-4288-85fe-d3e94667dc96
# ╠═dffce961-6844-4a44-beea-ab085c2f9f3f
# ╠═f708e6c0-cfac-4b4d-a3ed-69b98883294a
# ╠═d6cb95c1-a075-4544-9031-58aef65c7577
# ╠═1bb841e0-ddd2-4571-83a5-d929e0a8a69c
# ╠═7fbcfbde-0b5e-4bf2-9eda-9b15a4dd6bec
# ╟─ff77f7b4-52d1-4fc8-abfc-e623f7bcd423
# ╠═c7c8581c-e82c-428e-a599-3c003cc0151c
# ╟─69698f6a-110d-49ea-bd57-16af0de90844
# ╠═7d2d3618-7810-444d-90c0-0d592f8eba8c
# ╠═12968f5e-6e75-4429-bff8-0d1c644330a7
# ╟─2c22e86e-1139-4f66-8992-093d38f4c6cb
# ╠═60293fac-4d20-409b-bfd2-d5283f189320
# ╠═bef59a5b-6952-42aa-b700-4ad81d848fe4
# ╠═ce563ab5-6324-4a86-be61-7a107ff0e3b3
# ╟─4ccb52c2-a177-4b13-8ce5-2af6f286c737
# ╠═1f23193e-0bca-4bee-82b1-453100721ab6
# ╠═4a1fd4a3-71a0-4947-8b75-cc2100fb636b
# ╟─3496c181-c4c3-4b1b-a5e5-83df27182c99
# ╠═8f15ccf5-860a-47ca-b5c2-842c4e6f861a
# ╠═3758f272-297c-4325-9a7c-042b4f41d615
# ╠═a349099a-64aa-4e86-9bf2-b6157feff394
# ╠═c4032c7a-73f4-4342-a5a7-19fd14017402
# ╠═e8d198fa-11db-43b2-846e-311abe9f59aa
# ╠═2711b045-b879-4225-afae-13328a4e14b3
# ╟─b54dc329-7764-41e6-8716-ef20bef0b29b
# ╠═cc9d4ea3-4d1b-4f6e-8d49-77fec05e2804
# ╠═1d8be056-aa73-4491-8d6e-57502ccc48be
# ╠═5ef89a83-8c7c-4fcd-8498-9c2f452a13c8
# ╠═c9775558-94f2-4c8c-9d0e-ab948fa5ead4
# ╠═0a0c51ec-9443-4901-a3db-f205c4e94e99
# ╠═01f0e271-acb6-44f8-85c5-ada44f8d401b
# ╠═8ed4dff0-c0b5-4247-a779-59ef7aa500a1
# ╠═b5a91cd7-6e0b-4690-9dfa-36a2986ac8db
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

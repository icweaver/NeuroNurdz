### A Pluto.jl notebook ###
# v0.14.4

using Markdown
using InteractiveUtils

# ╔═╡ acc935c2-d23d-4e64-bd84-fb96cd67d1ed
begin
	using Markdown
	using PlutoUI
	using CairoMakie
	using StatsBase
	using BenchmarkTools
end

# ╔═╡ 3fc988b1-803f-4cb4-a8e3-0277c0785e83
md"""
In this notebook, we will take a look at quantifying the correlation between the signals fired from two separate neurons.

**The data:** A vector of time measurements from each neuron recording a firing event (in seconds):

``x\_i = [1, 2]\text{ s}``  \
``x\_j = [3]\text{ s}``

**The desired output:** A list of "circular lags" between the firing times between the two neurons. For the example above, we would expect the (sorted) output to be:

```julia
lags = [-3, -3, -2, -2, -1, -1, 0, 0] s
```

For example, a forward shift of one time interval would lead to the following lags:

```math
[2, 3]\text{ s} \ominus [3]\text{ s} = [-1, 0]\text{ s} \quad,
```

where:
```math
\mathbf{u} \ominus \mathbf{v} \equiv [u_i - v_j],
\text { for } u_i \in \mathbf{u},\ v_j \in \mathbf{v}
\iff |u_i - v_j| \le \epsilon
```
is the vector of all pair-wise differences in signal firing times between the two neurons within a threshold ``\epsilon``.

If we repeat this for all shifts in a single cycle (including a zero-shift), we would then measure the following lags between the timeseries measurements for each neuron:

```math
\begin{align}
\left([1, 2], [2, 3], [0, 3], [0, 1]\right)\text{ s} \ominus [3]\text{ s}
&= [-2, -1, -1, 0, -3, 0, -3, -2]\text{ s} \\
&\stackrel{\text{sort}}{=} [-3, -3, -2, -2, -1, -1, 0, 0] \quad.
\end{align}
```

Because order does not matter for this operation, we can accomplish this relatively straightfowardly using modulus arithmetic:
"""

# ╔═╡ 4a98643c-a0fa-40e8-89e8-19fdee5429e0
xᵢ = [1, 2]

# ╔═╡ f4241a88-e1dc-4a38-b05c-d59f2a1bca5b
xⱼ = [3]

# ╔═╡ 1e1d8f53-3e05-4646-9c34-f18e0585d787
begin
	@doc raw"""
		compute_lags(u::Vector{Real}, v::Vector{Real}; ϵ::Real)
	
	Returns ``\mathbf{u} \ominus \mathbf{v}`` for the given threshold
	``\epsilon``, where:
	
	```math
	\mathbf{u} \ominus \mathbf{v} \equiv [u_i - v_j],
	\text { for } u_i \in \mathbf{u},\ v_j \in \mathbf{v}
	\iff |u_i - v_j| \le \epsilon
	```
	"""
	function compute_lags(u, v; ϵ=10.0)
		lags = Float64[]
		for uᵢ ∈ u, vⱼ ∈ v
			lag!(lags, uᵢ, vⱼ; ϵ=ϵ)
		end
		return lags
	end
	
	# Compute within-threshold lag
	lag!(lags, uᵢ, vⱼ; ϵ=10) = abs(uᵢ - vⱼ) ≤ ϵ && push!(lags, uᵢ - vⱼ)
end

# ╔═╡ e8751abc-c980-4efe-a327-8905394b386e
function lag_vector(xᵢ_measurements, xⱼ_measurements)
	N_intervals = maximum((xᵢ, xⱼ))[1] + 1 # Latest time that others will wrap around
	us = [(xᵢ_measurements .+ i) .% N_intervals for i in 0:N_intervals-1]
	v = xⱼ_measurements
	return Base.Iterators.flatten(compute_lags.(us, v)) |> collect
end

# ╔═╡ 3608d193-1ead-45fb-aad0-8075793c00a5
lag_vector(xᵢ, xⱼ) |> sort

# ╔═╡ Cell order:
# ╟─3fc988b1-803f-4cb4-a8e3-0277c0785e83
# ╠═4a98643c-a0fa-40e8-89e8-19fdee5429e0
# ╠═f4241a88-e1dc-4a38-b05c-d59f2a1bca5b
# ╠═3608d193-1ead-45fb-aad0-8075793c00a5
# ╠═1e1d8f53-3e05-4646-9c34-f18e0585d787
# ╠═e8751abc-c980-4efe-a327-8905394b386e
# ╠═acc935c2-d23d-4e64-bd84-fb96cd67d1ed

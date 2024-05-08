### A Pluto.jl notebook ###
# v0.19.42

using Markdown
using InteractiveUtils

# ╔═╡ d5e4988a-09a0-11ef-2eb9-97277c169c65
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
	Pkg.status()
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    using WaveSim
	using Plots
	theme(:dracula)
end

# ╔═╡ 206b9281-0aab-435d-8592-6e5169449832
md"""
# Activate Local Environment
"""

# ╔═╡ b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
md"""
# Problem Setting
"""

# ╔═╡ 2b94be9f-5df9-42c3-bdd9-85515efd6c3e
begin
	# Space Discretization
	x_details = Dict("x_min" => 0, "x_max" => 100, "nx" => 101)

	# Time Discretization
	t_details = Dict("t_min" => 0, "t_max" => 1, "nt" => 1001)

	# Medium Velocity
	c = 334.

	# Source Details
	src_details = Dict("f0" => 10., "t0" => .1, "isrc" => 50.)
end;

# ╔═╡ c6d58808-b224-46f9-91f4-965c02dd32ec
begin
	# Space Discretization
	x_min = x_details["x_min"]
	x_max = x_details["x_max"]
	nx = x_details["nx"]
	x = range(start=x_min, stop=x_max, length=nx)
	
	# Time Discretization
	t_min = t_details["t_min"]
	t_max = t_details["t_max"]
	nt = t_details["nt"]
	t = range(start=t_min, stop=t_max, length=nt)
	
	# Source Function
	f0 = src_details["f0"] # Dominant Frequency
	t0 = src_details["t0"] # Source Time Shift
	
	src = -8. .* (t .- t0) .* f0 .* (exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2))
	
	# Plotting source Function
    p = plot(t, src, title="Source Function", label="Source Function", lw=2)
	xlabel!("t")
	ylabel!("Amplitude")
end

# ╔═╡ 2985d0e2-ae9b-41f3-996a-b8422c17ded8
md"""
# Finite Difference Simulation
"""

# ╔═╡ 8452844e-00cf-4fb2-9220-538ac9fdcbfc
# Solution
p_sols = WaveSim.wave_sim_1d(x_details, t_details, src_details, c);

# ╔═╡ f2e2b4b4-31ac-4d15-948c-539ddae6b5e4
begin
	wave_plot = plot(x, p_sols[100,:], ylims=(minimum(p_sols), maximum(p_sols)), label="Wave")
	xlabel!("x(meters)")
	ylabel!("u(x,t=$(100*(t_max-t_min)/(nt-1)))")
	title!("Wave at t=$(100*(t_max-t_min)/(nt-1)) secs")
	savefig("./docs/src/img/index/1d_wave.png")
end

# ╔═╡ a41a2651-35da-4c97-aaa1-74b3f35a6e67
wave_plot

# ╔═╡ 5f829cd4-1bca-4473-a8a7-aca6a91cbd91
# ╠═╡ disabled = true
#=╠═╡
begin
	# Animation
    anim = @animate for it in range(start=1, stop=t_details["nt"], length=min(t_details["nt"], 60*(t_details["t_max"]-t_details["t_min"])*10))
    plot(x, p_sols[Int(floor(it)),:], ylims=(minimum(p_sols), maximum(p_sols)))
    end
	
	mp4(anim, "anims/1D.mp4", fps=60)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─206b9281-0aab-435d-8592-6e5169449832
# ╠═d5e4988a-09a0-11ef-2eb9-97277c169c65
# ╟─b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
# ╠═2b94be9f-5df9-42c3-bdd9-85515efd6c3e
# ╟─c6d58808-b224-46f9-91f4-965c02dd32ec
# ╟─2985d0e2-ae9b-41f3-996a-b8422c17ded8
# ╠═8452844e-00cf-4fb2-9220-538ac9fdcbfc
# ╠═f2e2b4b4-31ac-4d15-948c-539ddae6b5e4
# ╟─a41a2651-35da-4c97-aaa1-74b3f35a6e67
# ╠═5f829cd4-1bca-4473-a8a7-aca6a91cbd91

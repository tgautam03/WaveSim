### A Pluto.jl notebook ###
# v0.19.41

using Markdown
using InteractiveUtils

# ╔═╡ d5e4988a-09a0-11ef-2eb9-97277c169c65
begin
    import Pkg
    # activate the shared project environment
    Pkg.activate(Base.current_project())
    # instantiate, i.e. make sure that all packages are downloaded
    Pkg.instantiate()

    import WaveSim
end

# ╔═╡ a2939c30-789f-4e22-842e-41477d733339
# Dependencies
using Plots

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
	x_min = 0
	x_max = 100
	nx = 101
	dx = (x_max-x_min)/(nx-1)
	x = range(start=x_min, stop=x_max, length=nx)
	c = 334

	# Time Discretization
	t_min = 0
	t_max = 1
	nt = 10001
	dt = (t_max-t_min)/(nt-1)
	t = range(start=t_min, stop=t_max, length=nt)
end;

# ╔═╡ fabd4596-76fb-4423-aa47-6cccf6d2eb77
c * dt / dx

# ╔═╡ 49fc12b4-bcf5-4fc7-8117-9937ed996d0d
begin
	# Source Function
	f0 = 20 # Dominant Frequency
	t0 = 4 / f0 # Source Time Shift
	isrc = Int(floor(nx/2)) # Source Location
	
	src = -8. .* (t .- t0) .* f0 .* (exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2))
end;

# ╔═╡ 56fde8df-bb92-4f9e-8c74-ad25cb67f03c
begin
	# Plotting Source Function
	plot(t, src, title="Source Function", label="Source Function", lw=2)
	xlabel!("t")
	ylabel!("Amplitude")
end

# ╔═╡ 2985d0e2-ae9b-41f3-996a-b8422c17ded8
md"""
# Finite Difference Simulation
"""

# ╔═╡ 6a19974a-c7c5-442b-9bf3-c450ade8d84a
begin
	# Pressure Field at p(x,tᵢ)
	p = zeros(nx)
	# Pressure Field at p(x,tᵢ-dt)
	p_prev = zeros(nx)
	# Pressure Field at p(x,tᵢ+dt)
	p_next = zeros(nx)

	# Solution at each step
	p_sols = zeros(nt,nx)
end;

# ╔═╡ 8452844e-00cf-4fb2-9220-538ac9fdcbfc
# # Looping over time
# for it in range(start=2, stop=nt)
# 	# Looping over Space to evaluate 2nd order derivative wrt to x
# 	d2p_dx2 = zeros(nx)
# 	for ix in range(start=2, stop=nx-1)
# 		d2p_dx2[ix] = p[ix+1] - 2*p[ix] + p[ix-1]
# 	end
	
# 	# Updating Solution
# 	p_next = (c*dt/dx)^2 * d2p_dx2 + 2*p - p_prev 
# 	p_next[isrc] += dt^2 * src[it]
	
# 	p_prev[:] = p[:]
# 	p[:] = p_next[:]

# 	# plot!(x, p)
# 	p_sols[it,:] = p
	
# end

# ╔═╡ 436698bb-c42d-4af1-bfcb-9dc298825a46
# Looping over time
for it in range(start=2, stop=nt)
	# Looping over Space
	for ix in range(start=2, stop=nx-1)
		# Evaluating 2nd derivative wrt x
		d2p_dx2 = p[ix+1] - 2*p[ix] + p[ix-1]
		# Updating Solution
		if ix == isrc
			p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix] + dt^2 * src[it]
		else
			p_next[ix] = (c*dt/dx)^2 * d2p_dx2 + 2*p[ix] - p_prev[ix]
		end
	end

	# Current Sol becomes Previous Sol
	p_prev[:] = p[:]
	# Next Sol becomes Current Sol
	p[:] = p_next[:]

	# Storing solutions at each time step
	p_sols[it,:] = p
	
end

# ╔═╡ 5f829cd4-1bca-4473-a8a7-aca6a91cbd91
anim = @animate for it in range(start=1, stop=nt)
	plot(x, p_sols[it,:], ylims=(minimum(p_sols), maximum(p_sols)))
end

# ╔═╡ 842ba3b3-1bad-4093-bdf0-66bf043927c8
mp4(anim, "1D.mp4", fps=60)

# ╔═╡ Cell order:
# ╟─206b9281-0aab-435d-8592-6e5169449832
# ╠═d5e4988a-09a0-11ef-2eb9-97277c169c65
# ╠═a2939c30-789f-4e22-842e-41477d733339
# ╟─b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
# ╠═2b94be9f-5df9-42c3-bdd9-85515efd6c3e
# ╠═fabd4596-76fb-4423-aa47-6cccf6d2eb77
# ╠═49fc12b4-bcf5-4fc7-8117-9937ed996d0d
# ╟─56fde8df-bb92-4f9e-8c74-ad25cb67f03c
# ╟─2985d0e2-ae9b-41f3-996a-b8422c17ded8
# ╠═6a19974a-c7c5-442b-9bf3-c450ade8d84a
# ╟─8452844e-00cf-4fb2-9220-538ac9fdcbfc
# ╠═436698bb-c42d-4af1-bfcb-9dc298825a46
# ╠═5f829cd4-1bca-4473-a8a7-aca6a91cbd91
# ╠═842ba3b3-1bad-4093-bdf0-66bf043927c8

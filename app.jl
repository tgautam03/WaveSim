### A Pluto.jl notebook ###
# v0.19.42

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
	using PlutoUI
end

# ╔═╡ 206b9281-0aab-435d-8592-6e5169449832
md"""
# Activate Local Environment
"""

# ╔═╡ ed6ee4ce-fcd8-403a-9061-d688737f87be
function select(options, UItype)
	
	return PlutoUI.combine() do Child
		
		inputs = [
			md""" $(name): $(
				Child(name, ui)
			)"""
			
			for (name, ui) in zip(options, UItype)
		]
		
		md"""
		$(inputs)
		"""
	end
end

# ╔═╡ b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
md"""
# Problem Setting
"""

# ╔═╡ 8ab95f04-2dad-4274-a7b5-69c23d5c2bb1
md"""
## Space and Time Discretization
"""

# ╔═╡ ce9660d5-a355-464c-9142-1bfd084fc658
@bind space select(["Length of the String (in meters)", "Number of points in x dimension"], [PlutoUI.Select(100:100:1000), PlutoUI.Select(101:100:1001)])

# ╔═╡ 241fc0e7-8b6c-4d61-b7c8-95f331469dc4
begin
	x_max, nx = space;

	md"""
	**Space Discretization Step Size (in meters): $(round((x_max)/(nx-1), digits=2))**
	"""
end

# ╔═╡ ba518ad0-391d-45ca-8617-73fcc80bc8f4
md"""
## Time Discretization
"""

# ╔═╡ 36b894da-0428-459c-92ff-0f3e34943769
@bind time select(["How long to run the simulation (in seconds)", "Number of points in t dimension"], [PlutoUI.Select(1:10), PlutoUI.Select(1001:1000:10001)])

# ╔═╡ 17fe1153-8514-4876-b5b7-7ba0bd2d1d89
begin
	t_max, nt = time;

	md"""
	**Time Discretization Step Size (in seconds): $(round((t_max)/(nt-1), digits=10))**
	"""
end

# ╔═╡ d9417b50-b232-47cc-8d78-506aca36fbca
md"""
## Source Function
"""

# ╔═╡ bd5d0e35-4bdd-488f-b8a5-bea04c27d764
@bind source select(["Pick a source function", "Source Frequency (in Hz)", "Source Initiation Time (in seconds)", "Plot Source Function"], [PlutoUI.Select(["Derivative of Gaussian", "Gaussian"]), PlutoUI.Select(5:5:30), PlutoUI.Slider(0.1:0.01:round(t_max/2, digits=2), show_value=true), PlutoUI.Select([true, false])])

# ╔═╡ 2985d0e2-ae9b-41f3-996a-b8422c17ded8
md"""
# Finite Difference Simulation
"""

# ╔═╡ c26a9e81-4dab-40f4-a9c8-d20ef3a1ae4a
md"""
Source Location from x=0 (in meters)
"""

# ╔═╡ 14b955b6-f4ea-40ed-8340-996d75973135
@bind isrc PlutoUI.Select([25,x_max/2,x_max-25])

# ╔═╡ c6d58808-b224-46f9-91f4-965c02dd32ec
begin
	func_name, f0, t0, plt = source
	
	# Space Discretization
	x = range(start=0, stop=x_max, length=nx)
	
	# Time Discretization
	t = range(start=0, stop=t_max, length=nt)
	
	# Source Function
	if func_name == "Derivative of Gaussian"
		src = -8. .* (t .- t0) .* f0 .* (exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2))
	elseif func_name == "Gaussian"
		src = exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2)
	end

	# Space Discretization
	x_details = Dict("x_min" => 0, "x_max" => x_max, "nx" => nx)

	# Time Discretization
	t_details = Dict("t_min" => 0, "t_max" => t_max, "nt" => nt)

	# Source Details
	src_details = Dict("f0" => f0, "t0" => t0, "isrc" => isrc, "src" => src)
	

	if plt
		# Plotting source Function
	    p = scatter(t, src, title="Source Function", label="Source Function", markersize=1)
		xlabel!("t")
		ylabel!("Amplitude")
	end
end

# ╔═╡ e1c436ce-f879-40f3-8d2e-bd8562f5abbd
md"""
Medium Velocity (in m/s)
"""

# ╔═╡ d2b91e55-c921-486d-92ce-9cd04c54f60b
@bind c PlutoUI.Select(100.:50:500.)

# ╔═╡ e090a86d-7ed8-4428-b777-06383d62f492
# Solution
p_sols = wave_sim_1d(x_details, t_details, src_details, c);

# ╔═╡ 2a602d52-d4f6-4894-b214-ba8a1c16fa89
@bind it PlutoUI.Slider(range(start=1, stop=nt, length=min(nt, 60*(t_max)*10)), show_value=false)

# ╔═╡ cf0e9ff0-dd7e-40e2-8544-e8d18ffa9e25
begin
	# Plotting
	t_val = round(it*(t_max)/(nt-1), digits=2)
	plot(x, p_sols[Int(floor(it)),:], ylims=(1.25*minimum(p_sols), 1.5*maximum(p_sols)), label="Wave (c=$(c) m/s)")
	xlabel!("x(meters)")
	ylabel!("u(x,t=$(t_val))")
	title!("Wave at t=$(t_val) secs")
end

# ╔═╡ 5f829cd4-1bca-4473-a8a7-aca6a91cbd91
# ╠═╡ disabled = true
#=╠═╡
begin
	# Animation
    anim = @animate for it in range(start=1, stop=nt, length=min(nt, 60*(t_max)*10))
		t_val = round(it*(t_max)/(nt-1), digits=2)
	    plot(x, p_sols[Int(floor(it)),:], ylims=(minimum(p_sols), maximum(p_sols)), label="Wave (c=$(c) m/s)")
		xlabel!("x(meters)")
		ylabel!("u(x,t=$(t_val))")
		title!("Wave at t=$(t_val) secs")
    end

	# gif(anim, "anims/1D.gif", fps=60)
	mp4(anim, "anims/1D.mp4", fps=60)
end
  ╠═╡ =#

# ╔═╡ Cell order:
# ╟─206b9281-0aab-435d-8592-6e5169449832
# ╠═d5e4988a-09a0-11ef-2eb9-97277c169c65
# ╟─ed6ee4ce-fcd8-403a-9061-d688737f87be
# ╟─b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
# ╟─8ab95f04-2dad-4274-a7b5-69c23d5c2bb1
# ╟─ce9660d5-a355-464c-9142-1bfd084fc658
# ╟─241fc0e7-8b6c-4d61-b7c8-95f331469dc4
# ╟─ba518ad0-391d-45ca-8617-73fcc80bc8f4
# ╟─36b894da-0428-459c-92ff-0f3e34943769
# ╟─17fe1153-8514-4876-b5b7-7ba0bd2d1d89
# ╟─d9417b50-b232-47cc-8d78-506aca36fbca
# ╟─bd5d0e35-4bdd-488f-b8a5-bea04c27d764
# ╟─c6d58808-b224-46f9-91f4-965c02dd32ec
# ╟─2985d0e2-ae9b-41f3-996a-b8422c17ded8
# ╟─c26a9e81-4dab-40f4-a9c8-d20ef3a1ae4a
# ╟─14b955b6-f4ea-40ed-8340-996d75973135
# ╟─e1c436ce-f879-40f3-8d2e-bd8562f5abbd
# ╟─d2b91e55-c921-486d-92ce-9cd04c54f60b
# ╟─e090a86d-7ed8-4428-b777-06383d62f492
# ╟─2a602d52-d4f6-4894-b214-ba8a1c16fa89
# ╟─cf0e9ff0-dd7e-40e2-8544-e8d18ffa9e25
# ╟─5f829cd4-1bca-4473-a8a7-aca6a91cbd91

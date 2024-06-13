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
end;

# ╔═╡ b56e0a7c-c0a2-43f3-ba46-7f5c38d1cf97
md"""
# Problem Setting
"""

# ╔═╡ 8ab95f04-2dad-4274-a7b5-69c23d5c2bb1
md"""
## Space and Time Discretization
"""

# ╔═╡ ce9660d5-a355-464c-9142-1bfd084fc658
@bind space select(["Spatial domain size (in meters)", "Number of points in x/z dimension"], [PlutoUI.Select(100:100:1000), PlutoUI.Select(101:100:1001)])

# ╔═╡ 241fc0e7-8b6c-4d61-b7c8-95f331469dc4
begin
	x_max, nx = space;
	dx = (x_max)/(nx-1)

	z_max = x_max 
	nz = nx
	dz = dx;

	md"""
	**Space discretization step size (in meters): $(round(dx, digits=2))**
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
	dt = (t_max)/(nt-1);

	md"""
	**Time Discretization Step Size (in seconds): $(round(dt, digits=10))**
	"""
end

# ╔═╡ d9417b50-b232-47cc-8d78-506aca36fbca
md"""
## Source Function
"""

# ╔═╡ bd5d0e35-4bdd-488f-b8a5-bea04c27d764
@bind source select(["Pick a source function", "Source Frequency (in Hz)", "Source Initiation Time (in seconds)", "Plot Source Function"], [PlutoUI.Select(["Derivative of Gaussian", "Gaussian"]), PlutoUI.Select(5:5:30), PlutoUI.Slider(0.05:0.01:round(t_max/2, digits=2), show_value=true), PlutoUI.Select([true, false])])

# ╔═╡ c6d58808-b224-46f9-91f4-965c02dd32ec
begin
	func_name, f0, t0, plt = source
	
	# Space Discretization
	x = range(start=0, stop=x_max, length=nx)
	
	# Time Discretization
	t = range(start=0, stop=t_max, length=nt)
	
	# Source Function
	if func_name == "Derivative of Gaussian"
		src = -8. .* (t .- t0) .* f0 .* (exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2)) / (dx*dz)
	elseif func_name == "Gaussian"
		src = exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2) / (dx*dz)
	end

	if plt
		# Plotting source Function
	    p = scatter(t, src, title="Source Function", label="Source Function", markersize=1, dpi=1000)
		xlabel!("t")
		ylabel!("Amplitude")
	end
end

# ╔═╡ 2985d0e2-ae9b-41f3-996a-b8422c17ded8
md"""
# Finite Difference Simulation
"""

# ╔═╡ f18d9481-13ae-49eb-b6fb-5501d00c2f4a
@bind fdm_details select(["Source location from x=0 (in meters)",
	"Source location from z=0 (in meters)",
	"Medium velocity 1 (in m/s)", 
	"Medium velocity 2 (in m/s)", 
	"Interface location from z=0 (in meters)",
	"FD scheme (3 point or 5 point)", 
	"Boundary condition"],[
	PlutoUI.NumberField(25:x_max-25, default=x_max/2),
	PlutoUI.NumberField(25:z_max-25, default=z_max/2),
	PlutoUI.Select(100.:50:500.), 
	PlutoUI.Select(100.:50:500.), 
	PlutoUI.NumberField(25:z_max-25, default=z_max/4),
	PlutoUI.Select([3, 5]),
	PlutoUI.Select(["zero", "neumann", "absorbing"]),])

# ╔═╡ b9dd9d24-bff2-4e84-b188-6dec8e865d60
begin
	isrc_z, isrc_x, c1, c2, interface_loc, op, boundary = fdm_details

	# Velocity Profile
	c = ones(nx,nz)
	c[1:Int(interface_loc/dz), : ] .= c2
	c[Int(interface_loc/dz):nz , :] .= c1

	# Space Discretization
	s_details = Dict("s_min" => 0, "s_max" => x_max, "ns" => nx)

	# Time Discretization
	t_details = Dict("t_min" => 0, "t_max" => t_max, "nt" => nt)

	# Source Details
	src_details = Dict("isrc_x" => isrc_x, "isrc_z" => isrc_z, "src" => src)

	# Solution
	p_sols = wave_sim_2d(s_details, t_details, src_details, c, op, boundary);

	md"""
	### CFL criteria: $(round((maximum(c)*dt/dx)+(maximum(c)*dt/dz), digits=3))
	"""
end

# ╔═╡ dcb8b5fb-93a8-496c-ac67-39c9197a370a
begin
	# Plotting
	heatmap(x, x, c, dpi=1000)
	xlabel!("x(meters)")
	ylabel!("z(meters)")
	title!("Velocity profile")
end

# ╔═╡ edbd96e2-9f8b-4919-882c-a41c6970e3a7
@bind sim_details select(["Time",
	"Camera Axis 1",
	"Camera Axis 2"],[
	PlutoUI.Slider(range(start=1, stop=nt, length=min(nt, 60*(t_max)*10)), show_value=false),
	PlutoUI.Slider(range(start=30, stop=150, length=150-30+1), show_value=true, default=45),
	PlutoUI.Slider(range(start=1, stop=89, length=89), show_value=true, default=30),])

# ╔═╡ cf0e9ff0-dd7e-40e2-8544-e8d18ffa9e25
begin
	# Plotting
	it, xcam, ycam = sim_details
	t_val = round(it*(t_max)/(nt-1), digits=2)
	plot(x, x, p_sols[Int(floor(it)),:,:], st=:surface, zshowaxis=false, clims=(1e-1*minimum(p_sols), 1e-1*maximum(p_sols)), zlims=(1*minimum(p_sols), 1*maximum(p_sols)), dpi=300, legend = :none, camera=(xcam, ycam), size=(800,800))
	xlabel!("x(meters)")
	ylabel!("z(meters)")
	title!("Wave at t=$(t_val) secs")
end

# ╔═╡ 8eb7e4f2-d39a-43ba-9538-00868c40969d
# ╠═╡ disabled = true
#=╠═╡
begin
	# Animation
    anim = @animate for it in range(start=1, stop=nt, length=min(nt, 60*(t_max)*10))
		t_val__ = round(it*(t_max)/(nt-1), digits=2)
		plot(x, x, p_sols[Int(floor(it)),:,:], st=:surface, clims=(1e-1*minimum(p_sols), 1e-1*maximum(p_sols)), zlims=(1*minimum(p_sols), 1*maximum(p_sols)), dpi=1000)
		xlabel!("x(meters)")
		ylabel!("z(meters)")
		title!("Wave at t=$(t_val__) secs")
    end

	# gif(anim, "anims/1D.gif", fps=60)
	mp4(anim, "anims/2D_tsunami.mp4", fps=60)
end
  ╠═╡ =#

# ╔═╡ 3a8559f5-ac3d-444a-82ec-0182743598ca
# ╠═╡ disabled = true
#=╠═╡
begin
	# Animation
    anim_ = @animate for it in range(start=1, stop=nt, length=min(nt, 60*(t_max)*10))
		t_val___ = round(it*(t_max)/(nt-1), digits=2)
		heatmap(x, x, p_sols[Int(floor(it)),:,:], clims=(10^(-1)*minimum(p_sols), 10^(-1)*maximum(p_sols)), dpi=1000)
		xlabel!("x(meters)")
		ylabel!("z(meters)")
		title!("Wave at t=$(t_val___) secs")
    end

	# gif(anim_, "anims/1D.gif", fps=60)
	mp4(anim_, "anims/2D_eq.mp4", fps=60)
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
# ╟─f18d9481-13ae-49eb-b6fb-5501d00c2f4a
# ╟─b9dd9d24-bff2-4e84-b188-6dec8e865d60
# ╟─dcb8b5fb-93a8-496c-ac67-39c9197a370a
# ╟─edbd96e2-9f8b-4919-882c-a41c6970e3a7
# ╟─cf0e9ff0-dd7e-40e2-8544-e8d18ffa9e25
# ╟─8eb7e4f2-d39a-43ba-9538-00868c40969d
# ╟─3a8559f5-ac3d-444a-82ec-0182743598ca

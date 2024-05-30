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

# ╔═╡ 4c551534-1ea6-11ef-003c-d99308c82076
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

# ╔═╡ bbae74a4-76fc-470f-b612-eb71a8c8926d
begin
	# Space discretization
	x_min = 0
	x_max = 100
	nx = 501
	dx = (x_max-x_min)/(nx-1) 
	x = range(start=x_min, stop=x_max, length=nx)

	z_min = x_min
	z_max = x_max
	nz = nx
	dz = dx
	z = range(start=z_min, stop=z_max, length=nz)

	# Time discretization
	t_min = 0
	t_max = 1
	nt = 1001
	dt = (t_max-t_min)/(nt-1)
	t = range(start=t_min, stop=t_max, length=nt)

	# Source function
	f0 = 10
	t0 = 0.1
	isrc = [floor(nx/2), floor(nz/2)]
	src = -8. .* (t .- t0) .* f0 .* (exp.(-1.0 .* (4*f0)^2 .* (t .- t0).^2)) / (dx*dz)
end;

# ╔═╡ 9e1efb9b-dd09-4223-9ee4-87ae7d484d12
begin
	plt = scatter(t, src, title="Source Function", label="Source Function", markersize=1, dpi=1000)
	xlabel!("t")
	ylabel!("Amplitude")
end

# ╔═╡ 2b9377e5-8b1e-410d-a161-c7051e34fc24
function wave_sim_2d(nx, nz, nt, dx, dt, isrc, src)
	boundary = "absorbing"
	c = 100
	
	# Pressure Field at p(x,z,t)
	p = zeros(nx, nz)
	# Pressure Field at p(x,z,t-dt)
	p_prev = zeros(nx, nz)
	# Pressure Field at p(x,z,t+dt)
	p_next = zeros(nx, nz)
	# Solution at each step
	p_sols = zeros(nt, nx, nz)

	# Looping over time
    for it in range(start=2, stop=nt)
        # Looping over Space
        for ix in range(start=2, stop=nx-1)
			for iz in range(start=2, stop=nz-1)
	            # Evaluating 2nd derivative wrt x
	            d2p_dx2 = p[ix+1,iz] - 2*p[ix,iz] + p[ix-1,iz]
				# Evaluating 2nd derivative wrt z
	            d2p_dz2 = p[ix,iz+1] - 2*p[ix,iz] + p[ix,iz-1]
	            # Updating Solution
	            if ix == isrc[1] && iz == isrc[2]
	                p_next[ix,iz] = (c*dt/dx)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz] + dt^2 * src[it]
	            else
	                p_next[ix,iz] = (c*dt/dx)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz]
	            end
	        end
		end

        # Boundary conditions
        if boundary == "zero"
            p_next[1,:] .= 0
            p_next[nx,:] .= 0
			p_next[:,1] .= 0
            p_next[:,nz] .= 0
        elseif boundary == "neumann"
            p_next[1,:] .= p_next[2,:]
            p_next[nx,:] .= p_next[nx-1,:]
			p_next[:,1] .= p_next[:,2]
            p_next[:,nz] .= p_next[:,nz-1]
        elseif boundary == "absorbing"
            p_next[1,:] = p[2,:] + (c*dt-dx)/(c*dt+dx) * (p_next[2,:]-p[1,:])
            p_next[nx,:] = p[nx-1,:] + (c*dt-dx)/(c*dt+dx) * (p_next[nx-1,:]-p[nx,:])
			p_next[:,1] = p[:,2] + (c*dt-dx)/(c*dt+dx) * (p_next[:,2]-p[:,1])
            p_next[:,nz] = p[:,nz-1] + (c*dt-dx)/(c*dt+dx) * (p_next[:,nz-1]-p[:,nz])
        end

        # Current Sol becomes Previous Sol
        p_prev[:,:] = p[:,:]
        # Next Sol becomes Current Sol
        p[:,:] = p_next[:,:]
    
        # Storing solutions at each time step
        p_sols[it,:,:] = p[:,:]
    end
	return p_sols
end

# ╔═╡ 3660afaa-3bc8-4f82-b95b-aa978c94dfbf
begin
	p_sols = wave_sim_2d(nx, nz, nt, dx, dt, isrc, src)
	@bind it PlutoUI.Slider(range(start=1, stop=nt, length=min(nt, 60*(t_max)*10)), show_value=false)
end

# ╔═╡ 53f2ed3f-a37c-4022-8563-20273cf77cf1
begin
	# Plotting
	t_val = round(it*(t_max)/(nt-1), digits=2)
	plot(x, z, p_sols[Int(floor(it)),:,:], st=:surface, zlims=(1*minimum(p_sols), 1*maximum(p_sols)), dpi=1000)
	xlabel!("x(meters)")
	ylabel!("z(meters)")
	title!("Wave at t=$(t_val) secs")
end

# ╔═╡ 3e865de5-567e-4269-a698-99be1d14d10a
minimum(p_sols)

# ╔═╡ Cell order:
# ╠═4c551534-1ea6-11ef-003c-d99308c82076
# ╠═bbae74a4-76fc-470f-b612-eb71a8c8926d
# ╟─9e1efb9b-dd09-4223-9ee4-87ae7d484d12
# ╠═2b9377e5-8b1e-410d-a161-c7051e34fc24
# ╟─3660afaa-3bc8-4f82-b95b-aa978c94dfbf
# ╠═53f2ed3f-a37c-4022-8563-20273cf77cf1
# ╠═3e865de5-567e-4269-a698-99be1d14d10a

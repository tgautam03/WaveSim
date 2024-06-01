function wave_sim_2d(s_details::Dict{String, Int64}, t_details::Dict{String, Int64}, src_details::Dict{String, Any}, c::Matrix, op::Int64, boundary::String)
	# Space Discretization
	x_min = s_details["s_min"]
	x_max = s_details["s_max"]
	nx = s_details["ns"]
	dx = (x_max-x_min)/(nx-1) 
	x = range(start=x_min, stop=x_max, length=nx)

    z_min = x_min
	z_max = x_max
	nz = nx
	dz = dx
	z = range(start=z_min, stop=z_max, length=nz)

	# Time Discretization
	t_min = t_details["t_min"]
	t_max = t_details["t_max"]
	nt = t_details["nt"]
	dt = (t_max-t_min)/(nt-1)
	t = range(start=t_min, stop=t_max, length=nt)
    
    # Source Function
	isrc_x = Int(floor(src_details["isrc_x"]/dx)+1) # Source Location in x
    isrc_z = Int(floor(src_details["isrc_z"]/dz)+1) # Source Location in z
    isrc = [isrc_x, isrc_z]
	src = src_details["src"] # Source Function
	
    ########################################################
    ############ Finite Difference Simulation ##############
    ########################################################
    # O(dt^2 + dx^2)
    if op == 3 
        p_sols = point_stencil_3(src, isrc, c, dx, nx, dz, nz, dt, nt, boundary)
    # O(dt^2 + dx^4)
    elseif op == 5 
        p_sols = point_stencil_5(src, isrc, c, dx, nx, dz, nz, dt, nt, boundary)
    end

    return p_sols
end

function point_stencil_3(src, isrc, c, dx, nx, dz, nz, dt, nt, boundary)
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
	            d2p_dx2 = (p[ix+1,iz] - 2*p[ix,iz] + p[ix-1,iz])/(dx^2)
				# Evaluating 2nd derivative wrt z
	            d2p_dz2 = (p[ix,iz+1] - 2*p[ix,iz] + p[ix,iz-1])/(dz^2)
	            # Updating Solution
	            if ix == isrc[1] && iz == isrc[2]
	                p_next[ix,iz] = (c[ix,iz]*dt)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz] + dt^2 * src[it]
	            else
	                p_next[ix,iz] = (c[ix,iz]*dt)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz]
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
            p_next[1,:] = p[2,:] + (c[1,:]*dt .- dx)/(c[1,:]*dt .+ dx) * (p_next[2,:]-p[1,:])
            p_next[nx,:] = p[nx-1,:] + (c[nx,:]*dt .- dx)/(c[nx,:]*dt .+ dx) * (p_next[nx-1,:]-p[nx,:])
			p_next[:,1] = p[:,2] + (c[:,1]*dt .- dx)/(c[:,1]*dt .+ dx) * (p_next[:,2]-p[:,1])
            p_next[:,nz] = p[:,nz-1] + (c[:,nz]*dt .- dx)/(c[:,nz]*dt .+ dx) * (p_next[:,nz-1]-p[:,nz])
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

function point_stencil_5(src, isrc, c, dx, nx, dz, nz, dt, nt, boundary)
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
        for ix in range(start=3, stop=nx-2)
			for iz in range(start=3, stop=nz-2)
	            # Evaluating 2nd derivative wrt x
	            d2p_dx2 = (-1/12 * p[ix+2,iz] + 4/3  * p[ix+1,iz] - 5/2 * p[ix,iz] +4/3  * p[ix-1,iz] - 1/12 * p[ix-2,iz])/(dx^2)
				# Evaluating 2nd derivative wrt z
	            d2p_dz2 = (-1/12 * p[ix, iz+2] + 4/3  * p[ix,iz+1] - 5/2 * p[ix,iz] +4/3  * p[ix,iz-1] - 1/12 * p[ix,iz-2])/(dz^2)
	            # Updating Solution
	            if ix == isrc[1] && iz == isrc[2]
	                p_next[ix,iz] = (c[ix,iz]*dt)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz] + dt^2 * src[it]
	            else
	                p_next[ix,iz] = (c[ix,iz]*dt)^2 * (d2p_dx2 + d2p_dz2) + 2*p[ix,iz] - p_prev[ix,iz]
	            end
	        end
		end

        # Boundary conditions
        if boundary == "zero"
            p_next[1:2,:] .= 0
            p_next[nx-1:nx,:] .= 0
			p_next[:,1:2] .= 0
            p_next[:,nz-1:nz] .= 0
        elseif boundary == "neumann"
            p_next[1:2,:] .= reshape(p_next[3,:], (1, size(p_next)[2]))
            p_next[nx-1:nx,:] .= reshape(p_next[nx-2,:], (1, size(p_next)[2]))
			p_next[:,1:2] .= reshape(p_next[:,3], (size(p_next)[2], 1))
            p_next[:,nz-1:nz] .= reshape(p_next[:,nz-2], (size(p_next)[2], 1))
        elseif boundary == "absorbing"
            p_next[1:2,:] = p[2:3,:] + (c[1:2,:]*dt .- dx)/(c[1:2,:]*dt .+ dx) * (p_next[2:3,:]-p[1:2,:])
            p_next[nx-1:nx,:] = p[nx-2:nx-1,:] + (c[nx-1:nx,:]*dt .- dx)/(c[nx-1:nx,:]*dt .+ dx) * (p_next[nx-2:nx-1,:]-p[nx-1:nx,:])
			p_next[:,1:2] = p[:,2:3] + (c[:,1:2]*dt .- dx)/(c[:,1:2]*dt .+ dx) * (p_next[:,2:3]-p[:,1:2])
            p_next[:,nz-1:nz] = p[:,nz-2:nz-1] + (c[:,nz-1:nz]*dt .- dx)/(c[:,nz-1:nz]*dt .+ dx) * (p_next[:,nz-2:nz-1]-p[:,nz-1:nz])
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
